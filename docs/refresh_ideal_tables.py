from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from statistics import mean

from vur_cfd.model import (
    AgeGroup,
    BBDProfile,
    BladderCapacityMethod,
    Sex,
    apply_initial_vur_grade,
    default_patient_from_literature,
    simulate_patient,
    with_bbd_modifiers,
    with_ureter_modifiers,
)
from vur_cfd.techniques import (
    BULKING_TECHNIQUES,
    TechniqueName,
    TechniquePlan,
    apply_technique,
    bulking_injection_layout,
)


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs"

AGE_GROUPS = [a.value for a in AgeGroup if a != AgeGroup.INFANT_0_12M]
SEXES = [Sex.FEMALE.value, Sex.MALE.value]
GRADES = [1, 2, 3, 4, 5]
TECHNIQUES = [t.value for t in BULKING_TECHNIQUES]

VOLUME_GRID = [0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
GRADE5_TORT_SWEEP_FINE = [1.5, 2.0, 2.5, 3.0, 3.5]
GRADE5_TORT_SWEEP_OPT = [1.5, 2.5, 3.5]
WORST_CASE_TORT_GRADE5 = 3.5
BALANCED_OBSTRUCTION_CEILING = 0.15

BBD_SCENARIOS = [
    ("none", BBDProfile.NONE.value, 0.0),
    ("moderate", BBDProfile.MIXED.value, 0.5),
    ("severe", BBDProfile.MIXED.value, 1.0),
]

SAMPLE_AGE_GROUPS = ["18_24m", "24_60m", "5_10y", "10_16y"]
SAMPLE_GRADES = [2, 3, 4]


@dataclass(frozen=True)
class SimPoint:
    sex: str
    age_group: str
    initial_grade: int
    technique: str
    deflux_volume_ml: float
    tortuosity_index: float
    bbd_profile: str
    bbd_severity: float


def _n_sites(technique: str) -> int:
    return len(bulking_injection_layout(TechniqueName(technique)))


def _simulate(point: SimPoint, cache: dict[tuple, dict]) -> dict:
    key = (
        point.sex,
        point.age_group,
        point.initial_grade,
        point.technique,
        point.deflux_volume_ml,
        point.tortuosity_index,
        point.bbd_profile,
        point.bbd_severity,
    )
    if key in cache:
        return cache[key]

    patient = default_patient_from_literature(
        age_group=AgeGroup(point.age_group),
        sex=Sex(point.sex),
        capacity_method=BladderCapacityMethod.KOFF,
        bladder_fill_fraction=1.0,
    )
    patient = apply_initial_vur_grade(patient, point.initial_grade)
    patient = with_ureter_modifiers(
        patient,
        tortuosity_index=point.tortuosity_index,
        compliance_factor=1.0,
        peristalsis_efficiency=1.0,
        bladder_fill_fraction=1.0,
    )
    patient = with_bbd_modifiers(
        patient,
        profile=BBDProfile(point.bbd_profile),
        severity=point.bbd_severity,
    )
    treated = apply_technique(
        patient,
        TechniquePlan(
            technique=TechniqueName(point.technique),
            placement_quality=0.75,
            deflux_volume_ml=point.deflux_volume_ml,
            tunnel_length_mm=0.0,
        ),
    )
    result = simulate_patient(treated, total_time_s=24.0, dt_s=0.05)

    out = {
        "pred_postop_reflux_fraction": float(result.reflux_fraction),
        "pred_postop_grade": int(result.vur_grade),
        "pred_postop_obstruction_index": float(result.obstruction_index),
        "pred_postop_severe_obstruction": int(bool(result.severe_obstruction)),
        "pred_postop_peak_pelvis_pa": float(result.peak_renal_pelvis_pressure_pa),
        # Absolute volumes are reported alongside reflux fraction because a
        # technique that suppresses antegrade flow can raise RF even when
        # absolute reflux falls (see generate_absolute_volume_table.py).
        "pred_postop_reflux_volume_ml": float(result.reflux_volume_ml),
        "pred_postop_antegrade_volume_ml": float(result.antegrade_volume_ml),
        "n_injection_sites": _n_sites(point.technique),
    }
    cache[key] = out
    return out


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _best_simple(rows: list[dict]) -> dict:
    return min(
        rows,
        key=lambda r: (
            int(r["pred_postop_severe_obstruction"]),
            int(r["pred_postop_grade"]),
            float(r["pred_postop_reflux_fraction"]),
            float(r["pred_postop_obstruction_index"]),
            int(r["n_injection_sites"]),
            str(r["technique"]),
        ),
    )


def _aggregate_candidate(
    age_group: str,
    initial_grade: int,
    bbd_profile: str,
    bbd_severity: float,
    technique: str,
    deflux_volume_ml: float,
    cache: dict[tuple, dict],
) -> dict:
    if initial_grade == 5:
        tort_sweep = GRADE5_TORT_SWEEP_OPT
        sweep_label = "1.5|2.5|3.5"
    else:
        tort_sweep = [1.0]
        sweep_label = "n/a"

    sims = []
    for sex in SEXES:
        for tort in tort_sweep:
            point = SimPoint(
                sex=sex,
                age_group=age_group,
                initial_grade=initial_grade,
                technique=technique,
                deflux_volume_ml=deflux_volume_ml,
                tortuosity_index=tort,
                bbd_profile=bbd_profile,
                bbd_severity=bbd_severity,
            )
            sims.append(_simulate(point, cache))

    grade_mean = mean(float(s["pred_postop_grade"]) for s in sims)
    row = {
        "age_group": age_group,
        "initial_grade": initial_grade,
        "bbd_profile": bbd_profile,
        "bbd_severity": bbd_severity,
        "technique": technique,
        "deflux_volume_ml": deflux_volume_ml,
        "n_injection_sites": _n_sites(technique),
        "pred_postop_grade_mean": grade_mean,
        "pred_postop_grade_rounded": int(round(grade_mean)),
        "pred_postop_reflux_fraction_mean": mean(
            float(s["pred_postop_reflux_fraction"]) for s in sims
        ),
        "pred_postop_obstruction_index_mean": mean(
            float(s["pred_postop_obstruction_index"]) for s in sims
        ),
        "pred_postop_peak_pelvis_pa_mean": mean(
            float(s["pred_postop_peak_pelvis_pa"]) for s in sims
        ),
        "pred_postop_reflux_volume_ml_mean": mean(
            float(s["pred_postop_reflux_volume_ml"]) for s in sims
        ),
        "pred_postop_antegrade_volume_ml_mean": mean(
            float(s["pred_postop_antegrade_volume_ml"]) for s in sims
        ),
        "severe_obstruction_any": int(
            any(int(s["pred_postop_severe_obstruction"]) == 1 for s in sims)
        ),
        "grade5_tortuosity_sweep": sweep_label,
    }
    return row


def _best_aggregated(rows: list[dict]) -> dict:
    return min(
        rows,
        key=lambda r: (
            int(r["severe_obstruction_any"]),
            float(r["pred_postop_obstruction_index_mean"]),
            float(r["pred_postop_reflux_fraction_mean"]),
            float(r["pred_postop_grade_mean"]),
            float(r["deflux_volume_ml"]),
            int(r["n_injection_sites"]),
            str(r["technique"]),
        ),
    )


def _best_balanced(rows: list[dict], obstruction_ceiling: float) -> dict:
    safe_rows = [
        r
        for r in rows
        if int(r["severe_obstruction_any"]) == 0
        and float(r["pred_postop_obstruction_index_mean"]) <= obstruction_ceiling
    ]

    if safe_rows:
        chosen = min(
            safe_rows,
            key=lambda r: (
                float(r["pred_postop_reflux_fraction_mean"]),
                float(r["pred_postop_grade_mean"]),
                float(r["pred_postop_obstruction_index_mean"]),
                float(r["deflux_volume_ml"]),
                int(r["n_injection_sites"]),
                str(r["technique"]),
            ),
        )
        return {
            **chosen,
            "balanced_mode": "within_ceiling",
            "balanced_obstruction_ceiling": obstruction_ceiling,
            "balanced_safe_candidates_n": len(safe_rows),
        }

    fallback = min(
        rows,
        key=lambda r: (
            int(r["severe_obstruction_any"]),
            float(r["pred_postop_obstruction_index_mean"]),
            float(r["pred_postop_reflux_fraction_mean"]),
            float(r["pred_postop_grade_mean"]),
            float(r["deflux_volume_ml"]),
            int(r["n_injection_sites"]),
            str(r["technique"]),
        ),
    )
    return {
        **fallback,
        "balanced_mode": "fallback_min_obstruction",
        "balanced_obstruction_ceiling": obstruction_ceiling,
        "balanced_safe_candidates_n": 0,
    }


def main() -> None:
    cache: dict[tuple, dict] = {}

    # 1) Dense grid by age/grade/technique (volume=1.0, tort=1.0, no BBD)
    grid_rows: list[dict] = []
    for sex in SEXES:
        for age_group in AGE_GROUPS:
            for initial_grade in GRADES:
                for technique in TECHNIQUES:
                    sim = _simulate(
                        SimPoint(
                            sex=sex,
                            age_group=age_group,
                            initial_grade=initial_grade,
                            technique=technique,
                            deflux_volume_ml=1.0,
                            tortuosity_index=1.0,
                            bbd_profile=BBDProfile.NONE.value,
                            bbd_severity=0.0,
                        ),
                        cache,
                    )
                    grid_rows.append(
                        {
                            "sex": sex,
                            "age_group": age_group,
                            "initial_grade": initial_grade,
                            "technique": technique,
                            **sim,
                        }
                    )
    _write_csv(
        OUT / "technique_grid_predictions_by_age_grade.csv",
        [
            "sex",
            "age_group",
            "initial_grade",
            "technique",
            "pred_postop_reflux_fraction",
            "pred_postop_grade",
            "pred_postop_obstruction_index",
            "pred_postop_severe_obstruction",
            "pred_postop_peak_pelvis_pa",
            "n_injection_sites",
        ],
        grid_rows,
    )

    # 2) Grade-5 volume+tortuosity sweep table (all techniques; no BBD)
    sweep_rows: list[dict] = []
    for sex in SEXES:
        for age_group in AGE_GROUPS:
            for initial_grade in GRADES:
                vol_list = VOLUME_GRID if initial_grade == 5 else [1.0]
                tort_list = GRADE5_TORT_SWEEP_FINE if initial_grade == 5 else [1.0]
                for technique in TECHNIQUES:
                    for volume in vol_list:
                        for tort in tort_list:
                            sim = _simulate(
                                SimPoint(
                                    sex=sex,
                                    age_group=age_group,
                                    initial_grade=initial_grade,
                                    technique=technique,
                                    deflux_volume_ml=volume,
                                    tortuosity_index=tort,
                                    bbd_profile=BBDProfile.NONE.value,
                                    bbd_severity=0.0,
                                ),
                                cache,
                            )
                            sweep_rows.append(
                                {
                                    "sex": sex,
                                    "age_group": age_group,
                                    "initial_grade": initial_grade,
                                    "technique": technique,
                                    "deflux_volume_ml": volume,
                                    "tortuosity_index": tort,
                                    **sim,
                                }
                            )
    _write_csv(
        OUT / "technique_grid_predictions_grade5_volume_tortuosity_sweep.csv",
        [
            "sex",
            "age_group",
            "initial_grade",
            "technique",
            "deflux_volume_ml",
            "tortuosity_index",
            "pred_postop_reflux_fraction",
            "pred_postop_grade",
            "pred_postop_obstruction_index",
            "pred_postop_severe_obstruction",
            "pred_postop_peak_pelvis_pa",
            "n_injection_sites",
        ],
        sweep_rows,
    )

    # 3) Best technique by sex/age/grade (no BBD; fixed volume 1.0, tort 1.0)
    best_simple_rows: list[dict] = []
    for sex in SEXES:
        for age_group in AGE_GROUPS:
            for initial_grade in GRADES:
                candidates = [
                    r
                    for r in grid_rows
                    if r["sex"] == sex
                    and r["age_group"] == age_group
                    and int(r["initial_grade"]) == initial_grade
                ]
                best_simple_rows.append(_best_simple(candidates))
    _write_csv(
        OUT / "ideal_technique_by_age_grade.csv",
        [
            "sex",
            "age_group",
            "initial_grade",
            "technique",
            "pred_postop_reflux_fraction",
            "pred_postop_grade",
            "pred_postop_obstruction_index",
            "pred_postop_severe_obstruction",
            "pred_postop_peak_pelvis_pa",
            "n_injection_sites",
        ],
        best_simple_rows,
    )
    _write_csv(
        OUT / "balanced_technique_by_age_grade.csv",
        [
            "sex",
            "age_group",
            "initial_grade",
            "technique",
            "pred_postop_reflux_fraction",
            "pred_postop_grade",
            "pred_postop_obstruction_index",
            "pred_postop_severe_obstruction",
            "pred_postop_peak_pelvis_pa",
            "n_injection_sites",
        ],
        best_simple_rows,
    )

    # 4) Best technique by sex/age/grade with grade-5 worst-case tortuosity
    # for scenario testing table.
    best_grade5_sweep_rows: list[dict] = []
    for sex in SEXES:
        for age_group in AGE_GROUPS:
            for initial_grade in GRADES:
                if initial_grade == 5:
                    candidates = [
                        r
                        for r in sweep_rows
                        if r["sex"] == sex
                        and r["age_group"] == age_group
                        and int(r["initial_grade"]) == 5
                        and float(r["tortuosity_index"]) == WORST_CASE_TORT_GRADE5
                    ]
                else:
                    candidates = [
                        r
                        for r in sweep_rows
                        if r["sex"] == sex
                        and r["age_group"] == age_group
                        and int(r["initial_grade"]) == initial_grade
                        and float(r["deflux_volume_ml"]) == 1.0
                        and float(r["tortuosity_index"]) == 1.0
                    ]
                best_grade5_sweep_rows.append(_best_simple(candidates))
    _write_csv(
        OUT / "ideal_technique_by_age_grade_grade5_sweep.csv",
        [
            "sex",
            "age_group",
            "initial_grade",
            "technique",
            "deflux_volume_ml",
            "tortuosity_index",
            "pred_postop_reflux_fraction",
            "pred_postop_grade",
            "pred_postop_obstruction_index",
            "pred_postop_severe_obstruction",
            "pred_postop_peak_pelvis_pa",
            "n_injection_sites",
        ],
        best_grade5_sweep_rows,
    )

    sample_rows = [
        r
        for r in best_simple_rows
        if r["age_group"] in SAMPLE_AGE_GROUPS and int(r["initial_grade"]) in SAMPLE_GRADES
    ]
    _write_csv(
        OUT / "ideal_technique_sample_patients.csv",
        [
            "sex",
            "age_group",
            "initial_grade",
            "technique",
            "pred_postop_reflux_fraction",
            "pred_postop_grade",
            "pred_postop_obstruction_index",
        ],
        sample_rows,
    )

    sample_grade5_rows = [
        r
        for r in best_grade5_sweep_rows
        if r["age_group"] in SAMPLE_AGE_GROUPS and int(r["initial_grade"]) in [2, 3, 4, 5]
    ]
    _write_csv(
        OUT / "ideal_technique_sample_patients_grade5_sweep.csv",
        [
            "sex",
            "age_group",
            "initial_grade",
            "technique",
            "deflux_volume_ml",
            "tortuosity_index",
            "pred_postop_reflux_fraction",
            "pred_postop_grade",
            "pred_postop_obstruction_index",
            "pred_postop_severe_obstruction",
            "pred_postop_peak_pelvis_pa",
            "n_injection_sites",
        ],
        sample_grade5_rows,
    )

    # 5) Main optimization table (age x grade x BBD), averaging over sex and
    # grade-5 tortuosity sweep.
    opt_rows: list[dict] = []
    opt_rows_balanced: list[dict] = []
    all_candidates_best_vol: list[dict] = []
    for age_group in AGE_GROUPS:
        for initial_grade in GRADES:
            for sev_label, bbd_profile, bbd_severity in BBD_SCENARIOS:
                candidates = []
                for technique in TECHNIQUES:
                    for volume in VOLUME_GRID:
                        cand = _aggregate_candidate(
                            age_group=age_group,
                            initial_grade=initial_grade,
                            bbd_profile=bbd_profile,
                            bbd_severity=bbd_severity,
                            technique=technique,
                            deflux_volume_ml=volume,
                            cache=cache,
                        )
                        cand["bbd_severity_label"] = sev_label
                        candidates.append(cand)
                # Per-technique best volume (for comparison table).
                by_tech: dict[str, list[dict]] = {}
                for c in candidates:
                    by_tech.setdefault(c["technique"], []).append(c)
                for tech_name, tech_cands in by_tech.items():
                    best_vol = min(
                        tech_cands,
                        key=lambda r: (
                            int(r["severe_obstruction_any"]),
                            float(r["pred_postop_obstruction_index_mean"]),
                            float(r["pred_postop_grade_mean"]),
                            float(r["pred_postop_reflux_fraction_mean"]),
                            float(r["deflux_volume_ml"]),
                        ),
                    )
                    all_candidates_best_vol.append(best_vol)

                opt_rows.append(_best_aggregated(candidates))
                opt_rows_balanced.append(
                    _best_balanced(candidates, BALANCED_OBSTRUCTION_CEILING)
                )

    _write_csv(
        OUT / "ideal_technique_by_age_grade_bbd_volume_grade5tortuosity.csv",
        [
            "age_group",
            "initial_grade",
            "bbd_profile",
            "bbd_severity",
            "bbd_severity_label",
            "technique",
            "deflux_volume_ml",
            "n_injection_sites",
            "pred_postop_grade_mean",
            "pred_postop_grade_rounded",
            "pred_postop_reflux_fraction_mean",
            "pred_postop_obstruction_index_mean",
            "pred_postop_peak_pelvis_pa_mean",
            "severe_obstruction_any",
            "grade5_tortuosity_sweep",
        ],
        opt_rows,
    )

    opt_rows_no_sweep = [
        {
            k: v
            for k, v in r.items()
            if k != "grade5_tortuosity_sweep"
        }
        for r in opt_rows
    ]
    _write_csv(
        OUT / "ideal_technique_by_age_grade_bbd_volume.csv",
        [
            "age_group",
            "initial_grade",
            "bbd_profile",
            "bbd_severity",
            "bbd_severity_label",
            "technique",
            "deflux_volume_ml",
            "n_injection_sites",
            "pred_postop_grade_mean",
            "pred_postop_grade_rounded",
            "pred_postop_reflux_fraction_mean",
            "pred_postop_obstruction_index_mean",
            "pred_postop_peak_pelvis_pa_mean",
            "severe_obstruction_any",
        ],
        opt_rows_no_sweep,
    )

    _write_csv(
        OUT / "ideal_technique_by_age_grade_bbd_volume_balanced.csv",
        [
            "age_group",
            "initial_grade",
            "bbd_profile",
            "bbd_severity",
            "bbd_severity_label",
            "technique",
            "deflux_volume_ml",
            "n_injection_sites",
            "pred_postop_grade_mean",
            "pred_postop_grade_rounded",
            "pred_postop_reflux_fraction_mean",
            "pred_postop_obstruction_index_mean",
            "pred_postop_peak_pelvis_pa_mean",
            "severe_obstruction_any",
            "balanced_mode",
            "balanced_obstruction_ceiling",
            "balanced_safe_candidates_n",
        ],
        opt_rows_balanced,
    )

    _write_csv(
        OUT / "ideal_technique_summary_clean.csv",
        [
            "age_group",
            "initial_grade",
            "bbd_severity_label",
            "technique",
            "deflux_volume_ml",
            "pred_postop_grade_rounded",
            "pred_postop_reflux_fraction_mean",
            "pred_postop_obstruction_index_mean",
        ],
        [
            {
                "age_group": r["age_group"],
                "initial_grade": r["initial_grade"],
                "bbd_severity_label": r["bbd_severity_label"],
                "technique": r["technique"],
                "deflux_volume_ml": r["deflux_volume_ml"],
                "pred_postop_grade_rounded": r["pred_postop_grade_rounded"],
                "pred_postop_reflux_fraction_mean": r["pred_postop_reflux_fraction_mean"],
                "pred_postop_obstruction_index_mean": r["pred_postop_obstruction_index_mean"],
            }
            for r in opt_rows
        ],
    )

    # 6) Manuscript tables.
    manuscript_all = [
        {
            "age_group": r["age_group"],
            "preop_vur_grade": int(r["initial_grade"]),
            "bbd_severity": r["bbd_severity_label"],
            "recommended_technique": r["technique"],
            "deflux_volume_ml": float(r["deflux_volume_ml"]),
            "predicted_reflux_fraction": f"{float(r['pred_postop_reflux_fraction_mean']):.3f}",
            "reflux_resolution_threshold_0p05": (
                "resolved"
                if float(r["pred_postop_reflux_fraction_mean"]) < 0.05
                else "persistent"
            ),
            "predicted_obstruction_index": f"{float(r['pred_postop_obstruction_index_mean']):.3f}",
            "predicted_postop_grade": int(r["pred_postop_grade_rounded"]),
        }
        for r in opt_rows
    ]
    _write_csv(
        OUT / "manuscript_table_ideal_all_ages.csv",
        [
            "age_group",
            "preop_vur_grade",
            "bbd_severity",
            "recommended_technique",
            "deflux_volume_ml",
            "predicted_reflux_fraction",
            "reflux_resolution_threshold_0p05",
            "predicted_obstruction_index",
            "predicted_postop_grade",
        ],
        manuscript_all,
    )

    manuscript_18_24 = [
        {
            "preop_vur_grade": int(r["initial_grade"]),
            "bbd_severity": r["bbd_severity_label"],
            "recommended_technique": r["technique"],
            "deflux_volume_ml": float(r["deflux_volume_ml"]),
            "predicted_reflux_fraction": f"{float(r['pred_postop_reflux_fraction_mean']):.3f}",
            "reflux_resolution_threshold_0p05": (
                "resolved"
                if float(r["pred_postop_reflux_fraction_mean"]) < 0.05
                else "persistent"
            ),
            "predicted_obstruction_index": f"{float(r['pred_postop_obstruction_index_mean']):.3f}",
            "predicted_postop_grade": int(r["pred_postop_grade_rounded"]),
        }
        for r in opt_rows
        if r["age_group"] == AgeGroup.TODDLER_18_24M.value
    ]
    manuscript_18_24.sort(key=lambda x: (x["preop_vur_grade"], x["bbd_severity"]))
    _write_csv(
        OUT / "manuscript_table_ideal_18_24m.csv",
        [
            "preop_vur_grade",
            "bbd_severity",
            "recommended_technique",
            "deflux_volume_ml",
            "predicted_reflux_fraction",
            "reflux_resolution_threshold_0p05",
            "predicted_obstruction_index",
            "predicted_postop_grade",
        ],
        manuscript_18_24,
    )

    tech_counts: dict[str, int] = {}
    vol_counts: dict[float, int] = {}
    for row in manuscript_all:
        tech = str(row["recommended_technique"])
        vol = float(row["deflux_volume_ml"])
        tech_counts[tech] = tech_counts.get(tech, 0) + 1
        vol_counts[vol] = vol_counts.get(vol, 0) + 1

    overview_rows = [{"metric": "n_total_strata", "value": len(manuscript_all)}]
    for tech in sorted(tech_counts):
        overview_rows.append(
            {
                "metric": f"technique_count_{tech}",
                "value": tech_counts[tech],
            }
        )
    for vol in sorted(vol_counts):
        overview_rows.append(
            {
                "metric": f"volume_count_{vol}_ml",
                "value": vol_counts[vol],
            }
        )
    _write_csv(
        OUT / "manuscript_table_optimization_overview.csv",
        ["metric", "value"],
        overview_rows,
    )

    manuscript_all_balanced = [
        {
            "age_group": r["age_group"],
            "preop_vur_grade": int(r["initial_grade"]),
            "bbd_severity": r["bbd_severity_label"],
            "recommended_technique": r["technique"],
            "deflux_volume_ml": float(r["deflux_volume_ml"]),
            "predicted_reflux_fraction": f"{float(r['pred_postop_reflux_fraction_mean']):.3f}",
            "reflux_resolution_threshold_0p05": (
                "resolved"
                if float(r["pred_postop_reflux_fraction_mean"]) < 0.05
                else "persistent"
            ),
            "predicted_obstruction_index": f"{float(r['pred_postop_obstruction_index_mean']):.3f}",
            "predicted_postop_grade": int(r["pred_postop_grade_rounded"]),
            "selection_mode": str(r["balanced_mode"]),
            "safe_candidates_n": int(r["balanced_safe_candidates_n"]),
        }
        for r in opt_rows_balanced
    ]
    _write_csv(
        OUT / "manuscript_table_ideal_all_ages_balanced.csv",
        [
            "age_group",
            "preop_vur_grade",
            "bbd_severity",
            "recommended_technique",
            "deflux_volume_ml",
            "predicted_reflux_fraction",
            "reflux_resolution_threshold_0p05",
            "predicted_obstruction_index",
            "predicted_postop_grade",
            "selection_mode",
            "safe_candidates_n",
        ],
        manuscript_all_balanced,
    )

    manuscript_18_24_balanced = [
        {
            "preop_vur_grade": int(r["initial_grade"]),
            "bbd_severity": r["bbd_severity_label"],
            "recommended_technique": r["technique"],
            "deflux_volume_ml": float(r["deflux_volume_ml"]),
            "predicted_reflux_fraction": f"{float(r['pred_postop_reflux_fraction_mean']):.3f}",
            "reflux_resolution_threshold_0p05": (
                "resolved"
                if float(r["pred_postop_reflux_fraction_mean"]) < 0.05
                else "persistent"
            ),
            "predicted_obstruction_index": f"{float(r['pred_postop_obstruction_index_mean']):.3f}",
            "predicted_postop_grade": int(r["pred_postop_grade_rounded"]),
            "selection_mode": str(r["balanced_mode"]),
            "safe_candidates_n": int(r["balanced_safe_candidates_n"]),
        }
        for r in opt_rows_balanced
        if r["age_group"] == AgeGroup.TODDLER_18_24M.value
    ]
    manuscript_18_24_balanced.sort(
        key=lambda x: (x["preop_vur_grade"], x["bbd_severity"])
    )
    _write_csv(
        OUT / "manuscript_table_ideal_18_24m_balanced.csv",
        [
            "preop_vur_grade",
            "bbd_severity",
            "recommended_technique",
            "deflux_volume_ml",
            "predicted_reflux_fraction",
            "reflux_resolution_threshold_0p05",
            "predicted_obstruction_index",
            "predicted_postop_grade",
            "selection_mode",
            "safe_candidates_n",
        ],
        manuscript_18_24_balanced,
    )

    tech_counts_balanced: dict[str, int] = {}
    vol_counts_balanced: dict[float, int] = {}
    fallback_n = 0
    for row in manuscript_all_balanced:
        tech = str(row["recommended_technique"])
        vol = float(row["deflux_volume_ml"])
        tech_counts_balanced[tech] = tech_counts_balanced.get(tech, 0) + 1
        vol_counts_balanced[vol] = vol_counts_balanced.get(vol, 0) + 1
        if row["selection_mode"] != "within_ceiling":
            fallback_n += 1

    overview_rows_balanced = [
        {"metric": "n_total_strata", "value": len(manuscript_all_balanced)},
        {
            "metric": "obstruction_ceiling",
            "value": BALANCED_OBSTRUCTION_CEILING,
        },
        {"metric": "fallback_strata_count", "value": fallback_n},
    ]
    for tech in sorted(tech_counts_balanced):
        overview_rows_balanced.append(
            {
                "metric": f"technique_count_{tech}",
                "value": tech_counts_balanced[tech],
            }
        )
    for vol in sorted(vol_counts_balanced):
        overview_rows_balanced.append(
            {
                "metric": f"volume_count_{vol}_ml",
                "value": vol_counts_balanced[vol],
            }
        )
    _write_csv(
        OUT / "manuscript_table_optimization_overview_balanced.csv",
        ["metric", "value"],
        overview_rows_balanced,
    )

    _write_csv(
        OUT / "technique_comparison_by_grade_bbd.csv",
        [
            "age_group",
            "initial_grade",
            "bbd_severity_label",
            "technique",
            "deflux_volume_ml",
            "pred_postop_reflux_fraction_mean",
            "pred_postop_obstruction_index_mean",
            "pred_postop_grade_rounded",
            "severe_obstruction_any",
        ],
        all_candidates_best_vol,
    )

    print("Refreshed ideal-technique tables from current physics.")
    print(f"- Cache points solved: {len(cache)}")
    print(f"- Wrote: {OUT}")


if __name__ == "__main__":
    main()
