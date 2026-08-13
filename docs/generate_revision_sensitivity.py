from __future__ import annotations

import csv
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = str(PROJECT_ROOT / "outputs" / ".mplconfig")

import vur_cfd.techniques as techniques
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
from vur_cfd.techniques import BULKING_TECHNIQUES, TechniqueName, TechniquePlan, apply_technique


OUT = PROJECT_ROOT / "outputs"
REVISION_OUT = OUT / "revision"
REVISION_OUT.mkdir(parents=True, exist_ok=True)

AGE_GROUPS = [a.value for a in AgeGroup if a != AgeGroup.INFANT_0_12M]
MOUND_AGE_GROUPS = [AgeGroup.TODDLER_18_24M.value]
GRADES = [1, 2, 3, 4, 5]
SEXES = [Sex.FEMALE.value, Sex.MALE.value]
MOUND_REFERENCE_SEXES = [Sex.FEMALE.value]
TECHNIQUES = [t.value for t in BULKING_TECHNIQUES]
PRIMARY_VOLUME_GRID = [0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
GRADE5_TORT_SWEEP = [1.5, 2.5, 3.5]
MOUND_GRADE5_TORT_SWEEP = [2.5]
OI_THRESHOLDS = [0.10, 0.15, 0.20, 0.25]
MOUND_V50_VALUES = [0.25, 0.35, 0.45]
MOUND_K_VALUES = [8.0, 12.0, 16.0]
BBD_SCENARIOS = [
    ("none", BBDProfile.NONE.value, 0.0),
    ("moderate", BBDProfile.MIXED.value, 0.5),
    ("severe", BBDProfile.MIXED.value, 1.0),
]


@dataclass(frozen=True)
class SimPoint:
    sex: str
    age_group: str
    initial_grade: int
    bbd_profile: str
    bbd_severity: float
    technique: str
    deflux_volume_ml: float
    tortuosity_index: float


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _n_sites(technique: str) -> int:
    return len(techniques.bulking_injection_layout(TechniqueName(technique)))


def _select_balanced(candidates: list[dict], obstruction_ceiling: float) -> dict:
    safe = [
        c
        for c in candidates
        if int(float(c.get("severe_obstruction_any", 0))) == 0
        and float(c["obstruction_index"]) <= obstruction_ceiling
    ]
    if safe:
        chosen = min(
            safe,
            key=lambda r: (
                float(r["reflux_fraction"]),
                float(r["obstruction_index"]),
                float(r["deflux_volume_ml"]),
                _n_sites(str(r["technique"])),
                str(r["technique"]),
            ),
        )
        return {**chosen, "selection_mode": "within_ceiling", "safe_candidates_n": len(safe)}

    fallback = min(
        candidates,
        key=lambda r: (
            float(r["obstruction_index"]),
            float(r["reflux_fraction"]),
            float(r["deflux_volume_ml"]),
            _n_sites(str(r["technique"])),
            str(r["technique"]),
        ),
    )
    return {**fallback, "selection_mode": "fallback_min_obstruction", "safe_candidates_n": 0}


def generate_oi_threshold_sensitivity() -> list[dict]:
    sweep_path = OUT / "optimal_obstruction_sweep.csv"
    rows = _read_csv(sweep_path)
    primary = [
        r
        for r in rows
        if float(r["deflux_volume_ml"]) in PRIMARY_VOLUME_GRID
        and int(r["initial_grade"]) in GRADES
    ]

    by_stratum: dict[tuple[str, int, str], list[dict]] = defaultdict(list)
    for row in primary:
        key = (row["age_group"], int(row["initial_grade"]), row["bbd_severity_label"])
        by_stratum[key].append(row)

    reference: dict[tuple[str, int, str], dict] = {}
    for key, candidates in by_stratum.items():
        reference[key] = _select_balanced(candidates, 0.15)

    out_rows: list[dict] = []
    for threshold in OI_THRESHOLDS:
        for key in sorted(by_stratum):
            age_group, grade, bbd = key
            chosen = _select_balanced(by_stratum[key], threshold)
            ref = reference[key]
            changed = (
                str(chosen["technique"]) != str(ref["technique"])
                or float(chosen["deflux_volume_ml"]) != float(ref["deflux_volume_ml"])
            )
            out_rows.append(
                {
                    "obstruction_ceiling": f"{threshold:.2f}",
                    "age_group": age_group,
                    "initial_grade": grade,
                    "bbd_severity_label": bbd,
                    "recommended_technique": chosen["technique"],
                    "deflux_volume_ml": f"{float(chosen['deflux_volume_ml']):.2f}",
                    "reflux_fraction": f"{float(chosen['reflux_fraction']):.4f}",
                    "obstruction_index": f"{float(chosen['obstruction_index']):.4f}",
                    "selection_mode": chosen["selection_mode"],
                    "safe_candidates_n": int(chosen["safe_candidates_n"]),
                    "reference_0p15_technique": ref["technique"],
                    "reference_0p15_volume_ml": f"{float(ref['deflux_volume_ml']):.2f}",
                    "changed_vs_0p15": int(changed),
                }
            )

    _write_csv(
        REVISION_OUT / "oi_threshold_sensitivity.csv",
        [
            "obstruction_ceiling",
            "age_group",
            "initial_grade",
            "bbd_severity_label",
            "recommended_technique",
            "deflux_volume_ml",
            "reflux_fraction",
            "obstruction_index",
            "selection_mode",
            "safe_candidates_n",
            "reference_0p15_technique",
            "reference_0p15_volume_ml",
            "changed_vs_0p15",
        ],
        out_rows,
    )
    return out_rows


def _simulate(point: SimPoint, cache: dict[tuple, dict]) -> dict:
    key = (
        techniques.MOUND_EFFICACY_MIDPOINT_ML,
        techniques.MOUND_EFFICACY_STEEPNESS,
        point,
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
        "reflux_fraction": float(result.reflux_fraction),
        "obstruction_index": float(result.obstruction_index),
        "severe_obstruction_any": int(bool(result.severe_obstruction)),
    }
    cache[key] = out
    return out


def _aggregate_candidate(
    *,
    age_group: str,
    initial_grade: int,
    bbd_profile: str,
    bbd_severity: float,
    technique: str,
    deflux_volume_ml: float,
    cache: dict[tuple, dict],
) -> dict:
    tort_values = MOUND_GRADE5_TORT_SWEEP if initial_grade == 5 else [1.0]
    sims: list[dict] = []
    for sex in MOUND_REFERENCE_SEXES:
        for tort in tort_values:
            sims.append(
                _simulate(
                    SimPoint(
                        sex=sex,
                        age_group=age_group,
                        initial_grade=initial_grade,
                        bbd_profile=bbd_profile,
                        bbd_severity=bbd_severity,
                        technique=technique,
                        deflux_volume_ml=deflux_volume_ml,
                        tortuosity_index=tort,
                    ),
                    cache,
                )
            )
    return {
        "age_group": age_group,
        "initial_grade": initial_grade,
        "bbd_profile": bbd_profile,
        "bbd_severity": bbd_severity,
        "technique": technique,
        "deflux_volume_ml": deflux_volume_ml,
        "reflux_fraction": mean(s["reflux_fraction"] for s in sims),
        "obstruction_index": mean(s["obstruction_index"] for s in sims),
        "severe_obstruction_any": int(any(s["severe_obstruction_any"] for s in sims)),
    }


def _select_all_strata_for_current_sigmoid(cache: dict[tuple, dict]) -> list[dict]:
    selections: list[dict] = []
    for age_group in MOUND_AGE_GROUPS:
        for grade in GRADES:
            for bbd_label, bbd_profile, bbd_severity in BBD_SCENARIOS:
                candidates = [
                    _aggregate_candidate(
                        age_group=age_group,
                        initial_grade=grade,
                        bbd_profile=bbd_profile,
                        bbd_severity=bbd_severity,
                        technique=technique,
                        deflux_volume_ml=volume,
                        cache=cache,
                    )
                    for technique in TECHNIQUES
                    for volume in PRIMARY_VOLUME_GRID
                ]
                chosen = _select_balanced(candidates, 0.15)
                selections.append(
                    {
                        "age_group": age_group,
                        "initial_grade": grade,
                        "bbd_severity_label": bbd_label,
                        "recommended_technique": chosen["technique"],
                        "deflux_volume_ml": float(chosen["deflux_volume_ml"]),
                        "reflux_fraction": float(chosen["reflux_fraction"]),
                        "obstruction_index": float(chosen["obstruction_index"]),
                        "selection_mode": chosen["selection_mode"],
                    }
                )
    return selections


def generate_mound_sigmoid_sensitivity() -> tuple[list[dict], list[dict]]:
    original_v50 = techniques.MOUND_EFFICACY_MIDPOINT_ML
    original_k = techniques.MOUND_EFFICACY_STEEPNESS
    cache: dict[tuple, dict] = {}

    all_selections: dict[tuple[float, float], list[dict]] = {}
    try:
        for v50 in MOUND_V50_VALUES:
            for k in MOUND_K_VALUES:
                techniques.MOUND_EFFICACY_MIDPOINT_ML = v50
                techniques.MOUND_EFFICACY_STEEPNESS = k
                all_selections[(v50, k)] = _select_all_strata_for_current_sigmoid(cache)
    finally:
        techniques.MOUND_EFFICACY_MIDPOINT_ML = original_v50
        techniques.MOUND_EFFICACY_STEEPNESS = original_k

    baseline_key = (0.35, 12.0)
    baseline = {
        (r["age_group"], r["initial_grade"], r["bbd_severity_label"]): r
        for r in all_selections[baseline_key]
    }

    summary_rows: list[dict] = []
    detail_rows: list[dict] = []
    for key in sorted(all_selections):
        v50, k = key
        rows = all_selections[key]
        tech_counts = Counter(str(r["recommended_technique"]) for r in rows)
        volume_counts = Counter(f"{float(r['deflux_volume_ml']):.1f}" for r in rows)
        changed = 0
        changed_technique = 0
        changed_volume = 0
        for row in rows:
            stratum_key = (row["age_group"], row["initial_grade"], row["bbd_severity_label"])
            ref = baseline[stratum_key]
            tech_changed = row["recommended_technique"] != ref["recommended_technique"]
            volume_changed = float(row["deflux_volume_ml"]) != float(ref["deflux_volume_ml"])
            any_changed = tech_changed or volume_changed
            changed += int(any_changed)
            changed_technique += int(tech_changed)
            changed_volume += int(volume_changed)
            detail_rows.append(
                {
                    "v50_ml": f"{v50:.2f}",
                    "k": f"{k:.1f}",
                    "age_group": row["age_group"],
                    "initial_grade": row["initial_grade"],
                    "bbd_severity_label": row["bbd_severity_label"],
                    "recommended_technique": row["recommended_technique"],
                    "deflux_volume_ml": f"{float(row['deflux_volume_ml']):.2f}",
                    "reflux_fraction": f"{float(row['reflux_fraction']):.4f}",
                    "obstruction_index": f"{float(row['obstruction_index']):.4f}",
                    "reference_technique": ref["recommended_technique"],
                    "reference_volume_ml": f"{float(ref['deflux_volume_ml']):.2f}",
                "changed_vs_reference": int(any_changed),
                    "sampling_protocol": "18_24m; female_reference; grade5_tortuosity_2.5",
                }
            )
        summary_rows.append(
            {
                "v50_ml": f"{v50:.2f}",
                "k": f"{k:.1f}",
                "n_strata": len(rows),
                "changed_strata_count": changed,
                "changed_technique_count": changed_technique,
                "changed_volume_count": changed_volume,
                "double_hit_count": tech_counts.get("double_hit", 0),
                "double_hit_plus_sting_count": tech_counts.get("double_hit_plus_sting", 0),
                "hit_count": tech_counts.get("hit", 0),
                "sting_count": tech_counts.get("sting", 0),
                "volume_count_0p6": volume_counts.get("0.6", 0),
                "volume_count_0p8": volume_counts.get("0.8", 0),
                "volume_count_1p0": volume_counts.get("1.0", 0),
                "volume_count_1p2": volume_counts.get("1.2", 0),
                "volume_count_1p5": volume_counts.get("1.5", 0),
                "volume_count_2p0": volume_counts.get("2.0", 0),
                "sampling_protocol": "18_24m; female_reference; grade5_tortuosity_2.5",
            }
        )

    _write_csv(
        REVISION_OUT / "mound_sigmoid_sensitivity.csv",
        [
            "v50_ml",
            "k",
            "n_strata",
            "changed_strata_count",
            "changed_technique_count",
            "changed_volume_count",
            "double_hit_count",
            "double_hit_plus_sting_count",
            "hit_count",
            "sting_count",
            "volume_count_0p6",
            "volume_count_0p8",
            "volume_count_1p0",
            "volume_count_1p2",
            "volume_count_1p5",
            "volume_count_2p0",
            "sampling_protocol",
        ],
        summary_rows,
    )
    _write_csv(
        REVISION_OUT / "mound_sigmoid_sensitivity_by_stratum.csv",
        [
            "v50_ml",
            "k",
            "age_group",
            "initial_grade",
            "bbd_severity_label",
            "recommended_technique",
            "deflux_volume_ml",
            "reflux_fraction",
            "obstruction_index",
            "reference_technique",
            "reference_volume_ml",
            "changed_vs_reference",
            "sampling_protocol",
        ],
        detail_rows,
    )
    return summary_rows, detail_rows


def generate_reviewer_evidence_table() -> list[dict]:
    rows = [
        {
            "reviewer_issue": "Model described as CFD despite reduced-order equations",
            "revision_action": "Reframe as reduced-order CFD-informed pressure-flow model",
            "evidence_or_artifact": "Revised manuscript Methods; Appendix parameter tables",
        },
        {
            "reviewer_issue": "Grade templates and parameters under-described",
            "revision_action": "Add explicit grade template, pressure, BBD, and technique coefficient tables",
            "evidence_or_artifact": "Revised Appendix Tables A1-A6",
        },
        {
            "reviewer_issue": "OI=0.15 not clinically validated",
            "revision_action": "State as modeling ceiling and add OI threshold sensitivity",
            "evidence_or_artifact": "outputs/revision/oi_threshold_sensitivity.csv",
        },
        {
            "reviewer_issue": "Sigmoid k and V50 not reported",
            "revision_action": "Report V50=0.35 mL/site and k=12.0; add sensitivity",
            "evidence_or_artifact": "outputs/revision/mound_sigmoid_sensitivity.csv",
        },
        {
            "reviewer_issue": "Bilateral effects appear numerically small",
            "revision_action": "Reframe bilateral result as relative independence under current architecture",
            "evidence_or_artifact": "Revised Results and Discussion",
        },
        {
            "reviewer_issue": "Sensitivity analysis only grade V",
            "revision_action": "Add grade-II/no-BBD tornado figure",
            "evidence_or_artifact": "outputs/revision/low_grade_sensitivity_tornado.*",
        },
        {
            "reviewer_issue": "STING grade-V curve appears ineffective and unsafe",
            "revision_action": "Explicitly acknowledge modeled grade-V STING failure and separate from broad clinical literature",
            "evidence_or_artifact": "Revised Results and response letter",
        },
        {
            "reviewer_issue": "No external validation",
            "revision_action": "Remove validation claims; label findings hypothesis-generating",
            "evidence_or_artifact": "Revised Abstract, Discussion, Limitations",
        },
    ]
    _write_csv(
        REVISION_OUT / "reviewer_response_evidence_table.csv",
        ["reviewer_issue", "revision_action", "evidence_or_artifact"],
        rows,
    )
    return rows


def main() -> None:
    oi_rows = generate_oi_threshold_sensitivity()
    mound_rows, mound_detail_rows = generate_mound_sigmoid_sensitivity()
    evidence_rows = generate_reviewer_evidence_table()
    print("Generated revision sensitivity outputs:")
    print(f"- {REVISION_OUT / 'oi_threshold_sensitivity.csv'} ({len(oi_rows)} rows)")
    print(f"- {REVISION_OUT / 'mound_sigmoid_sensitivity.csv'} ({len(mound_rows)} rows)")
    print(
        f"- {REVISION_OUT / 'mound_sigmoid_sensitivity_by_stratum.csv'} "
        f"({len(mound_detail_rows)} rows)"
    )
    print(f"- {REVISION_OUT / 'reviewer_response_evidence_table.csv'} ({len(evidence_rows)} rows)")


if __name__ == "__main__":
    main()
