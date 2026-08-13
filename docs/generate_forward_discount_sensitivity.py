from __future__ import annotations

import csv
import os
import sys
from collections import Counter
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
SEXES = [Sex.FEMALE.value, Sex.MALE.value]
GRADES = [1, 2, 3, 4, 5]
TECHNIQUES = [t.value for t in BULKING_TECHNIQUES]
VOLUME_GRID = [0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
GRADE5_TORT_SWEEP = [1.5, 2.5, 3.5]
OBSTRUCTION_CEILING = 0.15
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
    site_exponent: float


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


def _simulate(point: SimPoint, cache: dict[SimPoint, dict]) -> dict:
    if point in cache:
        return cache[point]

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
        "mean_forward_flow_ml_s": float(result.mean_forward_flow_ml_s),
    }
    cache[point] = out
    return out


def _aggregate_candidate(
    *,
    age_group: str,
    initial_grade: int,
    bbd_profile: str,
    bbd_severity: float,
    technique: str,
    deflux_volume_ml: float,
    site_exponent: float,
    cache: dict[SimPoint, dict],
) -> dict:
    tort_values = GRADE5_TORT_SWEEP if initial_grade == 5 else [1.0]
    sims: list[dict] = []
    for sex in SEXES:
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
                        site_exponent=site_exponent,
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
        "n_sites": _n_sites(technique),
        "reflux_fraction": mean(s["reflux_fraction"] for s in sims),
        "obstruction_index": mean(s["obstruction_index"] for s in sims),
        "severe_obstruction_any": int(any(s["severe_obstruction_any"] for s in sims)),
        "mean_forward_flow_ml_s": mean(s["mean_forward_flow_ml_s"] for s in sims),
    }


def _select_balanced(candidates: list[dict]) -> dict:
    safe = [
        c
        for c in candidates
        if int(c["severe_obstruction_any"]) == 0
        and float(c["obstruction_index"]) <= OBSTRUCTION_CEILING
    ]
    if safe:
        chosen = min(
            safe,
            key=lambda r: (
                float(r["reflux_fraction"]),
                float(r["obstruction_index"]),
                float(r["deflux_volume_ml"]),
                int(r["n_sites"]),
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
            int(r["n_sites"]),
            str(r["technique"]),
        ),
    )
    return {**fallback, "selection_mode": "fallback_min_obstruction", "safe_candidates_n": 0}


def _load_baseline_reference() -> dict[tuple[str, int, str], dict]:
    path = REVISION_OUT / "oi_threshold_sensitivity.csv"
    if not path.exists():
        raise FileNotFoundError(
            "Run generate_revision_sensitivity.py first so the OI=0.15 baseline exists."
        )
    rows = [
        row for row in _read_csv(path)
        if abs(float(row["obstruction_ceiling"]) - OBSTRUCTION_CEILING) < 1e-9
    ]
    return {
        (row["age_group"], int(row["initial_grade"]), row["bbd_severity_label"]): row
        for row in rows
    }


def generate_no_attenuation_recommendations() -> tuple[list[dict], list[dict]]:
    baseline = _load_baseline_reference()
    original_exponent = techniques.FORWARD_SITE_COUNT_ATTENUATION_EXPONENT
    no_discount_exponent = 0.0
    cache: dict[SimPoint, dict] = {}

    out_rows: list[dict] = []
    try:
        techniques.FORWARD_SITE_COUNT_ATTENUATION_EXPONENT = no_discount_exponent
        for age_group in AGE_GROUPS:
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
                            site_exponent=no_discount_exponent,
                            cache=cache,
                        )
                        for technique in TECHNIQUES
                        for volume in VOLUME_GRID
                    ]
                    chosen = _select_balanced(candidates)
                    ref = baseline[(age_group, grade, bbd_label)]
                    changed = (
                        str(chosen["technique"]) != str(ref["recommended_technique"])
                        or float(chosen["deflux_volume_ml"]) != float(ref["deflux_volume_ml"])
                    )
                    out_rows.append(
                        {
                            "age_group": age_group,
                            "initial_grade": grade,
                            "bbd_severity_label": bbd_label,
                            "baseline_site_exponent": f"{original_exponent:.2f}",
                            "no_discount_site_exponent": f"{no_discount_exponent:.2f}",
                            "baseline_technique": ref["recommended_technique"],
                            "baseline_volume_ml": ref["deflux_volume_ml"],
                            "baseline_reflux_fraction": ref["reflux_fraction"],
                            "baseline_obstruction_index": ref["obstruction_index"],
                            "no_discount_technique": chosen["technique"],
                            "no_discount_volume_ml": f"{float(chosen['deflux_volume_ml']):.2f}",
                            "no_discount_reflux_fraction": f"{float(chosen['reflux_fraction']):.4f}",
                            "no_discount_obstruction_index": f"{float(chosen['obstruction_index']):.4f}",
                            "no_discount_selection_mode": chosen["selection_mode"],
                            "no_discount_safe_candidates_n": int(chosen["safe_candidates_n"]),
                            "changed_vs_baseline": int(changed),
                            "sampling_protocol": "full 75-stratum balanced grid; sex averaged; grade5 tortuosity 1.5|2.5|3.5",
                        }
                    )
    finally:
        techniques.FORWARD_SITE_COUNT_ATTENUATION_EXPONENT = original_exponent

    baseline_tech_counts = Counter(row["baseline_technique"] for row in out_rows)
    no_discount_tech_counts = Counter(row["no_discount_technique"] for row in out_rows)
    summary_rows = [
        {
            "comparison": "baseline_0p5_vs_no_discount_0p0",
            "n_strata": len(out_rows),
            "changed_strata_count": sum(int(row["changed_vs_baseline"]) for row in out_rows),
            "baseline_double_hit_count": baseline_tech_counts.get("double_hit", 0),
            "baseline_double_hit_plus_sting_count": baseline_tech_counts.get("double_hit_plus_sting", 0),
            "baseline_hit_count": baseline_tech_counts.get("hit", 0),
            "baseline_sting_count": baseline_tech_counts.get("sting", 0),
            "no_discount_double_hit_count": no_discount_tech_counts.get("double_hit", 0),
            "no_discount_double_hit_plus_sting_count": no_discount_tech_counts.get("double_hit_plus_sting", 0),
            "no_discount_hit_count": no_discount_tech_counts.get("hit", 0),
            "no_discount_sting_count": no_discount_tech_counts.get("sting", 0),
            "sampling_protocol": "full 75-stratum balanced grid; OI<=0.15",
        }
    ]

    _write_csv(
        REVISION_OUT / "forward_site_discount_sensitivity.csv",
        [
            "age_group",
            "initial_grade",
            "bbd_severity_label",
            "baseline_site_exponent",
            "no_discount_site_exponent",
            "baseline_technique",
            "baseline_volume_ml",
            "baseline_reflux_fraction",
            "baseline_obstruction_index",
            "no_discount_technique",
            "no_discount_volume_ml",
            "no_discount_reflux_fraction",
            "no_discount_obstruction_index",
            "no_discount_selection_mode",
            "no_discount_safe_candidates_n",
            "changed_vs_baseline",
            "sampling_protocol",
        ],
        out_rows,
    )
    _write_csv(
        REVISION_OUT / "forward_site_discount_sensitivity_summary.csv",
        [
            "comparison",
            "n_strata",
            "changed_strata_count",
            "baseline_double_hit_count",
            "baseline_double_hit_plus_sting_count",
            "baseline_hit_count",
            "baseline_sting_count",
            "no_discount_double_hit_count",
            "no_discount_double_hit_plus_sting_count",
            "no_discount_hit_count",
            "no_discount_sting_count",
            "sampling_protocol",
        ],
        summary_rows,
    )
    return out_rows, summary_rows


def generate_matched_volume_demo() -> list[dict]:
    original_exponent = techniques.FORWARD_SITE_COUNT_ATTENUATION_EXPONENT
    rows: list[dict] = []
    cache: dict[SimPoint, dict] = {}
    try:
        for label, exponent in [("baseline_1_over_sqrt_n", original_exponent), ("no_site_discount", 0.0)]:
            techniques.FORWARD_SITE_COUNT_ATTENUATION_EXPONENT = exponent
            for technique in TECHNIQUES:
                point = SimPoint(
                    sex=Sex.FEMALE.value,
                    age_group=AgeGroup.TODDLER_18_24M.value,
                    initial_grade=5,
                    bbd_profile=BBDProfile.NONE.value,
                    bbd_severity=0.0,
                    technique=technique,
                    deflux_volume_ml=1.5,
                    tortuosity_index=2.5,
                    site_exponent=exponent,
                )
                sim = _simulate(point, cache)
                rows.append(
                    {
                        "scenario": label,
                        "site_exponent": f"{exponent:.2f}",
                        "age_group": point.age_group,
                        "sex": point.sex,
                        "initial_grade": point.initial_grade,
                        "bbd_profile": point.bbd_profile,
                        "tortuosity_index": f"{point.tortuosity_index:.1f}",
                        "technique": technique,
                        "n_sites": _n_sites(technique),
                        "deflux_volume_ml": f"{point.deflux_volume_ml:.2f}",
                        "reflux_fraction": f"{sim['reflux_fraction']:.4f}",
                        "obstruction_index": f"{sim['obstruction_index']:.4f}",
                        "mean_forward_flow_ml_s": f"{sim['mean_forward_flow_ml_s']:.4f}",
                    }
                )
    finally:
        techniques.FORWARD_SITE_COUNT_ATTENUATION_EXPONENT = original_exponent

    _write_csv(
        REVISION_OUT / "forward_site_discount_matched_volume_demo.csv",
        [
            "scenario",
            "site_exponent",
            "age_group",
            "sex",
            "initial_grade",
            "bbd_profile",
            "tortuosity_index",
            "technique",
            "n_sites",
            "deflux_volume_ml",
            "reflux_fraction",
            "obstruction_index",
            "mean_forward_flow_ml_s",
        ],
        rows,
    )
    return rows


def main() -> None:
    detail_rows, summary_rows = generate_no_attenuation_recommendations()
    demo_rows = generate_matched_volume_demo()
    print("Generated forward site-count discount sensitivity outputs:")
    print(f"- {REVISION_OUT / 'forward_site_discount_sensitivity.csv'} ({len(detail_rows)} rows)")
    print(f"- {REVISION_OUT / 'forward_site_discount_sensitivity_summary.csv'} ({len(summary_rows)} rows)")
    print(f"- {REVISION_OUT / 'forward_site_discount_matched_volume_demo.csv'} ({len(demo_rows)} rows)")


if __name__ == "__main__":
    main()
