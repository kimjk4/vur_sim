"""Fine-grained volume sweep to map the full dose-response curve and identify
the optimal obstruction index per stratum.

Produces three CSVs:
  1. optimal_obstruction_sweep.csv         — full sweep data
  2. optimal_obstruction_index_by_stratum.csv — one optimal point per stratum
  3. optimal_obstruction_index_summary.csv  — age-averaged summary per stratum
"""
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

VOLUME_SWEEP = [
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.55, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0,
]

GRADE5_TORT_SWEEP_OPT = [1.5, 2.5, 3.5]

BBD_SCENARIOS = [
    ("none", BBDProfile.NONE.value, 0.0),
    ("moderate", BBDProfile.MIXED.value, 0.5),
    ("severe", BBDProfile.MIXED.value, 1.0),
]

REFLUX_THRESHOLD = 0.05


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

    if point.deflux_volume_ml > 0.0:
        treated = apply_technique(
            patient,
            TechniquePlan(
                technique=TechniqueName(point.technique),
                placement_quality=0.75,
                deflux_volume_ml=point.deflux_volume_ml,
                tunnel_length_mm=0.0,
            ),
        )
    else:
        treated = patient

    result = simulate_patient(treated, total_time_s=24.0, dt_s=0.05)

    out = {
        "mean_forward_flow_ml_s": float(result.mean_forward_flow_ml_s),
        "reflux_fraction": float(result.reflux_fraction),
        "filling_reflux_fraction": float(result.filling_reflux_fraction),
        "voiding_reflux_fraction": float(result.voiding_reflux_fraction),
        "obstruction_index": float(result.obstruction_index),
        "peak_pelvis_pa": float(result.peak_renal_pelvis_pressure_pa),
        "severe_obstruction": int(bool(result.severe_obstruction)),
    }
    cache[key] = out
    return out


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _aggregate_sweep_point(
    age_group: str,
    initial_grade: int,
    bbd_profile: str,
    bbd_severity: float,
    technique: str,
    deflux_volume_ml: float,
    cache: dict[tuple, dict],
) -> dict:
    """Average simulation across sex and grade-5 tortuosity sweep."""
    tort_sweep = GRADE5_TORT_SWEEP_OPT if initial_grade == 5 else [1.0]

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

    return {
        "age_group": age_group,
        "initial_grade": initial_grade,
        "bbd_severity_label": "",  # filled by caller
        "technique": technique,
        "deflux_volume_ml": deflux_volume_ml,
        "mean_forward_flow_ml_s": mean(s["mean_forward_flow_ml_s"] for s in sims),
        "reflux_fraction": mean(s["reflux_fraction"] for s in sims),
        "filling_reflux_fraction": mean(s["filling_reflux_fraction"] for s in sims),
        "voiding_reflux_fraction": mean(s["voiding_reflux_fraction"] for s in sims),
        "obstruction_index": mean(s["obstruction_index"] for s in sims),
        "peak_pelvis_pa": mean(s["peak_pelvis_pa"] for s in sims),
        "severe_obstruction_any": int(
            any(int(s["severe_obstruction"]) == 1 for s in sims)
        ),
    }


def _find_optimal_point(rows_for_stratum: list[dict]) -> dict:
    """Find the minimum volume where reflux_fraction < REFLUX_THRESHOLD.

    If no volume achieves the threshold, fall back to the volume with minimum
    reflux fraction.
    """
    sorted_rows = sorted(rows_for_stratum, key=lambda r: float(r["deflux_volume_ml"]))

    for row in sorted_rows:
        if float(row["reflux_fraction"]) < REFLUX_THRESHOLD:
            return {
                "age_group": row["age_group"],
                "initial_grade": row["initial_grade"],
                "bbd_severity_label": row["bbd_severity_label"],
                "technique": row["technique"],
                "optimal_volume_ml": float(row["deflux_volume_ml"]),
                "optimal_obstruction_index": float(row["obstruction_index"]),
                "reflux_at_optimal": float(row["reflux_fraction"]),
                "forward_flow_at_optimal": float(row["mean_forward_flow_ml_s"]),
                "selection_mode": "within_threshold",
            }

    best = min(sorted_rows, key=lambda r: float(r["reflux_fraction"]))
    return {
        "age_group": best["age_group"],
        "initial_grade": best["initial_grade"],
        "bbd_severity_label": best["bbd_severity_label"],
        "technique": best["technique"],
        "optimal_volume_ml": float(best["deflux_volume_ml"]),
        "optimal_obstruction_index": float(best["obstruction_index"]),
        "reflux_at_optimal": float(best["reflux_fraction"]),
        "forward_flow_at_optimal": float(best["mean_forward_flow_ml_s"]),
        "selection_mode": "fallback_min_reflux",
    }


def main() -> None:
    cache: dict[tuple, dict] = {}

    # 1) Full volume sweep across all strata.
    sweep_rows: list[dict] = []
    print("Running volume sweep simulations...")

    total = len(AGE_GROUPS) * len(GRADES) * len(BBD_SCENARIOS) * len(TECHNIQUES) * len(VOLUME_SWEEP)
    done = 0

    for age_group in AGE_GROUPS:
        for initial_grade in GRADES:
            for sev_label, bbd_profile, bbd_severity in BBD_SCENARIOS:
                for technique in TECHNIQUES:
                    for volume in VOLUME_SWEEP:
                        row = _aggregate_sweep_point(
                            age_group=age_group,
                            initial_grade=initial_grade,
                            bbd_profile=bbd_profile,
                            bbd_severity=bbd_severity,
                            technique=technique,
                            deflux_volume_ml=volume,
                            cache=cache,
                        )
                        row["bbd_severity_label"] = sev_label
                        sweep_rows.append(row)
                        done += 1
                        if done % 500 == 0:
                            print(f"  {done}/{total} sweep points computed...")

    sweep_fields = [
        "age_group",
        "initial_grade",
        "bbd_severity_label",
        "technique",
        "deflux_volume_ml",
        "mean_forward_flow_ml_s",
        "reflux_fraction",
        "filling_reflux_fraction",
        "voiding_reflux_fraction",
        "obstruction_index",
        "peak_pelvis_pa",
        "severe_obstruction_any",
    ]
    _write_csv(OUT / "optimal_obstruction_sweep.csv", sweep_fields, sweep_rows)
    print(f"Wrote optimal_obstruction_sweep.csv ({len(sweep_rows)} rows)")

    # 2) Find optimal point per (age, grade, BBD, technique).
    optimal_rows: list[dict] = []
    from collections import defaultdict

    strata: dict[tuple, list[dict]] = defaultdict(list)
    for row in sweep_rows:
        key = (
            row["age_group"],
            row["initial_grade"],
            row["bbd_severity_label"],
            row["technique"],
        )
        strata[key].append(row)

    for _key, stratum_rows in strata.items():
        optimal_rows.append(_find_optimal_point(stratum_rows))

    optimal_fields = [
        "age_group",
        "initial_grade",
        "bbd_severity_label",
        "technique",
        "optimal_volume_ml",
        "optimal_obstruction_index",
        "reflux_at_optimal",
        "forward_flow_at_optimal",
        "selection_mode",
    ]
    _write_csv(
        OUT / "optimal_obstruction_index_by_stratum.csv",
        optimal_fields,
        optimal_rows,
    )
    print(f"Wrote optimal_obstruction_index_by_stratum.csv ({len(optimal_rows)} rows)")

    # 3) Age-averaged summary per (technique, grade, BBD).
    summary_strata: dict[tuple, list[dict]] = defaultdict(list)
    for row in optimal_rows:
        key = (row["technique"], row["initial_grade"], row["bbd_severity_label"])
        summary_strata[key].append(row)

    summary_rows: list[dict] = []
    for (technique, grade, bbd), group in summary_strata.items():
        summary_rows.append(
            {
                "technique": technique,
                "initial_grade": grade,
                "bbd_severity_label": bbd,
                "optimal_volume_ml": mean(r["optimal_volume_ml"] for r in group),
                "optimal_obstruction_index": mean(
                    r["optimal_obstruction_index"] for r in group
                ),
                "reflux_at_optimal": mean(r["reflux_at_optimal"] for r in group),
                "forward_flow_at_optimal": mean(
                    r["forward_flow_at_optimal"] for r in group
                ),
                "n_age_groups": len(group),
            }
        )

    summary_fields = [
        "technique",
        "initial_grade",
        "bbd_severity_label",
        "optimal_volume_ml",
        "optimal_obstruction_index",
        "reflux_at_optimal",
        "forward_flow_at_optimal",
        "n_age_groups",
    ]
    _write_csv(
        OUT / "optimal_obstruction_index_summary.csv",
        summary_fields,
        summary_rows,
    )
    print(f"Wrote optimal_obstruction_index_summary.csv ({len(summary_rows)} rows)")
    print(f"Cache points solved: {len(cache)}")


if __name__ == "__main__":
    main()
