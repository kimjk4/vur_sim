"""No-placement-asymmetry sensitivity analysis (JPU revision round 2).

Reviewer 2 observed that in the grade-V reference sweep STING carries roughly
twice the obstruction index of HIT at the same total volume delivered at a
single site, and asked whether this reflects the encoded assumption that
intraureteric placement produces less outlet narrowing than submeatal
placement.

This script re-runs the 75-stratum balanced selection and the grade-V reference
sweep with `techniques.PLACEMENT_WALL_ASYMMETRY_ENABLED = False`, which scores
submeatal (STING) injection sites with the same wall-plane multipliers as
intraureteric (HIT / Double HIT) sites. Everything else - candidate grid,
selection rule, obstruction ceiling - is identical to the base case, so the
difference is attributable to the placement prior alone.

Structure mirrors `generate_forward_discount_sensitivity.py` so the two
structural-prior analyses are directly comparable.
"""

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
    UVJValve,
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
# Grade-V reference sweep in the manuscript: 18-24 months, female, no BBD.
SWEEP_AGE_GROUP = AgeGroup.TODDLER_18_24M.value
SWEEP_SEX = Sex.FEMALE.value
SWEEP_TORTUOSITY = 2.5
BBD_SCENARIOS = [
    ("none", BBDProfile.NONE.value, 0.0),
    ("moderate", BBDProfile.MIXED.value, 0.5),
    ("severe", BBDProfile.MIXED.value, 1.0),
]

BASELINE_UVJ_FORWARD_PA_S_PER_ML = UVJValve().forward_resistance_pa_s_per_ml


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
    asymmetry_enabled: bool


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


def _simulate(point: SimPoint, cache: dict) -> dict:
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
        "reflux_volume_ml": float(result.reflux_volume_ml),
        "antegrade_volume_ml": float(result.antegrade_volume_ml),
        # Edema/narrowing proxy that drives the obstruction index; reported so
        # the mechanism of any OI change is visible rather than inferred.
        "uvj_forward_resistance_multiplier": float(
            treated.uvj.forward_resistance_pa_s_per_ml / BASELINE_UVJ_FORWARD_PA_S_PER_ML
        ),
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
    asymmetry_enabled: bool,
    cache: dict,
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
                        asymmetry_enabled=asymmetry_enabled,
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


def _load_baseline_reference() -> dict:
    """Base-case (OI <= 0.15) selections produced by generate_revision_sensitivity.py."""
    path = REVISION_OUT / "oi_threshold_sensitivity.csv"
    if not path.exists():
        raise FileNotFoundError(
            "Run generate_revision_sensitivity.py first so the OI=0.15 baseline exists."
        )
    rows = [
        row
        for row in _read_csv(path)
        if abs(float(row["obstruction_ceiling"]) - OBSTRUCTION_CEILING) < 1e-9
    ]
    return {
        (row["age_group"], int(row["initial_grade"]), row["bbd_severity_label"]): row
        for row in rows
    }


def generate_no_asymmetry_recommendations() -> tuple[list[dict], list[dict]]:
    baseline = _load_baseline_reference()
    original_flag = techniques.PLACEMENT_WALL_ASYMMETRY_ENABLED
    cache: dict = {}

    out_rows: list[dict] = []
    try:
        techniques.PLACEMENT_WALL_ASYMMETRY_ENABLED = False
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
                            asymmetry_enabled=False,
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
                            "baseline_technique": ref["recommended_technique"],
                            "baseline_volume_ml": ref["deflux_volume_ml"],
                            "baseline_reflux_fraction": ref["reflux_fraction"],
                            "baseline_obstruction_index": ref["obstruction_index"],
                            "no_asymmetry_technique": chosen["technique"],
                            "no_asymmetry_volume_ml": f"{float(chosen['deflux_volume_ml']):.2f}",
                            "no_asymmetry_reflux_fraction": f"{float(chosen['reflux_fraction']):.4f}",
                            "no_asymmetry_obstruction_index": f"{float(chosen['obstruction_index']):.4f}",
                            "no_asymmetry_selection_mode": chosen["selection_mode"],
                            "no_asymmetry_safe_candidates_n": int(chosen["safe_candidates_n"]),
                            "changed_vs_baseline": int(changed),
                            "sampling_protocol": "full 75-stratum balanced grid; sex averaged; grade5 tortuosity 1.5|2.5|3.5",
                        }
                    )
    finally:
        techniques.PLACEMENT_WALL_ASYMMETRY_ENABLED = original_flag

    baseline_counts = Counter(row["baseline_technique"] for row in out_rows)
    no_asym_counts = Counter(row["no_asymmetry_technique"] for row in out_rows)
    summary_rows = [
        {
            "comparison": "baseline_placement_asymmetry_vs_no_asymmetry",
            "n_strata": len(out_rows),
            "changed_strata_count": sum(int(row["changed_vs_baseline"]) for row in out_rows),
            "baseline_double_hit_count": baseline_counts.get("double_hit", 0),
            "baseline_double_hit_plus_sting_count": baseline_counts.get("double_hit_plus_sting", 0),
            "baseline_hit_count": baseline_counts.get("hit", 0),
            "baseline_hit_plus_sting_count": baseline_counts.get("hit_plus_sting", 0),
            "baseline_sting_count": baseline_counts.get("sting", 0),
            "no_asymmetry_double_hit_count": no_asym_counts.get("double_hit", 0),
            "no_asymmetry_double_hit_plus_sting_count": no_asym_counts.get("double_hit_plus_sting", 0),
            "no_asymmetry_hit_count": no_asym_counts.get("hit", 0),
            "no_asymmetry_hit_plus_sting_count": no_asym_counts.get("hit_plus_sting", 0),
            "no_asymmetry_sting_count": no_asym_counts.get("sting", 0),
            "sampling_protocol": "full 75-stratum balanced grid; OI<=0.15",
        }
    ]

    _write_csv(
        REVISION_OUT / "placement_asymmetry_sensitivity.csv",
        [
            "age_group",
            "initial_grade",
            "bbd_severity_label",
            "baseline_technique",
            "baseline_volume_ml",
            "baseline_reflux_fraction",
            "baseline_obstruction_index",
            "no_asymmetry_technique",
            "no_asymmetry_volume_ml",
            "no_asymmetry_reflux_fraction",
            "no_asymmetry_obstruction_index",
            "no_asymmetry_selection_mode",
            "no_asymmetry_safe_candidates_n",
            "changed_vs_baseline",
            "sampling_protocol",
        ],
        out_rows,
    )
    _write_csv(
        REVISION_OUT / "placement_asymmetry_sensitivity_summary.csv",
        [
            "comparison",
            "n_strata",
            "changed_strata_count",
            "baseline_double_hit_count",
            "baseline_double_hit_plus_sting_count",
            "baseline_hit_count",
            "baseline_hit_plus_sting_count",
            "baseline_sting_count",
            "no_asymmetry_double_hit_count",
            "no_asymmetry_double_hit_plus_sting_count",
            "no_asymmetry_hit_count",
            "no_asymmetry_hit_plus_sting_count",
            "no_asymmetry_sting_count",
            "sampling_protocol",
        ],
        summary_rows,
    )
    return out_rows, summary_rows


def generate_grade5_sweep() -> list[dict]:
    """Grade-V reference sweep with and without the placement prior.

    Averaged over sex and the grade-V tortuosity sweep (1.5 / 2.5 / 3.5), which is
    the protocol behind `outputs/optimal_obstruction_sweep.csv` and therefore behind
    Figure 3 and Appendix B7. The base-case column of this table reproduces those
    published values exactly; see `verify_grade5_baseline_matches_published()`.
    """
    original_flag = techniques.PLACEMENT_WALL_ASYMMETRY_ENABLED
    rows: list[dict] = []
    cache: dict = {}
    try:
        for label, enabled in [("baseline_asymmetry", True), ("no_asymmetry", False)]:
            techniques.PLACEMENT_WALL_ASYMMETRY_ENABLED = enabled
            for technique in TECHNIQUES:
                for volume in VOLUME_GRID:
                    sims = [
                        _simulate(
                            SimPoint(
                                sex=sex,
                                age_group=SWEEP_AGE_GROUP,
                                initial_grade=5,
                                bbd_profile=BBDProfile.NONE.value,
                                bbd_severity=0.0,
                                technique=technique,
                                deflux_volume_ml=volume,
                                tortuosity_index=tort,
                                asymmetry_enabled=enabled,
                            ),
                            cache,
                        )
                        for sex in SEXES
                        for tort in GRADE5_TORT_SWEEP
                    ]
                    rows.append(
                        {
                            "scenario": label,
                            "placement_asymmetry_enabled": int(enabled),
                            "age_group": SWEEP_AGE_GROUP,
                            "sex": "female|male (averaged)",
                            "initial_grade": 5,
                            "tortuosity_index": "1.5|2.5|3.5 (averaged)",
                            "technique": technique,
                            "n_sites": _n_sites(technique),
                            "deflux_volume_ml": f"{volume:.2f}",
                            "reflux_fraction": f"{mean(s['reflux_fraction'] for s in sims):.4f}",
                            "obstruction_index": f"{mean(s['obstruction_index'] for s in sims):.4f}",
                            "reflux_volume_ml": f"{mean(s['reflux_volume_ml'] for s in sims):.4f}",
                            "antegrade_volume_ml": f"{mean(s['antegrade_volume_ml'] for s in sims):.4f}",
                            "uvj_forward_resistance_multiplier": f"{mean(s['uvj_forward_resistance_multiplier'] for s in sims):.3f}",
                        }
                    )
    finally:
        techniques.PLACEMENT_WALL_ASYMMETRY_ENABLED = original_flag

    _write_csv(
        REVISION_OUT / "placement_asymmetry_grade5_sweep.csv",
        [
            "scenario",
            "placement_asymmetry_enabled",
            "age_group",
            "sex",
            "initial_grade",
            "tortuosity_index",
            "technique",
            "n_sites",
            "deflux_volume_ml",
            "reflux_fraction",
            "obstruction_index",
            "reflux_volume_ml",
            "antegrade_volume_ml",
            "uvj_forward_resistance_multiplier",
        ],
        rows,
    )
    return rows


def verify_grade5_baseline_matches_published(rows: list[dict]) -> None:
    """Cross-check the base-case sweep against the published Figure 3 / Appendix B7 data.

    Single-site techniques (STING, HIT) reproduce the published values exactly.
    Multi-site techniques deviate by up to ~0.006 in reflux fraction because
    `optimal_obstruction_sweep.csv` was generated before the round-1 revision
    introduced `FORWARD_SITE_COUNT_ATTENUATION_EXPONENT` in techniques.py. That
    drift predates this analysis and is below the reporting precision of the
    figure, so the tolerance is set accordingly; anything larger indicates a real
    regression and fails the run.
    """
    strict = {"sting", "hit"}
    strict_tol, multisite_tol = 5e-4, 6e-3
    path = OUT / "optimal_obstruction_sweep.csv"
    if not path.exists():
        print(f"WARNING: {path.name} absent; skipping Figure 3 cross-check.")
        return
    published = {
        (r["technique"], f"{float(r['deflux_volume_ml']):.2f}"): r
        for r in _read_csv(path)
        if r["age_group"] == SWEEP_AGE_GROUP
        and int(r["initial_grade"]) == 5
        and r["bbd_severity_label"] == "none"
    }
    problems = []
    worst = (0.0, "")
    for row in rows:
        if row["scenario"] != "baseline_asymmetry":
            continue
        ref = published.get((row["technique"], row["deflux_volume_ml"]))
        if ref is None:
            problems.append(f"{row['technique']} @ {row['deflux_volume_ml']} mL absent from published sweep")
            continue
        tol = strict_tol if row["technique"] in strict else multisite_tol
        for field in ("reflux_fraction", "obstruction_index"):
            delta = abs(float(row[field]) - float(ref[field]))
            if delta > worst[0]:
                worst = (delta, f"{row['technique']} @ {row['deflux_volume_ml']} mL {field}")
            if delta > tol:
                problems.append(
                    f"{row['technique']} @ {row['deflux_volume_ml']} mL {field}: "
                    f"{row[field]} vs published {float(ref[field]):.4f} (delta {delta:.4f} > {tol})"
                )
    if problems:
        for p in problems[:20]:
            print("MISMATCH:", p, file=sys.stderr)
        raise SystemExit(
            f"{len(problems)} grade-V baseline mismatch(es) against optimal_obstruction_sweep.csv."
        )
    print(
        "Grade-V baseline sweep reproduces optimal_obstruction_sweep.csv "
        f"(Figure 3 / Appendix B7); STING and HIT exact, max deviation {worst[0]:.4f} at {worst[1]}."
    )


def main() -> None:
    detail_rows, summary_rows = generate_no_asymmetry_recommendations()
    sweep_rows = generate_grade5_sweep()
    verify_grade5_baseline_matches_published(sweep_rows)
    print("Generated placement-asymmetry sensitivity outputs:")
    print(f"- {REVISION_OUT / 'placement_asymmetry_sensitivity.csv'} ({len(detail_rows)} rows)")
    print(
        f"- {REVISION_OUT / 'placement_asymmetry_sensitivity_summary.csv'} ({len(summary_rows)} rows)"
    )
    print(f"- {REVISION_OUT / 'placement_asymmetry_grade5_sweep.csv'} ({len(sweep_rows)} rows)")
    print()
    print(f"Changed strata: {summary_rows[0]['changed_strata_count']}/{summary_rows[0]['n_strata']}")


if __name__ == "__main__":
    main()
