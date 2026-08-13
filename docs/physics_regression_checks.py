from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from vur_cfd.model import (
    AgeGroup,
    BBDProfile,
    BladderCapacityMethod,
    Sex,
    apply_initial_vur_grade,
    default_patient_from_literature,
    simulate_patient,
    simulate_patient_with_trace,
    with_bbd_modifiers,
    with_ureter_modifiers,
)


def _make_patient(
    grade: int,
    *,
    bbd_profile: BBDProfile = BBDProfile.NONE,
    bbd_severity: float = 0.0,
    tortuosity_index: float = 1.0,
) -> object:
    patient = default_patient_from_literature(
        age_group=AgeGroup.TODDLER_18_24M,
        sex=Sex.FEMALE,
        capacity_method=BladderCapacityMethod.KOFF,
        bladder_fill_fraction=1.0,
    )
    patient = apply_initial_vur_grade(patient, grade)
    patient = with_bbd_modifiers(patient, profile=bbd_profile, severity=bbd_severity)
    patient = with_ureter_modifiers(
        patient,
        tortuosity_index=tortuosity_index,
        compliance_factor=1.0,
        peristalsis_efficiency=1.0,
        bladder_fill_fraction=1.0,
    )
    return patient


def check_grade0_low_reflux() -> None:
    result = simulate_patient(_make_patient(0), total_time_s=24.0, dt_s=0.05)
    assert result.reflux_fraction <= 0.01, (
        f"Grade-0 reflux too high: {result.reflux_fraction:.4f} (expected <= 0.01)."
    )


def check_no_hard_antegrade_floor() -> None:
    patient = _make_patient(2)
    _, trace = simulate_patient_with_trace(patient, total_time_s=24.0, dt_s=0.05)
    floor = patient.urine_production_ml_per_s * 0.1
    forward = np.asarray(trace.forward_flow_ml_s, dtype=float)
    at_floor_fraction = float(np.mean(np.isclose(forward, floor, atol=1e-9)))
    assert at_floor_fraction < 0.01, (
        f"Forward flow still appears hard-clamped at floor: {100*at_floor_fraction:.1f}% steps."
    )


def check_bbd_obstruction_decoupling() -> None:
    neutral = simulate_patient(_make_patient(4), total_time_s=24.0, dt_s=0.05)
    bbd = simulate_patient(
        _make_patient(4, bbd_profile=BBDProfile.MIXED, bbd_severity=1.0),
        total_time_s=24.0,
        dt_s=0.05,
    )
    assert abs(neutral.obstruction_index - bbd.obstruction_index) < 1e-12, (
        "Obstruction index changed under BBD despite decoupling setting."
    )
    assert neutral.severe_obstruction == bbd.severe_obstruction, (
        "Severe obstruction flag changed under BBD despite decoupling setting."
    )


def check_grade_severity_ordering() -> None:
    g0 = simulate_patient(_make_patient(0), total_time_s=24.0, dt_s=0.05)
    g2 = simulate_patient(_make_patient(2), total_time_s=24.0, dt_s=0.05)
    g5 = simulate_patient(_make_patient(5, tortuosity_index=3.0), total_time_s=24.0, dt_s=0.05)
    assert g0.reflux_fraction < g2.reflux_fraction < g5.reflux_fraction, (
        "Baseline reflux fraction is not ordered by increasing pre-op grade."
    )


def main() -> None:
    checks = [
        ("grade0_low_reflux", check_grade0_low_reflux),
        ("no_hard_antegrade_floor", check_no_hard_antegrade_floor),
        ("bbd_obstruction_decoupling", check_bbd_obstruction_decoupling),
        ("grade_severity_ordering", check_grade_severity_ordering),
    ]
    for name, fn in checks:
        fn()
        print(f"[OK] {name}")


if __name__ == "__main__":
    main()
