"""VUR simulation package (verification subset).

`model.py` and `techniques.py` are byte-identical to the versions used to produce
the published results. This ``__init__`` differs from the development package in
one respect: it no longer re-exports ``cohort`` (cohort-run helpers),
``visualization`` (GIF/figure export) or ``main`` (the CLI), because those
modules are not needed to reproduce any reported number and are omitted from this
archive. No solver or technique code is affected.
"""

from .model import (
    AgeGroup,
    BBDProfile,
    BladderCycle,
    BladderCapacityMethod,
    BilateralSimulation,
    FluidProperties,
    PatientModel,
    Peristalsis,
    Sex,
    SimulationResult,
    SimulationTrace,
    UVJValve,
    apply_initial_vur_grade,
    available_age_groups,
    available_bbd_profiles,
    available_bladder_capacity_methods,
    cmh2o_to_pa,
    default_patient_from_literature,
    estimate_full_bladder_capacity_ml,
    pa_to_cmh2o,
    simulate_coupled_bilateral_with_trace,
    simulate_patient,
    simulate_patient_with_trace,
    ureter_resistance_from_geometry,
    with_bbd_modifiers,
    with_ureter_modifiers,
)
from .techniques import TechniqueName, TechniquePlan, apply_technique, bulking_injection_layout

__all__ = [
    "AgeGroup",
    "BBDProfile",
    "BladderCycle",
    "BladderCapacityMethod",
    "BilateralSimulation",
    "FluidProperties",
    "PatientModel",
    "Peristalsis",
    "Sex",
    "SimulationResult",
    "SimulationTrace",
    "TechniqueName",
    "TechniquePlan",
    "UVJValve",
    "apply_initial_vur_grade",
    "apply_technique",
    "bulking_injection_layout",
    "available_age_groups",
    "available_bbd_profiles",
    "available_bladder_capacity_methods",
    "cmh2o_to_pa",
    "default_patient_from_literature",
    "estimate_full_bladder_capacity_ml",
    "pa_to_cmh2o",
    "simulate_coupled_bilateral_with_trace",
    "simulate_patient",
    "simulate_patient_with_trace",
    "ureter_resistance_from_geometry",
    "with_bbd_modifiers",
    "with_ureter_modifiers",
]
