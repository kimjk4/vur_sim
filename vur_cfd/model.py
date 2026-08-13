from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
from math import cos, exp, pi, sin
from typing import Optional, Tuple


CMH2O_TO_PA = 98.0665
OUNCE_TO_ML = 29.5735

# Obstruction index weighting calibration knobs.
# Weights are normalized to sum to 1.0 so the index stays in [0, 1] when
# all component risks are in [0, 1].  Relative proportions preserve the
# original ratios (16 : 4 : 3 : 1 : 6).
OBSTRUCTION_PRESSURE_WEIGHT = 0.533
OBSTRUCTION_LOW_FORWARD_WEIGHT = 0.133
OBSTRUCTION_LOW_VOIDING_WEIGHT = 0.100
OBSTRUCTION_TORTUOSITY_WEIGHT = 0.034
OBSTRUCTION_EDEMA_WEIGHT = 0.200

# If enabled, obstruction metrics are evaluated from a BBD-neutral reference
# simulation so bladder-bowel dysfunction does not drive obstruction scoring.
DECOUPLE_BBD_FROM_OBSTRUCTION = True


def cmh2o_to_pa(value_cmh2o: float) -> float:
    return value_cmh2o * CMH2O_TO_PA


def pa_to_cmh2o(value_pa: float) -> float:
    return value_pa / CMH2O_TO_PA


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _pressure_gated_passive_forward(
    forward_dp_pa: float,
    passive_floor_ml_s: float,
) -> float:
    """Pressure-aware basal antegrade transport.

    A small baseline trickle can persist near neutral/favorable gradients
    (gravity + elastic recoil), but it should collapse when the bladder-to-pelvis
    gradient is clearly adverse.
    """
    dp_mid_pa = cmh2o_to_pa(1.0)
    dp_scale_pa = cmh2o_to_pa(1.5)
    gate = 1.0 / (1.0 + exp(-(forward_dp_pa - dp_mid_pa) / max(dp_scale_pa, 1e-9)))
    return passive_floor_ml_s * gate


def _peristalsis_compliance_penalty(compliance_factor: float) -> float:
    """Asymmetric peristalsis penalty for ureter wall compliance.

    Low compliance (stiff/scarred wall) impairs contractile force more
    aggressively than high compliance (floppy/dilated wall) which loses
    bolus containment — a mechanistically different failure mode.
    """
    if compliance_factor >= 1.0:
        return exp(-0.60 * (compliance_factor - 1.0))
    return exp(-0.90 * (1.0 - compliance_factor))


def _phase_transition_factor(
    bladder_cycle: "BladderCycle",
    t_s: float,
    competence: float,
    ramp_duration_s: float = 1.0,
) -> float:
    """Smooth sigmoid ramp for reverse_phase_factor at filling→voiding transition.

    Instead of a step from 1.0 (filling) to (0.65 + 0.35 * competence) (voiding),
    we ramp over ``ramp_duration_s`` at each phase boundary.
    """
    voiding_target = 0.65 + 0.35 * competence
    local_t = t_s % bladder_cycle.cycle_duration_s
    if local_t < bladder_cycle.filling_duration_s:
        # In filling phase — check proximity to end of filling.
        time_to_voiding = bladder_cycle.filling_duration_s - local_t
        if time_to_voiding < 0.5 * ramp_duration_s:
            # Pre-ramp: begin transitioning before voiding onset.
            alpha = 0.5 - time_to_voiding / ramp_duration_s
            return 1.0 + alpha * (voiding_target - 1.0)
        return 1.0
    # In voiding phase.
    time_since_voiding = local_t - bladder_cycle.filling_duration_s
    if time_since_voiding < 0.5 * ramp_duration_s:
        alpha = 0.5 + time_since_voiding / ramp_duration_s
        return 1.0 + alpha * (voiding_target - 1.0)
    return voiding_target


def _obstruction_index_from_components(
    pressure_risk: float,
    low_forward_fraction: float,
    low_voiding_outflow_fraction: float,
    tortuosity_risk: float,
    edema_risk: float,
) -> float:
    return (
        OBSTRUCTION_PRESSURE_WEIGHT * pressure_risk
        + OBSTRUCTION_LOW_FORWARD_WEIGHT * low_forward_fraction
        + OBSTRUCTION_LOW_VOIDING_WEIGHT * low_voiding_outflow_fraction
        + OBSTRUCTION_TORTUOSITY_WEIGHT * tortuosity_risk
        + OBSTRUCTION_EDEMA_WEIGHT * edema_risk
    )


def _bbd_profile_value(profile: BBDProfile | str) -> str:
    return profile.value if isinstance(profile, BBDProfile) else str(profile)


def _has_active_bbd(patient: "PatientModel") -> bool:
    profile = _bbd_profile_value(patient.bbd_profile)
    return profile != BBDProfile.NONE.value and float(patient.bbd_severity) > 1e-6


def _without_bbd(patient: "PatientModel") -> "PatientModel":
    return replace(patient, bbd_profile=BBDProfile.NONE, bbd_severity=0.0)


class Sex(str, Enum):
    FEMALE = "female"
    MALE = "male"


class AgeGroup(str, Enum):
    INFANT_0_12M = "0_12m"
    TODDLER_12_18M = "12_18m"
    TODDLER_18_24M = "18_24m"
    EARLY_CHILD_24_60M = "24_60m"
    CHILD_5_10Y = "5_10y"
    ADOLESCENT_10_16Y = "10_16y"


def available_age_groups() -> list[AgeGroup]:
    # 0-12m is intentionally excluded from current simulation workflows.
    return [a for a in AgeGroup if a != AgeGroup.INFANT_0_12M]


class BladderCapacityMethod(str, Enum):
    KOFF = "koff"
    INFANT_FOCUSED = "infant_focused"
    KAEFER = "kaefer"
    HOLMDAHL = "holmdahl"


def available_bladder_capacity_methods() -> list[BladderCapacityMethod]:
    # Keep legacy enum values for backward compatibility with old payloads,
    # but expose only the unified Koff method in current tooling.
    return [BladderCapacityMethod.KOFF]


class BBDProfile(str, Enum):
    NONE = "none"
    OVERACTIVE = "overactive"
    DYSFUNCTIONAL_VOIDING = "dysfunctional_voiding"
    MIXED = "mixed"


def available_bbd_profiles() -> list[BBDProfile]:
    return list(BBDProfile)


@dataclass(frozen=True)
class BBDTuning:
    baseline_pressure_gain: float
    filling_pressure_gain: float
    voiding_pressure_gain: float
    filling_spike_pa: float
    filling_spike_hz: float
    outlet_resistance_gain: float
    min_residual_fraction: float


BBD_PROFILE_TUNING = {
    BBDProfile.NONE: BBDTuning(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.00),
    BBDProfile.OVERACTIVE: BBDTuning(0.08, 0.22, 0.12, cmh2o_to_pa(18.0), 0.45, 1.15, 0.08),
    BBDProfile.DYSFUNCTIONAL_VOIDING: BBDTuning(
        0.05,
        0.12,
        0.24,
        cmh2o_to_pa(9.0),
        0.30,
        1.80,
        0.22,
    ),
    BBDProfile.MIXED: BBDTuning(0.10, 0.28, 0.30, cmh2o_to_pa(24.0), 0.50, 2.10, 0.30),
}


AGE_GROUP_MID_YEARS = {
    AgeGroup.INFANT_0_12M: 0.50,
    AgeGroup.TODDLER_12_18M: 1.25,
    AgeGroup.TODDLER_18_24M: 1.75,
    AgeGroup.EARLY_CHILD_24_60M: 3.50,
    AgeGroup.CHILD_5_10Y: 7.50,
    AgeGroup.ADOLESCENT_10_16Y: 13.00,
}


AGE_GROUP_WEIGHT_KG = {
    AgeGroup.INFANT_0_12M: {Sex.FEMALE: 7.5, Sex.MALE: 8.0},
    AgeGroup.TODDLER_12_18M: {Sex.FEMALE: 9.6, Sex.MALE: 10.3},
    AgeGroup.TODDLER_18_24M: {Sex.FEMALE: 11.0, Sex.MALE: 11.5},
    AgeGroup.EARLY_CHILD_24_60M: {Sex.FEMALE: 15.0, Sex.MALE: 15.7},
    AgeGroup.CHILD_5_10Y: {Sex.FEMALE: 27.0, Sex.MALE: 28.0},
    AgeGroup.ADOLESCENT_10_16Y: {Sex.FEMALE: 50.0, Sex.MALE: 55.0},
}


def age_group_mid_years(age_group: AgeGroup) -> float:
    return AGE_GROUP_MID_YEARS[age_group]


def estimated_weight_kg(age_group: AgeGroup, sex: Sex) -> float:
    return AGE_GROUP_WEIGHT_KG[age_group][sex]


def estimate_full_bladder_capacity_ml(
    age_group: AgeGroup,
    sex: Sex,
    method: BladderCapacityMethod = BladderCapacityMethod.KOFF,
) -> float:
    """
    Estimate pediatric full bladder capacity using the unified Koff equation.

    Current simulation policy uses Koff for all supported ages:
    Vcap (mL) = (age_years + 2) * 30

    `method` and `sex` are retained in the signature for backward compatibility.
    """
    age_y = age_group_mid_years(age_group)
    _ = sex
    _ = method
    capacity_ml = (age_y + 2.0) * 30.0

    return _clamp(capacity_ml, 30.0, 700.0)


@dataclass(frozen=True)
class AgeAnatomyProfile:
    ureter_diameter_mm: float
    ureter_length_mm: float
    uvj_orifice_diameter_mm: float
    female_urethra_length_mm: float
    male_urethra_length_mm: float


@dataclass(frozen=True)
class AgePressureProfile:
    filling_baseline_cmh2o: float
    filling_end_cmh2o: float
    voiding_peak_female_cmh2o: float
    voiding_peak_male_cmh2o: float


# Anatomy defaults (see docs/citations.md C1-C5).
AGE_ANATOMY_PROFILES = {
    AgeGroup.INFANT_0_12M: AgeAnatomyProfile(3.20, 125.0, 1.20, 25.0, 75.0),
    AgeGroup.TODDLER_12_18M: AgeAnatomyProfile(3.25, 132.5, 1.30, 23.5, 94.0),
    AgeGroup.TODDLER_18_24M: AgeAnatomyProfile(3.30, 137.5, 1.40, 23.1, 97.0),
    AgeGroup.EARLY_CHILD_24_60M: AgeAnatomyProfile(3.50, 155.0, 1.60, 26.0, 106.0),
    AgeGroup.CHILD_5_10Y: AgeAnatomyProfile(4.00, 195.0, 1.90, 28.0, 128.0),
    AgeGroup.ADOLESCENT_10_16Y: AgeAnatomyProfile(4.70, 250.0, 2.30, 32.0, 158.0),
}


# Pressure defaults (see docs/citations.md C6-C9).
AGE_PRESSURE_PROFILES = {
    AgeGroup.INFANT_0_12M: AgePressureProfile(5.0, 10.0, 95.0, 105.0),
    AgeGroup.TODDLER_12_18M: AgePressureProfile(5.0, 10.0, 75.0, 90.0),
    AgeGroup.TODDLER_18_24M: AgePressureProfile(5.0, 10.0, 60.0, 75.0),
    AgeGroup.EARLY_CHILD_24_60M: AgePressureProfile(4.5, 10.0, 55.0, 70.0),
    AgeGroup.CHILD_5_10Y: AgePressureProfile(4.0, 10.0, 50.0, 65.0),
    AgeGroup.ADOLESCENT_10_16Y: AgePressureProfile(4.0, 10.0, 45.0, 60.0),
}


@dataclass(frozen=True)
class FluidProperties:
    # NOTE: These literature values are carried for documentation/reference
    # only.  The reduced-order model prescribes hydraulic resistances via
    # empirical calibration constants (e.g. REFERENCE_URETER_RESISTANCE) and
    # the r^4 diameter law; it does not compute resistance from mu or rho, so
    # neither field is read in any flow calculation.
    density_kg_m3: float = 1050.0
    viscosity_pa_s: float = 0.0010


@dataclass(frozen=True)
class BladderCycle:
    """Explicit two-phase cycle: filling then voiding."""

    baseline_pa: float
    mcbc_pa: float
    max_voiding_pa: float
    filling_duration_s: float = 8.0
    voiding_duration_s: float = 4.0

    @property
    def cycle_duration_s(self) -> float:
        return self.filling_duration_s + self.voiding_duration_s

    def phase_at(self, t_s: float) -> str:
        local_t = t_s % self.cycle_duration_s
        return "filling" if local_t < self.filling_duration_s else "voiding"

    def pressure_at(self, t_s: float) -> float:
        local_t = t_s % self.cycle_duration_s
        if local_t < self.filling_duration_s:
            alpha = local_t / max(self.filling_duration_s, 1e-9)
            return self.baseline_pa + alpha * (self.mcbc_pa - self.baseline_pa)

        alpha_v = (local_t - self.filling_duration_s) / max(self.voiding_duration_s, 1e-9)
        if alpha_v <= 0.25:
            ramp = alpha_v / 0.25
            return self.mcbc_pa + ramp * (self.max_voiding_pa - self.mcbc_pa)
        decay = (alpha_v - 0.25) / 0.75
        return self.max_voiding_pa + decay * (self.baseline_pa - self.max_voiding_pa)


@dataclass(frozen=True)
class Peristalsis:
    frequency_hz: float = 0.80
    amplitude_pa: float = 320.0
    duty_cycle: float = 0.35      # retained; controls peak sharpness via power-shaping
    baseline_fraction: float = 0.15  # tonic tone floor (fraction of amplitude)

    def pressure_at(self, t_s: float, phase_offset: float = 0.0) -> float:
        phase = ((t_s * self.frequency_hz) + phase_offset) % 1.0
        wave = 0.5 * (1.0 - cos(2.0 * pi * phase))  # raised cosine, 0→1
        sharpness = max(self.duty_cycle, 0.1) / 0.5
        shaped = wave ** (1.0 / max(sharpness, 0.2))
        return self.amplitude_pa * (
            self.baseline_fraction + (1.0 - self.baseline_fraction) * shaped
        )


@dataclass(frozen=True)
class UVJValve:
    forward_resistance_pa_s_per_ml: float = 120.0
    reverse_resistance_pa_s_per_ml: float = 850.0
    closure_barrier_pa: float = 450.0
    competence: float = 0.35

    def effective_reverse_resistance(self) -> float:
        return self.reverse_resistance_pa_s_per_ml * (1.0 + 2.0 * self.competence)

    def effective_barrier_pa(self) -> float:
        return self.closure_barrier_pa * self.competence


REFERENCE_URETER_RESISTANCE = 90.0
REFERENCE_URETER_LENGTH_MM = 137.5
REFERENCE_URETER_DIAMETER_MM = 3.3


def ureter_resistance_from_geometry(
    length_mm: float,
    diameter_mm: float,
    tortuosity_index: float,
    compliance_factor: float,
) -> float:
    length_term = max(length_mm, 1e-9) / REFERENCE_URETER_LENGTH_MM
    diameter_term = (REFERENCE_URETER_DIAMETER_MM / max(diameter_mm, 1e-9)) ** 4
    # Tortuosity exponent 1.25: empirical calibration parameter.  For a
    # tortuous tube the Dean-number correction gives friction-factor increase
    # ~1 + 0.033·De^0.5.  At pediatric Dean numbers (De ~ 10-50) this yields
    # effective path-length and secondary-flow penalties broadly consistent
    # with the 1.25 power law over the tortuosity range 1.0-2.5.
    tortuosity_term = max(tortuosity_index, 0.2) ** 1.25
    compliance_term = max(compliance_factor, 0.2) ** 0.55
    return REFERENCE_URETER_RESISTANCE * length_term * diameter_term * tortuosity_term / compliance_term


def _elliptical_resistance_term(major_radius_mm: float, minor_radius_mm: float, length_mm: float) -> float:
    # Laminar extension of Hagen-Poiseuille for elliptical cross-sections:
    # R ~ l * (a^2 + b^2) / (a^3 * b^3).
    a = max(major_radius_mm, 1e-6)
    b = max(minor_radius_mm, 1e-6)
    l = max(length_mm, 1e-6)
    return l * (a * a + b * b) / (a * a * a * b * b * b)


def _uvj_deformation_from_bladder_fill(
    patient: PatientModel,
    bladder_volume_ml: float,
    full_capacity_ml: float,
    reflux_flow_ml_s: float = 0.0,
) -> tuple[float, float, float]:
    """
    Deformation-driven UVJ closure model inspired by Kalayeh et al. (2020):
    - reference (zero-stress) state around 10% capacity,
    - bladder wall stretch changes UVJ length and eccentricity,
    - resistance follows elliptical Poiseuille scaling.

    When *reflux_flow_ml_s* > 0 the retrograde pressure distends the UVJ from
    the ureter side, widening the minor axis — a positive-feedback mechanism
    most relevant for grade V with massive reflux volumes.
    """
    volume_ratio = _clamp(bladder_volume_ml / max(full_capacity_ml, 1e-9), 0.10, 1.20)
    stretch = _clamp((volume_ratio - 0.10) / 0.90, 0.0, 1.15)

    # UVJ competence encodes effective intramural geometry tendency.
    base_length_to_diameter = _clamp(2.8 + 4.0 * patient.uvj.competence, 2.5, 7.0)
    closure_tendency = _clamp((base_length_to_diameter - 3.0) / 2.0, -0.35, 1.0)

    major0 = max(0.5 * patient.uvj_orifice_diameter_mm, 0.20)
    minor0 = max(0.55 * major0, 0.08)
    length0 = max(5.0 * patient.uvj_orifice_diameter_mm, 3.0)

    major_scale = _clamp(
        1.0 + stretch * (0.10 + 0.16 * max(closure_tendency, 0.0)),
        0.80,
        1.55,
    )
    minor_scale = _clamp(
        1.0 - stretch * (0.10 + 0.52 * closure_tendency),
        0.24,
        1.55,
    )
    # Reflux-driven distension: retrograde pressure opens UVJ from ureter side.
    if reflux_flow_ml_s > 0.0:
        reflux_pressure_pa = reflux_flow_ml_s * patient.uvj.effective_reverse_resistance()
        # Normalize to a ~0-1 scale; 2000 Pa ≈ 20 cmH2O is a strong reflux pressure.
        reflux_stretch = _clamp(reflux_pressure_pa / 2000.0, 0.0, 1.0)
        # Low competence means the UVJ is more susceptible to reflux distension.
        susceptibility = 0.15 * (1.0 - patient.uvj.competence)
        minor_scale = _clamp(
            minor_scale * (1.0 + susceptibility * reflux_stretch),
            0.24,
            1.55,
        )
    length_scale = _clamp(
        1.0 + stretch * (0.05 + 0.30 * closure_tendency),
        0.70,
        1.45,
    )

    base_term = _elliptical_resistance_term(major0, minor0, length0)
    deformed_term = _elliptical_resistance_term(
        major0 * major_scale,
        minor0 * minor_scale,
        length0 * length_scale,
    )
    forward_multiplier = _clamp(deformed_term / max(base_term, 1e-9), 0.15, 40.0)
    reverse_multiplier = _clamp(
        forward_multiplier
        * (1.0 + 0.35 * max(closure_tendency, 0.0) * stretch),
        0.20,
        65.0,
    )
    barrier_gain_pa = cmh2o_to_pa(2.5) * max(closure_tendency, 0.0) * stretch
    return forward_multiplier, reverse_multiplier, barrier_gain_pa


@dataclass(frozen=True)
class PatientModel:
    name: str = "baseline"
    sex: Sex = Sex.FEMALE
    age_group: AgeGroup = AgeGroup.TODDLER_18_24M
    bladder_capacity_method: BladderCapacityMethod = BladderCapacityMethod.KOFF

    fluid: FluidProperties = field(default_factory=FluidProperties)
    bladder: BladderCycle = field(
        default_factory=lambda: BladderCycle(
            baseline_pa=cmh2o_to_pa(5.0),
            mcbc_pa=cmh2o_to_pa(10.0),
            max_voiding_pa=cmh2o_to_pa(60.0),
        )
    )
    peristalsis: Peristalsis = field(default_factory=Peristalsis)
    uvj: UVJValve = field(default_factory=UVJValve)

    ureter_diameter_mm: float = 3.3
    ureter_length_mm: float = 137.5
    uvj_orifice_diameter_mm: float = 1.4
    ureter_tortuosity_index: float = 1.0
    ureter_compliance_factor: float = 1.0
    ureter_dilation_multiplier: float = 1.0
    ureter_resistance_pa_s_per_ml: float = 90.0
    initial_vur_grade: int = 0

    peristalsis_efficiency: float = 1.0

    urethra_length_mm: float = 23.1
    urethra_diameter_mm: float = 3.8
    urethra_resistance_scale: float = 3900.0

    bbd_profile: BBDProfile = BBDProfile.NONE
    bbd_severity: float = 0.0

    estimated_full_bladder_capacity_ml: float = 165.0
    bladder_fill_fraction: float = 1.0

    renal_compliance_ml_per_pa: float = 0.003
    baseline_renal_pelvis_pressure_pa: float = cmh2o_to_pa(7.0)
    baseline_renal_pelvis_volume_ml: float = 8.0
    urine_production_ml_per_s: float = 0.03

    # Nonlinear renal pelvis compliance parameters (Whitaker-type stiffening).
    # Beyond pelvis_stiffening_volume_threshold_ml the collecting system
    # stiffens at rate pelvis_stiffening_alpha.  The implemented law is
    # linear_term * (1 + alpha * excess), i.e. quadratic in volume above the
    # threshold (not a true exponential), used as a reduced-order surrogate for
    # the steep pressure-volume behavior seen in perfusion studies.
    pelvis_stiffening_alpha: float = 0.25
    pelvis_stiffening_volume_threshold_ml: float = 12.0

    def urethral_resistance_pa_s_per_ml(self) -> float:
        diameter = max(self.urethra_diameter_mm, 1e-9)
        return self.urethra_resistance_scale * self.urethra_length_mm / (diameter**4)

    def effective_collecting_compliance_ml_per_pa(self) -> float:
        return self.renal_compliance_ml_per_pa * _clamp(self.ureter_compliance_factor, 0.35, 3.5)

    def renal_pelvis_pressure_pa(self, volume_ml: float) -> float:
        """Nonlinear pelvis pressure with Whitaker-type stiffening at high volumes."""
        compliance = self.effective_collecting_compliance_ml_per_pa()
        dv = volume_ml - self.baseline_renal_pelvis_volume_ml
        linear_term = dv / max(compliance, 1e-12)
        excess = max(volume_ml - self.pelvis_stiffening_volume_threshold_ml, 0.0)
        stiffening = self.pelvis_stiffening_alpha * excess
        return self.baseline_renal_pelvis_pressure_pa + linear_term * (1.0 + stiffening)

    def initial_bladder_volume_ml(self) -> float:
        return _clamp(self.bladder_fill_fraction, 0.05, 1.25) * self.estimated_full_bladder_capacity_ml


def with_ureter_modifiers(
    patient: PatientModel,
    tortuosity_index: Optional[float] = None,
    compliance_factor: Optional[float] = None,
    peristalsis_efficiency: Optional[float] = None,
    peristalsis_frequency_hz: Optional[float] = None,
    peristalsis_amplitude_pa: Optional[float] = None,
    bladder_fill_fraction: Optional[float] = None,
) -> PatientModel:
    requested_tortuosity = (
        _clamp(tortuosity_index, 1.0, 3.5)
        if tortuosity_index is not None
        else patient.ureter_tortuosity_index
    )
    if patient.initial_vur_grade <= 3:
        # Tortuosity is modeled as negligible in grades 0-III.
        new_tortuosity = 1.0
    elif patient.initial_vur_grade == 4:
        # Grade IV: mild-moderate tortuosity only.
        new_tortuosity = _clamp(requested_tortuosity, 1.0, 1.35)
    else:
        # Grade V: severe tortuosity can be expressed.
        new_tortuosity = requested_tortuosity
    new_compliance = (
        _clamp(compliance_factor, 0.35, 2.5)
        if compliance_factor is not None
        else patient.ureter_compliance_factor
    )
    new_per_eff = (
        _clamp(peristalsis_efficiency, 0.05, 2.0)
        if peristalsis_efficiency is not None
        else patient.peristalsis_efficiency
    )
    new_fill = (
        _clamp(bladder_fill_fraction, 0.05, 1.25)
        if bladder_fill_fraction is not None
        else patient.bladder_fill_fraction
    )

    new_peristalsis = patient.peristalsis
    if peristalsis_frequency_hz is not None or peristalsis_amplitude_pa is not None:
        new_peristalsis = replace(
            patient.peristalsis,
            frequency_hz=(
                _clamp(peristalsis_frequency_hz, 0.1, 2.0)
                if peristalsis_frequency_hz is not None
                else patient.peristalsis.frequency_hz
            ),
            amplitude_pa=(
                _clamp(peristalsis_amplitude_pa, 20.0, 2500.0)
                if peristalsis_amplitude_pa is not None
                else patient.peristalsis.amplitude_pa
            ),
        )

    new_resistance = ureter_resistance_from_geometry(
        length_mm=patient.ureter_length_mm,
        diameter_mm=patient.ureter_diameter_mm,
        tortuosity_index=new_tortuosity,
        compliance_factor=new_compliance,
    )

    return replace(
        patient,
        ureter_tortuosity_index=new_tortuosity,
        ureter_compliance_factor=new_compliance,
        peristalsis_efficiency=new_per_eff,
        bladder_fill_fraction=new_fill,
        peristalsis=new_peristalsis,
        ureter_resistance_pa_s_per_ml=new_resistance,
    )


def _resolve_bbd_state(patient: PatientModel) -> tuple[float, BBDTuning]:
    severity = _clamp(patient.bbd_severity, 0.0, 1.0)
    profile = patient.bbd_profile
    if not isinstance(profile, BBDProfile):
        profile = BBDProfile(str(profile))
    if severity <= 1e-6 or profile == BBDProfile.NONE:
        return 0.0, BBD_PROFILE_TUNING[BBDProfile.NONE]
    return severity, BBD_PROFILE_TUNING[profile]


def with_bbd_modifiers(
    patient: PatientModel,
    profile: BBDProfile | str = BBDProfile.NONE,
    severity: float = 0.0,
) -> PatientModel:
    parsed = profile if isinstance(profile, BBDProfile) else BBDProfile(str(profile))
    clamped_severity = _clamp(severity, 0.0, 1.0)
    if parsed == BBDProfile.NONE or clamped_severity <= 1e-6:
        return replace(patient, bbd_profile=BBDProfile.NONE, bbd_severity=0.0)
    return replace(patient, bbd_profile=parsed, bbd_severity=clamped_severity)


@dataclass(frozen=True)
class GradeTemplate:
    competence: float
    reverse_resistance_multiplier: float
    barrier_multiplier: float
    voiding_pressure_multiplier: float
    dilation_multiplier: float
    tortuosity_index: float
    compliance_multiplier: float
    peristalsis_multiplier: float


# International reflux grade morphology trend:
# high grade corresponds to greater dilation; tortuosity is mainly grades IV-V.
VUR_GRADE_TEMPLATES = {
    0: GradeTemplate(0.998, 36.0, 6.0, 0.90, 1.00, 1.00, 0.95, 1.05),
    1: GradeTemplate(0.85, 12.0, 2.5, 1.00, 1.02, 1.00, 1.00, 1.00),
    2: GradeTemplate(0.68, 7.0, 1.9, 1.10, 1.05, 1.00, 1.08, 0.95),
    3: GradeTemplate(0.48, 4.2, 1.3, 1.30, 1.25, 1.00, 1.20, 0.85),
    4: GradeTemplate(0.30, 2.4, 0.90, 1.60, 1.55, 1.10, 1.35, 0.70),
    5: GradeTemplate(0.15, 1.2, 0.60, 2.00, 2.05, 1.75, 1.55, 0.50),
}


# Grade-dependent standard peristaltic settings used when initializing
# the baseline reflux state. Higher grades are modeled with weaker/slower waves.
GRADE_PERISTALSIS_STANDARD = {
    0: (0.90, 230.0),
    1: (0.85, 210.0),
    2: (0.80, 195.0),
    3: (0.72, 170.0),
    4: (0.62, 145.0),
    5: (0.52, 115.0),
}


@dataclass(frozen=True)
class SimulationResult:
    reflux_volume_ml: float
    antegrade_volume_ml: float
    reflux_fraction: float
    peak_renal_pelvis_pressure_pa: float
    mean_forward_flow_ml_s: float
    mean_reflux_flow_ml_s: float
    mean_effective_peristalsis_gain: float
    voiding_urethral_outflow_ml: float
    filling_reflux_volume_ml: float
    voiding_reflux_volume_ml: float
    filling_antegrade_volume_ml: float
    voiding_antegrade_volume_ml: float
    filling_reflux_fraction: float
    voiding_reflux_fraction: float
    low_forward_fraction: float
    estimated_full_bladder_capacity_ml: float
    initial_bladder_volume_ml: float
    final_bladder_volume_ml: float
    cumulative_urine_passed_ml: float
    ureter_tortuosity_index: float
    ureter_compliance_factor: float
    vur_grade: int
    obstruction_index: float
    severe_obstruction: bool


@dataclass(frozen=True)
class SimulationTrace:
    time_s: list[float]
    phase: list[str]
    bladder_pressure_pa: list[float]
    pelvis_pressure_pa: list[float]
    bladder_volume_ml: list[float]
    cumulative_urine_passed_ml: list[float]
    forward_flow_ml_s: list[float]
    reflux_flow_ml_s: list[float]
    urethral_outflow_ml_s: list[float]
    effective_peristalsis_gain: list[float]


@dataclass(frozen=True)
class BilateralSimulation:
    combined_result: SimulationResult
    combined_trace: SimulationTrace
    left_result: SimulationResult
    right_result: SimulationResult
    left_trace: SimulationTrace
    right_trace: SimulationTrace


def default_patient_from_literature(
    age_group: AgeGroup = AgeGroup.TODDLER_18_24M,
    sex: Sex = Sex.FEMALE,
    capacity_method: BladderCapacityMethod = BladderCapacityMethod.KOFF,
    bladder_fill_fraction: float = 1.0,
) -> PatientModel:
    """
    Build age/sex-specific defaults from literature ranges.

    Source IDs are listed in docs/citations.md.
    """
    anatomy = AGE_ANATOMY_PROFILES[age_group]
    pressure = AGE_PRESSURE_PROFILES[age_group]
    full_capacity_ml = estimate_full_bladder_capacity_ml(
        age_group=age_group,
        sex=sex,
        method=capacity_method,
    )

    voiding_peak_cmh2o = (
        pressure.voiding_peak_female_cmh2o
        if sex == Sex.FEMALE
        else pressure.voiding_peak_male_cmh2o
    )
    urethra_length = (
        anatomy.female_urethra_length_mm
        if sex == Sex.FEMALE
        else anatomy.male_urethra_length_mm
    )

    ureter_resistance = ureter_resistance_from_geometry(
        length_mm=anatomy.ureter_length_mm,
        diameter_mm=anatomy.ureter_diameter_mm,
        tortuosity_index=1.0,
        compliance_factor=1.0,
    )

    return PatientModel(
        sex=sex,
        age_group=age_group,
        bladder_capacity_method=capacity_method,
        bladder=BladderCycle(
            baseline_pa=cmh2o_to_pa(pressure.filling_baseline_cmh2o),
            mcbc_pa=cmh2o_to_pa(pressure.filling_end_cmh2o),
            max_voiding_pa=cmh2o_to_pa(voiding_peak_cmh2o),
        ),
        ureter_diameter_mm=anatomy.ureter_diameter_mm,
        ureter_length_mm=anatomy.ureter_length_mm,
        uvj_orifice_diameter_mm=anatomy.uvj_orifice_diameter_mm,
        ureter_resistance_pa_s_per_ml=ureter_resistance,
        urethra_length_mm=urethra_length,
        estimated_full_bladder_capacity_ml=full_capacity_ml,
        bladder_fill_fraction=_clamp(bladder_fill_fraction, 0.05, 1.25),
        baseline_renal_pelvis_pressure_pa=cmh2o_to_pa(7.0),
    )


def apply_initial_vur_grade(
    patient: PatientModel,
    grade: int,
    *,
    use_grade_voiding_pressure_multiplier: bool = False,
) -> PatientModel:
    """Apply a grade template to initialize a pre-op patient state."""
    if grade not in VUR_GRADE_TEMPLATES:
        raise ValueError("grade must be in [0, 1, 2, 3, 4, 5].")

    template = VUR_GRADE_TEMPLATES[grade]
    uvj = patient.uvj
    bladder = patient.bladder

    infant_grade_pressure_enabled = (
        use_grade_voiding_pressure_multiplier
        and patient.age_group == AgeGroup.INFANT_0_12M
    )
    voiding_pressure_multiplier = (
        template.voiding_pressure_multiplier if infant_grade_pressure_enabled else 1.0
    )
    new_bladder = replace(
        bladder,
        mcbc_pa=bladder.mcbc_pa,
        max_voiding_pa=bladder.max_voiding_pa * voiding_pressure_multiplier,
    )
    new_uvj = UVJValve(
        forward_resistance_pa_s_per_ml=uvj.forward_resistance_pa_s_per_ml,
        reverse_resistance_pa_s_per_ml=(
            uvj.reverse_resistance_pa_s_per_ml * template.reverse_resistance_multiplier
        ),
        closure_barrier_pa=uvj.closure_barrier_pa * template.barrier_multiplier,
        competence=template.competence,
    )

    new_diameter = patient.ureter_diameter_mm * template.dilation_multiplier
    # Higher-grade reflux is typically associated with a wider/distorted UVJ orifice.
    # Keep this tied to grade severity but slightly less aggressive than full ureter dilation.
    new_orifice_diameter = patient.uvj_orifice_diameter_mm * (
        1.0 + 0.80 * (template.dilation_multiplier - 1.0)
    )
    if grade <= 3:
        new_tortuosity = 1.0
    else:
        new_tortuosity = max(patient.ureter_tortuosity_index, template.tortuosity_index)
    new_length = patient.ureter_length_mm * (1.0 + 0.20 * (new_tortuosity - 1.0))
    new_compliance = _clamp(
        patient.ureter_compliance_factor * template.compliance_multiplier,
        0.35,
        2.5,
    )
    new_peristalsis_eff = _clamp(
        patient.peristalsis_efficiency * template.peristalsis_multiplier,
        0.05,
        2.0,
    )
    standard_frequency_hz, standard_amplitude_pa = GRADE_PERISTALSIS_STANDARD[grade]
    new_peristalsis = replace(
        patient.peristalsis,
        frequency_hz=standard_frequency_hz,
        amplitude_pa=standard_amplitude_pa,
    )
    new_resistance = ureter_resistance_from_geometry(
        length_mm=new_length,
        diameter_mm=new_diameter,
        tortuosity_index=new_tortuosity,
        compliance_factor=new_compliance,
    )

    return replace(
        patient,
        name=f"{patient.name}_grade{grade}",
        initial_vur_grade=grade,
        bladder=new_bladder,
        uvj=new_uvj,
        ureter_diameter_mm=new_diameter,
        uvj_orifice_diameter_mm=new_orifice_diameter,
        ureter_length_mm=new_length,
        ureter_dilation_multiplier=template.dilation_multiplier,
        ureter_tortuosity_index=new_tortuosity,
        ureter_compliance_factor=new_compliance,
        ureter_resistance_pa_s_per_ml=new_resistance,
        peristalsis_efficiency=new_peristalsis_eff,
        peristalsis=new_peristalsis,
    )


def _vur_grade_from_surrogate(
    reflux_fraction: float,
    peak_pressure_pa: float,
    filling_reflux_fraction: float,
) -> int:
    """Map continuous simulation metrics to a discrete VUR grade (0-5).

    CALIBRATION NOTE: The component weights (0.55/0.20/0.25) and bin
    thresholds (0.05/0.12/0.25/0.37/0.45) are the single most important
    calibration target for clinical validity.  They are currently hand-tuned
    to produce plausible grade distributions in silico and should be
    validated against a dataset of patients with known VCUG grades and
    corresponding urodynamic measurements when such data become available.
    """
    reflux_component = _clamp((reflux_fraction - 0.01) / 0.60, 0.0, 1.0)
    pressure_component = _clamp(
        (peak_pressure_pa - cmh2o_to_pa(35.0)) / cmh2o_to_pa(90.0), 0.0, 1.0
    )
    filling_component = _clamp(filling_reflux_fraction / 0.45, 0.0, 1.0)
    score = 0.55 * reflux_component + 0.20 * pressure_component + 0.25 * filling_component

    if score < 0.05:
        return 0
    if score < 0.12:
        return 1
    if score < 0.25:
        return 2
    if score < 0.37:
        return 3
    if score < 0.45:
        return 4
    return 5


def _run_simulation(
    patient: PatientModel,
    total_time_s: float,
    dt_s: float,
    obstruction_pressure_threshold_pa: float,
    severe_pressure_threshold_pa: float,
    record_trace: bool,
) -> Tuple[SimulationResult, Optional[SimulationTrace]]:
    steps = int(total_time_s / dt_s)
    if steps <= 0:
        raise ValueError("Simulation must contain at least one time step.")

    trace_time: list[float] = []
    trace_phase: list[str] = []
    trace_bladder: list[float] = []
    trace_pelvis: list[float] = []
    trace_bladder_volume: list[float] = []
    trace_cumulative_passed: list[float] = []
    trace_forward: list[float] = []
    trace_reflux: list[float] = []
    trace_urethral: list[float] = []
    trace_per_gain: list[float] = []

    volume_ml = patient.baseline_renal_pelvis_volume_ml
    full_capacity_ml = max(patient.estimated_full_bladder_capacity_ml, 1.0)
    bladder_volume_ml = patient.initial_bladder_volume_ml()
    initial_bladder_volume_ml = bladder_volume_ml
    bbd_severity, bbd_tuning = _resolve_bbd_state(patient)
    min_residual_ml = bbd_severity * bbd_tuning.min_residual_fraction * full_capacity_ml

    forward_total_ml = 0.0
    reflux_total_ml = 0.0
    urethral_outflow_total_ml = 0.0
    prev_reflux_flow_ml_s = 0.0

    filling_reflux_ml = 0.0
    voiding_reflux_ml = 0.0
    filling_forward_ml = 0.0
    voiding_forward_ml = 0.0

    forward_sum_ml_s = 0.0
    reflux_sum_ml_s = 0.0
    per_gain_sum = 0.0
    low_forward_time_s = 0.0
    voiding_time_s = 0.0
    low_voiding_outflow_time_s = 0.0
    peak_pressure_pa = patient.baseline_renal_pelvis_pressure_pa
    cumulative_passed_ml = 0.0

    for step in range(steps):
        t_s = step * dt_s
        phase = patient.bladder.phase_at(t_s)
        is_filling = phase == "filling"

        raw_bladder_pa = patient.bladder.pressure_at(t_s)
        stretch_fraction = _clamp(bladder_volume_ml / full_capacity_ml, 0.0, 1.5)
        phase_scale = (0.20 + 0.80 * stretch_fraction) if is_filling else (0.35 + 0.65 * stretch_fraction)
        if is_filling:
            phase_scale *= 1.0 + bbd_severity * bbd_tuning.filling_pressure_gain
        else:
            phase_scale *= 1.0 + bbd_severity * bbd_tuning.voiding_pressure_gain
        bladder_pa = patient.bladder.baseline_pa + phase_scale * (
            raw_bladder_pa - patient.bladder.baseline_pa
        )
        bladder_pa += (
            patient.bladder.baseline_pa
            * bbd_severity
            * bbd_tuning.baseline_pressure_gain
        )
        if is_filling and bbd_tuning.filling_spike_pa > 0.0:
            # Filling contractions in BBD are modeled as smooth undulations.
            spike_phase = (t_s * max(bbd_tuning.filling_spike_hz, 1e-6)) % 1.0
            spike_wave = 0.5 * (1.0 - cos(2.0 * pi * spike_phase))
            bladder_pa += 0.70 * bbd_severity * bbd_tuning.filling_spike_pa * spike_wave

        per_gain = patient.peristalsis_efficiency
        per_gain *= _peristalsis_compliance_penalty(patient.ureter_compliance_factor)
        # Tortuosity exponent 0.35: empirical calibration parameter capturing
        # how ureteral kinking impairs peristaltic bolus transport.  Unlike the
        # resistance exponent (1.25) this has no direct Dean-flow derivation;
        # it was chosen to give ~12% gain reduction at tortuosity 1.75.
        per_gain /= max(patient.ureter_tortuosity_index, 0.2) ** 0.35
        per_gain = _clamp(per_gain, 0.03, 1.6)
        peristaltic_pa = patient.peristalsis.pressure_at(t_s) * per_gain

        pelvis_pa = patient.renal_pelvis_pressure_pa(volume_ml)

        (
            uvj_forward_deformation_mult,
            uvj_reverse_deformation_mult,
            uvj_deformation_barrier_gain_pa,
        ) = _uvj_deformation_from_bladder_fill(
            patient,
            bladder_volume_ml=bladder_volume_ml,
            full_capacity_ml=full_capacity_ml,
            reflux_flow_ml_s=prev_reflux_flow_ml_s,
        )
        forward_dp_pa = pelvis_pa + peristaltic_pa - bladder_pa
        forward_resistance = (
            patient.ureter_resistance_pa_s_per_ml
            + patient.uvj.forward_resistance_pa_s_per_ml * uvj_forward_deformation_mult
        )
        passive_floor_ml_s = patient.urine_production_ml_per_s * 0.1
        raw_forward = forward_dp_pa / max(forward_resistance, 1e-12)
        passive_forward_ml_s = _pressure_gated_passive_forward(
            forward_dp_pa,
            passive_floor_ml_s,
        )
        forward_flow_ml_s = max(raw_forward + passive_forward_ml_s, 0.0)

        reflux_dp_pa = (
            bladder_pa
            - pelvis_pa
            - patient.uvj.effective_barrier_pa()
            - uvj_deformation_barrier_gain_pa
        )
        reverse_phase_factor = _phase_transition_factor(
            patient.bladder, t_s, patient.uvj.competence,
        )
        reverse_tortuosity_factor = 1.0 + 0.20 * max(patient.ureter_tortuosity_index - 1.0, 0.0)
        reverse_compliance_factor = 1.0 / (max(patient.ureter_compliance_factor, 0.2) ** 0.25)
        reverse_resistance = (
            patient.uvj.effective_reverse_resistance()
            * uvj_reverse_deformation_mult
            * reverse_phase_factor
            * reverse_tortuosity_factor
            * reverse_compliance_factor
        )
        reflux_flow_ml_s = max(reflux_dp_pa / max(reverse_resistance, 1e-12), 0.0)

        if is_filling:
            urethral_outflow_ml_s = 0.0
        else:
            urethral_resistance = (
                patient.urethral_resistance_pa_s_per_ml()
                * (1.0 + bbd_severity * (bbd_tuning.outlet_resistance_gain - 1.0))
            )
            outflow_scale = _clamp(bladder_volume_ml / full_capacity_ml, 0.05, 1.2)
            urethral_outflow_ml_s = max(
                bladder_pa / max(urethral_resistance, 1e-12),
                0.0,
            ) * outflow_scale
            urethral_outflow_ml_s = min(urethral_outflow_ml_s, bladder_volume_ml / max(dt_s, 1e-9))

        net_volume_change = (
            patient.urine_production_ml_per_s + reflux_flow_ml_s - forward_flow_ml_s
        ) * dt_s
        volume_ml = max(volume_ml + net_volume_change, 0.2)

        bladder_delta_ml = (forward_flow_ml_s - reflux_flow_ml_s - urethral_outflow_ml_s) * dt_s
        bladder_volume_ml = _clamp(
            bladder_volume_ml + bladder_delta_ml,
            min_residual_ml,
            1.4 * full_capacity_ml,
        )

        cumulative_passed_ml += urethral_outflow_ml_s * dt_s
        peak_pressure_pa = max(peak_pressure_pa, pelvis_pa)
        forward_total_ml += forward_flow_ml_s * dt_s
        reflux_total_ml += reflux_flow_ml_s * dt_s
        urethral_outflow_total_ml += urethral_outflow_ml_s * dt_s
        forward_sum_ml_s += forward_flow_ml_s
        reflux_sum_ml_s += reflux_flow_ml_s
        per_gain_sum += per_gain
        prev_reflux_flow_ml_s = reflux_flow_ml_s

        if is_filling:
            filling_reflux_ml += reflux_flow_ml_s * dt_s
            filling_forward_ml += forward_flow_ml_s * dt_s
        else:
            voiding_reflux_ml += reflux_flow_ml_s * dt_s
            voiding_forward_ml += forward_flow_ml_s * dt_s
            voiding_time_s += dt_s
            if urethral_outflow_ml_s < 0.8:
                low_voiding_outflow_time_s += dt_s

        if (
            forward_flow_ml_s < 0.5 * patient.urine_production_ml_per_s
            and pelvis_pa >= obstruction_pressure_threshold_pa
        ):
            low_forward_time_s += dt_s

        if record_trace:
            trace_time.append(t_s)
            trace_phase.append(phase)
            trace_bladder.append(bladder_pa)
            trace_pelvis.append(pelvis_pa)
            trace_bladder_volume.append(bladder_volume_ml)
            trace_cumulative_passed.append(cumulative_passed_ml)
            trace_forward.append(forward_flow_ml_s)
            trace_reflux.append(reflux_flow_ml_s)
            trace_urethral.append(urethral_outflow_ml_s)
            trace_per_gain.append(per_gain)

    total_transported_ml = max(forward_total_ml + reflux_total_ml, 1e-12)
    reflux_fraction = reflux_total_ml / total_transported_ml
    low_forward_fraction = low_forward_time_s / total_time_s

    filling_transport_ml = max(filling_reflux_ml + filling_forward_ml, 1e-12)
    voiding_transport_ml = max(voiding_reflux_ml + voiding_forward_ml, 1e-12)
    filling_reflux_fraction = filling_reflux_ml / filling_transport_ml
    voiding_reflux_fraction = voiding_reflux_ml / voiding_transport_ml

    vur_grade = _vur_grade_from_surrogate(
        reflux_fraction,
        peak_pressure_pa,
        filling_reflux_fraction,
    )

    pressure_risk = _clamp(
        (peak_pressure_pa - obstruction_pressure_threshold_pa)
        / max(severe_pressure_threshold_pa - obstruction_pressure_threshold_pa, 1e-12),
        0.0,
        2.0,
    )
    # Edema/narrowing proxy from UVJ forward resistance multiplier.
    baseline_uvj_forward = UVJValve().forward_resistance_pa_s_per_ml
    forward_resistance_multiplier = patient.uvj.forward_resistance_pa_s_per_ml / max(
        baseline_uvj_forward, 1e-9
    )
    edema_risk = _clamp((forward_resistance_multiplier - 1.0) / 4.0, 0.0, 2.0)
    low_voiding_outflow_fraction = (
        0.0 if voiding_time_s <= 0.0 else low_voiding_outflow_time_s / voiding_time_s
    )
    tortuosity_risk = _clamp((patient.ureter_tortuosity_index - 1.0) / 1.0, 0.0, 1.5)
    obstruction_index = _obstruction_index_from_components(
        pressure_risk=pressure_risk,
        low_forward_fraction=low_forward_fraction,
        low_voiding_outflow_fraction=low_voiding_outflow_fraction,
        tortuosity_risk=tortuosity_risk,
        edema_risk=edema_risk,
    )
    severe_obstruction = (
        peak_pressure_pa >= severe_pressure_threshold_pa or obstruction_index >= 1.0
    )

    result = SimulationResult(
        reflux_volume_ml=reflux_total_ml,
        antegrade_volume_ml=forward_total_ml,
        reflux_fraction=reflux_fraction,
        peak_renal_pelvis_pressure_pa=peak_pressure_pa,
        mean_forward_flow_ml_s=forward_sum_ml_s / steps,
        mean_reflux_flow_ml_s=reflux_sum_ml_s / steps,
        mean_effective_peristalsis_gain=per_gain_sum / steps,
        voiding_urethral_outflow_ml=urethral_outflow_total_ml,
        filling_reflux_volume_ml=filling_reflux_ml,
        voiding_reflux_volume_ml=voiding_reflux_ml,
        filling_antegrade_volume_ml=filling_forward_ml,
        voiding_antegrade_volume_ml=voiding_forward_ml,
        filling_reflux_fraction=filling_reflux_fraction,
        voiding_reflux_fraction=voiding_reflux_fraction,
        low_forward_fraction=low_forward_fraction,
        estimated_full_bladder_capacity_ml=full_capacity_ml,
        initial_bladder_volume_ml=initial_bladder_volume_ml,
        final_bladder_volume_ml=bladder_volume_ml,
        cumulative_urine_passed_ml=cumulative_passed_ml,
        ureter_tortuosity_index=patient.ureter_tortuosity_index,
        ureter_compliance_factor=patient.ureter_compliance_factor,
        vur_grade=vur_grade,
        obstruction_index=obstruction_index,
        severe_obstruction=severe_obstruction,
    )

    trace = None
    if record_trace:
        trace = SimulationTrace(
            time_s=trace_time,
            phase=trace_phase,
            bladder_pressure_pa=trace_bladder,
            pelvis_pressure_pa=trace_pelvis,
            bladder_volume_ml=trace_bladder_volume,
            cumulative_urine_passed_ml=trace_cumulative_passed,
            forward_flow_ml_s=trace_forward,
            reflux_flow_ml_s=trace_reflux,
            urethral_outflow_ml_s=trace_urethral,
            effective_peristalsis_gain=trace_per_gain,
        )
    return result, trace


def simulate_patient(
    patient: PatientModel,
    total_time_s: float = 120.0,
    dt_s: float = 0.02,
    obstruction_pressure_threshold_pa: float = cmh2o_to_pa(30.0),
    severe_pressure_threshold_pa: float = cmh2o_to_pa(45.0),
) -> SimulationResult:
    if total_time_s <= 0.0 or dt_s <= 0.0:
        raise ValueError("total_time_s and dt_s must be positive.")
    result, _ = _run_simulation(
        patient=patient,
        total_time_s=total_time_s,
        dt_s=dt_s,
        obstruction_pressure_threshold_pa=obstruction_pressure_threshold_pa,
        severe_pressure_threshold_pa=severe_pressure_threshold_pa,
        record_trace=False,
    )
    if DECOUPLE_BBD_FROM_OBSTRUCTION and _has_active_bbd(patient):
        reference_patient = _without_bbd(patient)
        reference_result, _ = _run_simulation(
            patient=reference_patient,
            total_time_s=total_time_s,
            dt_s=dt_s,
            obstruction_pressure_threshold_pa=obstruction_pressure_threshold_pa,
            severe_pressure_threshold_pa=severe_pressure_threshold_pa,
            record_trace=False,
        )
        result = replace(
            result,
            obstruction_index=reference_result.obstruction_index,
            severe_obstruction=reference_result.severe_obstruction,
        )
    return result


def simulate_patient_with_trace(
    patient: PatientModel,
    total_time_s: float = 120.0,
    dt_s: float = 0.02,
    obstruction_pressure_threshold_pa: float = cmh2o_to_pa(30.0),
    severe_pressure_threshold_pa: float = cmh2o_to_pa(45.0),
) -> Tuple[SimulationResult, SimulationTrace]:
    if total_time_s <= 0.0 or dt_s <= 0.0:
        raise ValueError("total_time_s and dt_s must be positive.")
    result, trace = _run_simulation(
        patient=patient,
        total_time_s=total_time_s,
        dt_s=dt_s,
        obstruction_pressure_threshold_pa=obstruction_pressure_threshold_pa,
        severe_pressure_threshold_pa=severe_pressure_threshold_pa,
        record_trace=True,
    )
    if DECOUPLE_BBD_FROM_OBSTRUCTION and _has_active_bbd(patient):
        reference_patient = _without_bbd(patient)
        reference_result, _ = _run_simulation(
            patient=reference_patient,
            total_time_s=total_time_s,
            dt_s=dt_s,
            obstruction_pressure_threshold_pa=obstruction_pressure_threshold_pa,
            severe_pressure_threshold_pa=severe_pressure_threshold_pa,
            record_trace=False,
        )
        result = replace(
            result,
            obstruction_index=reference_result.obstruction_index,
            severe_obstruction=reference_result.severe_obstruction,
        )
    assert trace is not None
    return result, trace


def _shared_bladder_cycle(left: PatientModel, right: PatientModel) -> BladderCycle:
    return BladderCycle(
        baseline_pa=0.5 * (left.bladder.baseline_pa + right.bladder.baseline_pa),
        mcbc_pa=max(left.bladder.mcbc_pa, right.bladder.mcbc_pa),
        max_voiding_pa=max(left.bladder.max_voiding_pa, right.bladder.max_voiding_pa),
        filling_duration_s=0.5 * (left.bladder.filling_duration_s + right.bladder.filling_duration_s),
        voiding_duration_s=0.5 * (left.bladder.voiding_duration_s + right.bladder.voiding_duration_s),
    )


def _new_side_state(patient: PatientModel) -> dict:
    return {
        "patient": patient,
        "volume_ml": patient.baseline_renal_pelvis_volume_ml,
        "peak_pressure_pa": patient.baseline_renal_pelvis_pressure_pa,
        "forward_total_ml": 0.0,
        "reflux_total_ml": 0.0,
        "filling_forward_ml": 0.0,
        "voiding_forward_ml": 0.0,
        "filling_reflux_ml": 0.0,
        "voiding_reflux_ml": 0.0,
        "forward_sum_ml_s": 0.0,
        "reflux_sum_ml_s": 0.0,
        "per_gain_sum": 0.0,
        "low_forward_time_s": 0.0,
        "prev_reflux_flow_ml_s": 0.0,
        "trace_pelvis": [],
        "trace_forward": [],
        "trace_reflux": [],
        "trace_per_gain": [],
    }


def _build_side_result(
    side_state: dict,
    *,
    steps: int,
    total_time_s: float,
    full_capacity_ml: float,
    initial_bladder_volume_ml: float,
    final_bladder_volume_ml: float,
    cumulative_urine_passed_ml: float,
    urethral_outflow_total_ml: float,
    low_voiding_outflow_fraction: float,
    obstruction_pressure_threshold_pa: float,
    severe_pressure_threshold_pa: float,
) -> SimulationResult:
    patient = side_state["patient"]
    forward_total_ml = side_state["forward_total_ml"]
    reflux_total_ml = side_state["reflux_total_ml"]
    filling_forward_ml = side_state["filling_forward_ml"]
    voiding_forward_ml = side_state["voiding_forward_ml"]
    filling_reflux_ml = side_state["filling_reflux_ml"]
    voiding_reflux_ml = side_state["voiding_reflux_ml"]
    peak_pressure_pa = side_state["peak_pressure_pa"]

    total_transported_ml = max(forward_total_ml + reflux_total_ml, 1e-12)
    reflux_fraction = reflux_total_ml / total_transported_ml

    filling_transport_ml = max(filling_forward_ml + filling_reflux_ml, 1e-12)
    voiding_transport_ml = max(voiding_forward_ml + voiding_reflux_ml, 1e-12)
    filling_reflux_fraction = filling_reflux_ml / filling_transport_ml
    voiding_reflux_fraction = voiding_reflux_ml / voiding_transport_ml

    low_forward_fraction = side_state["low_forward_time_s"] / total_time_s
    pressure_risk = _clamp(
        (peak_pressure_pa - obstruction_pressure_threshold_pa)
        / max(severe_pressure_threshold_pa - obstruction_pressure_threshold_pa, 1e-12),
        0.0,
        2.0,
    )
    baseline_uvj_forward = UVJValve().forward_resistance_pa_s_per_ml
    forward_resistance_multiplier = patient.uvj.forward_resistance_pa_s_per_ml / max(
        baseline_uvj_forward, 1e-9
    )
    edema_risk = _clamp((forward_resistance_multiplier - 1.0) / 4.0, 0.0, 2.0)
    tortuosity_risk = _clamp((patient.ureter_tortuosity_index - 1.0) / 1.0, 0.0, 1.5)
    obstruction_index = _obstruction_index_from_components(
        pressure_risk=pressure_risk,
        low_forward_fraction=low_forward_fraction,
        low_voiding_outflow_fraction=low_voiding_outflow_fraction,
        tortuosity_risk=tortuosity_risk,
        edema_risk=edema_risk,
    )
    severe_obstruction = (
        peak_pressure_pa >= severe_pressure_threshold_pa or obstruction_index >= 1.0
    )
    vur_grade = _vur_grade_from_surrogate(
        reflux_fraction,
        peak_pressure_pa,
        filling_reflux_fraction,
    )

    return SimulationResult(
        reflux_volume_ml=reflux_total_ml,
        antegrade_volume_ml=forward_total_ml,
        reflux_fraction=reflux_fraction,
        peak_renal_pelvis_pressure_pa=peak_pressure_pa,
        mean_forward_flow_ml_s=side_state["forward_sum_ml_s"] / steps,
        mean_reflux_flow_ml_s=side_state["reflux_sum_ml_s"] / steps,
        mean_effective_peristalsis_gain=side_state["per_gain_sum"] / steps,
        voiding_urethral_outflow_ml=urethral_outflow_total_ml,
        filling_reflux_volume_ml=filling_reflux_ml,
        voiding_reflux_volume_ml=voiding_reflux_ml,
        filling_antegrade_volume_ml=filling_forward_ml,
        voiding_antegrade_volume_ml=voiding_forward_ml,
        filling_reflux_fraction=filling_reflux_fraction,
        voiding_reflux_fraction=voiding_reflux_fraction,
        low_forward_fraction=low_forward_fraction,
        estimated_full_bladder_capacity_ml=full_capacity_ml,
        initial_bladder_volume_ml=initial_bladder_volume_ml,
        final_bladder_volume_ml=final_bladder_volume_ml,
        cumulative_urine_passed_ml=cumulative_urine_passed_ml,
        ureter_tortuosity_index=patient.ureter_tortuosity_index,
        ureter_compliance_factor=patient.ureter_compliance_factor,
        vur_grade=vur_grade,
        obstruction_index=obstruction_index,
        severe_obstruction=severe_obstruction,
    )


def _build_side_trace(
    side_state: dict,
    *,
    trace_time: list[float],
    trace_phase: list[str],
    trace_bladder: list[float],
    trace_bladder_volume: list[float],
    trace_cumulative_passed: list[float],
    trace_urethral: list[float],
) -> SimulationTrace:
    return SimulationTrace(
        time_s=trace_time,
        phase=trace_phase,
        bladder_pressure_pa=trace_bladder,
        pelvis_pressure_pa=side_state["trace_pelvis"],
        bladder_volume_ml=trace_bladder_volume,
        cumulative_urine_passed_ml=trace_cumulative_passed,
        forward_flow_ml_s=side_state["trace_forward"],
        reflux_flow_ml_s=side_state["trace_reflux"],
        urethral_outflow_ml_s=trace_urethral,
        effective_peristalsis_gain=side_state["trace_per_gain"],
    )


def simulate_coupled_bilateral_with_trace(
    left_patient: PatientModel,
    right_patient: PatientModel,
    total_time_s: float = 120.0,
    dt_s: float = 0.02,
    obstruction_pressure_threshold_pa: float = cmh2o_to_pa(30.0),
    severe_pressure_threshold_pa: float = cmh2o_to_pa(45.0),
) -> BilateralSimulation:
    """
    Coupled bilateral solver with a single shared bladder state.

    Left and right tracts are integrated simultaneously at each timestep, while
    bladder pressure/volume and urethral outflow are shared.
    """
    if total_time_s <= 0.0 or dt_s <= 0.0:
        raise ValueError("total_time_s and dt_s must be positive.")

    steps = int(total_time_s / dt_s)
    if steps <= 0:
        raise ValueError("Simulation must contain at least one time step.")

    shared_bladder = _shared_bladder_cycle(left_patient, right_patient)
    urethral_resistance_pa_s_per_ml = 0.5 * (
        left_patient.urethral_resistance_pa_s_per_ml()
        + right_patient.urethral_resistance_pa_s_per_ml()
    )
    full_capacity_ml = max(
        0.5
        * (
            left_patient.estimated_full_bladder_capacity_ml
            + right_patient.estimated_full_bladder_capacity_ml
        ),
        1.0,
    )
    bladder_volume_ml = 0.5 * (
        left_patient.initial_bladder_volume_ml() + right_patient.initial_bladder_volume_ml()
    )
    bladder_volume_ml = _clamp(bladder_volume_ml, 0.0, 1.4 * full_capacity_ml)
    initial_bladder_volume_ml = bladder_volume_ml

    left_bbd_severity, left_bbd_tuning = _resolve_bbd_state(left_patient)
    right_bbd_severity, right_bbd_tuning = _resolve_bbd_state(right_patient)
    if left_bbd_severity >= right_bbd_severity:
        bbd_severity = left_bbd_severity
        bbd_tuning = left_bbd_tuning
    else:
        bbd_severity = right_bbd_severity
        bbd_tuning = right_bbd_tuning
    left_profile = (
        left_patient.bbd_profile.value
        if isinstance(left_patient.bbd_profile, BBDProfile)
        else str(left_patient.bbd_profile)
    )
    right_profile = (
        right_patient.bbd_profile.value
        if isinstance(right_patient.bbd_profile, BBDProfile)
        else str(right_patient.bbd_profile)
    )
    if bbd_severity > 0.0 and left_profile != right_profile:
        bbd_tuning = BBD_PROFILE_TUNING[BBDProfile.MIXED]
    min_residual_ml = bbd_severity * bbd_tuning.min_residual_fraction * full_capacity_ml

    side_states = {
        "left": _new_side_state(left_patient),
        "right": _new_side_state(right_patient),
    }

    urethral_outflow_total_ml = 0.0
    cumulative_passed_ml = 0.0
    voiding_time_s = 0.0
    low_voiding_outflow_time_s = 0.0

    trace_time: list[float] = []
    trace_phase: list[str] = []
    trace_bladder: list[float] = []
    trace_bladder_volume: list[float] = []
    trace_cumulative_passed: list[float] = []
    trace_combined_pelvis: list[float] = []
    trace_combined_forward: list[float] = []
    trace_combined_reflux: list[float] = []
    trace_urethral: list[float] = []
    trace_combined_per_gain: list[float] = []

    for step in range(steps):
        t_s = step * dt_s
        phase = shared_bladder.phase_at(t_s)
        is_filling = phase == "filling"

        raw_bladder_pa = shared_bladder.pressure_at(t_s)
        stretch_fraction = _clamp(bladder_volume_ml / full_capacity_ml, 0.0, 1.5)
        phase_scale = (0.20 + 0.80 * stretch_fraction) if is_filling else (0.35 + 0.65 * stretch_fraction)
        if is_filling:
            phase_scale *= 1.0 + bbd_severity * bbd_tuning.filling_pressure_gain
        else:
            phase_scale *= 1.0 + bbd_severity * bbd_tuning.voiding_pressure_gain
        bladder_pa = shared_bladder.baseline_pa + phase_scale * (
            raw_bladder_pa - shared_bladder.baseline_pa
        )
        bladder_pa += (
            shared_bladder.baseline_pa
            * bbd_severity
            * bbd_tuning.baseline_pressure_gain
        )
        if is_filling and bbd_tuning.filling_spike_pa > 0.0:
            spike_phase = (t_s * max(bbd_tuning.filling_spike_hz, 1e-6)) % 1.0
            spike_wave = 0.5 * (1.0 - cos(2.0 * pi * spike_phase))
            bladder_pa += 0.70 * bbd_severity * bbd_tuning.filling_spike_pa * spike_wave

        combined_forward_flow_ml_s = 0.0
        combined_reflux_flow_ml_s = 0.0
        combined_peak_pelvis_pa = -1e12
        combined_per_gain = 0.0

        for side_name, side in side_states.items():
            patient = side["patient"]
            per_gain = patient.peristalsis_efficiency
            per_gain *= _peristalsis_compliance_penalty(patient.ureter_compliance_factor)
            # See tortuosity exponent note in _run_simulation.
            per_gain /= max(patient.ureter_tortuosity_index, 0.2) ** 0.35
            per_gain = _clamp(per_gain, 0.03, 1.6)
            side_phase_offset = 0.5 if side_name == "right" else 0.0
            peristaltic_pa = patient.peristalsis.pressure_at(t_s, phase_offset=side_phase_offset) * per_gain

            pelvis_pa = patient.renal_pelvis_pressure_pa(side["volume_ml"])

            (
                uvj_forward_deformation_mult,
                uvj_reverse_deformation_mult,
                uvj_deformation_barrier_gain_pa,
            ) = _uvj_deformation_from_bladder_fill(
                patient,
                bladder_volume_ml=bladder_volume_ml,
                full_capacity_ml=full_capacity_ml,
                reflux_flow_ml_s=side["prev_reflux_flow_ml_s"],
            )
            forward_dp_pa = pelvis_pa + peristaltic_pa - bladder_pa
            forward_resistance = (
                patient.ureter_resistance_pa_s_per_ml
                + patient.uvj.forward_resistance_pa_s_per_ml * uvj_forward_deformation_mult
            )
            passive_floor_ml_s = patient.urine_production_ml_per_s * 0.1
            raw_forward = forward_dp_pa / max(forward_resistance, 1e-12)
            passive_forward_ml_s = _pressure_gated_passive_forward(
                forward_dp_pa,
                passive_floor_ml_s,
            )
            forward_flow_ml_s = max(raw_forward + passive_forward_ml_s, 0.0)

            reflux_dp_pa = (
                bladder_pa
                - pelvis_pa
                - patient.uvj.effective_barrier_pa()
                - uvj_deformation_barrier_gain_pa
            )
            reverse_phase_factor = _phase_transition_factor(
                shared_bladder, t_s, patient.uvj.competence,
            )
            reverse_tortuosity_factor = 1.0 + 0.20 * max(patient.ureter_tortuosity_index - 1.0, 0.0)
            reverse_compliance_factor = 1.0 / (max(patient.ureter_compliance_factor, 0.2) ** 0.25)
            reverse_resistance = (
                patient.uvj.effective_reverse_resistance()
                * uvj_reverse_deformation_mult
                * reverse_phase_factor
                * reverse_tortuosity_factor
                * reverse_compliance_factor
            )
            reflux_flow_ml_s = max(reflux_dp_pa / max(reverse_resistance, 1e-12), 0.0)

            side["volume_ml"] = max(
                side["volume_ml"]
                + (
                    patient.urine_production_ml_per_s
                    + reflux_flow_ml_s
                    - forward_flow_ml_s
                )
                * dt_s,
                0.2,
            )
            side["peak_pressure_pa"] = max(side["peak_pressure_pa"], pelvis_pa)
            side["forward_total_ml"] += forward_flow_ml_s * dt_s
            side["reflux_total_ml"] += reflux_flow_ml_s * dt_s
            side["forward_sum_ml_s"] += forward_flow_ml_s
            side["reflux_sum_ml_s"] += reflux_flow_ml_s
            side["per_gain_sum"] += per_gain
            side["prev_reflux_flow_ml_s"] = reflux_flow_ml_s

            if is_filling:
                side["filling_forward_ml"] += forward_flow_ml_s * dt_s
                side["filling_reflux_ml"] += reflux_flow_ml_s * dt_s
            else:
                side["voiding_forward_ml"] += forward_flow_ml_s * dt_s
                side["voiding_reflux_ml"] += reflux_flow_ml_s * dt_s

            if (
                forward_flow_ml_s < 0.5 * patient.urine_production_ml_per_s
                and pelvis_pa >= obstruction_pressure_threshold_pa
            ):
                side["low_forward_time_s"] += dt_s

            side["trace_pelvis"].append(pelvis_pa)
            side["trace_forward"].append(forward_flow_ml_s)
            side["trace_reflux"].append(reflux_flow_ml_s)
            side["trace_per_gain"].append(per_gain)

            combined_forward_flow_ml_s += forward_flow_ml_s
            combined_reflux_flow_ml_s += reflux_flow_ml_s
            combined_peak_pelvis_pa = max(combined_peak_pelvis_pa, pelvis_pa)
            combined_per_gain += per_gain

        if is_filling:
            urethral_outflow_ml_s = 0.0
        else:
            urethral_resistance = (
                urethral_resistance_pa_s_per_ml
                * (1.0 + bbd_severity * (bbd_tuning.outlet_resistance_gain - 1.0))
            )
            outflow_scale = _clamp(bladder_volume_ml / full_capacity_ml, 0.05, 1.2)
            urethral_outflow_ml_s = max(
                bladder_pa / max(urethral_resistance, 1e-12),
                0.0,
            ) * outflow_scale
            urethral_outflow_ml_s = min(
                urethral_outflow_ml_s,
                bladder_volume_ml / max(dt_s, 1e-9),
            )
            voiding_time_s += dt_s
            if urethral_outflow_ml_s < 0.8:
                low_voiding_outflow_time_s += dt_s

        bladder_delta_ml = (
            combined_forward_flow_ml_s
            - combined_reflux_flow_ml_s
            - urethral_outflow_ml_s
        ) * dt_s
        bladder_volume_ml = _clamp(
            bladder_volume_ml + bladder_delta_ml,
            min_residual_ml,
            1.4 * full_capacity_ml,
        )

        urethral_outflow_total_ml += urethral_outflow_ml_s * dt_s
        cumulative_passed_ml += urethral_outflow_ml_s * dt_s

        trace_time.append(t_s)
        trace_phase.append(phase)
        trace_bladder.append(bladder_pa)
        trace_bladder_volume.append(bladder_volume_ml)
        trace_cumulative_passed.append(cumulative_passed_ml)
        trace_combined_pelvis.append(combined_peak_pelvis_pa)
        trace_combined_forward.append(combined_forward_flow_ml_s)
        trace_combined_reflux.append(combined_reflux_flow_ml_s)
        trace_urethral.append(urethral_outflow_ml_s)
        trace_combined_per_gain.append(0.5 * combined_per_gain)

    low_voiding_outflow_fraction = (
        0.0 if voiding_time_s <= 0.0 else low_voiding_outflow_time_s / voiding_time_s
    )
    final_bladder_volume_ml = bladder_volume_ml

    left_result = _build_side_result(
        side_states["left"],
        steps=steps,
        total_time_s=total_time_s,
        full_capacity_ml=full_capacity_ml,
        initial_bladder_volume_ml=initial_bladder_volume_ml,
        final_bladder_volume_ml=final_bladder_volume_ml,
        cumulative_urine_passed_ml=cumulative_passed_ml,
        urethral_outflow_total_ml=urethral_outflow_total_ml,
        low_voiding_outflow_fraction=low_voiding_outflow_fraction,
        obstruction_pressure_threshold_pa=obstruction_pressure_threshold_pa,
        severe_pressure_threshold_pa=severe_pressure_threshold_pa,
    )
    right_result = _build_side_result(
        side_states["right"],
        steps=steps,
        total_time_s=total_time_s,
        full_capacity_ml=full_capacity_ml,
        initial_bladder_volume_ml=initial_bladder_volume_ml,
        final_bladder_volume_ml=final_bladder_volume_ml,
        cumulative_urine_passed_ml=cumulative_passed_ml,
        urethral_outflow_total_ml=urethral_outflow_total_ml,
        low_voiding_outflow_fraction=low_voiding_outflow_fraction,
        obstruction_pressure_threshold_pa=obstruction_pressure_threshold_pa,
        severe_pressure_threshold_pa=severe_pressure_threshold_pa,
    )

    combined_reflux_ml = left_result.reflux_volume_ml + right_result.reflux_volume_ml
    combined_antegrade_ml = left_result.antegrade_volume_ml + right_result.antegrade_volume_ml
    combined_transport_ml = max(combined_reflux_ml + combined_antegrade_ml, 1e-12)
    combined_filling_reflux_ml = (
        left_result.filling_reflux_volume_ml + right_result.filling_reflux_volume_ml
    )
    combined_filling_antegrade_ml = (
        left_result.filling_antegrade_volume_ml + right_result.filling_antegrade_volume_ml
    )
    combined_voiding_reflux_ml = (
        left_result.voiding_reflux_volume_ml + right_result.voiding_reflux_volume_ml
    )
    combined_voiding_antegrade_ml = (
        left_result.voiding_antegrade_volume_ml + right_result.voiding_antegrade_volume_ml
    )
    combined_filling_transport = max(
        combined_filling_reflux_ml + combined_filling_antegrade_ml,
        1e-12,
    )
    combined_voiding_transport = max(
        combined_voiding_reflux_ml + combined_voiding_antegrade_ml,
        1e-12,
    )
    combined_peak_pressure_pa = max(
        left_result.peak_renal_pelvis_pressure_pa,
        right_result.peak_renal_pelvis_pressure_pa,
    )
    combined_filling_reflux_fraction = combined_filling_reflux_ml / combined_filling_transport
    combined_reflux_fraction = combined_reflux_ml / combined_transport_ml

    combined_result = SimulationResult(
        reflux_volume_ml=combined_reflux_ml,
        antegrade_volume_ml=combined_antegrade_ml,
        reflux_fraction=combined_reflux_fraction,
        peak_renal_pelvis_pressure_pa=combined_peak_pressure_pa,
        mean_forward_flow_ml_s=left_result.mean_forward_flow_ml_s + right_result.mean_forward_flow_ml_s,
        mean_reflux_flow_ml_s=left_result.mean_reflux_flow_ml_s + right_result.mean_reflux_flow_ml_s,
        mean_effective_peristalsis_gain=0.5
        * (
            left_result.mean_effective_peristalsis_gain
            + right_result.mean_effective_peristalsis_gain
        ),
        voiding_urethral_outflow_ml=urethral_outflow_total_ml,
        filling_reflux_volume_ml=combined_filling_reflux_ml,
        voiding_reflux_volume_ml=combined_voiding_reflux_ml,
        filling_antegrade_volume_ml=combined_filling_antegrade_ml,
        voiding_antegrade_volume_ml=combined_voiding_antegrade_ml,
        filling_reflux_fraction=combined_filling_reflux_fraction,
        voiding_reflux_fraction=combined_voiding_reflux_ml / combined_voiding_transport,
        low_forward_fraction=max(left_result.low_forward_fraction, right_result.low_forward_fraction),
        estimated_full_bladder_capacity_ml=full_capacity_ml,
        initial_bladder_volume_ml=initial_bladder_volume_ml,
        final_bladder_volume_ml=final_bladder_volume_ml,
        cumulative_urine_passed_ml=cumulative_passed_ml,
        ureter_tortuosity_index=max(
            left_result.ureter_tortuosity_index,
            right_result.ureter_tortuosity_index,
        ),
        ureter_compliance_factor=0.5
        * (
            left_result.ureter_compliance_factor + right_result.ureter_compliance_factor
        ),
        vur_grade=_vur_grade_from_surrogate(
            combined_reflux_fraction,
            combined_peak_pressure_pa,
            combined_filling_reflux_fraction,
        ),
        obstruction_index=max(left_result.obstruction_index, right_result.obstruction_index),
        severe_obstruction=bool(
            left_result.severe_obstruction or right_result.severe_obstruction
        ),
    )

    left_trace = _build_side_trace(
        side_states["left"],
        trace_time=trace_time,
        trace_phase=trace_phase,
        trace_bladder=trace_bladder,
        trace_bladder_volume=trace_bladder_volume,
        trace_cumulative_passed=trace_cumulative_passed,
        trace_urethral=trace_urethral,
    )
    right_trace = _build_side_trace(
        side_states["right"],
        trace_time=trace_time,
        trace_phase=trace_phase,
        trace_bladder=trace_bladder,
        trace_bladder_volume=trace_bladder_volume,
        trace_cumulative_passed=trace_cumulative_passed,
        trace_urethral=trace_urethral,
    )
    combined_trace = SimulationTrace(
        time_s=trace_time,
        phase=trace_phase,
        bladder_pressure_pa=trace_bladder,
        pelvis_pressure_pa=trace_combined_pelvis,
        bladder_volume_ml=trace_bladder_volume,
        cumulative_urine_passed_ml=trace_cumulative_passed,
        forward_flow_ml_s=trace_combined_forward,
        reflux_flow_ml_s=trace_combined_reflux,
        urethral_outflow_ml_s=trace_urethral,
        effective_peristalsis_gain=trace_combined_per_gain,
    )

    if DECOUPLE_BBD_FROM_OBSTRUCTION and (
        _has_active_bbd(left_patient) or _has_active_bbd(right_patient)
    ):
        reference_bilateral = simulate_coupled_bilateral_with_trace(
            left_patient=_without_bbd(left_patient),
            right_patient=_without_bbd(right_patient),
            total_time_s=total_time_s,
            dt_s=dt_s,
            obstruction_pressure_threshold_pa=obstruction_pressure_threshold_pa,
            severe_pressure_threshold_pa=severe_pressure_threshold_pa,
        )
        left_result = replace(
            left_result,
            obstruction_index=reference_bilateral.left_result.obstruction_index,
            severe_obstruction=reference_bilateral.left_result.severe_obstruction,
        )
        right_result = replace(
            right_result,
            obstruction_index=reference_bilateral.right_result.obstruction_index,
            severe_obstruction=reference_bilateral.right_result.severe_obstruction,
        )
        combined_result = replace(
            combined_result,
            obstruction_index=reference_bilateral.combined_result.obstruction_index,
            severe_obstruction=reference_bilateral.combined_result.severe_obstruction,
        )

    return BilateralSimulation(
        combined_result=combined_result,
        combined_trace=combined_trace,
        left_result=left_result,
        right_result=right_result,
        left_trace=left_trace,
        right_trace=right_trace,
    )
