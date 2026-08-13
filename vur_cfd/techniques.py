from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from math import exp, pi, sqrt

from .model import PatientModel, UVJValve


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


ML_TO_MM3 = 1000.0
BASE_INTRAMURAL_RATIO = 5.0
MIN_EFFECTIVE_RADIUS_MM = 0.10
MAX_COMPRESSION = 0.88

# Global calibration knobs for obstruction sensitivity.
# - `FORWARD_EDEMA_GLOBAL_SCALE` scales the per-site edema narrowing term.
#   This is a primary calibration target: it directly controls predicted
#   post-injection obstruction risk for all bulking techniques.  The current
#   value of 12.0 is empirically chosen and should be validated against
#   postoperative obstruction rates when clinical outcome data become available.
# - `FORWARD_GEOMETRIC_NARROWING_SCALE` scales instantaneous geometric narrowing
#   from the local Poiseuille radius transform.
# - `FORWARD_RESISTANCE_MAX_MULTIPLIER` limits maximal UVJ forward resistance.
FORWARD_EDEMA_GLOBAL_SCALE = 12.0
FORWARD_GEOMETRIC_NARROWING_SCALE = 0.65
FORWARD_RESISTANCE_MAX_MULTIPLIER = 12.0
# Multi-site injections distribute focal outlet narrowing across sites.
# Default exponent 0.5 gives the previous 1/sqrt(n_sites) attenuation. Set to
# 0.0 only in sensitivity analyses to remove this structural advantage.
FORWARD_SITE_COUNT_ATTENUATION_EXPONENT = 0.5
# Submeatal (STING) placement is assigned greater outlet narrowing and weaker
# anti-reflux coupling than intraureteric (HIT / Double HIT) placement. This is
# a structural prior, not a fitted estimate. Set to False only in sensitivity
# analyses to remove the submeatal/intraureteric asymmetry.
PLACEMENT_WALL_ASYMMETRY_ENABLED = True

# Mound efficacy sigmoid: models the clinical observation that sub-threshold
# Deflux volumes fail to form a mechanically competent submucosal mound.
# The sigmoid scales effective lumen displacement per injection site, so that
# very small per-site volumes produce negligible anti-reflux effect.
# - MOUND_EFFICACY_MIDPOINT_ML: per-site volume (mL) at 50% mound efficacy.
# - MOUND_EFFICACY_STEEPNESS: sigmoid steepness parameter.
MOUND_EFFICACY_MIDPOINT_ML = 0.35
MOUND_EFFICACY_STEEPNESS = 12.0
MOUND_EFFICACY_REFERENCE_DIAMETER_MM = 1.4


def _mound_efficacy(per_site_volume_ml: float, orifice_diameter_mm: float = 1.4) -> float:
    """Sigmoid scaling: sub-threshold per-site volumes produce negligible mound effect.

    The midpoint scales linearly with orifice diameter so that smaller (infant)
    orifices reach competent mound efficacy at lower per-site volumes while
    larger (adolescent) orifices require more volume.
    """
    scaled_midpoint = MOUND_EFFICACY_MIDPOINT_ML * (
        orifice_diameter_mm / MOUND_EFFICACY_REFERENCE_DIAMETER_MM
    )
    return 1.0 / (1.0 + exp(-MOUND_EFFICACY_STEEPNESS * (per_site_volume_ml - scaled_midpoint)))


class TechniqueName(str, Enum):
    STING = "sting"
    HIT = "hit"
    DOUBLE_HIT = "double_hit"
    HIT_PLUS_STING = "hit_plus_sting"
    DOUBLE_HIT_PLUS_STING = "double_hit_plus_sting"
    COHEN = "cohen"
    POLITANO_LEADBETTER = "politano_leadbetter"


BULKING_TECHNIQUES: tuple[TechniqueName, ...] = (
    TechniqueName.STING,
    TechniqueName.HIT,
    TechniqueName.DOUBLE_HIT,
    TechniqueName.HIT_PLUS_STING,
    TechniqueName.DOUBLE_HIT_PLUS_STING,
)


@dataclass(frozen=True)
class TechniquePlan:
    technique: TechniqueName
    placement_quality: float = 0.75

    # Endoscopic bulking controls.
    deflux_volume_ml: float = 1.0

    # Reimplantation control.
    tunnel_length_mm: float = 7.0


def _bulking_component_coeffs(component: TechniqueName) -> dict[str, float]:
    if component == TechniqueName.STING:
        return {
            # Distal submeatal mound: short affected segment, modest lumen displacement.
            "length_fraction": 0.45,
            "displacement_efficiency": 0.20,
            "reverse_dynamic_gain": 0.18,
            "barrier_gain_scale_pa": 1800.0,
            "competence_gain_scale": 0.22,
            "forward_edema_penalty": 0.08,
        }
    if component == TechniqueName.HIT:
        return {
            # Intraureteric mound: longer intramural coupling and stronger coaptation.
            "length_fraction": 0.68,
            "displacement_efficiency": 0.28,
            "reverse_dynamic_gain": 0.32,
            "barrier_gain_scale_pa": 2600.0,
            "competence_gain_scale": 0.31,
            "forward_edema_penalty": 0.12,
        }
    if component == TechniqueName.DOUBLE_HIT:
        return {
            # Two intraureteric mounds: broad segment coverage and strongest recoil support.
            "length_fraction": 0.88,
            "displacement_efficiency": 0.34,
            "reverse_dynamic_gain": 0.46,
            "barrier_gain_scale_pa": 3300.0,
            "competence_gain_scale": 0.40,
            "forward_edema_penalty": 0.18,
        }
    raise ValueError(f"Unsupported bulking component: {component}")


def bulking_injection_layout(technique: TechniqueName) -> list[dict]:
    if technique == TechniqueName.STING:
        return [
            {
                "component": TechniqueName.STING.value,
                "fraction": 1.0,
                "label": "Distal mound (6 o'clock)",
                "location": "submeatal at ureteral orifice (6 o'clock)",
                "clockface_hour": 6.0,
                "axial_fraction": 0.08,
                "wall_plane": "submeatal",
            }
        ]
    if technique == TechniqueName.HIT:
        return [
            {
                "component": TechniqueName.HIT.value,
                "fraction": 1.0,
                "label": "HIT intramural tunnel mound (6 o'clock)",
                "location": "intraureteric within intramural detrusor tunnel",
                "clockface_hour": 6.0,
                "axial_fraction": 0.34,
                "wall_plane": "detrusor_tunnel",
            }
        ]
    if technique == TechniqueName.DOUBLE_HIT:
        return [
            {
                "component": TechniqueName.DOUBLE_HIT.value,
                "fraction": 0.5,
                "label": "Double HIT superior mound (6 o'clock)",
                "location": "more proximal intraureteric site within detrusor tunnel",
                "clockface_hour": 6.0,
                "axial_fraction": 0.56,
                "wall_plane": "detrusor_tunnel",
            },
            {
                "component": TechniqueName.DOUBLE_HIT.value,
                "fraction": 0.5,
                "label": "Double HIT distal intramural mound (6 o'clock)",
                "location": "distal intraureteric site within detrusor tunnel near meatus",
                "clockface_hour": 6.0,
                "axial_fraction": 0.30,
                "wall_plane": "detrusor_tunnel",
            },
        ]
    if technique == TechniqueName.HIT_PLUS_STING:
        return [
            {
                "component": TechniqueName.STING.value,
                "fraction": 0.5,
                "label": "STING distal mound (6 o'clock)",
                "location": "submeatal at ureteral orifice (6 o'clock)",
                "clockface_hour": 6.0,
                "axial_fraction": 0.08,
                "wall_plane": "submeatal",
            },
            {
                "component": TechniqueName.HIT.value,
                "fraction": 0.5,
                "label": "HIT intramural tunnel mound (6 o'clock)",
                "location": "intraureteric within intramural detrusor tunnel",
                "clockface_hour": 6.0,
                "axial_fraction": 0.34,
                "wall_plane": "detrusor_tunnel",
            },
        ]
    if technique == TechniqueName.DOUBLE_HIT_PLUS_STING:
        one_third = 1.0 / 3.0
        return [
            {
                "component": TechniqueName.STING.value,
                "fraction": one_third,
                "label": "STING distal mound (6 o'clock)",
                "location": "submeatal at ureteral orifice (6 o'clock)",
                "clockface_hour": 6.0,
                "axial_fraction": 0.08,
                "wall_plane": "submeatal",
            },
            {
                "component": TechniqueName.DOUBLE_HIT.value,
                "fraction": one_third,
                "label": "Double HIT superior mound (6 o'clock)",
                "location": "more proximal intraureteric site within detrusor tunnel",
                "clockface_hour": 6.0,
                "axial_fraction": 0.56,
                "wall_plane": "detrusor_tunnel",
            },
            {
                "component": TechniqueName.DOUBLE_HIT.value,
                "fraction": one_third,
                "label": "Double HIT distal intramural mound (6 o'clock)",
                "location": "distal intraureteric site within detrusor tunnel near meatus",
                "clockface_hour": 6.0,
                "axial_fraction": 0.30,
                "wall_plane": "detrusor_tunnel",
            },
        ]
    raise ValueError(f"Unsupported bulking layout technique: {technique}")


def _poiseuille_resistance_ratio(
    baseline_radius_mm: float,
    effective_radius_mm: float,
    length_ratio: float,
) -> float:
    r0 = max(baseline_radius_mm, MIN_EFFECTIVE_RADIUS_MM)
    r1 = max(effective_radius_mm, MIN_EFFECTIVE_RADIUS_MM)
    return max(length_ratio, 0.01) * (r0 / r1) ** 4


def _site_placement_multipliers(site: dict, quality: float) -> dict[str, float]:
    """
    Placement-aware modifiers for reduced-order bulking effects.

    - `axial_fraction`: 0.0 distal/orifice side, 1.0 proximal intramural.
    - `clockface_hour`: expected injection clock-face location around the lumen.
    - `wall_plane`: submeatal vs intramural detrusor-tunnel context.
    """
    axial = _clamp(float(site.get("axial_fraction", 0.30)), 0.0, 1.0)
    hour = float(site.get("clockface_hour", 6.0))
    wall_plane = str(site.get("wall_plane", "intraureteric")).strip().lower()

    # 6 o'clock is the intended anti-reflux injection position.
    hour_delta = abs(hour - 6.0)
    hour_delta = min(hour_delta, 12.0 - hour_delta)
    circumferential_focus = _clamp(1.0 - 0.08 * hour_delta, 0.45, 1.05)

    distal = exp(-((axial - 0.12) / 0.22) ** 2)
    mid = exp(-((axial - 0.40) / 0.20) ** 2)
    proximal = exp(-((axial - 0.74) / 0.24) ** 2)
    # Intramural tunnel occupancy proxy (higher means deeper in detrusor tunnel).
    tunnel_confinement = exp(-((axial - 0.46) / 0.28) ** 2)
    intratunnel_bias = _clamp(0.50 + 0.50 * tunnel_confinement, 0.50, 1.00)
    if not PLACEMENT_WALL_ASYMMETRY_ENABLED and wall_plane == "submeatal":
        # No-asymmetry sensitivity analysis: score submeatal sites with the same
        # wall-plane multipliers as intraureteric sites.
        wall_plane = "detrusor_tunnel"
    if wall_plane == "submeatal":
        displacement_wall = 0.82
        barrier_wall = 1.16
        reverse_wall = 0.92
        forward_wall = 1.18
        length_wall = 0.90
    elif wall_plane in {"detrusor_tunnel", "intraureteric"}:
        # HIT/Double-HIT mounds are intraureteric inside the detrusor tunnel
        # and therefore show stronger anti-reflux coupling for similar volume.
        displacement_wall = 1.00 + 0.20 * intratunnel_bias
        barrier_wall = 1.08 + 0.34 * intratunnel_bias
        reverse_wall = 1.10 + 0.48 * intratunnel_bias
        # Intramural intraureteric placement is more coaptive and less outlet narrowing.
        forward_wall = 0.90 + 0.15 * (1.0 - intratunnel_bias)
        length_wall = 1.02 + 0.28 * intratunnel_bias
    else:
        # Fallback for unknown placement labels: keep near intramural behavior.
        displacement_wall = 0.98 + 0.18 * intratunnel_bias
        barrier_wall = 1.05 + 0.28 * intratunnel_bias
        reverse_wall = 1.06 + 0.42 * intratunnel_bias
        forward_wall = 0.92 + 0.16 * (1.0 - intratunnel_bias)
        length_wall = 1.00 + 0.24 * intratunnel_bias

    return {
        "length_multiplier": _clamp(
            length_wall * (0.74 + 0.26 * mid + 0.50 * proximal),
            0.55,
            1.45,
        ),
        "displacement_multiplier": _clamp(
            displacement_wall
            * circumferential_focus
            * (0.78 + 0.32 * mid + 0.20 * distal + 0.18 * intratunnel_bias),
            0.45,
            1.45,
        ),
        "barrier_multiplier": _clamp(
            barrier_wall
            * circumferential_focus
            * (
                0.66
                + 0.48 * distal
                + 0.34 * proximal
                + 0.18 * mid
                + 0.24 * quality
            ),
            0.45,
            2.10,
        ),
        "reverse_multiplier": _clamp(
            reverse_wall
            * circumferential_focus
            * (
                0.70
                + 0.20 * distal
                + 0.44 * proximal
                + 0.22 * mid
                + 0.18 * intratunnel_bias
            ),
            0.50,
            2.35,
        ),
        "forward_multiplier": _clamp(
            forward_wall
            * (
                0.66
                + 0.46 * distal
                + 0.12 * mid
                + 0.14 * (1.0 - quality)
                - 0.12 * intratunnel_bias
            ),
            0.40,
            2.00,
        ),
    }


def _apply_bulking(patient: PatientModel, plan: TechniquePlan) -> PatientModel:
    quality = _clamp(plan.placement_quality, 0.0, 1.0)
    volume = max(plan.deflux_volume_ml, 0.0)
    site_layout = bulking_injection_layout(plan.technique)
    total_fraction = sum(float(site["fraction"]) for site in site_layout)
    if total_fraction <= 0.0:
        raise ValueError("Bulking site fractions must sum to a positive value.")
    # Distributed multi-site injections reduce focal narrowing accumulation.
    site_count_forward_norm = float(len(site_layout)) ** (
        -FORWARD_SITE_COUNT_ATTENUATION_EXPONENT
    )

    baseline_orifice_radius_mm = max(0.5 * patient.uvj_orifice_diameter_mm, 0.20)
    baseline_orifice_area_mm2 = pi * baseline_orifice_radius_mm * baseline_orifice_radius_mm
    intramural_length_mm = max(BASE_INTRAMURAL_RATIO * patient.uvj_orifice_diameter_mm, 3.0)

    # Resistances are accumulated as incremental ratios over the intramural segment.
    forward_multiplier_increment = 0.0
    reverse_multiplier_increment = 0.0
    barrier_gain_pa = 0.0
    competence_gain = 0.0

    for site in site_layout:
        component = TechniqueName(str(site["component"]))
        normalized_fraction = float(site["fraction"]) / total_fraction
        local_volume = volume * normalized_fraction
        coeffs = _bulking_component_coeffs(component)
        placement = _site_placement_multipliers(site, quality)

        affected_length_mm = max(
            coeffs["length_fraction"] * intramural_length_mm * placement["length_multiplier"],
            0.25,
        )
        wall_plane = str(site.get("wall_plane", "")).strip().lower()
        if wall_plane == "submeatal":
            # STING effects are local and concentrated near the distal meatus.
            affected_length_mm = min(affected_length_mm, 0.72 * intramural_length_mm)
        if wall_plane in {"detrusor_tunnel", "intraureteric"}:
            # HIT/Double-HIT effects are confined to the intramural detrusor tunnel.
            affected_length_mm = min(affected_length_mm, 0.98 * intramural_length_mm)
        length_ratio = affected_length_mm / intramural_length_mm
        efficacy = _mound_efficacy(local_volume, orifice_diameter_mm=baseline_orifice_radius_mm * 2)
        displaced_volume_mm3 = (
            local_volume
            * ML_TO_MM3
            * coeffs["displacement_efficiency"]
            * quality
            * placement["displacement_multiplier"]
            * efficacy
        )
        area_reduction_mm2 = displaced_volume_mm3 / max(affected_length_mm, 1e-9)
        raw_compression = area_reduction_mm2 / max(baseline_orifice_area_mm2, 1e-9)
        local_compression = _clamp(raw_compression, 0.0, MAX_COMPRESSION)

        # Better quality improves anti-reflux coaptation; poorer quality shifts toward outlet penalty.
        anti_reflux_compression = (
            local_compression
            * (0.65 + 0.35 * quality)
            * placement["barrier_multiplier"]
        )
        forward_penalty_compression = (
            local_compression
            * (0.55 + 0.65 * (1.0 - quality))
            * placement["forward_multiplier"]
        )
        reverse_collapse = _clamp(
            anti_reflux_compression
            * (1.0 + coeffs["reverse_dynamic_gain"] * quality)
            * placement["reverse_multiplier"],
            0.0,
            0.94,
        )

        forward_radius_mm = baseline_orifice_radius_mm * sqrt(
            max(1.0 - forward_penalty_compression, 0.02)
        )
        reverse_radius_mm = baseline_orifice_radius_mm * sqrt(
            max(1.0 - reverse_collapse, 0.01)
        )

        local_forward_ratio = _poiseuille_resistance_ratio(
            baseline_orifice_radius_mm,
            forward_radius_mm,
            length_ratio,
        )
        local_reverse_ratio = _poiseuille_resistance_ratio(
            baseline_orifice_radius_mm,
            reverse_radius_mm,
            length_ratio,
        )

        forward_multiplier_increment += (
            site_count_forward_norm
            * FORWARD_GEOMETRIC_NARROWING_SCALE
            * max(local_forward_ratio - 1.0, 0.0)
        )
        effective_edema_volume = local_volume * (0.10 + 0.90 * efficacy)
        forward_multiplier_increment += (
            site_count_forward_norm
            *
            FORWARD_EDEMA_GLOBAL_SCALE
            *
            coeffs["forward_edema_penalty"]
            * effective_edema_volume
            * (0.40 + 0.60 * (1.0 - quality))
        )
        reverse_multiplier_increment += max(local_reverse_ratio - 1.0, 0.0)

        barrier_gain_pa += (
            coeffs["barrier_gain_scale_pa"]
            * anti_reflux_compression
            * length_ratio
            * (0.50 + 0.50 * quality)
        )
        competence_gain += (
            coeffs["competence_gain_scale"]
            * anti_reflux_compression
            * (0.40 + 0.60 * quality)
            * (0.78 + 0.22 * placement["reverse_multiplier"])
        )

    uvj = patient.uvj
    forward_multiplier = _clamp(
        1.0 + forward_multiplier_increment,
        1.0,
        FORWARD_RESISTANCE_MAX_MULTIPLIER,
    )
    reverse_multiplier = _clamp(1.0 + reverse_multiplier_increment, 1.0, 16.0)
    new_uvj = UVJValve(
        forward_resistance_pa_s_per_ml=uvj.forward_resistance_pa_s_per_ml
        * forward_multiplier,
        reverse_resistance_pa_s_per_ml=uvj.reverse_resistance_pa_s_per_ml
        * reverse_multiplier,
        closure_barrier_pa=uvj.closure_barrier_pa + barrier_gain_pa,
        competence=_clamp(uvj.competence + competence_gain, 0.0, 1.0),
    )

    return replace(patient, uvj=new_uvj, name=f"{patient.name}+{plan.technique.value}")


def _apply_reimplant(patient: PatientModel, plan: TechniquePlan) -> PatientModel:
    quality = _clamp(plan.placement_quality, 0.0, 1.0)
    # Anti-reflux behavior is governed by submucosal tunnel to intravesical
    # ureteral orifice diameter ratio (classic target around 5:1).
    tunnel_ratio = plan.tunnel_length_mm / max(patient.uvj_orifice_diameter_mm, 1e-9)

    if plan.technique == TechniqueName.COHEN:
        reverse_factor = 0.95
        forward_penalty_factor = 0.85
        competence_gain_base = 0.35
    else:
        reverse_factor = 1.10
        forward_penalty_factor = 1.00
        competence_gain_base = 0.38

    ratio_term = max(tunnel_ratio - 3.0, 0.0)
    reverse_multiplier = 1.0 + 1.80 * reverse_factor * ratio_term * quality
    barrier_gain_pa = 700.0 + 550.0 * reverse_factor * ratio_term * quality

    # Penalize only when tunnel ratio exceeds common anti-reflux targets (around 5:1).
    excess_ratio = max(tunnel_ratio - 5.0, 0.0)
    narrowing_quality_factor = 0.45 + 0.55 * (1.0 - quality)
    forward_multiplier = (
        1.0 + 0.24 * forward_penalty_factor * excess_ratio * narrowing_quality_factor
    )

    uvj = patient.uvj
    new_uvj = UVJValve(
        forward_resistance_pa_s_per_ml=uvj.forward_resistance_pa_s_per_ml * forward_multiplier,
        reverse_resistance_pa_s_per_ml=uvj.reverse_resistance_pa_s_per_ml * reverse_multiplier,
        closure_barrier_pa=uvj.closure_barrier_pa + barrier_gain_pa,
        competence=_clamp(uvj.competence + competence_gain_base * quality, 0.0, 1.0),
    )
    return replace(patient, uvj=new_uvj, name=f"{patient.name}+{plan.technique.value}")


def apply_technique(patient: PatientModel, plan: TechniquePlan) -> PatientModel:
    if plan.technique in BULKING_TECHNIQUES:
        return _apply_bulking(patient, plan)
    if plan.technique in (TechniqueName.COHEN, TechniqueName.POLITANO_LEADBETTER):
        return _apply_reimplant(patient, plan)
    raise ValueError(f"Unsupported technique: {plan.technique}")
