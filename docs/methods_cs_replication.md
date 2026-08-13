# Technical Methods for Computer Science Replication

## 1. Scope
This document specifies how to reproduce the current VUR simulator outputs (single-case and cohort mode) from the codebase in:
this repository.

Model class:
- reduced-order, explicit time-marching pressure-flow solver,
- unilateral mode,
- coupled bilateral mode with one shared bladder state,
- deformation-aware UVJ resistance update (elliptical Poiseuille scaling),
- BBD-aware bladder/outlet perturbation module.

Primary modules:
- `vur_cfd/model.py`
- `vur_cfd/techniques.py`
- `vur_cfd/interactive_server.py`

## 2. State and Time Discretization
For each timestep `k` with step `dt`:
- pelvis volume `V_pelvis[k]`
- bladder volume `V_bladder[k]`
- bladder pressure `P_bladder[k]`
- pelvis pressure `P_pelvis[k]`
- peristaltic pressure `P_per[k]`
- antegrade flow `Q_fwd[k]`
- reflux flow `Q_ref[k]`
- urethral outflow `Q_ur[k]`.

Cycle timing (`BladderCycle`):
- filling duration = `8.0 s`
- voiding duration = `4.0 s`
- cycle = `12.0 s`.

Pressure unit conversion:
- `1 cmH2O = 98.0665 Pa`.

## 3. Unilateral Solver Equations
### 3.1 Bladder pressure with volume-dependent scaling and BBD
`P_bladder_raw = bladder_cycle.pressure_at(t)`.

Let `stretch = clamp(V_bladder/full_capacity, 0, 1.5)`.

Phase scaling:
- filling: `phase_scale = 0.20 + 0.80*stretch`
- voiding: `phase_scale = 0.35 + 0.65*stretch`.

Then:
`P_bladder = P_baseline + phase_scale*(P_bladder_raw - P_baseline)`.

If BBD is enabled (`bbd_severity in [0,1]`), phase-scale is further multiplied by profile-dependent gains:
- filling gain (overactivity-dominant),
- voiding gain (dysfunctional-voiding dominant),
and baseline offset is increased by:
`P_bladder += P_baseline * bbd_severity * baseline_pressure_gain`.

For filling, BBD can add smooth detrusor contractions:
`spike_wave = 0.5 * (1 - cos(2*pi*f*t))`.
`P_spike = 0.70 * bbd_severity * filling_spike_pa * spike_wave`.
`P_bladder += P_spike`.

The raised-cosine waveform replaces a half-wave rectified sine, avoiding sharp zero-crossings. The 0.70 factor compensates for the higher mean of the smoother waveform.

### 3.2 Collecting-system pressure (nonlinear compliance)
`C_eff = renal_compliance_ml_per_pa * clamp(ureter_compliance_factor, 0.35, 3.5)`.

Linear term:
`linear = (V_pelvis - V0_pelvis) / C_eff`.

Stiffening term (Whitaker-type):
`excess = max(V_pelvis - V_stiffening_threshold, 0)`.
`stiffening = alpha * excess`.

Combined:
`P_pelvis = P0_pelvis + linear * (1 + stiffening)`.

Defaults:
- `pelvis_stiffening_alpha = 0.25`
- `pelvis_stiffening_volume_threshold_ml = 12.0`.

This nonlinear extension prevents underestimation of peak renal pelvis pressure in high-reflux scenarios (grades IV-V) where pelvis distension causes progressive wall stiffening, consistent with Whitaker perfusion study observations.

### 3.3 Peristalsis waveform and gain

#### 3.3.1 Waveform
Peristalsis is modeled as a continuous raised-cosine waveform with power-shaping and a tonic baseline floor:
```
phase = ((t * frequency_hz) + phase_offset) % 1.0
wave = 0.5 * (1 - cos(2*pi*phase))              # raised cosine, 0→1
sharpness = max(duty_cycle, 0.1) / 0.5
shaped = wave ^ (1 / max(sharpness, 0.2))
P_wave = amplitude * (baseline_fraction + (1 - baseline_fraction) * shaped)
```

Parameters:
- `baseline_fraction = 0.15` (tonic ureteral tone floor, fraction of amplitude),
- `duty_cycle = 0.35` (retained for backward compatibility; controls peak sharpness via the power exponent),
- `phase_offset = 0.0` for unilateral; `0.0` for left and `0.5` for right in bilateral mode (half-cycle offset models physiologic asynchrony of bilateral peristalsis).

The waveform smoothly oscillates and never hard-zeros, producing the rhythmic undulations seen in physiologic ureteral peristalsis rather than spike-then-silence impulses.

#### 3.3.2 Compliance penalty (asymmetric)
Compliance penalty uses an asymmetric exponential to reflect distinct impairment mechanisms for stiff versus floppy ureters:
```
if compliance >= 1.0:
    penalty = exp(-0.60 * (compliance - 1.0))    # floppy/dilated wall: poor tone
else:
    penalty = exp(-0.90 * (1.0 - compliance))    # stiff/scarred wall: poor contraction
```

Full gain:
`gain = clamp(per_eff * penalty / tortuosity^0.35, 0.03, 1.6)`.

`P_per = P_wave * gain`.

Note on tortuosity exponent: the `0.35` exponent on peristalsis gain is a calibration parameter representing reduced propulsive efficiency in tortuous ureters. It is empirically chosen and should be documented as such. A theoretical validation path exists via Dean-number secondary flow corrections for typical pediatric ureteral geometries.

### 3.4 Flow equations with deformation-aware UVJ
Antegrade:
`passive_floor = urine_production_rate * 0.1`.
`DeltaP_forward = P_pelvis + P_per - P_bladder`.
`Q_fwd_raw = DeltaP_forward/(R_ureter + R_uvj_fwd*F_uvj_fwd)`.
`Q_passive = passive_floor * sigmoid((DeltaP_forward - P_mid)/P_scale)`.
`Q_fwd = max(Q_fwd_raw + Q_passive, 0)`.

`P_mid = 1 cmH2O` and `P_scale = 1.5 cmH2O` in the implementation.
This pressure-gated passive term captures gravity/elastance-driven trickle near neutral gradients, but suppresses basal flow when the bladder-to-pelvis gradient is clearly adverse.

Reflux:
`Q_ref = max((P_bladder - P_pelvis - Barrier_eff - Barrier_deform)/R_rev_eff, 0)`.

Where:
- `Barrier_eff = closure_barrier * competence`
- `R_rev_eff = uvj.effective_reverse_resistance * F_uvj_rev * phase_factor * tortuosity_factor * compliance_factor`
- `phase_factor` = smooth sigmoid ramp between `1.0` (filling) and `0.65 + 0.35*competence` (voiding), with a `1.0 s` ramp centered on the filling→voiding boundary (see §3.6)
- `tortuosity_factor = 1 + 0.20*max(tortuosity-1, 0)`
- `compliance_factor = 1/(max(compliance, 0.2)^0.25)`.

`F_uvj_fwd`, `F_uvj_rev`, and `Barrier_deform` are computed from bladder fill-dependent UVJ deformation:
- reference state at ~10% capacity,
- UVJ major/minor axes and intramural length scale with stretch,
- resistance term uses elliptical Poiseuille extension:
`R ~ l*(a^2+b^2)/(a^3*b^3)`.

Reflux feedback on UVJ geometry (§3.7): retrograde flow pressure also contributes to UVJ distension from the ureteral side, creating a positive feedback loop in high-grade reflux.

Urethral outflow:
- filling: `Q_ur = 0`
- voiding:
`Q_ur = max(P_bladder / R_urethra, 0) * outflow_scale`,
`outflow_scale = clamp(V_bladder/full_capacity, 0.05, 1.2)`,
with cap `Q_ur <= V_bladder/dt`.

If BBD is enabled, voiding outlet resistance is scaled:
`R_urethra_eff = R_urethra * (1 + bbd_severity*(outlet_resistance_gain - 1))`.

### 3.5 Volume updates
Pelvis:
`V_pelvis[k+1] = max(V_pelvis + (urine_prod + Q_ref - Q_fwd)*dt, 0.2)`.

Bladder:
`V_bladder[k+1] = clamp(V_bladder + (Q_fwd - Q_ref - Q_ur)*dt, V_residual_min, 1.4*full_capacity)`.

For BBD, `V_residual_min = bbd_severity * min_residual_fraction * full_capacity`.

### 3.6 Phase transition smoothing
The reverse_phase_factor transitions smoothly between filling (`1.0`) and voiding (`0.65 + 0.35*competence`) over a `1.0 s` ramp centered on the filling→voiding boundary.

Let `voiding_target = 0.65 + 0.35 * competence`.
Let `local_t = t mod cycle_duration`.

- If filling and `time_to_voiding < 0.5 * ramp_duration`:
  `alpha = 0.5 - time_to_voiding / ramp_duration`
  `factor = 1.0 + alpha * (voiding_target - 1.0)`.
- If voiding and `time_since_voiding < 0.5 * ramp_duration`:
  `alpha = 0.5 + time_since_voiding / ramp_duration`
  `factor = 1.0 + alpha * (voiding_target - 1.0)`.
- Otherwise: `1.0` (filling) or `voiding_target` (voiding).

This eliminates the step discontinuity in effective reflux resistance at the phase boundary and produces smoother pressure-flow traces.

### 3.7 Reflux feedback on UVJ geometry
In addition to bladder-fill-driven UVJ deformation, retrograde flow exerts pressure on the UVJ from the ureteral side, potentially distending it further. This creates a positive feedback loop: reflux → UVJ distension → reduced resistance → more reflux.

Implementation uses the previous timestep's reflux flow:
```
reflux_pressure = Q_reflux_prev * R_uvj_reverse
reflux_stretch = clamp(reflux_pressure / 2000.0, 0, 1)
susceptibility = 0.15 * (1 - competence)
minor_scale *= (1 + susceptibility * reflux_stretch)
```

Clamped to `[0.24, 1.55]`. The effect is strongest in low-competence (high-grade) junctions and negligible at grade 0-I where the UVJ is nearly fully competent.

## 4. Coupled Bilateral Solver
Implemented in `simulate_coupled_bilateral_with_trace(...)`.

Left and right sides are integrated simultaneously with:
- shared `P_bladder(t)` and `V_bladder(t)`,
- side-specific pelvis states and side-specific UVJ/peristalsis/geometry,
- left and right peristaltic waves phase-offset by half a cycle (`phase_offset=0.0` for left, `0.5` for right), modeling the physiologic asynchrony of bilateral ureteral contractions.

Shared bladder construction:
- baseline pressure = mean(left, right),
- end-filling and peak-voiding pressure = max(left, right),
- urethral resistance = mean(left, right),
- full capacity = mean(left, right),
- initial volume = mean(left, right initial volumes).

Per-step bilateral update:
1. compute `Q_fwd,left`, `Q_ref,left`, `Q_fwd,right`, `Q_ref,right`;
2. sum flows:
`Q_fwd,sum = Q_fwd,left + Q_fwd,right`,
`Q_ref,sum = Q_ref,left + Q_ref,right`;
3. compute one shared `Q_ur`;
4. update shared bladder using summed flows;
5. update each pelvis state with its side-specific equation.

Outputs:
- side results (`left_result`, `right_result`),
- side traces (`left_trace`, `right_trace`),
- combined result/trace.

Combined result rules:
- reflux and antegrade volumes are summed,
- peak pelvis pressure = `max(left_peak, right_peak)`,
- secondary surrogate grade from combined reflux fractions + combined peak pressure,
- obstruction index = `max(left_index, right_index)`,
- severe obstruction = `left_severe OR right_severe`.

## 5. Parameterization and Defaults
### 5.1 Anatomy and pressure profiles
Age/sex tables are hard-coded in:
- `AGE_ANATOMY_PROFILES`
- `AGE_PRESSURE_PROFILES`.

### 5.2 Ureter resistance model
`R_ureter = R_ref * (L/L_ref) * (D_ref/D)^4 * tortuosity^1.25 / compliance^0.55`.

Constants:
- `R_ref = 90.0`
- `L_ref = 137.5 mm`
- `D_ref = 3.3 mm`.

Note on tortuosity exponent: the `1.25` exponent on resistance is a calibration parameter approximating combined effects of increased effective path length and secondary (Dean-type) flows in tortuous ureters. For typical pediatric Dean numbers, a theoretical friction-factor correction of `1 + 0.033*De^0.5` supports an effective exponent in the range 1.1-1.4.

### 5.3 Grade initialization
Pre-op grade template (`0..5`) controls:
- competence,
- reverse resistance multiplier,
- barrier multiplier,
- dilation multiplier,
- baseline tortuosity/compliance/peristalsis multipliers.

Optional infant-only switch:
- `use_grade_voiding_pressure_multiplier` is a legacy switch for <12 month runs; current simulator age groups start at `12_18m`, so this is inactive by default and in typical use.

Tortuosity constraints:
- grade `<=3`: forced `1.0`
- grade `4`: clamped `1.0-1.35`
- grade `5`: clamped `1.0-3.5`.

Grade-standard peristalsis defaults (smooth raised-cosine waveform, `baseline_fraction=0.15`):
- g0 `(0.90 Hz, 230 Pa)`
- g1 `(0.85 Hz, 210 Pa)`
- g2 `(0.80 Hz, 195 Pa)`
- g3 `(0.72 Hz, 170 Pa)`
- g4 `(0.62 Hz, 145 Pa)`
- g5 `(0.52 Hz, 115 Pa)`.

Amplitudes are calibrated for the smooth raised-cosine waveform (higher mean than the previous duty-cycled pulse).

## 6. Deflux Technique Transform
Interactive mode currently allows only bulking techniques:
- `sting`
- `hit`
- `double_hit`
- `hit_plus_sting`
- `double_hit_plus_sting`.

Site layout:
- STING: 1 site
- HIT: 1 site
- Double HIT: 2 sites (0.5/0.5)
- HIT+STING: 2 sites (0.5/0.5)
- Double HIT+STING: 3 sites (1/3 each).

Current site-placement map (hard-coded in `bulking_injection_layout`):
- STING: `wall_plane=submeatal`, `clockface=6`, `axial_fraction=0.08`.
- HIT: `wall_plane=detrusor_tunnel`, `clockface=6`, `axial_fraction=0.34`.
- Double HIT superior: `wall_plane=detrusor_tunnel`, `clockface=6`, `axial_fraction=0.56`.
- Double HIT distal intramural: `wall_plane=detrusor_tunnel`, `clockface=6`, `axial_fraction=0.30`.
- Combination techniques reuse these sites and split volume by `fraction`.

Placement-aware local-length capping:
- submeatal STING effects are capped to `0.72 * intramural_length`,
- detrusor-tunnel HIT/Double-HIT effects are capped to `0.98 * intramural_length`.

Component coefficients (`techniques.py`) are physics-oriented and now encode:
- affected intramural-length fraction,
- volume-to-lumen displacement efficiency,
- reverse dynamic gain,
- barrier and competence gains,
- forward edema penalty.

Calibrated global narrowing knobs (current build):
- `FORWARD_EDEMA_GLOBAL_SCALE = 12.0`
- `FORWARD_GEOMETRIC_NARROWING_SCALE = 0.65`
- `FORWARD_RESISTANCE_MAX_MULTIPLIER = 12.0`

Current per-component magnitudes are implemented directly in code for:
- STING,
- HIT,
- DOUBLE_HIT,
and reused via site-layout splitting for combination techniques.

For each site, with local volume and placement quality, additive gains are applied to UVJ parameters. In bilateral mode, left and right techniques are applied independently before bilateral solving.

## 7. Derived Outcomes
### 7.1 Primary efficacy endpoint: reflux metrics
`reflux_fraction = reflux_volume / (reflux_volume + antegrade_volume)`.

Also computed separately for filling and voiding phases.

### 7.2 Secondary classifier: surrogate grade
Score:
- reflux component: `clamp((reflux_fraction - 0.01)/0.60, 0, 1)`
- pressure component: `clamp((peak_pressure - 35cmH2O)/90cmH2O, 0, 1)`
- filling component: `clamp(filling_reflux_fraction/0.45, 0, 1)`.

Weighted sum:
`score = 0.55*reflux + 0.20*pressure + 0.25*filling`.

Bins:
- `<0.05 -> 0`
- `<0.12 -> 1`
- `<0.25 -> 2`
- `<0.37 -> 3`
- `<0.45 -> 4`
- else `5`.

### 7.3 Co-primary safety endpoint: obstruction index (normalized weights)
`obstruction_index = 0.533*pressure_risk + 0.133*low_forward_fraction + 0.100*low_voiding_outflow_fraction + 0.034*tortuosity_risk + 0.200*edema_risk`.

Weights are normalized to sum to 1.0 (relative ratios 16:4:3:1:6 preserved from original 0.80:0.20:0.15:0.05:0.30).

`edema_risk` is computed from post-technique UVJ forward-resistance multiplier:
- `edema_risk = clamp((forward_resistance_multiplier - 1.0)/4.0, 0.0, 2.0)`.

BBD decoupling rule:
- if BBD is active, `obstruction_index` and `severe_obstruction` are replaced by values from a BBD-neutral reference simulation (`bbd_profile=none`, `bbd_severity=0`) using identical geometry/technique settings.

Thresholds:
- obstruction pressure anchor: `30 cmH2O`
- severe pressure anchor: `45 cmH2O`
- severe flag if `peak_pressure >= severe_threshold` or `obstruction_index >= 1.0`.

## 8. API Input Schema and Bilateral Fields
`interactive_server._run_simulation(payload)` accepts:
- common: `age_group`, `sex`, `laterality`, `capacity_method`, `sim_time`, `dt`, `max_points`, `bladder_fill_fraction`
- left-side defaults: `initial_grade`, `technique`, `deflux_volume`, `placement_quality`, `tortuosity_index`, `ureter_compliance`, `peristalsis_efficiency`
- side-explicit:
`left_initial_grade`, `right_initial_grade`,
`right_technique`, `right_deflux_volume`, `right_placement_quality`,
`left_tortuosity_index`, `right_tortuosity_index`,
`left_ureter_compliance`, `right_ureter_compliance`,
`left_peristalsis_efficiency`, `right_peristalsis_efficiency`.

Response:
- unilateral: `patient/result/trace`
- bilateral: `patient/result/trace` plus `sides.left.*` and `sides.right.*`.

## 9. Cohort Mapping Conventions
For Deflux database tables used in this project:
- sex coding: `1 -> male`, `2 -> female`
- technique coding:
  - `1 -> hit`
  - `2 -> double_hit`
  - `3 -> sting`
  - `4 -> hit_plus_sting`
  - `5 -> double_hit_plus_sting`
- pre-op grade mapped from `VUR_Grade` and clamped to `[0,5]`
- deflux volume clamped to `[0,4] mL`.

Missing-outcome handling is analysis-dependent:
- primary tables may use non-missing denominators,
- sensitivity runs may recode missing outcomes to zero (no-event), as specified.

## 10. Reproducible Calls
### 10.1 CLI smoke test
```bash
python3 -m vur_cfd.main \
  --technique hit_plus_sting \
  --age-group 18_24m \
  --sex female \
  --initial-grade 4 \
  --sim-time 24 \
  --dt 0.05
```

### 10.2 Bilateral API call
```python
from vur_cfd.interactive_server import _run_simulation

out = _run_simulation({
    "age_group": "18_24m",
    "sex": "female",
    "laterality": "bilateral",
    "capacity_method": "koff",
    "left_initial_grade": 4,
    "right_initial_grade": 3,
    "technique": "double_hit_plus_sting",
    "right_technique": "hit",
    "deflux_volume": 1.2,
    "right_deflux_volume": 0.9,
    "placement_quality": 0.80,
    "right_placement_quality": 0.70,
    "left_tortuosity_index": 2.2,
    "right_tortuosity_index": 1.0,
    "left_ureter_compliance": 1.1,
    "right_ureter_compliance": 0.9,
    "left_peristalsis_efficiency": 0.9,
    "right_peristalsis_efficiency": 1.0,
    "bladder_fill_fraction": 1.0,
    "sim_time": 24.0,
    "dt": 0.05,
    "max_points": 700,
})
```

## 11. Generated Artifacts
Key output tables are under:
`outputs/`

including:
- calibration tables,
- per-patient predictions,
- grade/tortuosity/volume sweep outputs.
- publication figures in `outputs/figures/`, including:
  - `Figure3_deflux_physics.*` (submeatal vs intramural tunnel placement),
  - `Figure4_3D_render_bilateral.*` (upright kidneys/ureters, transparent bladder wall, intramural UVJ tunnel),
  - `Figure4_3D_render_rotation.gif`.

## 12. Limits
- reduced-order model, not voxel/mesh-resolved FSI,
- technique coefficients are calibration priors, requiring future validation against institutional postoperative outcomes when available,
- event prediction quality improves when combined with non-hydrodynamic covariates (e.g., BBD, recurrent UTI history, bilateral burden),
- quasi-steady flow equations (no inertial `L·dQ/dt` terms); acceptable at pediatric ureteral Womersley numbers but may slightly overpredict instantaneous reflux onset during rapid voiding pressure rise,
- constant urine production rate (`0.03 mL/s`); acceptable over the short simulation window but does not capture GFR changes under renal backpressure,
- surrogate grade bin thresholds and forward-edema global scale (`12.0`) are hand-tuned for in silico plausibility and should be treated as calibration targets for future validation,
- tortuosity exponents (resistance `1.25`, peristalsis `0.35`) are empirically chosen calibration parameters.

## 13. Citation Linkage
Source tracking for parameter priors is maintained in:
`docs/citations.md`
