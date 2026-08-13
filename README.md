# VUR-Sim — verification code

Reduced-order, CFD-informed pressure-flow simulation of endoscopic injection
therapy in pediatric vesicoureteral reflux.

This repository is the archival source code for:

> Kim JK, Roth J, Whittam B, Batra N, Chua M, Misseri R. Reduced-order,
> computational fluid dynamics-informed pressure-flow simulation of endoscopic
> injection therapy in pediatric vesicoureteral reflux: an in-silico evaluation
> of technique-dependent outcomes. *Journal of Pediatric Urology* (in press).

## No patient data

No patient data were used in the study and none are contained in this
repository. Every file here is source code, literature-derived parameter
documentation, or an output file produced by running that code. Nothing here
requires a data-use agreement.

## Verify everything with one command

```bash
git clone https://github.com/kimjk4/vur-sim.git && cd vur-sim
pip install numpy==1.24.4          # the only third-party dependency
python3 verify.py
```

`verify.py` copies the repository to a temporary directory, regenerates every
reported result from the model, compares each regenerated file byte-for-byte
against the deposited one, and prints the manuscript's headline numbers
recomputed from what it just produced. It exits non-zero if any check fails. The
deposited files in your clone are read but never overwritten.

Full run: about 7 minutes on a 2023 laptop, over half of it in the grade-V
dose-response sweep. `python3 verify.py --quick` skips the two slowest
regenerations and finishes in about 3 minutes.

The model contains no random number generation, so results are deterministic:
a correct run reproduces the published values exactly.

## What each check confirms

| Manuscript claim | Script | Checked against |
| --- | --- | --- |
| Solver behaviour invariants | `docs/physics_regression_checks.py` | internal assertions |
| Base-case selection over 75 strata — Double HIT 65/75, Double HIT + STING 9/75, HIT 1/75 (Results 4, Appendix B6, Supplementary Table S1) | `docs/refresh_ideal_tables.py` | `outputs/manuscript_table_ideal_all_ages_balanced.csv` |
| Absolute reflux and antegrade volumes per stratum (Supplementary Table S1b) | `docs/generate_absolute_volume_table.py` | `outputs/revision/balanced_selection_absolute_volumes.csv` |
| OI-threshold sensitivity — 55/75 strata change at a 0.10 ceiling (Appendix B8, Results 6) | `docs/generate_revision_sensitivity.py` | `outputs/revision/oi_threshold_sensitivity.csv` |
| Mound-sigmoid sensitivity (Appendix B9) | `docs/generate_revision_sensitivity.py` | `outputs/revision/mound_sigmoid_sensitivity*.csv` |
| No-site-discount sensitivity — Double HIT falls 65 → 27, 50/75 strata change (Appendix B10, Results 4 and 6, Abstract) | `docs/generate_forward_discount_sensitivity.py` | `outputs/revision/forward_site_discount_sensitivity*.csv` |
| No-placement-asymmetry sensitivity — 16/75 strata change, Double HIT 65 → 52, Double HIT + STING 9 → 20 (Appendix B11) | `docs/generate_placement_asymmetry_sensitivity.py` | `outputs/revision/placement_asymmetry_*.csv` |
| Grade-V dose-response sweep (Figure 3 data, Appendix B7) | `docs/optimal_obstruction_sweep.py` | `outputs/optimal_obstruction_sweep.csv` |

Two of these scripts also verify themselves and fail rather than write divergent
numbers: `generate_absolute_volume_table.py` asserts that all 75 regenerated
strata match the deposited base-case table on technique, volume, reflux fraction
and obstruction index; `generate_placement_asymmetry_sensitivity.py`
cross-checks its base-case grade-V sweep against `optimal_obstruction_sweep.csv`.

To run any single script on its own:

```bash
PYTHONPATH=. python3 docs/generate_forward_discount_sensitivity.py
```

Note that scripts write into `outputs/`, overwriting the deposited copy of
whatever they regenerate. `verify.py` avoids this by working on a temporary copy.
Steps that depend on another script's output are `generate_forward_discount_sensitivity.py`
and `generate_placement_asymmetry_sensitivity.py`, which read the OI = 0.15
baseline from `outputs/revision/oi_threshold_sensitivity.csv` — deposited here, so
they run on a fresh clone without any prerequisite.

## The model

| File | Contents | Dependencies |
| --- | --- | --- |
| `vur_cfd/model.py` | Pressure-flow solver, bladder/pelvis/ureter/UVJ compartments, grade templates, BBD profiles, obstruction index | standard library only |
| `vur_cfd/techniques.py` | Bulking-agent techniques, injection layouts, mound sigmoid, placement multipliers, structural-prior switches | standard library only |

Both are byte-identical to the versions used to produce the published results.

Two module-level switches in `vur_cfd/techniques.py` encode the structural priors
examined in the sensitivity analyses. Both default to the base-case behaviour and
are toggled only by the sensitivity scripts:

- `FORWARD_SITE_COUNT_ATTENUATION_EXPONENT` (default `0.5`, the 1/sqrt(n_sites)
  forward-narrowing attenuation; set to `0.0` to remove it — Appendix B10)
- `PLACEMENT_WALL_ASYMMETRY_ENABLED` (default `True`; set to `False` to score
  submeatal injection sites with the same wall-plane multipliers as
  intraureteric sites — Appendix B11)

## Documentation

| File | Contents |
| --- | --- |
| `docs/methods_cs_replication.md` | Full solver specification: state variables, time discretization, and governing equations, sufficient to re-implement the model independently |
| `docs/methods_pediatric_urology_parameters.md` | Every default numeric input with its literature source |
| `docs/citations.md` | Citation ledger keyed to the parameter IDs used above |

## Known discrepancy

`outputs/optimal_obstruction_sweep.csv`, the source of Figure 3 and Appendix B7,
was generated before `FORWARD_SITE_COUNT_ATTENUATION_EXPONENT` was introduced in
`vur_cfd/techniques.py`, so it is the one deposited file that does not regenerate
byte-for-byte. It is deposited as-is, because it is the file the published figure
was drawn from.

Across all 7,125 rows the largest deviations are 0.0153 in reflux fraction
(12-18 months, grade V, severe BBD, HIT + STING at 2.0 mL) and 0.0098 in
obstruction index. Almost none of those rows are tabulated anywhere. Restricted
to the single stratum Appendix B7 does report — 18-24 months, female reference,
no BBD, grade V — the largest deviations are 0.0058 and 0.0009, which move four
values in the published B7 table by 0.001 in the third decimal:

| Technique | Published RF range | Regenerated | Published OI range | Regenerated |
| --- | --- | --- | --- | --- |
| STING | 0.233-0.249 | 0.234-0.248 | 0.269-0.315 | unchanged |
| HIT | 0.047-0.059 | 0.047-0.058 | 0.136-0.205 | unchanged |
| Double HIT | 0.021-0.000 | unchanged | 0.084-0.244 | unchanged |
| Double HIT + STING | 0.010-0.000 | unchanged | 0.079-0.298 | 0.079-0.299 |

No interpretation in that table changes, and no other reported result depends on
this file. `verify.py` therefore compares it numerically against a 0.02 tolerance
rather than byte-for-byte, and prints both the whole-file and the
Appendix-B7-stratum deviations, so the size of the gap is visible rather than
something to take on trust.

## Deliberately not included

- **Figure and table rendering scripts.** They reproduce presentation, not
  results. The data behind every figure and table is deposited here, and the
  numbers those figures depict are what `verify.py` checks.
- **The browser demonstration** at <https://vur-sim.vercel.app>, and its source.
  The deployment is a convenience demonstration only; **this repository, not the
  deployment, is the archival source.** It reimplements no physics — it renders
  values produced by `vur_cfd/model.py`.
- **`vur_cfd/cohort.py`, `main.py` (CLI), `visualization.py`,
  `anatomy_drawing.py`, `interactive_server.py`** — cohort-run helpers, the
  command-line entry point, and animation/rendering code. None is used by any
  analysis reported in the paper. Because they are absent, `vur_cfd/__init__.py`
  no longer re-exports them; that is the only difference between this package and
  the development version, and it touches no solver or technique code.

## Important limits

- This is **not** a clinical decision tool.
- The surrogate VUR grade is derived from hemodynamic outputs; it is not a VCUG
  image classifier.
- Technique-effect coefficients are explicit, tunable calibration priors. They
  have not been validated against patient-level postoperative outcomes.
- Tortuosity, compliance, and peristalsis models are phenomenologic.

## Citation

See `CITATION.cff`. Please cite both the article and this software record.

## License

MIT. See `LICENSE`.
