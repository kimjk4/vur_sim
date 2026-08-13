#!/usr/bin/env python3
"""One-command verification of the published results.

Regenerates every quantitative claim in the manuscript from the model and
compares the result with the deposited output files.

    python3 verify.py            # all checks
    python3 verify.py --quick    # skip the two slow regenerations

The work is done in a temporary directory, so the deposited files in this
repository are read but never overwritten. Exit status is 0 if every check
passes, 1 otherwise.
"""
from __future__ import annotations

import argparse
import collections
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent

# Keep progress visible when stdout is redirected to a file or a pipe.
sys.stdout.reconfigure(line_buffering=True)

# (script, manuscript location, outputs compared byte-for-byte, slow?)
CHECKS = [
    (
        "docs/physics_regression_checks.py",
        "Methods - solver behaviour; Appendix A",
        [],
        False,
    ),
    (
        "docs/generate_absolute_volume_table.py",
        "Supplementary Table S1b; re-derives the 75-stratum base case and "
        "asserts it matches the deposited table",
        ["outputs/revision/balanced_selection_absolute_volumes.csv"],
        False,
    ),
    (
        "docs/generate_revision_sensitivity.py",
        "Appendix B8 (OI-threshold) and B9 (mound sigmoid)",
        [
            "outputs/revision/oi_threshold_sensitivity.csv",
            "outputs/revision/mound_sigmoid_sensitivity.csv",
            "outputs/revision/mound_sigmoid_sensitivity_by_stratum.csv",
            "outputs/revision/reviewer_response_evidence_table.csv",
        ],
        False,
    ),
    (
        "docs/generate_forward_discount_sensitivity.py",
        "Appendix B10 (no-site-discount); Results 4 and 6; Abstract",
        [
            "outputs/revision/forward_site_discount_sensitivity.csv",
            "outputs/revision/forward_site_discount_sensitivity_summary.csv",
            "outputs/revision/forward_site_discount_matched_volume_demo.csv",
        ],
        False,
    ),
    (
        "docs/generate_placement_asymmetry_sensitivity.py",
        "Appendix B11 (no-placement-asymmetry)",
        [
            "outputs/revision/placement_asymmetry_sensitivity.csv",
            "outputs/revision/placement_asymmetry_sensitivity_summary.csv",
            "outputs/revision/placement_asymmetry_grade5_sweep.csv",
        ],
        False,
    ),
    (
        "docs/refresh_ideal_tables.py",
        "Results 4 base-case selection; Appendix B6; Supplementary Table S1",
        ["outputs/manuscript_table_ideal_all_ages_balanced.csv"],
        True,
    ),
]

# Regenerating this one is expected to differ slightly; see README, "Known
# discrepancy". Compared numerically against a tolerance instead of byte-for-byte.
SWEEP_SCRIPT = "docs/optimal_obstruction_sweep.py"
SWEEP_CSV = "outputs/optimal_obstruction_sweep.csv"
SWEEP_TOLERANCE = 0.02
SWEEP_KEYS = ["age_group", "initial_grade", "bbd_severity_label", "technique", "deflux_volume_ml"]
SWEEP_NUMERIC = ["reflux_fraction", "obstruction_index"]
# The only stratum this file is tabulated for in the manuscript (Appendix B7).
SWEEP_REPORTED = {"age_group": "18_24m", "initial_grade": "5", "bbd_severity_label": "none"}


def rule(char: str = "-") -> None:
    print(char * 78)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def headline_numbers(work: Path) -> None:
    """Print the manuscript's headline counts, recomputed from regenerated files."""
    rule()
    print("HEADLINE NUMBERS, recomputed from the regenerated output files")
    rule()

    base = read_csv(work / "outputs/manuscript_table_ideal_all_ages_balanced.csv")
    counts = collections.Counter(row["recommended_technique"] for row in base)
    print(f"Results 4 / Appendix B6 - base-case selection over {len(base)} strata:")
    for technique, n in counts.most_common():
        print(f"    {technique:<24} {n:>3}/{len(base)}")

    oi = read_csv(work / "outputs/revision/oi_threshold_sensitivity.csv")
    per_ceiling: dict[str, list[int]] = collections.defaultdict(list)
    for row in oi:
        per_ceiling[row["obstruction_ceiling"]].append(int(row["changed_vs_0p15"]))
    print("\nAppendix B8 - strata whose selection changes vs the 0.15 ceiling:")
    for ceiling in sorted(per_ceiling):
        flags = per_ceiling[ceiling]
        if ceiling == "0.15":
            continue
        print(f"    OI ceiling {ceiling}      {sum(flags):>3}/{len(flags)}")

    fwd = read_csv(work / "outputs/revision/forward_site_discount_sensitivity_summary.csv")[0]
    print("\nAppendix B10 - removing the 1/sqrt(n) forward-narrowing attenuation:")
    print(f"    strata changed           {fwd['changed_strata_count']:>3}/{fwd['n_strata']}")
    print(
        f"    Double HIT selections    {fwd['baseline_double_hit_count']:>3}"
        f" -> {fwd['no_discount_double_hit_count']}"
    )

    pla = read_csv(work / "outputs/revision/placement_asymmetry_sensitivity_summary.csv")[0]
    print("\nAppendix B11 - removing the placement wall-asymmetry prior:")
    print(f"    strata changed           {pla['changed_strata_count']:>3}/{pla['n_strata']}")
    print(
        f"    Double HIT selections    {pla['baseline_double_hit_count']:>3}"
        f" -> {pla['no_asymmetry_double_hit_count']}"
    )
    print(
        f"    Double HIT + STING       {pla['baseline_double_hit_plus_sting_count']:>3}"
        f" -> {pla['no_asymmetry_double_hit_plus_sting_count']}"
    )


def compare_sweep(work: Path) -> bool:
    """Numeric comparison for the one file with a documented discrepancy."""
    deposited = {
        tuple(row[k] for k in SWEEP_KEYS): row for row in read_csv(REPO / SWEEP_CSV)
    }
    regenerated = read_csv(work / SWEEP_CSV)
    if len(deposited) != len(regenerated):
        print(f"  FAIL  row count {len(regenerated)} vs deposited {len(deposited)}")
        return False

    worst = {column: (0.0, None) for column in SWEEP_NUMERIC}
    worst_reported = {column: (0.0, None) for column in SWEEP_NUMERIC}
    for row in regenerated:
        key = tuple(row[k] for k in SWEEP_KEYS)
        old = deposited.get(key)
        if old is None:
            print(f"  FAIL  regenerated row absent from deposited file: {key}")
            return False
        in_reported = all(row[k] == v for k, v in SWEEP_REPORTED.items())
        for column in SWEEP_NUMERIC:
            delta = abs(float(row[column]) - float(old[column]))
            if delta > worst[column][0]:
                worst[column] = (delta, key)
            if in_reported and delta > worst_reported[column][0]:
                worst_reported[column] = (delta, key)

    print("  across all 7,125 rows, most of which the manuscript never tabulates:")
    ok = True
    for column, (delta, key) in worst.items():
        status = "ok" if delta <= SWEEP_TOLERANCE else "EXCEEDS TOLERANCE"
        print(f"    max |delta {column}| = {delta:.4f}  ({status})")
        if key is not None and delta > 0:
            print(f"        largest at {', '.join(f'{k}={v}' for k, v in zip(SWEEP_KEYS, key))}")
        ok = ok and delta <= SWEEP_TOLERANCE

    reported = ", ".join(f"{k}={v}" for k, v in SWEEP_REPORTED.items())
    print(f"  restricted to the stratum tabulated in Appendix B7 ({reported}):")
    for column, (delta, _) in worst_reported.items():
        print(f"    max |delta {column}| = {delta:.4f}")

    if ok:
        print(f"  PASS  within the documented tolerance of {SWEEP_TOLERANCE}")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="skip the two slow regenerations (base-case table and grade-V sweep)",
    )
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="vur-sim-verify-") as tmp:
        work = Path(tmp) / "repo"
        shutil.copytree(REPO, work, ignore=shutil.ignore_patterns(".git", "__pycache__"))

        failures: list[str] = []
        print(f"Python {sys.version.split()[0]}")
        print(f"Working copy: {work}")

        for script, location, compared, slow in CHECKS:
            if slow and args.quick:
                print(f"\nSKIPPED (--quick)  {script}")
                continue
            rule()
            print(f"{script}\n  {location}")
            started = time.time()
            result = subprocess.run(
                [sys.executable, script],
                cwd=work,
                env={**os.environ, "PYTHONPATH": str(work)},
                capture_output=True,
                text=True,
            )
            elapsed = time.time() - started
            if result.returncode != 0:
                print(f"  FAIL  exited {result.returncode} after {elapsed:.0f}s")
                print("  " + "\n  ".join(result.stdout.strip().splitlines()[-8:]))
                print("  " + "\n  ".join(result.stderr.strip().splitlines()[-8:]))
                failures.append(script)
                continue
            print(f"  ran in {elapsed:.0f}s")
            for line in result.stdout.splitlines():
                if "Verified" in line or "[OK]" in line or "Changed strata" in line:
                    print(f"  > {line.strip()}")
            for relative in compared:
                same = (work / relative).read_bytes() == (REPO / relative).read_bytes()
                print(f"  {'IDENTICAL' if same else 'DIFFERS  '}  {relative}")
                if not same:
                    failures.append(relative)

        if not args.quick:
            rule()
            print(f"{SWEEP_SCRIPT}\n  Figure 3 / Appendix B7 grade-V dose-response sweep")
            started = time.time()
            result = subprocess.run(
                [sys.executable, SWEEP_SCRIPT],
                cwd=work,
                env={**os.environ, "PYTHONPATH": str(work)},
                capture_output=True,
                text=True,
            )
            print(f"  ran in {time.time() - started:.0f}s")
            if result.returncode != 0:
                print(f"  FAIL  exited {result.returncode}")
                print("  " + "\n  ".join(result.stderr.strip().splitlines()[-8:]))
                failures.append(SWEEP_SCRIPT)
            elif not compare_sweep(work):
                failures.append(SWEEP_CSV)

        if not args.quick:
            headline_numbers(work)

        rule("=")
        if failures:
            print(f"RESULT: {len(failures)} check(s) failed")
            for item in failures:
                print(f"  - {item}")
            return 1
        print("RESULT: all checks passed")
        return 0


if __name__ == "__main__":
    sys.exit(main())
