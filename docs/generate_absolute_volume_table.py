"""Per-stratum balanced selections with absolute reflux and antegrade volumes.

Reviewer 2 (round 2) noted that the manuscript correctly cautions that reflux
fraction can rise when antegrade flow is suppressed, and should therefore be
read alongside absolute volumes - yet Figure 4 and Supplementary Table S1
report only RF and OI. This script regenerates the 75-stratum balanced
selection and additionally reports mean absolute reflux and antegrade volumes
for each selected strategy.

Selection logic is imported verbatim from `refresh_ideal_tables.py` rather than
reimplemented, and the result is asserted to match the published table
`outputs/manuscript_table_ideal_all_ages_balanced.csv` on technique, volume,
reflux fraction and obstruction index for all 75 strata. Any mismatch is a hard
error: this script must add columns, never change numbers.
"""

from __future__ import annotations

import csv
import importlib.util
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = str(PROJECT_ROOT / "outputs" / ".mplconfig")

OUT = PROJECT_ROOT / "outputs"
REVISION_OUT = OUT / "revision"
REVISION_OUT.mkdir(parents=True, exist_ok=True)

PUBLISHED_TABLE = OUT / "manuscript_table_ideal_all_ages_balanced.csv"
TARGET = REVISION_OUT / "balanced_selection_absolute_volumes.csv"


def _load_refresh_module():
    """Import refresh_ideal_tables.py so its selection logic is reused verbatim."""
    path = Path(__file__).resolve().parent / "refresh_ideal_tables.py"
    spec = importlib.util.spec_from_file_location("refresh_ideal_tables", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["refresh_ideal_tables"] = module
    spec.loader.exec_module(module)
    return module


R = _load_refresh_module()


def _read_csv(path: Path) -> list:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, fieldnames: list, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_rows() -> list:
    cache: dict = {}
    rows: list = []
    for age_group in R.AGE_GROUPS:
        for initial_grade in R.GRADES:
            for sev_label, bbd_profile, bbd_severity in R.BBD_SCENARIOS:
                candidates = []
                for technique in R.TECHNIQUES:
                    for volume in R.VOLUME_GRID:
                        cand = R._aggregate_candidate(
                            age_group=age_group,
                            initial_grade=initial_grade,
                            bbd_profile=bbd_profile,
                            bbd_severity=bbd_severity,
                            technique=technique,
                            deflux_volume_ml=volume,
                            cache=cache,
                        )
                        cand["bbd_severity_label"] = sev_label
                        candidates.append(cand)
                chosen = R._best_balanced(candidates, R.BALANCED_OBSTRUCTION_CEILING)
                rows.append(
                    {
                        "age_group": age_group,
                        "preop_vur_grade": initial_grade,
                        "bbd_severity": sev_label,
                        "recommended_technique": chosen["technique"],
                        "deflux_volume_ml": f"{float(chosen['deflux_volume_ml']):.1f}",
                        "predicted_reflux_fraction": f"{float(chosen['pred_postop_reflux_fraction_mean']):.3f}",
                        "predicted_reflux_volume_ml": f"{float(chosen['pred_postop_reflux_volume_ml_mean']):.3f}",
                        "predicted_antegrade_volume_ml": f"{float(chosen['pred_postop_antegrade_volume_ml_mean']):.3f}",
                        "predicted_obstruction_index": f"{float(chosen['pred_postop_obstruction_index_mean']):.3f}",
                        "predicted_postop_grade": int(round(float(chosen["pred_postop_grade_mean"]))),
                        "selection_mode": str(chosen["balanced_mode"]),
                    }
                )
    return rows


def verify_against_published(rows: list) -> None:
    published = _read_csv(PUBLISHED_TABLE)
    key = lambda r: (r["age_group"], str(r["preop_vur_grade"]), r["bbd_severity"])
    pub = {key(r): r for r in published}

    if len(rows) != len(pub):
        raise SystemExit(
            f"Row-count mismatch: regenerated {len(rows)} vs published {len(pub)}."
        )

    problems = []
    for row in rows:
        ref = pub.get(key(row))
        if ref is None:
            problems.append(f"{key(row)}: stratum absent from published table")
            continue
        for field in (
            "recommended_technique",
            "deflux_volume_ml",
            "predicted_reflux_fraction",
            "predicted_obstruction_index",
        ):
            got, want = str(row[field]), str(ref[field])
            if field == "recommended_technique":
                same = got == want
            else:
                same = abs(float(got) - float(want)) < 5e-4
            if not same:
                problems.append(f"{key(row)} {field}: regenerated {got} vs published {want}")

    if problems:
        for p in problems[:20]:
            print("MISMATCH:", p, file=sys.stderr)
        raise SystemExit(
            f"{len(problems)} mismatch(es) against {PUBLISHED_TABLE.name}. "
            "Refusing to write: this script must add columns, not change numbers."
        )
    print(f"Verified {len(rows)}/{len(pub)} strata match {PUBLISHED_TABLE.name} exactly.")


def main() -> None:
    rows = build_rows()
    verify_against_published(rows)
    _write_csv(
        TARGET,
        [
            "age_group",
            "preop_vur_grade",
            "bbd_severity",
            "recommended_technique",
            "deflux_volume_ml",
            "predicted_reflux_fraction",
            "predicted_reflux_volume_ml",
            "predicted_antegrade_volume_ml",
            "predicted_obstruction_index",
            "predicted_postop_grade",
            "selection_mode",
        ],
        rows,
    )
    print(f"Wrote {TARGET} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
