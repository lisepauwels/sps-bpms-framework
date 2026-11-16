"""
Build data/<md_date>/<group>/timestamp_table.csv for one group of one MD's config
(config/md_<date>.py, GROUPS). Same row schema as lossmap_framework: fields match
config.PARAMETERS exactly, `name` from naming.canonical_name(meta, keys=group["name_keys"]), every
row checked against group["validity_rules"], beam_present/crates_saved left empty for
fetch_data.py to fill in.

Run (from bpm_framework/):
    python scripts/build_timestamp_table.py --md-date 20260401 --group chirps
    python scripts/build_timestamp_table.py --md-date 20260401 --group references
"""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bpm_framework/

from bpm_core import io_data, naming  # noqa: E402


def _drop_excluded(index: pd.DataFrame, exclude: list[tuple[str, str]]) -> pd.DataFrame:
    if not exclude:
        return index
    mask = index[["study", "id"]].apply(tuple, axis=1).isin(set(exclude))
    if mask.sum() != len(exclude):
        found = set(index[mask][["study", "id"]].apply(tuple, axis=1))
        raise ValueError(f"exclude_cycles entries not found in nxcals_index.csv: {set(exclude) - found}")
    return index[~mask]


def _build_20260401(md_config, group: dict) -> pd.DataFrame:
    """From lossmaps_crab_cavities/nxcalsdata/nxcals_index.csv, rows whose plane is in
    group["planes"] (H/V for chirps, no_chirp for references)."""
    index = pd.read_csv(md_config.NXCALS_INDEX_CSV)

    assert (index["TCSM_gap"].astype(str) + "sigma" == md_config.TCSM_GAP).all(), (
        "TCSM_gap is no longer constant across nxcals_index.csv -- revisit config.TCSM_GAP."
    )
    assert (index["energy"] == md_config.ENERGY).all(), (
        "energy is no longer constant across nxcals_index.csv -- revisit config.ENERGY."
    )
    assert (index["md_user"] == md_config.SELECTOR).all(), (
        "md_user (selector) is no longer constant across nxcals_index.csv -- revisit config.SELECTOR."
    )

    index = index[index["plane"].isin(group["planes"])]
    index = _drop_excluded(index, group["exclude_cycles"])

    rows = []
    for _, r in index.iterrows():
        meta = {
            "crab_phase": r["phase"],
            "lm_type": str(r["LM_type"]).lower(),
            "plane": r["plane"],
            "coll_gap": md_config.TCSM_GAP,
            "energy": r["energy"],
        }
        row = {
            "name": naming.canonical_name(meta, keys=group["name_keys"]),
            "timestamp": r["cycle_stamp_utc"],
            "selector": md_config.SELECTOR,
            "cycle_length_s": md_config.CYCLE_LENGTH_S,
            **meta,
        }
        naming.validate(row, group["validity_rules"])
        row["beam_present"] = None
        row["crates_saved"] = None
        rows.append(row)

    df = pd.DataFrame(rows)
    dup = df["timestamp"].duplicated(keep=False)
    if dup.any():
        raise ValueError(f"Duplicate cycle stamps in table:\n{df[dup].to_string()}")
    return df


BUILDERS = {
    "20260401": _build_20260401,
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--md-date", required=True, choices=sorted(BUILDERS))
    parser.add_argument("--group", required=True)
    args = parser.parse_args()

    md_config = importlib.import_module(f"config.md_{args.md_date}")
    if args.group not in md_config.GROUPS:
        parser.error(f"--group must be one of {sorted(md_config.GROUPS)}")
    df = BUILDERS[args.md_date](md_config, md_config.GROUPS[args.group])

    out_csv = io_data.DATA_ROOT / args.md_date / args.group / "timestamp_table.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(df.to_string(index=False))
    print(f"\n{len(df)} rows, {df['name'].nunique()} names saved to:\n  {out_csv}")


if __name__ == "__main__":
    main()
