"""
Fetch fresh NXCALS data for rows of a timestamp table and fill in their quality flags -- same
shape as lossmap_framework's fetch_data.py. md_date and group are inferred from --table's path
(data/<md_date>/<group>/timestamp_table.csv).

For each row:
  - BCT: reuse data/<md_date>/<group>/<name>/bct_<ts>.parquet if present (e.g. imported by
    import_lossmap_bct.py), else query + save it.
  - Orbit: query all 6 ALPS crates (io_data.query_orbit) + the per-cycle injection anchor
    (io_data.query_injection_offset_ms), format into bpms_H/bpms_V frames with timestamp + the
    MD's metadata fields (config.PARAMETERS minus selector/cycle_length_s), save both.
  - Assess quality on the in-memory frames (data_treatment.assess_quality) and write
    beam_present/crates_saved back into the table -- one row at a time, so an interrupted run
    doesn't lose earlier progress.

Skips a row entirely if beam_present is already filled; skips the orbit re-query (but still runs
the quality check, on the saved files) if both bpms parquet already exist. --force redoes both.

Needs pysparkx/NXCALS credentials already exported in your shell:
    export PYSPARKX_RBAC_USER=lpauwels
    export PYSPARKX_RBAC_PASSWORD='your-password'
Run this yourself, in that shell -- not delegated to an agent.

Usage (from bpm_framework/):
    python scripts/fetch_data.py --table ../data/20260401/references/timestamp_table.csv
    python scripts/fetch_data.py --table ../data/20260401/chirps/timestamp_table.csv --name 0_betatron_H
    python scripts/fetch_data.py --table ../data/20260401/references/timestamp_table.csv --name=-120_no_chirp --timestamp "2026-04-01 14:00:54.135000+00:00"
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bpm_framework/

from bpm_core import data_treatment, io_data  # noqa: E402

NOT_METADATA = {"selector", "cycle_length_s"}  # describe how to fetch a cycle, not what it is

BPM_POSITIONS_JSON = io_data.ROOT / "bpm_positions.json"
CRATE_BPM_MAP_JSON = io_data.ROOT / "crate_bpm_map.json"


def fetch_one(
    md_date: str, group: str, row: pd.Series, bpm_positions: dict, crate_bpm_map: dict,
    force: bool = False,
) -> dict:
    """Query, format, assess and save one row (BCT only if not already on disk) -- or, if the
    orbit is already saved and not forced, just reload+assess. Returns the quality flags."""
    name, timestamp, selector = row["name"], row["timestamp"], row["selector"]
    cycle_length_s = float(row["cycle_length_s"])
    local_ts = io_data.local_timestamp_from_utc(timestamp)

    md_config = importlib.import_module(f"config.md_{md_date}")
    meta = {k: row[k] for k in md_config.PARAMETERS if k not in NOT_METADATA}

    bct_path = io_data.bct_path(md_date, group, name, local_ts)
    if force or not bct_path.exists():
        times_ms, intensity = io_data.query_bct(timestamp, selector)
        bct_df = io_data.format_bct(times_ms, intensity, meta)
        io_data.save_bct(bct_df, bct_path)
    else:
        bct_df = io_data.load_bct(md_date, group, name, local_ts)

    paths = {p: io_data.bpms_path(md_date, group, name, local_ts, p) for p in io_data.PLANES}
    if force or not all(path.exists() for path in paths.values()):
        orbit = io_data.query_orbit(timestamp, selector, cycle_length_s, crate_bpm_map)
        inj_offset_ms = io_data.query_injection_offset_ms(timestamp, selector)
        bpm_dfs = io_data.format_bpms(orbit, inj_offset_ms, bpm_positions, meta)
        for plane, path in paths.items():
            io_data.save_bpms(bpm_dfs[plane], path)
    else:
        print(f"[{name} / {timestamp}] orbit already saved on disk, skipping NXCALS re-query.")
        bpm_dfs = {p: io_data.load_bpms(md_date, group, name, local_ts, p) for p in io_data.PLANES}

    quality = data_treatment.assess_quality(bct_df, list(bpm_dfs.values()))
    if quality["bpm_crates_down"]:
        print(f"[{name} / {timestamp}] crates down: {quality['bpm_crates_down']}")
    return {"beam_present": not quality["noisy_fill"], "crates_saved": quality["n_bpm_crates_down"] == 0}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=Path, required=True, help="Path to a timestamp_table.csv.")
    parser.add_argument("--name", help="Only process this canonical case name.")
    parser.add_argument("--timestamp", help="Only process this exact timestamp (requires --name).")
    parser.add_argument("--force", action="store_true",
                        help="Re-fetch and overwrite rows/files even if already saved/processed.")
    args = parser.parse_args()
    if args.timestamp and not args.name:
        parser.error("--timestamp requires --name.")

    table = args.table.resolve()
    group, md_date = table.parent.name, table.parent.parent.name

    with open(BPM_POSITIONS_JSON) as f:
        bpm_positions = json.load(f)
    with open(CRATE_BPM_MAP_JSON) as f:
        crate_bpm_map = json.load(f)

    df = pd.read_csv(table, dtype={"beam_present": "object", "crates_saved": "object"})
    rows = df.index
    if args.name:
        rows = df.index[df["name"] == args.name]
        if args.timestamp:
            rows = df.index[(df["name"] == args.name) & (df["timestamp"] == args.timestamp)]
    if len(rows) == 0:
        parser.error(f"No matching rows for name={args.name!r} timestamp={args.timestamp!r}.")

    for i in rows:
        row = df.loc[i]
        if not args.force and pd.notna(row["beam_present"]):
            print(f"[{row['name']} / {row['timestamp']}] already processed, skipping (use --force to redo).")
            continue
        print(f"[{row['name']} / {row['timestamp']}] fetching...")
        try:
            flags = fetch_one(md_date, group, row, bpm_positions, crate_bpm_map, force=args.force)
        except Exception as exc:
            print(f"[{row['name']} / {row['timestamp']}] ERROR: {exc}")
            continue
        df.at[i, "beam_present"] = flags["beam_present"]
        df.at[i, "crates_saved"] = flags["crates_saved"]
        df.to_csv(table, index=False)
        print(f"[{row['name']} / {row['timestamp']}] beam_present={flags['beam_present']} "
              f"crates_saved={flags['crates_saved']}")


if __name__ == "__main__":
    main()
