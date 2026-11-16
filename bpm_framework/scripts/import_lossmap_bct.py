"""
Copy BCT parquet already fetched by lossmap_framework into this framework's data/ -- so cycles
shared with a loss-map table (e.g. 20260401 chirps) don't get re-queried. Files are copied, not
symlinked (this framework stays standalone), and keep their exact content/format (io_data's BCT
query/format is a copy of lossmap_framework's).

For each row of --table: looks for
    <lossmap-root>/data/<md_date>/<name>/bct_<local_ts>.parquet
and copies it to
    data/<md_date>/<group>/<name>/bct_<local_ts>.parquet
after checking the file's own embedded `timestamp` matches the row. Rows with no lossmap file
(e.g. every 20260401 reference) are just reported -- fetch_data.py queries those fresh.
Existing destination files are left alone unless --force.

Usage (from bpm_framework/):
    python scripts/import_lossmap_bct.py --table ../data/20260401/chirps/timestamp_table.csv
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bpm_framework/

from bpm_core import io_data  # noqa: E402

LOSSMAP_ROOT = io_data.ROOT.parents[1] / "sps-measurements" / "lossmap_framework"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=Path, required=True, help="Path to a timestamp_table.csv.")
    parser.add_argument("--lossmap-root", type=Path, default=LOSSMAP_ROOT)
    parser.add_argument("--force", action="store_true", help="Overwrite existing destination files.")
    args = parser.parse_args()

    table = args.table.resolve()
    group, md_date = table.parent.name, table.parent.parent.name
    df = pd.read_csv(table)

    copied = skipped = 0
    not_found = []
    for _, row in df.iterrows():
        name, timestamp = row["name"], row["timestamp"]
        local_ts = io_data.local_timestamp_from_utc(timestamp)
        src = args.lossmap_root / "data" / md_date / name / f"bct_{local_ts}.parquet"
        dst = io_data.bct_path(md_date, group, name, local_ts)
        if not src.exists():
            not_found.append(f"{name} / {timestamp}")
            continue
        if dst.exists() and not args.force:
            skipped += 1
            continue
        embedded = pd.read_parquet(src, columns=["timestamp"])["timestamp"].iloc[0]
        if pd.Timestamp(embedded) != pd.Timestamp(timestamp):
            raise ValueError(f"{src}: embedded timestamp {embedded} != table timestamp {timestamp}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        copied += 1

    print(f"{copied} copied, {skipped} already present, {len(not_found)} not in lossmap_framework.")
    for item in not_found:
        print(f"  not found: {item}")


if __name__ == "__main__":
    main()
