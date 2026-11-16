# SPS BPM orbit framework

Retrieve, store and plot the SPS closed orbit from the ALPS BPMs (NXCALS), cycle by cycle, for
any MD — the BPM counterpart of the loss-map framework (`sps-measurements/lossmap_framework`):
same timestamp-table format, same case names, same time axis and s convention as the BLM data.

- **`bpm_framework/`** — the code. How it works: [`bpm_framework/README.md`](bpm_framework/README.md).
- **`data/<md_date>/`** — per-MD timestamp tables, bad / excluded BPM lists. The orbit and BCT
  parquet files are not tracked; regenerate them with `scripts/fetch_data.py` (needs NXCALS access).
- **`md_20260401_crab_cavities.md`** — the 2026-04-01 crab-cavity phase scan: what was retrieved,
  the reference / chirp tables, exclusions, overview results.
- `figures/` — plot output (not tracked).

Quick start (from `bpm_framework/`):

```
python scripts/fetch_data.py --table ../data/20260401/references/timestamp_table.csv
python scripts/plot_data.py --table ../data/20260401/references/timestamp_table.csv --sigma-clip 3
python scripts/plot_overview.py --table ../data/20260401/references/timestamp_table.csv --sigma-clip 3
```

Dependencies: numpy, pandas, pyarrow, matplotlib; `pysparkx` (NXCALS) for fetching;
`xtrack` + `apertls` for the aperture overlay. The aperture and lattice paths
(`bpm_framework/bpm_core/lattice.py`) point to `sps-xsuite-model` and `sps-crab-cavities`
checked out next to this repository's parent directories.
