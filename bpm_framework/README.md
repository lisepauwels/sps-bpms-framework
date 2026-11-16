# BPM framework

Retrieve, store and plot the SPS closed orbit (ALPS BPMs) for any MD or time window, cycle by
cycle. Built to mirror `../../sps-measurements/lossmap_framework` (BLMs), so the two can be used
side by side: same table format, same case names, same time axis, same s convention.

Instrument background (what ALPS measures, timing, pitfalls): `../sps_bpm_orbit_notes_for_lise.md`.
MD-specific write-up for the 2026-04-01 crab-cavity MD: `../md_20260401_crab_cavities.md`.

---

## 1. Layout

```
data/bpms/
  bpm_framework/                  code -- run every script from here
    bpm_core/                     the package: pure, MD-agnostic logic, no dependency on config/
      naming.py                   canonical case names + validation (copied from lossmap_core)
      io_data.py                  NXCALS queries, formatting, saving, loading, paths
      data_treatment.py           quality checks, sigma clipping, mean orbit
      lattice.py                  s positions from sps.seq, aperture from an xsuite line (apertls)
      plotting.py                 all figures (pure functions, no table looping)
    config/md_<date>.py           per-MD case schema (see 4)
    scripts/                      runnable drivers (see 5)
    crate_bpm_map.json            BPMALPS_N -> 40 channel names in matrix order (incl. "SPARE")
    bpm_positions.json            channel name -> s [m]
  data/<md_date>/
    bad_bpms.json                 channel -> reason ("faulty", "offset", ...)
    <group>/timestamp_table.csv   one table per group of cycles (see 3)
    <group>/<name>/bct_<ts>.parquet, bpms_H_<ts>.parquet, bpms_V_<ts>.parquet
    <group>/excluded_bpms_<variant>.json   written by plot_overview.py, reusable on other groups
  figures/<md_date>/<group>/...   all plots (see 6)
```

`bpm_core` is usable on its own, e.g. in a notebook:

```python
import sys; sys.path.insert(0, "/Users/lisepauwels/phd/data/bpms/bpm_framework")
from bpm_core import io_data, data_treatment, plotting
h = io_data.load_bpms("20260401", "references", "-90_no_chirp", "20260401_141445", "H")
bct = io_data.load_bct("20260401", "references", "-90_no_chirp", "20260401_141445")
```

---

## 2. The instrument in five lines

- **ALPS**: 6 VME crates `BPMALPS_1..6`, one per sextant, each publishing a 40-row orbit matrix
  per cycle (`Orbit` property: `positions` [mm], `measStamp` [ms]). 227 named channels around the
  ring (113 H / 114 V), the rest are `SPARE`.
- Arc monitors are single-plane (BPH = H, BPV = V); the LSS couplers `BPCE` publish both planes at
  the same s. The plane is the suffix of the channel name: `BPH.60608.H`, `BPCE.61705.V`.
- One orbit sample every 10 ms (each a 1 ms average), `measStamp` = 1, 11, 21, ... ms counted
  **from first injection**. Sample counts differ between crates for the same cycle, so each crate
  is indexed with its own `measStamp`.
- The orbit stops when the beam does: samples after beam loss / dump are NaN.
- Readings below ~10 % of nominal intensity are not trustworthy (log-detector resolution).

---

## 3. Timestamp tables

One CSV per group of cycles, `data/<md_date>/<group>/timestamp_table.csv`, same columns as the
loss-map tables:

| column | meaning |
|---|---|
| `name` | canonical case name, `naming.canonical_name(meta, keys=NAME_KEYS)` |
| `timestamp` | exact UTC cycle stamp -- the NXCALS query key |
| `selector` | e.g. `SPS.USER.LHCMD1` (interleaved users exist, always filtered on) |
| `cycle_length_s` | sets the query window |
| *case parameters* | whatever the MD config defines (`crab_phase`, `lm_type`, `plane`, ...) |
| `beam_present` | filled by `fetch_data.py`: BCT peak above 1e10 charges |
| `crates_saved` | filled by `fetch_data.py`: all 6 crates returned a finite orbit |

A **group** is any set of cycles you want to treat together (e.g. reference cycles vs chirp
cycles). Each group gets its own directory, table, data and figures; the case schema can differ
per group (name keys, validity rules, excluded cycles).

---

## 4. Per-MD config: `config/md_<date>.py`

```python
SELECTOR, ENERGY, ..., CYCLE_LENGTH_S       # constants of the MD
PARAMETERS = [...]                           # every column of a table row (except the 2 flags)
GROUPS = {
    "<group>": {
        "name_keys":      [...],             # PARAMETERS that make up the case name
        "validity_rules": {...},             # cross-field checks (naming.validate)
        "exclude_cycles": [...],             # rows dropped when building the table
        ...                                  # anything the table builder needs
    },
}
CHANNEL_MAP_WINDOW = (t0, t1)                # few-minute window to read the channel names
OVERVIEW_EXCLUDE = [(name, utc_ts), ...]     # kept in table/per-cycle plots, not in overviews
ORBIT_REFERENCE_PHASE = ...                  # case subtracted in the delta-orbit overview
```

To add an MD: write `config/md_<date>.py`, add a `_build_<date>` function to
`scripts/build_timestamp_table.py` (reads whatever source lists that MD's cycles) and register it
in `BUILDERS`.

---

## 5. Pipeline

Run from `bpm_framework/`. Scripts marked **NXCALS** need credentials in the shell:

```
export PYSPARKX_RBAC_USER=lpauwels
export PYSPARKX_RBAC_PASSWORD='your-password'
```

| step | script | NXCALS | what it does |
|---|---|---|---|
| 1 | `build_bpm_map.py --md-date D` | yes | reads `channelNames` of the 6 crates over `CHANNEL_MAP_WINDOW`, matches every name to `sps.seq`, writes `crate_bpm_map.json` + `bpm_positions.json`. Fails on any name not in the lattice (add it to `LATTICE_ALIASES`). Once per channel-map change. |
| 2 | `build_timestamp_table.py --md-date D --group G` | no | builds `data/D/G/timestamp_table.csv` from the config |
| 3 | `import_lossmap_bct.py --table T` | no | optional: copies BCT files already fetched by lossmap_framework for the same cycles, so they are not re-queried |
| 4 | `fetch_data.py --table T [--name N] [--timestamp TS] [--force]` | yes | per row: BCT (if not on disk), orbit of the 6 crates, injection anchor -> saves parquet, fills `beam_present`/`crates_saved`. Resumable: done rows are skipped. |
| 5 | `plot_data.py --table T ...` | no | per-cycle orbit plots + per-case stats plot (see 6) |
| 6 | `plot_overview.py --table T ...` | no | cross-case overview: mean orbit per case, delta to reference (see 6) |

A name starting with `-` must be passed as `--name=-90_no_chirp`.

### What one fetch queries (per cycle)

| data | NXCALS device / property | fields |
|---|---|---|
| orbit, per crate | `BPMALPS_N` / `Orbit` | `cyclestamp, positions, measStamp` |
| injection anchor | `BPMALPS_6` / `LoggingAcquisition` | `acqStampFirstInj` (selector and `> 0` filtered: other users carry a `-1` placeholder) |
| BCT | `SPS.BCTDC24.51454` / `Acquisition` | `totalIntensity, measStamp` (lossmap_framework's query, verbatim) |

Every query is filtered on the selector and matched to the exact cycle stamp.

### Saved files (per cycle)

`bpms_H_<ts>.parquet` and `bpms_V_<ts>.parquet`, long format (one row per BPM per sample):

| column | meaning |
|---|---|
| `channel_name` | `BPH.60608.H` |
| `crate` | `BPMALPS_6` |
| `sample_time_ms` | ms **from the cycle stamp** -- same axis as the BLM files' `sample_time_ms` |
| `meas_ms` | raw ALPS `measStamp`, ms from first injection |
| `position_mm` | orbit [mm], NaN once the beam is gone |
| `s_m` | `bpm_positions.json` |
| `inj_offset_ms` | `acqStampFirstInj - cyclestamp` for this cycle (~1018 ms), so `sample_time_ms = meas_ms + inj_offset_ms` |
| case metadata | the MD's PARAMETERS (e.g. `timestamp, crab_phase, lm_type, plane, coll_gap, energy`) -- `plane` is the case's plane, not the BPM's |

`bct_<ts>.parquet`: `sample_time_ms, intensity` (x1e10 charges) + metadata, identical to
lossmap_framework's. Its time uses lossmap_framework's fixed 1015 ms offset, so it is ~3 ms off
the BPM axis (below the 10 ms orbit sampling).

Files hold **raw data only**: no masking, no clipping, no gating. `<ts>` in file names is local
(Europe/Zurich) time `YYYYMMDD_HHMMSS`; the table's `timestamp` is UTC. SPARE slots are dropped.
About 3.5 MB per cycle.

### Conventions

- **s**: lattice `at` from `acc-models-sps/sps.seq` (`refer = centre`), same as
  `blm_positions.json`. ALPS names that differ from the lattice go through `LATTICE_ALIASES` in
  `build_bpm_map.py` (so far: `BPH.41608` -> `BPH.41607`).
- **Aperture**: `apertls` on the xsuite line
  `code/sps-crab-cavities/studies/ccs_scans2/data/sps_crab_scan_matched.json`. That line's s is
  **not** the sps.seq s (it is 3.14 m shorter, up to 3.3 m off), so `lattice.compute_aperture`
  maps it onto sps.seq s by interpolating through every element present in both.
- **Time**: everything plotted against `sample_time_ms` (from cycle stamp), in s.

---

## 6. Treatment and plots

### BPM selection (applied at plot time, never saved)

1. `data/<md_date>/bad_bpms.json` -- known broken channels, dropped unless `--keep-bad`.
2. `--sigma-clip N` -- per cycle and plane, each BPM's median orbit over the cycle; iterative
   N-sigma clipping across the ring (`data_treatment.clip_mask`, the CCC orbit-correction
   behaviour, as in `sps-crab-cavities/studies/twiss_analysis.ipynb`). One in/out decision per BPM
   per cycle. Clipped BPMs are printed.
3. `--exclude-bpms FILE.json` -- drop a fixed list, typically the `excluded_bpms_<variant>.json`
   written by `plot_overview.py` for another group, so both groups use exactly the same BPMs.

### `plot_data.py` -- per cycle and per case

```
python scripts/plot_data.py --table T [--name N] [--timestamp TS]
       [--sigma-clip N | --exclude-bpms FILE] [--keep-bad]
       [--aperture [--aperture-line LINE.json]]
       [--s-range S0 S1] [--time-range T0_MS T1_MS] [--step-ms 100]
       [--cmap plasma] [--fmt png|pdf]
```

Only rows with `beam_present` and `crates_saved` True. Per cycle:

- `orbit_<name>_<ts>`: H (top) / V (bottom) orbit along s, one line per sample (every `--step-ms`,
  default 100 ms), colour = time in cycle. `--aperture` overlays the aperture, `--s-range` zooms.

Per case (all cycles of a name):

- `orbit_stats_<name>`: per sample, mean orbit across BPMs (point) and RMS = std across BPMs
  (error bar) vs time in cycle; one line per cycle, legend = local start time. A sample needs at
  least half the BPMs finite to get a point.

### `plot_overview.py` -- across cases

```
python scripts/plot_overview.py --table T [--sigma-clip N] [--keep-bad]
       [--s-range S0 S1] [--time-range T0_MS T1_MS] [--ref-phase P] [--cmap coolwarm]
```

- Drops cycles in `config.OVERVIEW_EXCLUDE`.
- Excluded BPMs = bad + the **union** of BPMs clipped in any cycle, so every case is averaged over
  the same BPMs. Written to `data/<md_date>/<group>/excluded_bpms_<variant>.json`.
- Per case and BPM (`data_treatment.mean_orbit`): **mean orbit** = mean over every sample of every
  cycle; **RMS** = std of those samples. **Delta** = mean orbit minus the mean orbit of the
  reference case (`ORBIT_REFERENCE_PHASE` / `--ref-phase`).

| figure | content |
|---|---|
| `mean_orbit_by_phase` | per case: mean over BPMs of the mean orbit, RMS over BPMs as error bar |
| `mean_orbit_vs_s` | mean orbit along s, one line per case, colour = phase |
| `rms_vs_s` | per-BPM RMS along s, colour = phase |
| `delta_orbit_vs_s` | delta orbit along s, colour = phase |
| `delta_orbit_by_phase` | per case: mean over BPMs of the delta, RMS over BPMs as error bar |
| `delta_rms_by_phase` | per case: sqrt(mean(delta^2)) over BPMs |
| `mean_orbit_summary.csv` | all numbers behind them, one row per (case, plane, BPM) |

The overview is written for cases that differ by `crab_phase`; another MD would swap the colour /
x-axis column.

### Figure directories

`figures/<md_date>/<group>/<name>/<variant>/` and `figures/<md_date>/<group>/overview/<variant>/`.
`<variant>` encodes the options: `default`, or tags joined by `_` in a fixed order -- time range
(`t1000-6000ms`), s range (`s5790-6920m`), `aperture`, `withbad`, `clip3sigma`, `step10ms`,
`excl-<group>-<variant>`. Same options always land in the same directory; no titles on figures.

---

## 7. Not done yet

- Gating orbit samples on the BCT (drop samples once intensity falls below a fraction of the
  peak). Needed for anything quantitative on loss cycles.
- Data-driven screening of crates 1-5 for `bad_bpms.json` (only crate-6 channels from the notes so far).
- A generic "time window" table builder (start/end/selector) for studies without a cycle index.
- Turn-by-turn (`CaptureAcquisition`) -- deliberately deferred, heavy queries.
