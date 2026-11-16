# BPM orbit — 2026-04-01 crab-cavity phase scan

What was retrieved, how it is organised, what was excluded and why, and what the plots show.
Framework documentation (scripts, file formats, options): `bpm_framework/README.md`.
Instrument background: `sps_bpm_orbit_notes_for_lise.md`.

## 1. The MD

- **2026-04-01**, `SPS.USER.LHCMD1`, 26 GeV, Q20, single bunch (~1.8e11 p), TCSM at 5 sigma.
- One crab cavity (`ACFCAH.61739`, s = 6313.6 m, LSS6) at ~1 MV; its phase setpoint `cavPhA` was
  scanned in 30 deg steps from -179 to +179 deg, plus reference cycles (`LossMapRef`).
- At each phase: a few cycles **without excitation** (`no_chirp`) and a few with a **tune chirp**
  in H or V that scrapes ~80 % of the beam within ~300 ms (the loss maps).
- Scan from 12:02 to 15:27 UTC. The cycle index of the whole scan is
  `../sps-measurements/lossmaps_crab_cavities/nxcalsdata/nxcals_index.csv` (138 cycles:
  44 H, 42 V, 52 no_chirp) -- the single source both the BLM and the BPM tables are built from.
- Beam circulates ~20 s per cycle (injection ~1.02 s after the cycle stamp, dump at ~21 s).

## 2. Two groups, two tables

The `no_chirp` and the chirp cycles answer different questions, so they are separate groups
(`config/md_20260401.py`, `GROUPS`), each with its own table, data and figures:

| | `references` | `chirps` |
|---|---|---|
| table | `data/20260401/references/timestamp_table.csv` | `data/20260401/chirps/timestamp_table.csv` |
| source rows | `plane == "no_chirp"` | `plane in (H, V)` |
| rows / cases | **50 / 13** | **82 / 28** |
| case name | `<phase>_no_chirp`, e.g. `-90_no_chirp` | `<phase>_betatron_<plane>`, e.g. `-90_betatron_H`, `ref_betatron_V` |
| what the orbit is | the stable closed orbit at that crab phase, for the whole cycle | the orbit before and during the chirp, until the beam is lost |
| use | crab-induced orbit vs phase, the reference orbit, BPM health | orbit next to the loss maps |

**references** — 4 cycles per phase, except -120 and -60 (3). Two `no_chirp` rows of
`nxcals_index.csv` were dropped because they are the **same physical cycle** as a V-chirp row
(identical cycle stamp); lossmap_framework treats both as V chirps:

- `Phase_-120_no_chirp id1` = `Phase_-150_V id3` (14:00:21.735 UTC)
- `Phase_-60_no_chirp id1` = `Phase_-120_V id3` (14:05:13.335 UTC)

The `ref` cycles (`LossMapRef_H`, `LossMapRefV_V`) are **not** in the references: per
`lossmaps_crab_cavities/saveddata/DATA_NOTES.md` they are chirped loss maps (the H / V chirp
references), so they sit in `chirps` as `ref_betatron_H` / `ref_betatron_V`.

**chirps** — exactly lossmap_framework's 20260401 table (same 82 rows, same names, same
exclusions, checked column by column): 3 cycles per case, except `-90_betatron_V` (1). Excluded,
as on the BLM side: `LossMapRef_H` id1 (discarded) and id4 (bad data), `LossMapRefV_V` id1
(discarded), `Phase_-90_V` id1 (same cycle as `Phase_-90_H` id3).

## 3. Retrieval

All done, from `bpm_framework/`:

```
python scripts/build_bpm_map.py --md-date 20260401                       # NXCALS, once
python scripts/build_timestamp_table.py --md-date 20260401 --group references
python scripts/build_timestamp_table.py --md-date 20260401 --group chirps
python scripts/import_lossmap_bct.py --table ../data/20260401/chirps/timestamp_table.csv
python scripts/fetch_data.py --table ../data/20260401/references/timestamp_table.csv   # NXCALS
python scripts/fetch_data.py --table ../data/20260401/chirps/timestamp_table.csv       # NXCALS
```

- **Channel map** read over 14:00-14:06 UTC (inside the scan): 227 named channels, 113 H / 114 V,
  per-crate counts identical to the BPM notes. One name differs from the lattice:
  ALPS `BPH.41608` = lattice `BPH.41607` (alias in `build_bpm_map.py`).
- **BCT**: chirps copied from lossmap_framework (82/82, not re-queried); references fetched here.
- **Status**: all 132 cycles fetched, `beam_present` and `crates_saved` True on every row.
  ~170 MB (references) + ~200 MB (chirps).
- **Time**: measured injection offset 1018.05 ms on every cycle checked (BPM notes: 1018.05).
  Files carry `sample_time_ms` from the cycle stamp, directly comparable with the BLM
  `sample_time_ms`.

Per cycle: `data/20260401/<group>/<name>/bct_<ts>.parquet`, `bpms_H_<ts>.parquet`,
`bpms_V_<ts>.parquet` (`<ts>` = local time, e.g. `20260401_141445` = 12:14:45 UTC).

## 4. BPM selection

| list | channels | why |
|---|---|---|
| `data/20260401/bad_bpms.json` | `BPH.62008.H` (faulty), `BPH.63008.H` (faulty), `BPH.63608.H` (-53 mm offset) | crate 6, from the BPM notes |
| 3-sigma clip on the references | `BPCE.21706.H` (+10 mm), `BPCE.21803.H` (+16 mm), `BPH.31808.H`, `BPH.41208.H`, `BPCN.22708.V` (-16 mm) | outliers of the cycle-median orbit in at least one reference cycle |

The two BPCE in LSS2 are most likely the extraction bump (real orbit, 16-19 mm per the notes),
not faults; the CCC clip removes them too. The union of both lists is saved as
`data/20260401/references/excluded_bpms_clip3sigma.json` and reused for the chirps, so every plot
uses the same 106 H / 113 V BPMs.

**Open:** `BPV.31908.V` (s = 2910 m) is noisy -- its sample-to-sample RMS is 0.67 mm against
~0.1 mm for the other V BPMs, and it swings by up to 1 mm between phases. It was only clipped in
the excluded 0-deg cycle (below), so it is currently kept. Candidate for `bad_bpms.json`.

**Cycle left out of the overview** (`config.OVERVIEW_EXCLUDE`): `0_no_chirp` 14:16:01.335 UTC.
The orbit steps by ~2.2 mm RMS in H and ~0.6 mm in V between 3 and 5 s into the cycle and stays
there until the dump; the other three 0-deg cycles are normal. V moves too, so it is not the
crab. Cause unknown. The cycle stays in the table and has its own per-cycle plot.

## 5. Plots

### Per cycle / per case (`plot_data.py`)

```
T=../data/20260401/references/timestamp_table.csv
python scripts/plot_data.py --table $T                                  # default/
python scripts/plot_data.py --table $T --sigma-clip 3                   # clip3sigma/
python scripts/plot_data.py --table $T --aperture                       # aperture/
python scripts/plot_data.py --table $T --s-range 5790 6920 --aperture   # s5790-6920m_aperture/ (sextant 6, the crab)

python scripts/plot_data.py --table ../data/20260401/chirps/timestamp_table.csv --aperture \
    --exclude-bpms ../data/20260401/references/excluded_bpms_clip3sigma.json
```

- `orbit_<name>_<ts>`: orbit along s, H top / V bottom, colour = time in cycle.
- `orbit_stats_<name>`: mean across BPMs vs time, RMS across BPMs as error bar, one line per cycle.

On the references (e.g. `-90_no_chirp`): after the first second the orbit is stable for the whole
cycle; the 4 cycles of a phase overlap (mean x ~ -0.5 mm, mean y ~ +0.8 mm, RMS across BPMs
~2.6 mm H / ~1.4 mm V, dominated by the static orbit).

On the chirps the orbit stops when the beam is lost (e.g. ~5.2 s into the cycle on
`-90_betatron_H`). The last samples before that are read by only a few BPMs at low intensity and
are not meaningful (see 7).

### Overview across phases (`plot_overview.py`)

```
python scripts/plot_overview.py --table ../data/20260401/references/timestamp_table.csv --sigma-clip 3
```

Output: `figures/20260401/references/overview/clip3sigma/`. 49 reference cycles (one excluded).
Per phase and BPM: **mean orbit** = mean over all samples of all cycles of that phase; **RMS** =
std of those samples. **Delta** = mean orbit minus the **0 deg** mean orbit
(`ORBIT_REFERENCE_PHASE = 0`; change with `--ref-phase`).

| figure | what it shows |
|---|---|
| `mean_orbit_vs_s` | H orbit pattern changes with phase (blue vs red); V identical for all phases |
| `mean_orbit_by_phase` | flat: mean over BPMs -0.45 to -0.55 mm (H), 0.81-0.83 mm (V); the error bar is the static orbit (2.2-2.8 mm H, 1.4 mm V) |
| `rms_vs_s` | per-BPM stability: ~0.2-0.5 mm H, ~0.1 mm V |
| `delta_orbit_vs_s` | the crab orbit: a clean betatron-like wave around the whole ring in H, opposite sign for negative / positive phases; V flat except `BPV.31908.V` |
| `delta_orbit_by_phase` | mean of the delta ~0 (the wave averages out), RMS over BPMs as error bar |
| `delta_rms_by_phase` | **size of the crab orbit change vs phase** (below) |
| `mean_orbit_summary.csv` | every number, per phase / plane / BPM |

RMS over BPMs of the orbit change relative to 0 deg [mm]:

| phase [deg] | -179 | -150 | -120 | -90 | -60 | -30 | 0 | 30 | 60 | 90 | 120 | 150 | 179 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H | 1.19 | 1.60 | 1.74 | **2.38** | 1.40 | 0.73 | 0 | 0.19 | 0.51 | 0.48 | 0.18 | 0.44 | 0.99 |
| V | 0.09 | 0.08 | 0.13 | 0.20 | 0.11 | 0.07 | 0 | 0.10 | 0.04 | 0.06 | 0.08 | 0.06 | 0.06 |

H: largest at -90 deg, small from 0 to +150 deg. V stays at the 0.05-0.2 mm level (the crab kicks
horizontally only -- V is the null test); its largest values follow H, which could be coupling or
the noisy `BPV.31908.V`.

## 6. Points to keep in mind

- **Phase convention**: the `crab_phase` here is the `cavPhA` setpoint (the BPM notes' "your knob").
  The notes' afternoon phase jumps used a different knob (`antPh1`) on other cycles; the notes'
  reference -35 deg is on `cavPhA`.
- **What ALPS reports for a crabbed bunch** is not the bunch centroid but a 200 MHz-weighted
  centroid, ~2.5-3x larger than the centroid shift (notes, section 2). Orbit numbers here are
  ALPS readings, not beam offsets at the collimator.
- **Time axes**: BPM `sample_time_ms` uses the measured 1018 ms injection offset, the BCT the
  fixed 1015 ms of lossmap_framework -- ~3 ms apart, below the 10 ms orbit sampling.
- **Aperture**: from the xsuite line, mapped onto the sps.seq s axis (up to 3.3 m correction).

## 7. Next

- Gate orbit samples on the BCT (drop samples once intensity is below a fraction of the peak),
  needed before reading anything from the chirp orbits during the loss.
- Decide on `BPV.31908.V` (bad list or not).
- Chirp overview: orbit change at chirp onset per phase, next to the loss maps.
