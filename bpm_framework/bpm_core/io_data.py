"""
Everything data-related for one cycle: query it from NXCALS, format it, save it to parquet, load
it back -- each a separate step (see scripts/fetch_data.py, which composes them: query -> format
-> data_treatment.assess_quality on the in-memory result -> save -> update the timestamp table).
Same shape as lossmap_framework/lossmap_core/io_data.py; query_bct/format_bct/save_bct/load_bct
are copied from there unchanged, so BCT files imported from lossmap_framework
(scripts/import_lossmap_bct.py) and BCT files fetched here are interchangeable.

Layout, under DATA_ROOT (= data/bpms/data, next to bpm_framework/):
    <md_date>/<group>/timestamp_table.csv
    <md_date>/<group>/<name>/bct_<timestamp>.parquet
    <md_date>/<group>/<name>/bpms_H_<timestamp>.parquet
    <md_date>/<group>/<name>/bpms_V_<timestamp>.parquet
`<group>` is one timestamp table's directory (e.g. "references" / "chirps" for 20260401), `<name>`
the canonical name from naming.canonical_name(). `<timestamp>` on disk is local (Europe/Zurich)
time, "YYYYMMDD_HHMMSS" -- see local_timestamp_from_utc (a timestamp table's own "timestamp"
column is the UTC cycle stamp itself, never this derived label).

ALPS orbit: BPMALPS_1..6 / Orbit, one row per cycle, `positions` a (40, n_samples) matrix in mm,
channel-major in the order of crate_bpm_map.json (SPARE slots included, so list index == matrix
row), `measStamp` the in-cycle time of each sample in ms FROM FIRST INJECTION (grid 1 + 10i ms).
Sample counts differ between crates for the same cycle, so each crate is indexed with its own
measStamp, never zipped positionally. Short/MD3 cycles log an empty (40, 0) matrix.

Time axis: saved BPM files carry both `meas_ms` (raw, from first injection) and `sample_time_ms`
(= meas_ms + inj_offset_ms, from the CYCLE STAMP -- the same axis as the BLM files'
sample_time_ms, which carry no injection constant). inj_offset_ms is measured per cycle from
BPMALPS_6/LoggingAcquisition (acqStampFirstInj - cyclestamp, ~1018.05 ms on 20260401), not a
constant. NB the copied BCT query still uses the fixed INJECTION_OFFSET_MS = 1015 of
lossmap_framework, so BCT and BPM times differ by ~3 ms -- below the 10 ms orbit sampling.

Saved BPM parquet carries only raw data + s_m (static geometry). Bad-channel masking, BCT gating
and reference subtraction are analysis-time concerns (data/<md_date>/bad_bpms.json,
data_treatment.py), never baked into saved files. SPARE rows are dropped at format time (always
NaN). No BPM-plane column: the files are split by plane, and the channel name carries it as a
suffix (".H"/".V") -- a "plane" column would also clash with the case metadata "plane" (the chirp
plane).
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent  # bpm_framework/ (this file lives in bpm_core/)
DATA_ROOT = ROOT.parent / "data"
FIGURES_ROOT = ROOT.parent / "figures"

TS_RE = re.compile(r"(\d{8}_\d{6})")

BPM_CRATES: list[str] = [f"BPMALPS_{i}" for i in range(1, 7)]
SPARE = "SPARE"
PLANES = ("H", "V")

ANCHOR_DEVICE = "BPMALPS_6"  # any crate carries LoggingAcquisition; one is enough per cycle
ANCHOR_SANITY_MS = (1000.0, 1050.0)  # acqStampFirstInj - cyclestamp outside this -> warn


def local_timestamp_from_utc(cycle_stamp_utc: str) -> str:
    """Europe/Zurich "YYYYMMDD_HHMMSS" label from a UTC cyclestamp -- for naming saved files only."""
    return pd.Timestamp(cycle_stamp_utc).tz_convert("Europe/Zurich").strftime("%Y%m%d_%H%M%S")


def bpm_plane(channel_name: str) -> str:
    """"H"/"V" from an ALPS channel name's suffix, e.g. "BPCE.61705.V" -> "V"."""
    plane = channel_name.rsplit(".", 1)[-1]
    if plane not in PLANES:
        raise ValueError(f"Cannot read a plane off channel name {channel_name!r}.")
    return plane


INJECTION_OFFSET_MS = 1015.0  # BCT measStamp is 0-based from injection; BLM times are 0-based
                               # from SPS cycle start. Injection happens ~1015 ms into the cycle.
                               # (copied from lossmap_framework -- BCT only, see module docstring)

_LOGGED_IN = False


def _ensure_login() -> None:
    """Import pysparkx and log in, lazily -- only when a function actually queries NXCALS.
    Reads PYSPARKX_RBAC_USER/PYSPARKX_RBAC_PASSWORD from the environment. Keeps every non-query
    helper usable without pysparkx installed."""
    global _LOGGED_IN
    if not _LOGGED_IN:
        import pysparkx
        pysparkx.login()
        _LOGGED_IN = True


def _vec(cell) -> np.ndarray:
    """Unpack a VECTOR cell (NUMERIC or STRING) into a 1-D numpy array."""
    if cell is None:
        return np.array([])
    if hasattr(cell, "asDict"):
        cell = cell.asDict()
    if isinstance(cell, dict):
        return np.array(cell.get("elements") or [])
    return np.array(list(cell))


def _mat(cell, n_rows: int | None = None) -> np.ndarray:
    """Unpack a MATRIX_NUMERIC cell into a 2-D float array, using its own `dimensions`.
    Falls back to reshaping a flat array to n_rows rows if dimensions are missing."""
    if cell is None:
        return np.empty((0, 0))
    if hasattr(cell, "asDict"):
        cell = cell.asDict()
    if isinstance(cell, dict):
        dims = list(cell.get("dimensions") or [])
        flat = np.array(cell.get("elements") or [], dtype=float)
    else:
        dims, flat = [], np.array(cell, dtype=float)
    if flat.ndim == 2:
        return flat
    if len(dims) == 2 and flat.size == dims[0] * dims[1]:
        return flat.reshape(dims[0], dims[1])
    if n_rows and flat.size % n_rows == 0:
        return flat.reshape(n_rows, -1)
    if flat.size == 0:
        return np.empty((n_rows or 0, 0))
    raise ValueError(f"Cannot reshape matrix: dims={dims}, {flat.size} elements.")


def _nearest_row(df: pd.DataFrame, cycle_ns: int) -> pd.Series:
    exact = df[df["cyclestamp"] == cycle_ns]
    if not exact.empty:
        return exact.iloc[0]
    df = df.copy()
    df["_d"] = (df["cyclestamp"] - cycle_ns).abs()
    return df.sort_values("_d").iloc[0]


# ---------------------------------------------------------------------------
# Query (NXCALS)
# ---------------------------------------------------------------------------

def query_bct(cycle_stamp_utc: str, selector: str) -> tuple[np.ndarray, np.ndarray]:
    """Return (times_ms, intensity) from SPS.BCTDC24.51454/Acquisition."""
    _ensure_login()
    from pysparkx.nxcals import DataQuery
    ts = pd.Timestamp(cycle_stamp_utc)
    df = (
        DataQuery(device="SPS.BCTDC24.51454", property="Acquisition",
                  start_time=(ts - pd.Timedelta(seconds=2)).isoformat(),
                  end_time=(ts + pd.Timedelta(seconds=45)).isoformat())
        .build()
        .where(f"selector = '{selector}'")
        .selectExpr("cyclestamp", "totalIntensity", "measStamp")
        .toPandas()
    )
    if df.empty:
        raise RuntimeError(f"No BCT data for {cycle_stamp_utc} / {selector}")
    row = _nearest_row(df, ts.value)
    times = _vec(row["measStamp"]).astype(float) + INJECTION_OFFSET_MS
    return times, _vec(row["totalIntensity"]).astype(float)


def query_injection_offset_ms(cycle_stamp_utc: str, selector: str) -> float:
    """
    acqStampFirstInj - cyclestamp [ms] for this cycle, from BPMALPS_6/LoggingAcquisition -- the
    offset between ALPS's measStamp axis (from first injection) and the cycle stamp.
    Other users interleaved in the same supercycle (e.g. MD3) carry an INTEGER -1 sentinel in
    acqStampFirstInj (not null -- .notna() doesn't filter it), hence the selector AND positivity
    filters below.
    """
    _ensure_login()
    from pysparkx.nxcals import DataQuery
    ts = pd.Timestamp(cycle_stamp_utc)
    df = (
        DataQuery(device=ANCHOR_DEVICE, property="LoggingAcquisition",
                  start_time=(ts - pd.Timedelta(seconds=2)).isoformat(),
                  end_time=(ts + pd.Timedelta(seconds=45)).isoformat())
        .build()
        .where(f"selector = '{selector}'")
        .selectExpr("cyclestamp", "acqStampFirstInj")
        .toPandas()
    )
    df = df[df["acqStampFirstInj"].astype("int64") > 0]
    if df.empty:
        raise RuntimeError(f"No LoggingAcquisition anchor for {cycle_stamp_utc} / {selector}")
    row = _nearest_row(df, ts.value)
    if int(row["cyclestamp"]) != ts.value:
        raise RuntimeError(f"LoggingAcquisition: no row at exactly cyclestamp {cycle_stamp_utc}")
    offset_ms = (int(row["acqStampFirstInj"]) - ts.value) / 1e6
    lo, hi = ANCHOR_SANITY_MS
    if not lo <= offset_ms <= hi:
        warnings.warn(f"{cycle_stamp_utc}: injection offset {offset_ms:.2f} ms outside [{lo}, {hi}].")
    return offset_ms


def query_orbit_crate(crate: str, cycle_stamp_utc: str, selector: str,
                      cycle_length_s: float,
                      crate_bpm_map: dict[str, list[str]]) -> pd.DataFrame:
    """
    Query positions + measStamp from one ALPS crate's Orbit property.
    Channel names come from crate_bpm_map (matrix row order, SPARE included).

    Returns long-format DataFrame: channel_name, crate, meas_ms, position_mm -- SPARE rows
    dropped -- or an empty DataFrame on failure / empty cycle.
    """
    _ensure_login()
    from pysparkx.nxcals import DataQuery
    ts = pd.Timestamp(cycle_stamp_utc)
    t0 = (ts - pd.Timedelta(seconds=2)).isoformat()
    t1 = (ts + pd.Timedelta(seconds=cycle_length_s + 5)).isoformat()

    ch_names = list(crate_bpm_map.get(crate, []))
    if not ch_names:
        warnings.warn(f"{crate}: not found in crate_bpm_map -- skipping.")
        return pd.DataFrame()

    try:
        df = (
            DataQuery(device=crate, property="Orbit", start_time=t0, end_time=t1)
            .build()
            .where(f"selector = '{selector}'")
            .selectExpr("cyclestamp", "positions", "measStamp")
            .toPandas()
        )
    except Exception as exc:
        warnings.warn(f"{crate}: query failed -- {exc}")
        return pd.DataFrame()

    if df.empty:
        warnings.warn(f"{crate}: no data in NXCALS for this cycle.")
        return pd.DataFrame()

    row = _nearest_row(df, ts.value)
    if int(row["cyclestamp"]) != ts.value:
        warnings.warn(f"{crate}: no row at exactly cyclestamp {cycle_stamp_utc} -- skipping.")
        return pd.DataFrame()

    times = _vec(row["measStamp"]).astype(float)
    raw = _mat(row["positions"], n_rows=len(ch_names))
    if raw.size == 0 or times.size == 0:
        warnings.warn(f"{crate}: empty orbit matrix for this cycle -- skipping.")
        return pd.DataFrame()

    n_ch, n_samp = raw.shape
    if n_ch != len(ch_names):
        warnings.warn(f"{crate}: matrix has {n_ch} rows but map has {len(ch_names)} names -- skipping.")
        return pd.DataFrame()
    if n_samp != len(times):
        warnings.warn(f"{crate}: matrix has {n_samp} samples but measStamp has "
                      f"{len(times)} entries -- truncating to min.")
        n_t = min(n_samp, len(times))
        raw, times = raw[:, :n_t], times[:n_t]

    keep = [i for i, name in enumerate(ch_names) if name != SPARE]
    names = [ch_names[i] for i in keep]
    raw = raw[keep, :]
    return pd.DataFrame({
        "channel_name": np.repeat(names, len(times)),
        "crate":        crate,
        "meas_ms":      np.tile(times, len(names)),
        "position_mm":  raw.ravel(),
    })


def query_orbit(cycle_stamp_utc: str, selector: str, cycle_length_s: float,
                crate_bpm_map: dict[str, list[str]]) -> pd.DataFrame:
    """Query all 6 ALPS crates and concatenate into one long-format DataFrame (both planes)."""
    frames: list[pd.DataFrame] = []
    for crate in BPM_CRATES:
        sub = query_orbit_crate(crate, cycle_stamp_utc, selector, cycle_length_s, crate_bpm_map)
        if not sub.empty:
            frames.append(sub)
    if not frames:
        raise RuntimeError("All ALPS crate queries returned empty data.")
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Format -- turn a raw query result into the DataFrame(s) that get both assessed and saved.
# ---------------------------------------------------------------------------

def format_bct(times_ms: np.ndarray, intensity: np.ndarray, meta: dict) -> pd.DataFrame:
    df = pd.DataFrame({"sample_time_ms": times_ms, "intensity": intensity})
    for k, v in meta.items():
        df[k] = v
    return df


def format_bpms(orbit: pd.DataFrame, inj_offset_ms: float, bpm_positions: dict[str, float],
                meta: dict) -> dict[str, pd.DataFrame]:
    """
    Long-format orbit (query_orbit) -> {"H": df, "V": df}, each with
    channel_name, crate, sample_time_ms (from cycle stamp), meas_ms (from first injection),
    position_mm, s_m, inj_offset_ms + meta columns.
    """
    df = orbit.copy()
    df["sample_time_ms"] = df["meas_ms"] + inj_offset_ms
    df["s_m"] = df["channel_name"].map(bpm_positions)
    df["inj_offset_ms"] = inj_offset_ms
    df = df[["channel_name", "crate", "sample_time_ms", "meas_ms", "position_mm", "s_m",
             "inj_offset_ms"]]
    for k, v in meta.items():
        df[k] = v
    planes = df["channel_name"].map(bpm_plane)
    return {p: df[planes == p].reset_index(drop=True) for p in PLANES}


# ---------------------------------------------------------------------------
# Save -- validate the expected shape, then write.
# ---------------------------------------------------------------------------

BCT_REQUIRED_COLUMNS = {"sample_time_ms", "intensity"}
BPM_REQUIRED_COLUMNS = {"channel_name", "crate", "sample_time_ms", "meas_ms", "position_mm",
                        "s_m", "inj_offset_ms"}


def _validate(df: pd.DataFrame, required_columns: set[str], kind: str) -> None:
    if df.empty:
        raise ValueError(f"Refusing to save an empty {kind} DataFrame.")
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"{kind} DataFrame missing required column(s): {missing}")


def save_bct(df: pd.DataFrame, output: Path) -> None:
    _validate(df, BCT_REQUIRED_COLUMNS, "bct")
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output, index=False)


def save_bpms(df: pd.DataFrame, output: Path) -> None:
    _validate(df, BPM_REQUIRED_COLUMNS, "bpms")
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output, index=False)


# ---------------------------------------------------------------------------
# Paths / load
# ---------------------------------------------------------------------------

def data_dir(md_date: str, group: str, name: str) -> Path:
    return DATA_ROOT / md_date / group / name


def bct_path(md_date: str, group: str, name: str, timestamp: str) -> Path:
    return data_dir(md_date, group, name) / f"bct_{timestamp}.parquet"


def bpms_path(md_date: str, group: str, name: str, timestamp: str, plane: str) -> Path:
    if plane not in PLANES:
        raise ValueError(f"plane must be one of {PLANES}, got {plane!r}")
    return data_dir(md_date, group, name) / f"bpms_{plane}_{timestamp}.parquet"


def available_timestamps(md_date: str, group: str, name: str) -> list[str]:
    """Sorted local timestamps with a bct and both bpms_H/bpms_V parquet present."""
    d = data_dir(md_date, group, name)
    sets = [
        {TS_RE.search(f.name).group(1) for f in d.glob(pattern) if TS_RE.search(f.name)}
        for pattern in ("bct_*.parquet", "bpms_H_*.parquet", "bpms_V_*.parquet")
    ]
    return sorted(set.intersection(*sets))


def load_bct(md_date: str, group: str, name: str, timestamp: str) -> pd.DataFrame:
    return pd.read_parquet(bct_path(md_date, group, name, timestamp))


def load_bpms(md_date: str, group: str, name: str, timestamp: str, plane: str) -> pd.DataFrame:
    return pd.read_parquet(bpms_path(md_date, group, name, timestamp, plane))
