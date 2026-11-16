"""
Per-cycle quality checks: is there beam (BCT), did every ALPS crate deliver an orbit. Operates on
plain bct/bpm DataFrames (io_data.load_*/format_*), doesn't care whether they came from disk or a
fresh query. flag_noisy_fill is copied from lossmap_framework/lossmap_core/data_treatment.py.

Analysis-time treatment (bad-channel masking from data/<md_date>/bad_bpms.json, BCT gating of
orbit samples after beam loss, stable means, reference orbit, deltas) comes later -- see README.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import io_data

NOISY_MAX_INTENSITY_THRESHOLD = 1.0  # same value as lossmap_framework (verified on 20260826)


def flag_noisy_fill(bct: pd.DataFrame, threshold: float = NOISY_MAX_INTENSITY_THRESHOLD) -> bool:
    """True if peak BCT intensity never exceeds `threshold` -- no beam this cycle."""
    return bool(bct["intensity"].max() < threshold)


def check_bpm_crates(bpms: list[pd.DataFrame]) -> dict:
    """
    A crate is "down" for this cycle when none of its rows carries a finite position (absent
    from the query, or present but all-NaN). `bpms` = the cycle's per-plane frames.
    """
    df = pd.concat(bpms, ignore_index=True)
    finite = df[np.isfinite(df["position_mm"])]
    present = set(finite["crate"].unique())
    down = sorted(set(io_data.BPM_CRATES) - present)
    return {"n_bpm_crates_down": len(down), "bpm_crates_down": down}


def assess_quality(
    bct: pd.DataFrame, bpms: list[pd.DataFrame],
    noisy_threshold: float = NOISY_MAX_INTENSITY_THRESHOLD,
) -> dict:
    noisy_fill = flag_noisy_fill(bct, noisy_threshold)
    crate_info = check_bpm_crates(bpms)
    return {"noisy_fill": noisy_fill, **crate_info}


# ---------------------------------------------------------------------------
# Loss window: (win_start_ms, win_end_ms, peak) for one cycle's BCT trace. Copied verbatim from
# lossmap_framework/lossmap_core/data_treatment.py, so BPM and BLM use the same window.
#
# get_window() is the single entry point -- automatic detection by default, or an exact override
# if you already know the window you want (e.g. 20260629's hand-tuned cases, or a future "slow
# losses" window). Same signature either way, so callers don't need to branch.
# ---------------------------------------------------------------------------

INJECTION_BLANK_MS  = 100.0    # ignore drops within this window after the peak -- fast
                                # post-injection settling, not a real loss
RATE_WINDOW_MS       = 50.0    # trailing window used to estimate the local loss rate
RATE_SMOOTH_SAMPLES  = 9       # centered moving average on the rate curve before thresholding
MAX_SLOW_RATE         = 0.00015  # fraction of peak lost per ms that still counts as "slow"
MIN_WINDOW_MS         = 200.0

ONSET_BACKOFF_SAMPLES = 2      # the rate-based onset search lags the true onset by ~1-2 samples
LOW_FRAC               = 0.03  # "settled" -- 3% of peak (above the true noise floor, below any
                                # genuine residual plateau seen in past campaigns)
SMOOTH_SAMPLES         = 9     # rolling-median window used only to decide when settled
PERSIST_SAMPLES        = 10    # smoothed signal must stay below the floor this many consecutive
                                # samples to count as settled -- one noisy/late sample can't move it


def _find_settled_end_idx(I: np.ndarray, start_idx: int, low: float) -> int:
    smoothed = pd.Series(I).rolling(SMOOTH_SAMPLES, center=True, min_periods=1).median().to_numpy()
    n = len(I)
    for i in range(start_idx, n):
        tail = smoothed[i:i + PERSIST_SAMPLES]
        if smoothed[i] < low and np.all(tail < low):
            return i
    return n - 1


def detect_window(times_ms: np.ndarray, intensity: np.ndarray) -> tuple[float, float, float]:
    """Automatic (win_start_ms, win_end_ms, peak). Only meaningful for a cycle with one real,
    confirmed loss event -- not a gradual whole-cycle loss with no single sharp onset."""
    times_ms = np.asarray(times_ms, float)
    t = times_ms[2:]                     # skip first 2 samples (acquisition artefacts)
    I = np.asarray(intensity, float)[2:]
    peak = float(np.max(I))
    norm = I / peak

    peak_idx = int(np.argmax(I))
    t_peak = float(t[peak_idx])
    search_start = t_peak + INJECTION_BLANK_MS + RATE_WINDOW_MS
    step_ms = float(np.median(np.diff(t)))
    lag = max(1, round(RATE_WINDOW_MS / step_ms))

    rate_raw = np.full(len(norm), np.nan)
    rate_raw[lag:] = (norm[:-lag] - norm[lag:]) / RATE_WINDOW_MS
    rate = pd.Series(rate_raw).rolling(RATE_SMOOTH_SAMPLES, center=True, min_periods=1).mean().to_numpy()

    i_search_start = int(np.searchsorted(t, search_start))
    global_peak_idx = i_search_start + int(np.nanargmax(rate[i_search_start:]))

    if rate[global_peak_idx] >= MAX_SLOW_RATE:
        onset_idx = global_peak_idx
        while onset_idx > i_search_start and rate[onset_idx - 1] >= MAX_SLOW_RATE:
            onset_idx -= 1
    else:
        onset_idx = int(np.argmax(t >= search_start))
    onset_idx = max(i_search_start, onset_idx - ONSET_BACKOFF_SAMPLES)

    low = LOW_FRAC * peak
    end_idx = _find_settled_end_idx(I, onset_idx, low)

    win_end = max(float(t[end_idx]), float(t[onset_idx]) + MIN_WINDOW_MS)
    return float(t[onset_idx]), win_end, peak


def get_window(times_ms: np.ndarray, intensity: np.ndarray,
               override_ms: tuple[float, float] | None = None) -> tuple[float, float, float]:
    """
    (win_start_ms, win_end_ms, peak) for one cycle. Automatic detection by default; pass
    override_ms=(start_ms, end_ms) to use an exact window instead (hand-tuned cases, a future
    slow-losses window, etc.) -- same return shape either way.
    """
    if override_ms is not None:
        ws, we = override_ms
        peak = float(np.max(np.asarray(intensity, float)[2:]))
        return float(ws), float(we), peak
    return detect_window(times_ms, intensity)


def gate_on_bct(bpms: pd.DataFrame, bct: pd.DataFrame, min_frac: float) -> pd.DataFrame:
    """
    Drop orbit samples taken while the BCT is below min_frac of its peak (BCT interpolated onto
    each sample's sample_time_ms; peak skips the first 2 BCT samples, as get_window does). ALPS
    readings at low intensity are not trustworthy (resolution degrades with bunch charge).
    """
    t = bct["sample_time_ms"].to_numpy(dtype=float)
    intensity = bct["intensity"].to_numpy(dtype=float)
    peak = float(np.max(intensity[2:]))
    frac = np.interp(bpms["sample_time_ms"].to_numpy(dtype=float), t, intensity) / peak
    return bpms[frac >= min_frac]


def relative_to_window(bpms: pd.DataFrame, t0_ms: float, t1_ms: float) -> pd.DataFrame:
    """position_mm minus each BPM's mean position over [t0_ms, t1_ms) (sample_time_ms), e.g. the
    second before a loss window -- leaves only the orbit change, without the static orbit."""
    pre = bpms[(bpms["sample_time_ms"] >= t0_ms) & (bpms["sample_time_ms"] < t1_ms)]
    ref = pre.groupby("channel_name")["position_mm"].mean()
    out = bpms.copy()
    out["position_mm"] = out["position_mm"] - out["channel_name"].map(ref)
    return out


def clip_mask(arr: np.ndarray, nsigma: float = 3.0) -> np.ndarray:
    """Iterative sigma clipping (CCC orbit-correction behaviour, as in
    sps-crab-cavities/studies/twiss_analysis.ipynb): drop points further than nsigma*std from the
    mean of the kept points, repeat until nothing changes. Returns a mask, True = kept. NaN
    entries are never kept and don't enter mean/std."""
    arr = np.asarray(arr, dtype=float)
    kept = np.isfinite(arr)
    while kept.any():
        sub = arr[kept]
        mean, std = sub.mean(), sub.std()
        new_kept = kept.copy()
        new_kept[kept] = np.abs(sub - mean) < nsigma * std
        if np.array_equal(new_kept, kept):
            break
        kept = new_kept
    return kept


def sigma_clipped_bpms(bpms: pd.DataFrame, nsigma: float = 3.0,
                       exclude: dict[str, str] | None = None) -> dict[str, float]:
    """
    BPMs of one plane whose cycle-median orbit is an nsigma outlier across the ring (clip_mask on
    the per-BPM medians over all finite samples). One decision per BPM per cycle, so every sample
    is drawn/analysed with the same BPM set. `exclude` channels (e.g. bad_bpms) are removed before
    clipping, so they don't inflate the std. Returns {channel_name: median_mm} of clipped BPMs.
    """
    df = bpms if not exclude else bpms[~bpms["channel_name"].isin(exclude)]
    medians = df.groupby("channel_name")["position_mm"].median()
    kept = clip_mask(medians.to_numpy(), nsigma)
    clipped = medians[~kept & np.isfinite(medians.to_numpy())]
    return clipped.to_dict()


def mean_orbit(
    frames: list[pd.DataFrame],
    exclude: dict[str, str] | None = None,
    time_range_ms: tuple[float, float] | None = None,
) -> pd.DataFrame:
    """
    Per BPM, over every finite sample of every frame (e.g. all cycles of one case, one plane):
    mean_mm = mean orbit, rms_mm = standard deviation of the samples around it (how much that BPM
    moves over time/cycles), n_samples. `exclude` channels dropped; time_range_ms filters on
    sample_time_ms (from cycle stamp). Returns channel_name, s_m, mean_mm, rms_mm, n_samples,
    sorted by s.
    """
    df = pd.concat(frames, ignore_index=True)
    if exclude:
        df = df[~df["channel_name"].isin(exclude)]
    if time_range_ms is not None:
        df = df[df["sample_time_ms"].between(*time_range_ms)]
    df = df[np.isfinite(df["position_mm"])]
    out = df.groupby("channel_name").agg(
        s_m=("s_m", "first"), mean_mm=("position_mm", "mean"),
        rms_mm=("position_mm", "std"), n_samples=("position_mm", "size"),
    )
    return out.reset_index().sort_values(["s_m", "channel_name"], ignore_index=True)


def load_bad_bpms(md_date: str) -> dict[str, str]:
    """data/<md_date>/bad_bpms.json -> {channel_name: reason} ("faulty"/"offset"/...)."""
    with open(io_data.DATA_ROOT / md_date / "bad_bpms.json") as f:
        return json.load(f)
