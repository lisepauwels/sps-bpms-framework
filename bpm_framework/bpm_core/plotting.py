"""
Orbit plotting: pure functions, data in (DataFrames), figures out -- no table looping here (that's
scripts/plot_data.py), same split as lossmap_framework's plotting.py. No titles on the figures
(name/timestamp are in the filename).

plot_orbit_vs_s: one line per orbit sample, BPMs along s, coloured by time in the cycle (s from
the cycle stamp, i.e. sample_time_ms -- same axis as the BLM files). Samples after the beam is
gone are all-NaN and simply not drawn. Optional aperture overlay (lattice.compute_aperture).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

plt.rcParams.update({"font.family": "serif", "font.size": 17})

PLANE_YLABEL = {"H": "x [mm]", "V": "y [mm]"}
DEFAULT_STEP_MS = 100.0
APERTURE_COLUMNS = {"H": ("x_min_mm", "x_max_mm"), "V": ("y_min_mm", "y_max_mm")}


def variant_name(
    time_range_ms: tuple[float, float] | None = None,
    s_range: tuple[float, float] | None = None,
    aperture: bool = False,
    keep_bad: bool = False,
    sigma_clip: float | None = None,
    step_ms: float | None = None,
    exclude_list: str | None = None,
    extra: list[str] | None = None,
) -> str:
    """Figure subdirectory name for one combination of plot options: "default" with none, else
    the option tags joined by "_" in a fixed order (step_ms: pass None for the default step;
    exclude_list: a short label for a reused excluded-BPM list, e.g. "references-clip3sigma")."""
    tags = []
    if time_range_ms:
        tags.append(f"t{time_range_ms[0]:.0f}-{time_range_ms[1]:.0f}ms")
    if s_range:
        tags.append(f"s{s_range[0]:.0f}-{s_range[1]:.0f}m")
    if aperture:
        tags.append("aperture")
    if keep_bad:
        tags.append("withbad")
    if sigma_clip:
        tags.append(f"clip{sigma_clip:g}sigma")
    if step_ms:
        tags.append(f"step{step_ms:g}ms")
    if exclude_list:
        tags.append(f"excl-{exclude_list}")
    tags.extend(extra or [])
    return "_".join(tags) or "default"


def orbit_matrix(
    bpms: pd.DataFrame,
    bad_bpms: dict[str, str] | None = None,
    time_range_ms: tuple[float, float] | None = None,
    step_ms: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """
    Long-format bpms frame (one plane) -> (s_m [n_bpm], t_ms [n_t], positions [n_t, n_bpm],
    channel names), BPMs sorted by s. Drops bad_bpms channels, keeps samples inside time_range_ms
    (sample_time_ms, from cycle stamp), thins to one sample every step_ms, and drops samples with
    no finite BPM (after the beam is gone).
    """
    df = bpms
    if bad_bpms:
        df = df[~df["channel_name"].isin(bad_bpms)]
    wide = df.pivot_table(index="sample_time_ms", columns="channel_name", values="position_mm",
                          aggfunc="first", dropna=False)
    s_by_channel = df.groupby("channel_name")["s_m"].first()
    channels = sorted(wide.columns, key=lambda c: (s_by_channel[c], c))
    wide = wide[channels]

    t = wide.index.to_numpy(dtype=float)
    keep = np.ones(len(t), dtype=bool)
    if time_range_ms is not None:
        keep &= (t >= time_range_ms[0]) & (t <= time_range_ms[1])
    if step_ms:
        keep &= np.isclose(np.mod(t - t[0], step_ms), 0.0) | np.isclose(np.mod(t - t[0], step_ms), step_ms)
    pos = wide.to_numpy(dtype=float)
    keep &= np.isfinite(pos).any(axis=1)
    return s_by_channel[channels].to_numpy(dtype=float), t[keep], pos[keep], channels


def plot_orbit_vs_s(
    ax: plt.Axes,
    bpms: pd.DataFrame,
    plane: str,
    bad_bpms: dict[str, str] | None = None,
    time_range_ms: tuple[float, float] | None = None,
    step_ms: float | None = DEFAULT_STEP_MS,
    norm: Normalize | None = None,
    cmap: str = "plasma",
    aperture: pd.DataFrame | None = None,
    s_range: tuple[float, float] | None = None,
    ylabel: str | None = None,
) -> ScalarMappable:
    """Draw one plane's orbit vs s on ax, one line per sample coloured by time [s]. Returns the
    ScalarMappable for a colorbar."""
    s, t_ms, pos, _ = orbit_matrix(bpms, bad_bpms, time_range_ms, step_ms)
    t_s = t_ms / 1000.0
    if norm is None:
        norm = Normalize(vmin=t_s.min(), vmax=t_s.max())
    segments = []
    for row in pos:
        finite = np.isfinite(row)
        segments.append(np.column_stack([s[finite], row[finite]]))
    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=0.8, alpha=0.8)
    lc.set_array(t_s)
    ax.add_collection(lc)

    if aperture is not None:
        lo, hi = APERTURE_COLUMNS[plane]
        ap = aperture.sort_values("s_m")
        ax.plot(ap["s_m"], ap[lo], color="k", lw=1.0, label="Aperture")
        ax.plot(ap["s_m"], ap[hi], color="k", lw=1.0)

    if s_range is not None:
        ax.set_xlim(*s_range)
        in_range = (s >= s_range[0]) & (s <= s_range[1])
        y = pos[:, in_range]
    else:
        ax.set_xlim(0, max(s.max(), 6911.5))
        y = pos
    if aperture is None and np.isfinite(y).any():
        lo_y, hi_y = np.nanmin(y), np.nanmax(y)
        pad = 0.05 * (hi_y - lo_y or 1.0)
        ax.set_ylim(lo_y - pad, hi_y + pad)
    elif aperture is not None:
        ap_lo, ap_hi = APERTURE_COLUMNS[plane]
        ap = aperture
        if s_range is not None:
            ap = ap[(ap["s_m"] >= s_range[0]) & (ap["s_m"] <= s_range[1])]
        ax.set_ylim(1.1 * ap[ap_lo].min(), 1.1 * ap[ap_hi].max())

    ax.set_ylabel(ylabel or PLANE_YLABEL[plane])
    ax.grid(alpha=0.3)
    return ScalarMappable(norm=norm, cmap=cmap)


def plot_orbit_both_planes(
    bpms_by_plane: dict[str, pd.DataFrame],
    output: Path,
    bad_bpms: dict[str, str] | None = None,
    time_range_ms: tuple[float, float] | None = None,
    step_ms: float | None = DEFAULT_STEP_MS,
    cmap: str = "plasma",
    aperture: pd.DataFrame | None = None,
    s_range: tuple[float, float] | None = None,
    ylabels: dict[str, str] | None = None,
) -> None:
    """H (top) and V (bottom) orbit vs s, shared s axis, one shared time colorbar. Saves to output."""
    t_all = np.concatenate([
        orbit_matrix(df, bad_bpms, time_range_ms, step_ms)[1] for df in bpms_by_plane.values()
    ]) / 1000.0
    norm = Normalize(vmin=t_all.min(), vmax=t_all.max())

    fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
    for ax, plane in zip(axes, ("H", "V")):
        mappable = plot_orbit_vs_s(
            ax, bpms_by_plane[plane], plane, bad_bpms=bad_bpms, time_range_ms=time_range_ms,
            step_ms=step_ms, norm=norm, cmap=cmap, aperture=aperture, s_range=s_range,
            ylabel=(ylabels or {}).get(plane),
        )
    axes[-1].set_xlabel("s [m]")
    fig.tight_layout()
    cbar = fig.colorbar(mappable, ax=axes, pad=0.01)
    cbar.set_label("Time in cycle [s]")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Orbit statistics across BPMs vs time
# ---------------------------------------------------------------------------

MIN_FINITE_FRAC = 0.5  # a sample needs at least this fraction of its plane's BPMs finite to get a
                        # mean/RMS -- as the beam disappears the few BPMs still reading dominate a
                        # plain nanmean (cf. the >=120-of-224 guard in the BPM notes, 7.6)


def orbit_stats(
    bpms: pd.DataFrame,
    bad_bpms: dict[str, str] | None = None,
    time_range_ms: tuple[float, float] | None = None,
    step_ms: float | None = DEFAULT_STEP_MS,
    s_range: tuple[float, float] | None = None,
    min_finite_frac: float = MIN_FINITE_FRAC,
) -> pd.DataFrame:
    """
    Per orbit sample, across the BPMs of one plane (bad_bpms dropped, optionally only those in
    s_range): mean and RMS = standard deviation around that mean [mm], plus the number of finite
    BPMs. Samples with fewer than min_finite_frac of the BPMs finite are dropped.
    Returns sample_time_ms, mean_mm, rms_mm, n_bpms.
    """
    s, t, pos, _ = orbit_matrix(bpms, bad_bpms, time_range_ms, step_ms)
    if s_range is not None:
        pos = pos[:, (s >= s_range[0]) & (s <= s_range[1])]
    n = np.isfinite(pos).sum(axis=1)
    keep = n >= min_finite_frac * pos.shape[1]
    with np.errstate(all="ignore"):
        mean = np.nanmean(pos[keep], axis=1)
        rms = np.nanstd(pos[keep], axis=1)
    return pd.DataFrame({"sample_time_ms": t[keep], "mean_mm": mean, "rms_mm": rms, "n_bpms": n[keep]})


def plot_orbit_stats(
    stats_by_plane: dict[str, dict[str, pd.DataFrame]], output: Path,
    ylabels: dict[str, str] | None = None, xlabel: str = "Time in cycle [s]",
) -> None:
    """
    H (top) / V (bottom): mean orbit across BPMs vs time in cycle, RMS across BPMs as error bar,
    one line per cycle (stats_by_plane[plane][label] from orbit_stats; label = legend entry).
    """
    fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
    colors = plt.get_cmap("tab10")
    for ax, plane in zip(axes, ("H", "V")):
        for i, (label, st) in enumerate(stats_by_plane[plane].items()):
            ax.errorbar(st["sample_time_ms"] / 1000.0, st["mean_mm"], yerr=st["rms_mm"],
                        color=colors(i % 10), marker="o", markersize=3, lw=1.0, elinewidth=0.6,
                        capsize=2, alpha=0.8, label=label)
        ax.set_ylabel((ylabels or {}).get(plane, f"Mean {PLANE_YLABEL[plane]}"))
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=13, ncol=2)
    axes[-1].set_xlabel(xlabel)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Overview across cases (e.g. crab phases): one mean orbit per case
# ---------------------------------------------------------------------------

PHASE_CMAP = "coolwarm"  # diverging, blue (negative phase) -> red (positive phase)


def plot_mean_orbit_by_phase(
    summary: pd.DataFrame, output: Path, column: str = "mean_mm",
    ylabel: dict[str, str] | None = None,
) -> None:
    """
    summary: one row per (phase, plane, BPM) with crab_phase, plane and `column` (e.g. mean_mm
    from data_treatment.mean_orbit, or delta_mm). Per phase, the mean over BPMs of `column`, with
    the RMS (std) over BPMs as error bar; H top, V bottom.
    """
    ylabel = ylabel or {p: f"Mean {PLANE_YLABEL[p]}" for p in ("H", "V")}
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    for ax, plane in zip(axes, ("H", "V")):
        g = summary[summary["plane"] == plane].groupby("crab_phase")[column]
        st = pd.DataFrame({"mean": g.mean(), "rms": g.std()}).sort_index()
        ax.errorbar(st.index, st["mean"], yerr=st["rms"], fmt="o-", color="k", ms=6,
                    capsize=4, lw=1.2)
        ax.set_ylabel(ylabel[plane])
        ax.grid(alpha=0.3)
    _phase_xaxis(axes[-1], summary)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rms_by_phase(
    summary: pd.DataFrame, output: Path, column: str = "delta_mm",
    ylabel: dict[str, str] | None = None,
) -> None:
    """Per phase, RMS over BPMs of `column` = sqrt(mean(column**2)) (e.g. the RMS of the orbit
    difference to the reference phase); H top, V bottom."""
    ylabel = ylabel or {p: f"RMS {PLANE_YLABEL[p]}" for p in ("H", "V")}
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    for ax, plane in zip(axes, ("H", "V")):
        sub = summary[summary["plane"] == plane]
        rms = sub.groupby("crab_phase")[column].apply(lambda v: np.sqrt(np.nanmean(v ** 2))).sort_index()
        ax.plot(rms.index, rms.to_numpy(), "o-", color="k", ms=6, lw=1.2)
        ax.set_ylabel(ylabel[plane])
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.3)
    _phase_xaxis(axes[-1], summary)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _phase_xaxis(ax: plt.Axes, summary: pd.DataFrame) -> None:
    ax.set_xlabel("Crab phase [deg]")
    ax.set_xticks(sorted(summary["crab_phase"].unique()))
    ax.tick_params(axis="x", labelrotation=45)


def plot_orbit_vs_s_by_phase(
    summary: pd.DataFrame, output: Path, column: str = "mean_mm",
    ylabel: dict[str, str] | None = None, s_range: tuple[float, float] | None = None,
    cmap: str = PHASE_CMAP,
) -> None:
    """
    summary as in plot_mean_orbit_by_phase. One line per phase along s (BPMs), coloured by phase,
    H top / V bottom. column="mean_mm" for the mean orbit, "rms_mm" for the per-BPM RMS.
    """
    ylabel = ylabel or {p: PLANE_YLABEL[p] for p in ("H", "V")}
    phases = sorted(summary["crab_phase"].unique())
    lim = max(abs(min(phases)), abs(max(phases)))
    norm = Normalize(vmin=-lim, vmax=lim)
    cm = plt.get_cmap(cmap)

    fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
    for ax, plane in zip(axes, ("H", "V")):
        sub = summary[summary["plane"] == plane]
        for phase in phases:
            ph = sub[sub["crab_phase"] == phase].sort_values("s_m")
            ax.plot(ph["s_m"], ph[column], color=cm(norm(phase)), lw=1.0, marker="o", ms=2.5)
        ax.set_ylabel(ylabel[plane])
        ax.grid(alpha=0.3)
        if s_range is not None:
            ax.set_xlim(*s_range)
        else:
            ax.set_xlim(0, 6911.5)
    axes[-1].set_xlabel("s [m]")
    fig.tight_layout()
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), ax=axes, pad=0.01)
    cbar.set_label("Crab phase [deg]")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
