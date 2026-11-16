"""
Plot the orbit along the ring (H and V, coloured by time in cycle) for every qualifying cycle in a
timestamp table -- same shape as fetch_data.py; md_date/group inferred from --table's path.
"Qualifying" = beam_present and crates_saved both True (filled in by fetch_data.py).

Bad BPMs (data/<md_date>/bad_bpms.json) are left out unless --keep-bad. --sigma-clip N also
drops, per cycle and plane, BPMs whose cycle-median orbit is an N-sigma outlier across the ring
(iterative clipping, as done in the CCC for orbit correction -- data_treatment.sigma_clipped_bpms);
the dropped BPMs are printed. --exclude-bpms drops a fixed list instead (e.g. the references'
clipped set, so loss cycles use exactly the same BPMs). --aperture overlays the
machine aperture (apertls on an xsuite line, mapped to the sps.seq s frame -- see
bpm_core/lattice.py); computed once per run, ~2 s.

Loss cycles:
  --loss-window     only the orbit samples inside each cycle's loss window, found on its BCT with
                    data_treatment.get_window (the lossmap_framework detector, same window as the
                    loss maps), from --pre-ms before the window start (default 100) to its end;
                    every 10 ms sample unless --step-ms is given. The stats plot's time axis is
                    then relative to the window start. Cycles whose window has no orbit samples
                    (detector latched onto the dump) are skipped.
  --relative        plot each BPM's change from its mean over the second before the window
                    (data_treatment.relative_to_window) instead of the absolute orbit -- removes
                    the static orbit. Needs --loss-window. Aperture is not drawn (absolute).
  --min-intensity F drop samples taken while the BCT is below F of its peak
                    (data_treatment.gate_on_bct) -- ALPS readings are unreliable at low intensity
                    (BPM notes: below ~10 %).

Output: figures/<md_date>/<group>/<name>/<variant>/orbit_<name>_<local_ts>.<fmt>
        figures/<md_date>/<group>/<name>/<variant>/orbit_stats_<name>.<fmt>
The second is one figure per case: mean orbit across BPMs vs time (RMS across BPMs as error bar),
one line per cycle, computed on the same BPM set as the orbit plots (bad + --sigma-clip dropped,
--s-range/--time-range/--step-ms applied).
<variant> is one subdirectory per combination of plot options: "default" with none, otherwise the
option tags joined by "_" in a fixed order (e.g. "clip3sigma", "s5790-6920m_aperture").

Usage (from bpm_framework/):
    python scripts/plot_data.py --table ../data/20260401/references/timestamp_table.csv --name=-90_no_chirp
    python scripts/plot_data.py --table ../data/20260401/references/timestamp_table.csv --aperture
    python scripts/plot_data.py --table ../data/20260401/references/timestamp_table.csv --sigma-clip 3
    python scripts/plot_data.py --table ../data/20260401/chirps/timestamp_table.csv --aperture \
        --exclude-bpms ../data/20260401/references/excluded_bpms_clip3sigma.json
    python scripts/plot_data.py --table ../data/20260401/chirps/timestamp_table.csv --loss-window --relative \
        --min-intensity 0.1 --exclude-bpms ../data/20260401/references/excluded_bpms_clip3sigma.json
    python scripts/plot_data.py --table ../data/20260401/references/timestamp_table.csv --s-range 5800 6900 --aperture
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bpm_framework/

from bpm_core import data_treatment, io_data, lattice, plotting  # noqa: E402


LOSS_WINDOW_STEP_MS = 10.0   # every orbit sample inside a loss window (it lasts ~0.2-1 s)
DEFAULT_PRE_MS = 100.0       # orbit samples shown before the loss-window start
RELATIVE_REF_MS = 1000.0     # --relative: reference orbit = mean over this long before the window


def _default_step(args) -> float:
    return LOSS_WINDOW_STEP_MS if args.loss_window else plotting.DEFAULT_STEP_MS


def _variant(args) -> str:
    """Subdirectory name for this combination of plot options."""
    extra = []
    if args.loss_window:
        extra.append("losswindow" if args.pre_ms == DEFAULT_PRE_MS else f"losswindow-pre{args.pre_ms:g}ms")
    if args.relative:
        extra.append("relative")
    if args.min_intensity:
        extra.append(f"minI{args.min_intensity:g}")
    return plotting.variant_name(
        time_range_ms=args.time_range, s_range=args.s_range,
        aperture=args.aperture and not args.relative,
        keep_bad=args.keep_bad, sigma_clip=args.sigma_clip,
        step_ms=None if args.step_ms in (None, _default_step(args)) else args.step_ms,
        exclude_list=_exclude_label(args.exclude_bpms) if args.exclude_bpms else None,
        extra=extra,
    )


def _exclude_label(path: Path) -> str:
    """data/<date>/references/excluded_bpms_clip3sigma.json -> "references-clip3sigma"."""
    return f"{path.resolve().parent.name}-{path.stem.removeprefix('excluded_bpms_')}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=Path, required=True, help="Path to a timestamp_table.csv.")
    parser.add_argument("--name", help="Only plot this canonical case name.")
    parser.add_argument("--timestamp", help="Only plot this exact timestamp (requires --name).")
    parser.add_argument("--time-range", nargs=2, type=float, metavar=("START_MS", "END_MS"),
                        help="Only samples with sample_time_ms (from cycle stamp) in this range.")
    parser.add_argument("--step-ms", type=float,
                        help="Draw one orbit every STEP_MS (default 100, or 10 = every sample with "
                             "--loss-window).")
    parser.add_argument("--loss-window", action="store_true",
                        help="Only samples inside each cycle's BCT loss window (see docstring).")
    parser.add_argument("--pre-ms", type=float, default=DEFAULT_PRE_MS,
                        help="With --loss-window: also show this long before the window start.")
    parser.add_argument("--relative", action="store_true",
                        help="Plot the change from the orbit over the second before the window "
                             "(needs --loss-window).")
    parser.add_argument("--min-intensity", type=float, metavar="FRAC",
                        help="Drop samples with BCT below FRAC of its peak (e.g. 0.1).")
    parser.add_argument("--s-range", nargs=2, type=float, metavar=("S_MIN", "S_MAX"),
                        help="Zoom on this s range [m].")
    parser.add_argument("--aperture", action="store_true", help="Overlay the machine aperture.")
    parser.add_argument("--aperture-line", type=Path, default=lattice.APERTURE_LINE,
                        help="xsuite line with apertures (default: sps-crab-cavities matched line).")
    parser.add_argument("--keep-bad", action="store_true", help="Don't drop bad_bpms.json channels.")
    parser.add_argument("--sigma-clip", type=float, metavar="N",
                        help="Also drop BPMs whose cycle-median orbit is an N-sigma outlier "
                             "(iterative, per plane, per cycle). E.g. 3. Default: off.")
    parser.add_argument("--exclude-bpms", type=Path, metavar="JSON",
                        help="Also drop every BPM in this {channel: reason} file, e.g. the list "
                             "written by plot_overview.py for the references "
                             "(data/<date>/references/excluded_bpms_clip3sigma.json) -- same BPMs "
                             "for every cycle, no per-cycle clipping.")
    parser.add_argument("--cmap", default="plasma")
    parser.add_argument("--fmt", default="png", choices=["png", "pdf"])
    args = parser.parse_args()
    if args.timestamp and not args.name:
        parser.error("--timestamp requires --name.")
    if args.relative and not args.loss_window:
        parser.error("--relative requires --loss-window.")
    if args.loss_window and args.time_range:
        parser.error("--loss-window and --time-range are exclusive.")
    if args.relative and args.aperture:
        print("--relative: aperture not drawn (it is absolute).")
    step_ms = args.step_ms or _default_step(args)

    table = args.table.resolve()
    group, md_date = table.parent.name, table.parent.parent.name
    df = pd.read_csv(table)
    df = df[(df["beam_present"] == True) & (df["crates_saved"] == True)]  # noqa: E712
    if args.name:
        df = df[df["name"] == args.name]
    if args.timestamp:
        df = df[df["timestamp"] == args.timestamp]
    if df.empty:
        parser.error("No qualifying rows (beam_present & crates_saved) match the selection.")

    bad_bpms = {} if args.keep_bad else data_treatment.load_bad_bpms(md_date)
    if args.exclude_bpms:
        with open(args.exclude_bpms) as f:
            bad_bpms = {**bad_bpms, **json.load(f)}
        print(f"Excluding {len(bad_bpms)} BPMs (bad_bpms.json + {args.exclude_bpms.name})")
    aperture = lattice.compute_aperture(args.aperture_line) if args.aperture and not args.relative else None
    ylabels = {"H": "x - x(pre) [mm]", "V": "y - y(pre) [mm]"} if args.relative else None
    stats_ylabels = {"H": "Mean x - x(pre) [mm]", "V": "Mean y - y(pre) [mm]"} if args.relative else None
    variant = _variant(args)

    stats: dict[str, dict[str, dict[str, pd.DataFrame]]] = {}  # name -> plane -> label -> stats
    for _, row in df.iterrows():
        name = row["name"]
        local_ts = io_data.local_timestamp_from_utc(row["timestamp"])
        bpms = {p: io_data.load_bpms(md_date, group, name, local_ts, p) for p in io_data.PLANES}
        time_range, t0_ms = args.time_range, 0.0
        if args.loss_window or args.min_intensity:
            bct = io_data.load_bct(md_date, group, name, local_ts)
        if args.loss_window:
            ws, we, _ = data_treatment.get_window(bct["sample_time_ms"].to_numpy(float),
                                                  bct["intensity"].to_numpy(float))
            time_range, t0_ms = (ws - args.pre_ms, we), ws
            print(f"[{name} / {local_ts}] loss window {ws:.0f}-{we:.0f} ms")
            if args.relative:
                bpms = {p: data_treatment.relative_to_window(d, ws - RELATIVE_REF_MS, ws)
                        for p, d in bpms.items()}
        excluded = dict(bad_bpms or {})
        if args.sigma_clip:
            for plane, plane_df in bpms.items():
                clipped = data_treatment.sigma_clipped_bpms(plane_df, args.sigma_clip, exclude=bad_bpms)
                excluded.update({ch: f"{args.sigma_clip:g}sigma" for ch in clipped})
                listing = ", ".join(f"{ch} ({med:+.1f} mm)" for ch, med in clipped.items()) or "none"
                print(f"[{name} / {local_ts}] {plane} clipped: {listing}")
        if args.min_intensity:
            bpms = {p: data_treatment.gate_on_bct(d, bct, args.min_intensity) for p, d in bpms.items()}
        in_range = [d if time_range is None else d[d["sample_time_ms"].between(*time_range)]
                    for d in bpms.values()]
        if not any(np.isfinite(d["position_mm"]).any() for d in in_range):
            print(f"[{name} / {local_ts}] no orbit samples left (beam gone / below intensity cut) -- skipped")
            continue
        output = (io_data.FIGURES_ROOT / md_date / group / name / variant
                  / f"orbit_{name}_{local_ts}.{args.fmt}")
        plotting.plot_orbit_both_planes(
            bpms, output, bad_bpms=excluded, time_range_ms=time_range, step_ms=step_ms,
            cmap=args.cmap, aperture=aperture, s_range=args.s_range, ylabels=ylabels,
        )
        print(f"[{name} / {local_ts}] -> {output}")

        label = pd.Timestamp(row["timestamp"]).tz_convert("Europe/Zurich").strftime("%H:%M:%S")
        for plane in io_data.PLANES:
            st = plotting.orbit_stats(
                bpms[plane], bad_bpms=excluded, time_range_ms=time_range,
                step_ms=step_ms, s_range=args.s_range,
            )
            st["sample_time_ms"] -= t0_ms
            stats.setdefault(name, {}).setdefault(plane, {})[label] = st

    for name, stats_by_plane in stats.items():
        output = io_data.FIGURES_ROOT / md_date / group / name / variant / f"orbit_stats_{name}.{args.fmt}"
        plotting.plot_orbit_stats(
            stats_by_plane, output, ylabels=stats_ylabels,
            xlabel="Time from loss-window start [s]" if args.loss_window else "Time in cycle [s]",
        )
        print(f"[{name}] stats -> {output}")


if __name__ == "__main__":
    main()
