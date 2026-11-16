"""
Overview across the cases of one timestamp table (e.g. the 20260401 references: one case per crab
phase). For each case and plane, a "mean orbit": per BPM, the mean over every sample of every
qualifying cycle of that case, plus the per-BPM RMS (std of those samples) -- see
data_treatment.mean_orbit.

BPM selection: bad_bpms.json (unless --keep-bad) plus, with --sigma-clip N, every BPM clipped in
ANY qualifying cycle of the table (data_treatment.sigma_clipped_bpms, per cycle and plane) -- the
union, so every case is averaged over the same BPMs and the per-case numbers stay comparable.

Cycles listed in config.OVERVIEW_EXCLUDE are left out entirely (also from the clipping).

Delta orbit: per phase, mean orbit minus the mean orbit of the reference phase
(config.ORBIT_REFERENCE_PHASE, or --ref-phase), per BPM -- removes the static orbit, leaving the
crab-induced orbit change.

The excluded BPM list (bad + clipped union) is also written to
data/<md_date>/<group>/excluded_bpms_<variant>.json, for reuse on other tables
(scripts/plot_data.py --exclude-bpms), e.g. the chirp cycles, which are then cut with exactly the
same BPMs as the references.

Figures, in figures/<md_date>/<group>/overview/<variant>/:
    mean_orbit_by_phase.<fmt>   per phase: mean over BPMs of the mean orbit, RMS over BPMs as
                                error bar (H top, V bottom)
    mean_orbit_vs_s.<fmt>       mean orbit along s, one line per phase, coloured by phase
    rms_vs_s.<fmt>              per-BPM RMS along s, one line per phase, coloured by phase
    delta_orbit_vs_s.<fmt>      mean orbit minus reference-phase orbit along s, coloured by phase
    delta_orbit_by_phase.<fmt>  per phase: mean over BPMs of the delta, RMS over BPMs as error bar
    delta_rms_by_phase.<fmt>    per phase: RMS over BPMs of the delta, sqrt(mean(delta**2))
    mean_orbit_summary.csv      the numbers behind them: crab_phase, name, plane, channel_name,
                                s_m, mean_mm, rms_mm, n_samples, n_cycles, delta_mm
--s-range restricts everything (figures and the per-phase means) to BPMs in that range.

Usage (from bpm_framework/):
    python scripts/plot_overview.py --table ../data/20260401/references/timestamp_table.csv --sigma-clip 3
    python scripts/plot_overview.py --table ../data/20260401/references/timestamp_table.csv --sigma-clip 3 --s-range 5790 6920
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bpm_framework/

from bpm_core import data_treatment, io_data, plotting  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=Path, required=True, help="Path to a timestamp_table.csv.")
    parser.add_argument("--sigma-clip", type=float, metavar="N",
                        help="Drop BPMs that are an N-sigma outlier in any cycle (e.g. 3). Default: off.")
    parser.add_argument("--keep-bad", action="store_true", help="Don't drop bad_bpms.json channels.")
    parser.add_argument("--time-range", nargs=2, type=float, metavar=("START_MS", "END_MS"),
                        help="Only samples with sample_time_ms (from cycle stamp) in this range.")
    parser.add_argument("--s-range", nargs=2, type=float, metavar=("S_MIN", "S_MAX"),
                        help="Only BPMs in this s range [m].")
    parser.add_argument("--ref-phase", type=float,
                        help="Crab phase subtracted in the delta plots (default: config.ORBIT_REFERENCE_PHASE).")
    parser.add_argument("--cmap", default=plotting.PHASE_CMAP)
    parser.add_argument("--fmt", default="png", choices=["png", "pdf"])
    args = parser.parse_args()

    table = args.table.resolve()
    group, md_date = table.parent.name, table.parent.parent.name
    df = pd.read_csv(table)
    df = df[(df["beam_present"] == True) & (df["crates_saved"] == True)]  # noqa: E712
    if df.empty:
        parser.error("No qualifying rows (beam_present & crates_saved) in the table.")
    df = df.assign(crab_phase=pd.to_numeric(df["crab_phase"]))

    md_config = importlib.import_module(f"config.md_{md_date}")
    skip = set(getattr(md_config, "OVERVIEW_EXCLUDE", []))
    skipped = df[[(n, t) in skip for n, t in zip(df["name"], df["timestamp"])]]
    for _, row in skipped.iterrows():
        print(f"Left out (config.OVERVIEW_EXCLUDE): {row['name']} / {row['timestamp']}")
    df = df.drop(skipped.index)
    ref_phase = args.ref_phase if args.ref_phase is not None else md_config.ORBIT_REFERENCE_PHASE
    if ref_phase not in set(df["crab_phase"]):
        parser.error(f"Reference phase {ref_phase} not among {sorted(df['crab_phase'].unique())}.")

    bad_bpms = {} if args.keep_bad else data_treatment.load_bad_bpms(md_date)

    # load every cycle once, collect the union of clipped BPMs
    frames: dict[str, dict[str, list[pd.DataFrame]]] = {}  # name -> plane -> frames
    excluded = dict(bad_bpms)
    for _, row in df.iterrows():
        local_ts = io_data.local_timestamp_from_utc(row["timestamp"])
        for plane in io_data.PLANES:
            bpms = io_data.load_bpms(md_date, group, row["name"], local_ts, plane)
            frames.setdefault(row["name"], {}).setdefault(plane, []).append(bpms)
            if args.sigma_clip:
                clipped = data_treatment.sigma_clipped_bpms(bpms, args.sigma_clip, exclude=bad_bpms)
                excluded.update({ch: f"{args.sigma_clip:g}sigma" for ch in clipped})
    clipped_all = sorted(ch for ch, why in excluded.items() if why != bad_bpms.get(ch))
    print(f"Excluded: {len(bad_bpms)} bad + {len(clipped_all)} clipped in at least one cycle: "
          f"{', '.join(clipped_all) or 'none'}")

    phase_by_name = df.groupby("name")["crab_phase"].first()
    rows = []
    for name, by_plane in frames.items():
        for plane, plane_frames in by_plane.items():
            mo = data_treatment.mean_orbit(plane_frames, exclude=excluded, time_range_ms=args.time_range)
            mo.insert(0, "plane", plane)
            mo.insert(0, "name", name)
            mo.insert(0, "crab_phase", phase_by_name[name])
            mo["n_cycles"] = len(plane_frames)
            rows.append(mo)
    summary = pd.concat(rows, ignore_index=True).sort_values(["crab_phase", "plane", "s_m"],
                                                             ignore_index=True)
    ref = summary[summary["crab_phase"] == ref_phase].set_index(["plane", "channel_name"])["mean_mm"]
    summary["delta_mm"] = summary["mean_mm"] - pd.Series(
        list(zip(summary["plane"], summary["channel_name"]))).map(ref).to_numpy()
    if args.s_range:
        summary = summary[summary["s_m"].between(*args.s_range)]

    variant = plotting.variant_name(time_range_ms=args.time_range, s_range=args.s_range,
                                    keep_bad=args.keep_bad, sigma_clip=args.sigma_clip)
    out_dir = io_data.FIGURES_ROOT / md_date / group / "overview" / variant
    out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_dir / "mean_orbit_summary.csv", index=False)
    excluded_json = io_data.DATA_ROOT / md_date / group / f"excluded_bpms_{variant}.json"
    with open(excluded_json, "w") as f:
        json.dump(dict(sorted(excluded.items())), f, indent=4)
    print(f"Excluded BPM list -> {excluded_json}")

    plotting.plot_mean_orbit_by_phase(summary, out_dir / f"mean_orbit_by_phase.{args.fmt}")
    plotting.plot_orbit_vs_s_by_phase(summary, out_dir / f"mean_orbit_vs_s.{args.fmt}",
                                      column="mean_mm", s_range=args.s_range, cmap=args.cmap,
                                      ylabel={"H": "Mean x [mm]", "V": "Mean y [mm]"})
    plotting.plot_orbit_vs_s_by_phase(summary, out_dir / f"rms_vs_s.{args.fmt}",
                                      column="rms_mm", s_range=args.s_range, cmap=args.cmap,
                                      ylabel={"H": "RMS x [mm]", "V": "RMS y [mm]"})
    ref_label = f"{ref_phase:g}"
    plotting.plot_orbit_vs_s_by_phase(summary, out_dir / f"delta_orbit_vs_s.{args.fmt}",
                                      column="delta_mm", s_range=args.s_range, cmap=args.cmap,
                                      ylabel={"H": f"x - x({ref_label}°) [mm]",
                                              "V": f"y - y({ref_label}°) [mm]"})
    plotting.plot_mean_orbit_by_phase(summary, out_dir / f"delta_orbit_by_phase.{args.fmt}",
                                      column="delta_mm",
                                      ylabel={"H": f"Mean x - x({ref_label}°) [mm]",
                                              "V": f"Mean y - y({ref_label}°) [mm]"})
    plotting.plot_rms_by_phase(summary, out_dir / f"delta_rms_by_phase.{args.fmt}",
                               column="delta_mm",
                               ylabel={"H": f"RMS x - x({ref_label}°) [mm]",
                                       "V": f"RMS y - y({ref_label}°) [mm]"})
    print(f"Reference phase {ref_label} deg. Saved to {out_dir}")


if __name__ == "__main__":
    main()
