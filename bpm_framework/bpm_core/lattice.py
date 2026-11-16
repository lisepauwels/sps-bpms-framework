"""
Static geometry: s positions from acc-models-sps/sps.seq (refer = centre -- the convention of both
bpm_positions.json and lossmap_framework's blm_positions.json), and the machine aperture from an
xsuite line via apertls, mapped onto that same s frame.

The xsuite lines in sps-crab-cavities are NOT on the sps.seq s frame: sps_crab_scan_matched.json
is 3.14 m shorter (6908.37 vs 6911.52 m), the difference building up through sextants 1-3 and
constant at ~3.27 m from s ~2975 m on (checked 2026-11-13 on 194 BPMs present in both).
line_s_to_seq_s() interpolates line s -> sps.seq s through every element present in both, so
the aperture lands on the same axis as BPM/BLM positions.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from . import io_data

SPS_SEQ = (
    io_data.ROOT.parents[2] / "code" / "sps-xsuite-model" / "external" / "acc-models-sps" / "sps.seq"
)
APERTURE_LINE = (
    io_data.ROOT.parents[2] / "code" / "sps-crab-cavities" / "studies" / "ccs_scans2" / "data"
    / "sps_crab_scan_matched.json"
)

SEQ_RE = re.compile(r"^\s*([A-Z][A-Z0-9]*\.[0-9A-Z.]+)\s*:\s*\S+\s*,\s*at\s*=\s*([0-9.eE+-]+)")
SEQ_LENGTH_RE = re.compile(r"SEQUENCE\s*,.*?\bL\s*=\s*([0-9.eE+-]+)")


def lattice_positions(seq_path: Path = SPS_SEQ) -> dict[str, float]:
    """element name -> `at` [m] for every element in sps.seq with an explicit `at`."""
    text = seq_path.read_text()
    if "refer = centre" not in text:
        raise ValueError(f"{seq_path}: expected 'refer = centre' -- s convention would differ from BLMs.")
    out = {}
    for line in text.splitlines():
        m = SEQ_RE.match(line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def lattice_length(seq_path: Path = SPS_SEQ) -> float:
    m = SEQ_LENGTH_RE.search(seq_path.read_text())
    if not m:
        raise ValueError(f"{seq_path}: no SEQUENCE length found.")
    return float(m.group(1))


def line_s_to_seq_s(line, seq_path: Path = SPS_SEQ):
    """Return f(s_line) -> s_seq, piecewise-linear through every element (centre) present in both
    the xsuite line and sps.seq, anchored at 0 and at both machine lengths."""
    tt = line.get_table()
    seq = lattice_positions(seq_path)
    s_line, s_seq = [0.0], [0.0]
    for i, name in enumerate(tt.name[:-1]):
        at = seq.get(str(name).upper())
        if at is None:
            continue
        s_line.append(0.5 * (tt.s[i] + tt.s[i + 1]))
        s_seq.append(at)
    s_line.append(line.get_length())
    s_seq.append(lattice_length(seq_path))
    s_line, s_seq = np.asarray(s_line), np.asarray(s_seq)
    order = np.argsort(s_line)
    s_line, s_seq = s_line[order], s_seq[order]
    if np.any(np.diff(s_seq) < -1.0):
        raise ValueError("line -> sps.seq s mapping is not monotonic -- element order differs.")
    return lambda s: np.interp(s, s_line, s_seq)


def compute_aperture(line_path: Path = APERTURE_LINE, seq_path: Path = SPS_SEQ) -> pd.DataFrame:
    """
    Aperture extent at every Limit* element of the xsuite line (apertls, as in
    sps-crab-cavities/studies/twiss_analysis.ipynb), on the sps.seq s frame.
    Returns s_m, x_min_mm, x_max_mm, y_min_mm, y_max_mm. Needs xtrack + apertls (~2 s).
    """
    import apertls
    import xtrack as xt

    line = xt.load(str(line_path))
    tt = line.get_table()
    ap_calc = apertls.ApertureCalculator(line)
    x_ext = ap_calc.compute_x_extent()
    y_ext = ap_calc.compute_y_extent()
    mask = np.array([line[nn].__class__.__name__.startswith("Limit") for nn in tt.name[:-1]])
    s_line = tt.s[:-1][mask]
    to_seq = line_s_to_seq_s(line, seq_path)
    return pd.DataFrame({
        "s_m": to_seq(s_line),
        "x_min_mm": x_ext[:, 0] * 1e3, "x_max_mm": x_ext[:, 1] * 1e3,
        "y_min_mm": y_ext[:, 0] * 1e3, "y_max_mm": y_ext[:, 1] * 1e3,
    })
