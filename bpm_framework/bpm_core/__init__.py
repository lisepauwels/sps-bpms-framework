"""
bpm_core -- reusable, MD-agnostic ALPS orbit pipeline logic, mirroring lossmap_framework's
lossmap_core. No dependency on config/; every function takes plain DataFrames/paths/dicts.

    naming          canonical case naming + cross-field validation (copied from lossmap_core)
    io_data         NXCALS query (orbit, injection anchor, BCT) / format / save / load
    data_treatment  per-cycle quality checks (beam present, ALPS crates complete)
    lattice         sps.seq s positions, apertls aperture mapped onto the sps.seq s frame
    plotting        orbit vs s coloured by time (pure functions, no table looping)

Usable standalone, e.g. from a notebook:
    from bpm_core import io_data
    orbit_h = io_data.load_bpms("20260401", "references", "-120_no_chirp", "20260401_160021", "H")
"""

from . import data_treatment, io_data, lattice, naming, plotting

__all__ = ["data_treatment", "io_data", "lattice", "naming", "plotting"]
