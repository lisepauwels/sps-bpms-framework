"""
Ground-truth constants for the 2026-04-01 SPS crab-cavity phase-scan MD (BPM side).

Every cycle's metadata comes from ../../sps-measurements/lossmaps_crab_cavities/nxcalsdata/
nxcals_index.csv (same source as lossmap_framework's config/md_20260401.py, whose constants are
copied here so case names/timestamps match the BLM side exactly). Constant across the whole
campaign (asserted in scripts/build_timestamp_table.py): selector, energy, TCSM gap.

Two groups, one timestamp table each (data/20260401/<group>/timestamp_table.csv):
    chirps      plane H/V rows -- exactly lossmap_framework's 20260401 table (82 rows, same
                names, same EXCLUDE_CYCLES).
    references  plane == "no_chirp" rows -- no excitation, the clean stable-orbit cycles.

Same 3-part schema as lossmap_framework (PARAMETERS / NAME_KEYS / VALIDITY_RULES), except
NAME_KEYS/VALIDITY_RULES/EXCLUDE_CYCLES are per group.
"""

from pathlib import Path

MD_DATE = "20260401"
SELECTOR = "SPS.USER.LHCMD1"
ENERGY = "26GeV"
TCSM_GAP = "5sigma"
CYCLE_LENGTH_S = 26.72

NXCALS_INDEX_CSV = (
    Path(__file__).resolve().parents[3]
    / "sps-measurements" / "lossmaps_crab_cavities" / "nxcalsdata" / "nxcals_index.csv"
)

# 6-minute window inside the scan for reading BPMALPS_N:Orbit:channelNames
# (scripts/build_bpm_map.py) -- channelNames is republished every cycle, longer windows blow the
# gRPC message limit.
CHANNEL_MAP_WINDOW = ("2026-04-01T14:00:00Z", "2026-04-01T14:06:00Z")

PARAMETERS = [
    "timestamp", "selector", "cycle_length_s",
    "crab_phase", "lm_type", "plane", "coll_gap", "energy",
]

GROUPS = {
    "chirps": {
        "planes": ["H", "V"],
        "name_keys": ["crab_phase", "lm_type", "plane"],
        "validity_rules": {"plane": ("lm_type", {"betatron": ["H", "V"]})},
        # copied from lossmap_framework/config/md_20260401.py -- see the reasons there
        "exclude_cycles": [
            ("LossMapRef_H", "id1"),
            ("LossMapRef_H", "id4"),
            ("LossMapRefV_V", "id1"),
            ("Phase_-90_V", "id1"),
        ],
    },
    "references": {
        "planes": ["no_chirp"],
        "name_keys": ["crab_phase", "plane"],
        "validity_rules": {"plane": ("lm_type", {"betatron": ["no_chirp"]})},
        # Same physical cycle as a V-chirp row (identical cycle_stamp_utc in nxcals_index.csv,
        # checked 2026-11-13); lossmap_framework treats them as V chirps:
        #   Phase_-120_no_chirp id1 == Phase_-150_V id3  (14:00:21.735 UTC)
        #   Phase_-60_no_chirp  id1 == Phase_-120_V id3  (14:05:13.335 UTC)
        "exclude_cycles": [
            ("Phase_-120_no_chirp", "id1"),
            ("Phase_-60_no_chirp", "id1"),
        ],
    },
}

# (name, UTC timestamp) cycles kept in the table and in per-cycle plots but left out of the
# cross-case overview (scripts/plot_overview.py). Checked 2026-11-13:
#   - 0_no_chirp 14:16:01.335: orbit steps by ~2.2 mm RMS in H and ~0.6 mm in V between 3 and 5 s
#     into the cycle and stays there (other three 0_no_chirp cycles normal). V moves too, so not a
#     crab effect -- cause unknown.
OVERVIEW_EXCLUDE = [
    ("0_no_chirp", "2026-04-01 14:16:01.335000+00:00"),
]

# Crab phase whose mean orbit is subtracted in the overview's delta-orbit plots
# (scripts/plot_overview.py --ref-phase overrides).
ORBIT_REFERENCE_PHASE = 0
