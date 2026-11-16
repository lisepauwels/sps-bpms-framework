"""
Canonical study/case naming -- the single identifier used both for a timestamp table's row groups
and for the data/<md_date>/<group>/<name>/ and figures/ directory names. Copied verbatim from
lossmap_framework/lossmap_core/naming.py, so case names match the BLM side exactly.

Also holds validate(), for checking a case's metadata against a per-MD config.VALIDITY_RULES --
lives here because it operates on the same "case metadata dict" as canonical_name().
"""


def canonical_name(study_meta: dict, keys: list[str] | None = None) -> str:
    """
    Joins the selected metadata fields with "_", in the given order.

    keys=None (default): every field in study_meta except "note", in the dict's own order --
    these fields ARE the characteristics of one case, so the name is built generically from
    whatever a given MD's config supplies (e.g. lm_type/plane/energy/collimator_gap/coast_or_cycle
    for 20260826), not a hardcoded subset of field names that would silently drop fields another
    MD's config defines differently.

    keys=[...]: only join these fields, in this order -- for when some fields are constant across
    an entire MD (e.g. crab-cavities' energy/coll_gap are always 26GeV/5sigma) and would just be
    redundant noise repeated in every name; pass e.g. ["crab_phase", "lm_type", "plane"] there.

    Raises if any selected field is still unconfirmed (None) -- such a study isn't identified well
    enough to name yet; the caller should be filtering these out before it ever needs a name.
    """
    if keys is None:
        keys = [k for k in study_meta if k != "note"]
    parts = []
    for key in keys:
        value = study_meta.get(key)
        if value is None:
            raise ValueError(
                f"canonical_name requires every selected field to be confirmed -- '{key}' is "
                "None for this study, which isn't identified well enough to name yet."
            )
        parts.append(str(value).replace(" ", ""))
    return "_".join(parts)


def validate(meta: dict, validity_rules: dict) -> None:
    """
    Checks meta's fields against validity_rules -- a per-MD config.VALIDITY_RULES dict of the
    shape {field: (depends_on_field, {value_of_depends_on: [allowed values of field]})}, e.g.
    {"plane": ("lm_type", {"betatron": ["H", "V"], "offmom": ["DPpos", "DPneg"]})} encodes
    "a betatron case's plane can only be H or V".

    Raises ValueError if a field's actual value isn't in the allowed set for its depends_on
    field's current value. Silently skips a rule if either field is missing from meta, or if
    depends_on's current value isn't a key in the mapping (nothing declared to check against for
    that value -- validity_rules should stay exhaustive over every value depends_on can actually
    take, per MD, but this function doesn't enforce that itself).
    """
    for field, (depends_on, allowed_by_value) in validity_rules.items():
        if field not in meta or depends_on not in meta:
            continue
        allowed = allowed_by_value.get(meta[depends_on])
        if allowed is None:
            continue
        if meta[field] not in allowed:
            raise ValueError(
                f"Invalid {field}={meta[field]!r} for {depends_on}={meta[depends_on]!r} "
                f"-- expected one of {allowed}."
            )
