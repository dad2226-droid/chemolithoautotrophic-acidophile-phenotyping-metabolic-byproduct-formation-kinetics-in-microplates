"""Condition repeatability."""
from __future__ import annotations
import hashlib
import json
import numpy as np
import pandas as pd
from .plate_comparisons import CONDITION_KEYS


METRICS = (("AUC_to_own_max", "auc_to_max"),
           ("time_to_observed_max", "time_to_max_d"))


DESIGN_FIELDS = ("medium_id", "preparation_id", "readout", "initial_OD600_final",
                 "final_well_volume_uL", "nominal_initial_pH_json",
                 "initial_pH_status", "initial_pH_measured", "initial_pH_scope",
                 "incubation_temperature_C", "additive_stock_pH")


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _scalar(value):
    if value is None or (not isinstance(value, (dict, list)) and pd.isna(value)):
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, str) and value.strip() in ("null", "nan", ""):
        return None
    return value


def _components(value):
    parsed = json.loads(value) if isinstance(value, str) else value
    return {str(k): float(v) for k, v in (parsed or {}).items() if float(v) != 0}


def _joined(values):
    return ";".join(sorted({part for value in values if pd.notna(value)
                            for part in str(value).split(";") if part}))


def _identity(fmt, definition):
    encoded = _json(definition)
    return fmt + "_" + hashlib.sha256(encoded.encode()).hexdigest()[:16], encoded


def _label(components, fmt):
    if not components:
        return "0 mM Ga" if fmt == "flask" else "no additive"
    return "; ".join(f"{value:g} mM {name.replace('_', ' ')}"
                     for name, value in sorted(components.items()))


def _plate_records(tables):
    metrics = tables.get("validation_plate_peak_metrics", pd.DataFrame())
    if metrics.empty:
        return []
    coverage = tables["af_iron_condition_coverage"]
    metadata = {tuple(r[c] for c in CONDITION_KEYS): r
                for r in coverage.to_dict("records")}
    if len(metadata) != len(coverage):
        raise ValueError("Condition repeatability requires unique plate condition metadata")
    records = []
    for number, r in enumerate(metrics.to_dict("records"), 1):
        context = tuple(r[c] for c in CONDITION_KEYS)
        design = metadata[context]
        components = _components(design["component_final_mM_json"])
        bases = _components_basis(design.get("component_basis_json", "{}"), components)
        definition = {field: _scalar(design.get(field)) for field in DESIGN_FIELDS}
        # Live producer frames carry this as JSON text; CSV reloads can infer
        # numbers. Both representations must identify the same design.
        nominal_ph = definition["nominal_initial_pH_json"]
        if isinstance(nominal_ph, str):
            definition["nominal_initial_pH_json"] = json.loads(nominal_ph)
        definition.update(species=r["species"], signal=r["signal"],
                          components_mM=components, concentration_bases=bases)
        condition, encoded = _identity("plate", definition)
        for metric, field in METRICS:
            value = float(r[field]) if pd.notna(r[field]) else np.nan
            if metric == "time_to_observed_max" and not r["positive_rise"]:
                value = np.nan
            records.append(dict(format="plate", metric=metric, condition_id=condition,
                condition_label=_label(components, "plate"), condition_definition_json=encoded,
                components_mM_json=_json(components), no_addition=not components,
                Ga_mM=components.get("gallium", np.nan) if bases.get("gallium") == "Ga"
                      else 0. if not components else np.nan,
                batch_id=r["batch_id"], unit_id=r["plate_id"],
                # Equal chemistry in different columns remains technical support
                # inside the same physical plate, never additional plate units.
                physical_location=str(r["column"]), raw_value=value,
                native_unit="AU d" if metric == "AUC_to_own_max" else "d",
                culture_ids=r.get("culture_ids", ""),
                source_context_ids="/".join(map(str, context)),
                source_record_ids=f"validation_plate_peak_metrics:{number}",
                source_status=str(r.get("status" if metric == "AUC_to_own_max"
                                        else "maximum_status", "unspecified")),
                final_read_maximum=bool(r["maximum_at_final_read"])))
    return records


def _components_basis(value, components):
    bases = json.loads(value) if isinstance(value, str) else (value or {})
    return {k: bases.get(k, "unspecified") for k in components}


def _flask_records(tables):
    metrics = tables.get("validation_flask_peak_metrics", pd.DataFrame())
    if metrics.empty:
        return []
    points = tables["validation_flask_points"]
    records = []
    metadata = {}
    for flask, g in points.groupby("flask_ID", sort=False):
        fields = [c for c in ("kind", "organism", "medium", "incubation_temperature_C",
                              "inoculum_culture_id") if c in g]
        unique = g[fields].drop_duplicates()
        if len(unique) != 1:
            raise ValueError(f"Conflicting fixed flask metadata for {flask}")
        metadata[flask] = unique.iloc[0].to_dict()
    for number, r in enumerate(metrics.to_dict("records"), 1):
        design = metadata[r["flask_ID"]]
        if design["kind"] != "biotic":
            continue
        components = _components(r["quantities_json"])
        # Actual initial pH/Fe readings vary between batches and belong to the
        # repeatability being described, not a new nominal formulation identity.
        definition = {"species": _scalar(design.get("organism")),
                      "medium_id": _scalar(design.get("medium")),
                      "incubation_temperature_C": _scalar(design.get("incubation_temperature_C")),
                      "signal": r["signal"], "components_mM": components}
        condition, encoded = _identity("flask", definition)
        for metric, field in METRICS:
            value = float(r[field]) if pd.notna(r[field]) else np.nan
            if metric == "time_to_observed_max" and not r["positive_rise"]:
                value = np.nan
            records.append(dict(format="flask", metric=metric, condition_id=condition,
                condition_label=_label(components, "flask"), condition_definition_json=encoded,
                components_mM_json=_json(components), no_addition=not components,
                Ga_mM=components.get("Ga", 0.) if set(components).issubset({"Ga"}) else np.nan,
                batch_id=r["batch"], unit_id=r["flask_ID"], physical_location=r["flask_ID"],
                raw_value=value, native_unit="mM d" if metric == "AUC_to_own_max" else "d",
                culture_ids=design.get("inoculum_culture_id", ""),
                source_context_ids=r["condition_ID"],
                source_record_ids=f"validation_flask_peak_metrics:{number}",
                source_status=str(r.get("status" if metric == "AUC_to_own_max"
                                        else "maximum_status", "unspecified")),
                final_read_maximum=bool(r["maximum_at_final_read"])))
    return records


CONDITION_COLUMNS = ["format", "metric", "condition_id", "condition_label",
                     "condition_definition_json", "components_mM_json", "no_addition",
                     "Ga_mM", "native_unit"]


def summarize_records(records):
    """Deduplicate physical locations, then average physical units and batches."""
    source = pd.DataFrame(records)
    if source.empty:
        return {"condition_repeatability_" + name: pd.DataFrame()
                for name in ("units", "batches", "summary")}
    identity = ["format", "condition_id", "metric"]
    location_keys = identity + ["batch_id", "unit_id", "physical_location"]
    locations = []
    for _, g in source.groupby(location_keys, sort=True, dropna=False):
        values = g.raw_value.to_numpy(float)
        if not np.allclose(values, values[0], rtol=0, atol=0, equal_nan=True):
            raise ValueError("Conflicting metrics for the same physical location and condition")
        r = g.iloc[0].to_dict()
        for field in ("culture_ids", "source_context_ids", "source_record_ids", "source_status"):
            r[field] = _joined(g[field])
        r["source_membership_n"] = len(g)
        locations.append(r)
    locations = pd.DataFrame(locations)
    units = []
    for _, g in locations.groupby(identity + ["batch_id", "unit_id"], sort=True):
        r = {c: g.iloc[0][c] for c in CONDITION_COLUMNS + ["batch_id", "unit_id"]}
        r.update(raw_value=g.raw_value.mean(), physical_location_n=len(g),
                 eligible_location_n=int(g.raw_value.notna().sum()),
                 source_membership_n=int(g.source_membership_n.sum()),
                 final_read_maximum=bool(g.final_read_maximum.any()))
        for field in ("culture_ids", "source_context_ids", "source_record_ids", "source_status"):
            r[field] = _joined(g[field])
        r["physical_locations"] = _joined(g.physical_location)
        r["raw_status"] = "available" if np.isfinite(r["raw_value"]) else "unavailable_native_metric"
        units.append(r)
    units = pd.DataFrame(units)
    batches = []
    for _, g in units.groupby(identity + ["batch_id"], sort=True):
        r = {c: g.iloc[0][c] for c in CONDITION_COLUMNS + ["batch_id"]}
        r.update(raw_value=g.raw_value.mean(), physical_unit_n=len(g),
                 eligible_unit_n=int(g.raw_value.notna().sum()),
                 final_read_maximum_n=int(g.final_read_maximum.sum()))
        r["unit_ids"] = _joined(g.unit_id)
        for field in ("culture_ids", "source_context_ids", "source_record_ids", "source_status"):
            r[field] = _joined(g[field])
        batches.append(r)
    batches = pd.DataFrame(batches)
    summaries = []
    for _, g in batches.groupby(identity, sort=True):
        values = g.raw_value.to_numpy(float)
        finite = values[np.isfinite(values)]
        mean = float(np.mean(finite)) if len(finite) else np.nan
        sd = float(np.std(finite, ddof=1)) if len(finite) > 1 else np.nan
        status = ("fewer_than_two_usable_batches" if len(finite) < 2 else
                  "nonpositive_condition_mean" if mean <= 0 else "estimated")
        r = {c: g.iloc[0][c] for c in CONDITION_COLUMNS}
        r.update(condition_mean=mean, raw_sd=sd, candidate_batch_n=len(g),
                 eligible_batch_n=len(finite), physical_unit_n=int(g.physical_unit_n.sum()),
                 eligible_unit_n=int(g.eligible_unit_n.sum()),
                 negative_batch_value_n=int(np.sum(finite < 0)),
                 final_read_maximum_n=int(g.final_read_maximum_n.sum()),
                 sd_percent=100 * sd / mean if status == "estimated" else np.nan,
                 eligible=status == "estimated", status=status)
        r["batch_ids"] = _joined(g.batch_id)
        for field in ("unit_ids", "culture_ids", "source_context_ids", "source_record_ids", "source_status"):
            r[field] = _joined(g[field])
        summaries.append(r)
    summary = pd.DataFrame(summaries)
    support = identity + ["condition_mean", "candidate_batch_n", "eligible_batch_n", "eligible", "status"]
    for frame in (units, batches):
        frame_support = frame.merge(summary[support], on=identity, validate="many_to_one")
        usable = frame_support.eligible & np.isfinite(frame_support.raw_value)
        frame_support["deviation_percent"] = np.nan
        frame_support.loc[usable, "deviation_percent"] = 100 * (
            frame_support.loc[usable, "raw_value"] / frame_support.loc[usable, "condition_mean"] - 1)
        frame_support.loc[~np.isfinite(frame_support.raw_value), "status"] = "unavailable_native_metric"
        frame_support["eligible"] = usable
        if frame is units:
            units_result = frame_support
        else:
            batches_result = frame_support
    return dict(condition_repeatability_units=units_result,
                condition_repeatability_batches=batches_result,
                condition_repeatability_summary=summary)


def build_tables(tables):
    return summarize_records(_plate_records(tables) + _flask_records(tables))
