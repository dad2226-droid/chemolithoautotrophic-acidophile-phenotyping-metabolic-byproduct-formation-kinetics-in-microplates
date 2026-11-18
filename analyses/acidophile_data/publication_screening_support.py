"""Screening support."""
from __future__ import annotations
import numpy as np
import pandas as pd


CONTEXT = ["species", "experiment_id", "group_id", "column"]


PHYSICAL = ["species", "batch_id", "plate_id", "column"]


SCREEN_KEYS = ["species", "stratum"]


OUTCOME_CATEGORIES = ("all_positive", "all_negative", "identical_mixed",
                      "discordant", "unevaluable")


def _finite(values):
    values = np.asarray(values, float)
    return values[np.isfinite(values)]


def _quantile(values, quantile):
    values = _finite(values)
    return float(np.quantile(values, quantile)) if len(values) else np.nan


def _percent(numerator, denominator):
    return 100. * numerator / denominator if denominator else np.nan


def _augment(original, additions, keys):
    overlap = set(original).intersection(additions) - set(keys)
    if overlap:
        raise ValueError("Screen support would replace existing fields: " + ", ".join(sorted(overlap)))
    result = original.merge(additions, on=keys, how="left", validate="one_to_one", sort=False)
    if len(result) != len(original):
        raise ValueError("Screen support changed the result population")
    return result


def _acquisition_margins(tables, metadata, units):
    """Use the producer's physical identities and exact duplicate rule."""
    source = pd.concat([tables["af_iron_plate_trajectories"],
                        tables["validation_lf_plate_trajectories"]], ignore_index=True)
    source = source[source.signal.eq("raw_A380")]
    fields = CONTEXT + ["batch_id", "plate_id", "ordinal", "days_mean",
                        "biotic_change", "abiotic_change", "net_change", "culture_ids"]
    # Full chemical definitions remain in the original unit-level support;
    # repeating those JSON objects at every acquisition obscures the gate data.
    identities = metadata[CONTEXT + ["condition_id", "stratum", "medium_id"]]
    source = source[fields].merge(identities, on=CONTEXT, validate="many_to_one")
    keys = PHYSICAL + ["ordinal"]
    for key, group in source[source.duplicated(keys, keep=False)].groupby(keys, sort=False):
        for field in ("days_mean", "biotic_change", "abiotic_change", "net_change"):
            values = group[field].to_numpy(float)
            if not np.allclose(values, values[0], rtol=0, atol=0, equal_nan=True):
                raise ValueError(f"Conflicting screening margin observation: {key}, {field}")
        if group.condition_id.nunique() != 1 or group.stratum.nunique() != 1:
            raise ValueError(f"Conflicting screening margin identity: {key}")
    source = source.drop_duplicates(keys).copy()
    unit_fields = PHYSICAL + ["decision", "decision_status", "complete_recorded_window"]
    source = source.merge(units[unit_fields], on=PHYSICAL, validate="many_to_one")
    source["is_postbaseline"] = source.days_mean.gt(0)
    finite_biotic = np.isfinite(source.biotic_change)
    finite_abiotic = np.isfinite(source.abiotic_change)
    source["biotic_margin_AU"] = np.where(finite_biotic, source.biotic_change - 1., np.nan)
    source["abiotic_margin_AU"] = np.where(finite_abiotic, 1. - source.abiotic_change, np.nan)
    source["matched_pair_eligible"] = finite_biotic & finite_abiotic
    source["biotic_gate_pass"] = finite_biotic & source.biotic_margin_AU.ge(0)
    source["abiotic_gate_pass"] = finite_abiotic & source.abiotic_margin_AU.gt(0)
    source["operational_event_at_acquisition"] = (source.is_postbaseline &
        source.matched_pair_eligible & source.biotic_gate_pass & source.abiotic_gate_pass)
    source["joint_gate_margin_AU"] = np.where(source.matched_pair_eligible,
        np.minimum(source.biotic_margin_AU, source.abiotic_margin_AU), np.nan)
    source["pair_status"] = np.select(
        [~finite_biotic & ~finite_abiotic, ~finite_biotic, ~finite_abiotic],
        ["both_signals_unavailable", "biotic_signal_unavailable", "abiotic_signal_unavailable"],
        default="matched_finite_pair")
    source["biotic_gate_rule"] = "biotic_margin_AU >= 0"
    source["abiotic_gate_rule"] = "abiotic_margin_AU > 0"
    return source.sort_values(keys).reset_index(drop=True)


def _selected_margins(acquisitions, units):
    rows = []
    groups = {key: group for key, group in acquisitions.groupby(PHYSICAL, sort=False)}
    for unit in units.to_dict("records"):
        key = tuple(unit[field] for field in PHYSICAL)
        group = groups[key]
        post = group[group.is_postbaseline]
        passing = post[post.operational_event_at_acquisition].sort_values("ordinal")
        observed = bool(len(passing))
        if observed != bool(unit["event"]):
            raise ValueError(f"Screen margin events disagree with original decisions: {key}")
        selected = None
        if unit["decision"] == 1:
            selected = passing.iloc[0]
            selection = "first_recorded_event"
        elif unit["decision"] == 0:
            if not unit["complete_recorded_window"] or not len(post) or not post.matched_pair_eligible.all():
                raise ValueError(f"Incomplete trajectory classified as no-event: {key}")
            selected = post.sort_values(["joint_gate_margin_AU", "ordinal"],
                                         ascending=[False, True]).iloc[0]
            selection = "nearest_joint_gate_approach"
        else:
            selection = "unavailable_unevaluable_decision"
        identity = {field: unit[field] for field in PHYSICAL +
                    ["condition_id", "stratum", "decision", "decision_status", "complete_recorded_window"]}
        row = identity | dict(selection=selection, selected_acquisition_n=int(selected is not None),
            selection_rule=("earliest_passing_postbaseline_ordinal" if selection == "first_recorded_event"
                            else "max_min_of_paired_gate_margins_then_earliest_ordinal" if selected is not None
                            else "no_selection_for_unevaluable_decision"))
        for field in ("ordinal", "days_mean", "biotic_change", "abiotic_change", "biotic_margin_AU",
                      "abiotic_margin_AU", "joint_gate_margin_AU", "biotic_gate_pass", "abiotic_gate_pass",
                      "operational_event_at_acquisition"):
            row[field] = selected[field] if selected is not None else np.nan
        if selected is not None and selection == "first_recorded_event":
            if not np.isclose(selected.days_mean, unit["first_event_time_d"], rtol=0, atol=0,
                              equal_nan=True) or not np.isclose(
                    selected.biotic_change-selected.abiotic_change, unit["first_event_net_AU"], rtol=0, atol=0):
                raise ValueError(f"Selected first event differs from the screening classification: {key}")
        rows.append(row)
    return pd.DataFrame(rows)


def add_decision_support(tables, metadata, units, summary):
    """Return augmented prevalence/margin summary and two detailed support frames.

    A no-event trajectory's nearest joint gate approach maximizes the smaller
    signed gate margin across matched postbaseline acquisitions, using earliest
    ordinal to break ties. Both reported margins belong to that one acquisition.
    A zero minimum margin can still fail the strict abiotic gate; eligibility
    and gate booleans, rather than the minimum alone, define the original event.
    Unevaluable trajectories have no selected nearest approach.
    """
    acquisitions = _acquisition_margins(tables, metadata, units)
    selected = _selected_margins(acquisitions, units)
    rows = []
    for key, group in selected.groupby(SCREEN_KEYS, sort=True):
        row = dict(zip(SCREEN_KEYS, key))
        counts = [int(group.decision.eq(1).sum()), int(group.decision.eq(0).sum()),
                  int(group.decision.isna().sum())]
        if sum(counts) != len(group):
            raise ValueError("Unrecognized operational decision")
        original = summary
        for field, value in row.items():
            original = original[original[field].eq(value)]
        if len(original) != 1:
            raise ValueError("Missing or duplicate decision summary")
        expected = original.iloc[0]
        if counts != [expected.event_n, expected.no_event_n, expected.unevaluable_decision_n] or len(group) != expected.candidate_plate_condition_n:
            raise ValueError("Screen prevalence support does not match the decision summary")
        row.update(complete_no_event_n=counts[1], event_percent_all_candidates=_percent(counts[0], len(group)),
                   complete_no_event_percent_evaluable=_percent(counts[1], counts[0]+counts[1]),
                   unevaluable_percent_all_candidates=_percent(counts[2], len(group)))
        for selection, prefix in (("first_recorded_event", "first_event"),
                                  ("nearest_joint_gate_approach", "nearest_joint_gate_approach")):
            subset = group[group.selection.eq(selection)]
            row[prefix + "_margin_unit_n"] = len(subset)
            for margin in ("biotic", "abiotic", "joint_gate"):
                for suffix, quantile in (("min", 0), ("q25", .25), ("median", .5), ("q75", .75), ("max", 1)):
                    row[f"{prefix}_{margin}_margin_{suffix}_AU"] = _quantile(subset[margin + "_margin_AU"], quantile)
        rows.append(row)
    augmented = _augment(summary, pd.DataFrame(rows), SCREEN_KEYS)
    return augmented, {"screen_acquisition_gate_margins": acquisitions,
                       "screen_selected_gate_margins": selected}


def add_agreement_support(summary, records):
    """Decompose identical event fractions without changing their denominator."""
    details = records.copy()
    categories = []
    for row in details.itertuples():
        if not row.complete:
            category = "unevaluable"
        elif row.fraction_range != 0:
            category = "discordant"
        elif row.event_fraction == 1:
            category = "all_positive"
        elif row.event_fraction == 0:
            category = "all_negative"
        elif 0 < row.event_fraction < 1 and row.level == "between batches":
            category = "identical_mixed"
        else:
            raise ValueError("Invalid complete repeated-screen outcome")
        categories.append(category)
    details["outcome_category"] = categories
    keys = SCREEN_KEYS + ["comparison_level"]
    rows = []
    for original in summary.to_dict("records"):
        group = details[details.species.eq(original["species"]) & details.stratum.eq(original["stratum"]) &
                        details.level.eq(original["comparison_level"])]
        row = {field: original[field] for field in keys}
        for category in OUTCOME_CATEGORIES:
            row[category + "_comparison_n"] = int(group.outcome_category.eq(category).sum())
        identical = sum(row[category + "_comparison_n"] for category in OUTCOME_CATEGORIES[:3])
        if (sum(row[category + "_comparison_n"] for category in OUTCOME_CATEGORIES) != original["candidate_comparison_n"]
                or identical != original["identical_event_fraction_n"]
                or row["unevaluable_comparison_n"] != original["incomplete_comparison_n"]
                or identical + row["discordant_comparison_n"] != original["evaluable_comparison_n"]):
            raise ValueError("Repeated-screen decomposition does not match the agreement summary")
        for category in OUTCOME_CATEGORIES[:-1]:
            row[category + "_percent_evaluable"] = _percent(row[category + "_comparison_n"], original["evaluable_comparison_n"])
        rows.append(row)
    return _augment(summary, pd.DataFrame(rows), keys), {"screen_agreement_outcome_mix": details}
