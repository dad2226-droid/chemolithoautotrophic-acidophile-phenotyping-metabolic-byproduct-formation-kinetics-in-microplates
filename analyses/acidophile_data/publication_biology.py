"""Biology."""
from __future__ import annotations
import hashlib
import json
from .plate_comparisons import CONDITION_KEYS


CLASSES = ("Metal sulfates", "Other salts", "Sulfur compounds", "Oxyanions", "Organic additions")


def control_plate_means(curves, conditions):
    meta = conditions[CONDITION_KEYS + ['no_addition_control','control_cohort','control_cohort_definition']]
    f = curves.merge(meta, on=CONDITION_KEYS, validate='many_to_one')
    f = f[f.no_addition_control].drop_duplicates(
        ['control_cohort','signal','batch_id','plate_id','column','ordinal'])
    out = f.groupby(['control_cohort','signal','batch_id','plate_id','ordinal'], as_index=False).agg(
        days_mean=('days_mean','mean'), days_min=('days_min','min'), days_max=('days_max','max'),
        biotic_change=('biotic_change','mean'), abiotic_change=('abiotic_change','mean'),
        net_change=('net_change','mean'), control_column_n=('column','nunique'),
        culture_ids=('culture_ids',lambda x:';'.join(sorted(set(';'.join(x).split(';'))))))
    return out.merge(meta[['control_cohort','control_cohort_definition']].drop_duplicates(),
                     on='control_cohort', validate='many_to_one')


def _classify(identifier: str) -> str:
    name = identifier.casefold()
    if any(word in name for word in ("citric", "citrate", "gluconate", "sucrose", "formate")):
        return CLASSES[4]
    if any(word in name for word in ("tetrathionate", "thiosulfate", "sulfite", "thiocyanate")):
        return CLASSES[2]
    if any(word in name for word in ("phosphate", "nitrate", "iodate", "chromate", "perrhenate", "molybdate", "tungstate", "selenite", "selenate")):
        return CLASSES[3]
    if "sulfate" in name and not any(word in name for word in ("sodium_", "potassium_", "ammonium_")):
        return CLASSES[0]
    return CLASSES[1]


def _condition_fields(row) -> dict:
    values = json.loads(row.component_final_mM_json)
    principal = "citric_acid" if "citric_acid" in values else next(iter(values))
    name = row.experiment_id.removeprefix("AF_AFM1_").replace("_OD0p001", "")
    # Use the measured reagent identity for the display label.
    if name == "sodium_selenate" and "sodium_selenite" in values:
        name = "sodium_selenite"
    label = name.replace("_", " ").replace("iii", "III").replace("ii", "II").replace("iv", "IV")
    if "citric_acid" in values and "gallium" in values:
        label = f"citrate + {values['gallium']:g} mM Ga"
    zero = all(float(value) == 0 for value in values.values())
    cohort_values = {key: str(row.get(key, "")) for key in (
        "medium_id", "preparation_id", "initial_OD600_final", "final_well_volume_uL",
        "initial_pH_status", "initial_pH_measured", "nominal_initial_pH_json", "incubation_temperature_C")}
    cohort = hashlib.sha256(json.dumps(cohort_values, sort_keys=True).encode()).hexdigest()[:12]
    return dict(additive_class=_classify(row.experiment_id), additive_label=label,
                principal_component=principal, principal_final_mM=float(values[principal]),
                no_addition_control=zero, control_cohort=cohort,
                control_cohort_definition=json.dumps(cohort_values, sort_keys=True))
