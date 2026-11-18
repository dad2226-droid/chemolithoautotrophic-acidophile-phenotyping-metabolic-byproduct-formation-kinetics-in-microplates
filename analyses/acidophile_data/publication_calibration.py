"""Calibration."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .metadata import load_metadata
from scipy.stats import t as t_distribution
from .optical import load_optical_registry, load_iron_references, final_well_concentration
from .plate_reader import read_plate


DENOMINATOR_GUARD_AU = 1e-6


CEILING_AU = 4.0


PH_SIGNALS = ("ratio_raw", "ratio_corrected")


IRON_SIGNALS = ("a380", "a380_corrected", "a428", "a428_corrected")


def _selected(row):
    ex = row["experiment_id"]
    a = json.loads(row["acquisition_metadata_json"])
    if ex.startswith("final_indicator_"):
        return a.get("read_index") == 0 or a.get("reading_label") == "spectrum"
    if ex in {"pH_UV_20uM_stock_initial_and_spectra", "pH_UV_200uM_stock_initial_spectra"}:
        return True
    if ex == "pH_configuration_23_C_UV" or (ex.startswith("iron_configuration_") and ex.endswith("_23_C_UV")):
        return a.get("reading_label") == "initial"
    return False


def _channel(plate, well, wavelength):
    if wavelength not in plate.wavelengths_nm:
        return np.nan, "missing_channel"
    if plate.is_parse_error(well, wavelength):
        return np.nan, "parse_error"
    if plate.is_censored_high(well, wavelength):
        return np.nan, "instrument_overrange"
    value = plate.value(well, wavelength)
    return value, "missing_value" if not np.isfinite(value) else "numeric_ceiling" if value > CEILING_AU else "eligible"


def _read_sources(root):
    registry = load_optical_registry(root)
    endpoints, spectra, coverage = [], [], []
    for entry in registry.to_dict("records"):
        if not _selected(entry):
            continue
        reference = json.loads(entry["reference_metadata_json"])[entry["reference_set_id"]]
        experiment = json.loads(entry["experiment_metadata_json"])
        acquisition = json.loads(entry["acquisition_metadata_json"])
        prep = experiment["preparation"]
        replication = experiment["replication"]
        plate = read_plate(root / entry["source_path"])
        common = {key: entry[key] for key in ("source_path", "physical_plate_id", "experiment_id", "reference_set_id")}
        common.update(run_number=acquisition["run_number"],
                      nominal_dye_uM=prep.get("nominal_dyed_well_thymol_blue_uM", 10),
                      dye_stock_uM=prep.get("thymol_blue_stock_uM", 200),
                      acid_M=prep.get("declared_h2so4_M", prep.get("declared_h2so4_mM", np.nan) / 1000),
                      independent_preparations=replication.get("independent_standard_preparations", 1),
                      series_ID=reference.get("axis_details", {}).get("series_ID", ""),
                      timestamp=plate.saved_at.isoformat() if plate.saved_at else "",
                      reference_pH_remeasured_for_scan=False,
                      physical_plate_identity_status=acquisition.get("plate_association_status", "declared_metadata"))
        if plate.mode == "Spectrum":
            for wavelength in (380, 428, 430, 485, 545, 750):
                present = wavelength in plate.wavelengths_nm
                coverage.append({**common, "wavelength_nm": wavelength, "measured": present,
                                 "status": "measured" if present else "outside_scan_range" if wavelength > max(plate.wavelengths_nm) else "not_on_native_grid",
                                 "scan_min_nm": min(plate.wavelengths_nm), "scan_max_nm": max(plate.wavelengths_nm), "interpolation_used": False})
        for group in reference["row_groups"]:
            concentrations = group["reported_concentration"].get("values_by_column")
            phs = group.get("reference_pH_by_column")
            unit = group["reported_concentration"].get("unit")
            for index, column in enumerate(group["plate_columns"]):
                ph = phs[index] if phs and phs[index] is not None else np.nan
                concentration = concentrations[index] if concentrations and concentrations[index] is not None else np.nan
                if np.isfinite(concentration):
                    concentration *= {"M": 1, "mM": .001, "uM": 1e-6}.get(unit, 1)
                    concentration = final_well_concentration(concentration, group["final_well_volume_uL"], reference_volume_uL=group["standard_volume_uL"])
                shared = {**common, "condition": group["condition"], "row_block": "".join(group["rows"]), "plate_column": column,
                          "reference_pH": ph, "reference_pH_status": group.get("reference_pH_status", ""),
                          "reference_concentration_M": concentration,
                          "reference_quantity": group["reported_concentration"].get("quantity", ""),
                          "final_well_volume_uL": group["final_well_volume_uL"],
                          "condition_id": f"{entry['reference_set_id']}:{''.join(group['rows'])}:{column}"}
                if plate.mode == "Spectrum":
                    for wavelength in plate.wavelengths_nm:
                        readings = [_channel(plate, f"{row}{column}", wavelength) for row in group["rows"]]
                        valid = [v for v, status in readings if status == "eligible"]
                        complete = len(valid) == len(group["rows"]) == 4
                        spectra.append({**shared, "wavelength_nm": wavelength, "mean_absorbance": np.mean(valid) if complete else np.nan,
                                        "available_mean_absorbance": np.mean(valid) if valid else np.nan,
                                        "sd_absorbance": np.std(valid, ddof=1) if complete else np.nan,
                                        "n_expected": 4, "n_eligible": len(valid), "complete": complete})
                else:
                    for row in group["rows"]:
                        well = f"{row}{column}"
                        record = {**shared, "well": well}
                        for wavelength in (380, 428, 430, 485, 545, 750):
                            record[f"a{wavelength}"], record[f"a{wavelength}_status"] = _channel(plate, well, wavelength)
                        endpoints.append(record)
    return pd.DataFrame(endpoints), pd.DataFrame(spectra), pd.DataFrame(coverage)


def signal_value(record, signal):
    """Return derived signal and explicit eligibility without altering raw values."""
    if signal.startswith("ratio"):
        wavelengths = (545, 485, 750) if signal.endswith("corrected") else (545, 485)
    else:
        wavelength = int(signal.split("_")[0][1:])
        wavelengths = (wavelength, 750) if signal.endswith("corrected") else (wavelength,)
    flags = [record[f"a{w}_status"] for w in wavelengths if record[f"a{w}_status"] != "eligible"]
    if flags:
        return np.nan, "|".join(sorted(set(flags)))
    if signal.startswith("ratio"):
        sub = record["a750"] if signal.endswith("corrected") else 0
        denominator = record["a485"] - sub
        if denominator <= DENOMINATOR_GUARD_AU:
            return np.nan, "nonpositive_or_near_zero_denominator"
        return (record["a545"] - sub) / denominator, "eligible"
    return record[f"a{wavelength}"] - (record["a750"] if signal.endswith("corrected") else 0), "eligible"


def _level_means(wells):
    records = []
    for (_, _, _), group in wells.groupby(["physical_plate_id", "condition_id", "condition"], sort=True):
        base = group.iloc[0].drop([c for c in group if c.startswith("a") and c != "acid_M"] + ["well"]).to_dict()
        is_ph = group.experiment_id.iloc[0].startswith(("pH", "final_indicator_acid"))
        signals = PH_SIGNALS if is_ph else IRON_SIGNALS + PH_SIGNALS
        for signal in signals:
            values = [signal_value(row, signal) for row in group.to_dict("records")]
            good = [v for v, status in values if status == "eligible"]
            records.append({**base, "signal": signal, "assay": "pH" if signal in PH_SIGNALS else "iron",
                            "axis_value": base["reference_pH"] if signal in PH_SIGNALS else base["reference_concentration_M"] * 1000,
                            "axis_unit": "pH" if signal in PH_SIGNALS else "mM_total_Fe",
                            "mean": np.mean(good) if len(good) == len(group) == 4 else np.nan,
                            "available_mean": np.mean(good) if good else np.nan,
                            "technical_sd": np.std(good, ddof=1) if len(good) > 1 else np.nan,
                            "n_expected": 4, "n_eligible": len(good), "n_observed": len(group),
                            "eligibility": "complete" if len(good) == len(group) == 4 else "incomplete_technical_wells",
                            "invalid_reasons": "|".join(sorted({s for _, s in values if s != "eligible"}))})
    return pd.DataFrame(records)


def _core(levels):
    ph = levels.experiment_id.isin(["pH_configuration_23_C_UV", "pH_UV_20uM_stock_initial_and_spectra"])
    iron = levels.experiment_id.str.startswith("iron_configuration_") & levels.signal.isin(("a380", "a380_corrected"))
    return levels[(ph | iron) & levels.condition.eq("thymol_blue") | (iron & levels.condition.eq("no_thymol_blue"))].copy()


def _chemical_references(root):
    raw = load_iron_references(root)
    descriptor = Path(root) / "data/01_indicator_calibration/02_iron_III_indicator/iron_III_indicator.yaml"
    provenance = load_metadata(descriptor).get("reference_concentration_methods", {})
    if (provenance.get("concentration_levels_individually_measured") is not True
            or not provenance.get("measurement_basis")):
        raise ValueError("Reference concentration levels must have explicit measurement provenance")
    if provenance.get("reported_total_Fe_records_per_reference_sample") != 1:
        raise ValueError("Reference total-Fe uncertainty handling requires one measured value per reference sample")
    rows = []
    for sample, group in raw.groupby("reference_sample_ID", sort=True):
        total = group[group.analyte.eq("total_Fe")].measured_concentration_M.iloc[0]
        ferrous = group[group.analyte.eq("Fe_II")].measured_concentration_M
        mean, sd, n = ferrous.mean(), ferrous.std(ddof=1), len(ferrous)
        half = t_distribution.ppf(.975, n - 1) * sd / np.sqrt(n) if n > 1 else np.nan
        rows.append({"reference_sample_ID": sample, "series_ID": group.series_ID.iloc[0], "plate_column": group.plate_column.iloc[0],
                     "nominal_H2SO4_concentration_M": group.nominal_H2SO4_concentration_M.iloc[0],
                     "total_Fe_M": total, "total_Fe_reference_record_n": int(group.analyte.eq("total_Fe").sum()),
                     "concentration_level_individually_measured": True,
                     "reference_provenance": provenance["measurement_basis"],
                     # The deposit supplies one total-Fe value per sample, without
                     # its analytical uncertainty. This is an analysis limitation,
                     # not an assumption that the measurement has zero uncertainty.
                     "total_Fe_uncertainty_scope": "conditional_on_single_deposited_total_Fe_value; analytical_uncertainty_not_supplied",
                     "total_Fe_uncertainty_propagated": False, "FeII_M": mean, "FeII_sd_M": sd, "FeII_n": n,
                     "FeIII_M": total - mean, "FeIII_conditional_sd_M": sd, "FeIII_conditional_ci_low_M": total - mean - half,
                     "FeIII_conditional_ci_high_M": total - mean + half, "FeIII_positive_all_titrations": bool((total - ferrous > 0).all()),
                     "FeIII_negative_mean": total - mean < 0, "reference_basis": "before_final_addition", "independent_standard_preparations": 1,
                     "uncertainty_scope": "FeII_repeatability_only_totalFe_treated_fixed"})
    return pd.DataFrame(rows)
