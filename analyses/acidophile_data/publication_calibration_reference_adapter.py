"""Calibration reference adapter."""
from __future__ import annotations
import pandas as pd


def canonical_ph_adapter(wells, signal_value):
    selected = wells[(wells.experiment_id.eq("pH_configuration_23_C_UV") | wells.experiment_id.str.startswith("iron_configuration_")) & wells.condition.eq("thymol_blue")]
    rows = []
    for row in selected.to_dict("records"):
        is_acid = row["experiment_id"] == "pH_configuration_23_C_UV"
        for signal in ("ratio_raw", "ratio_corrected"):
            response, status = signal_value(row, signal)
            rows.append({"dataset": "acid_calibration" if is_acid else "ferric_matrix",
                         "plate_id": row["physical_plate_id"], "layout_group": row["row_block"], "column": row["plate_column"],
                         "well": row["well"], "acid_m": 0.0 if is_acid else row["acid_M"],
                         "reference_ph_after_dye": row["reference_pH"], "reference_pH_measurement_stage": row["reference_pH_status"],
                         "source_file": row["source_path"], "signal": signal, "response": response,
                         "eligibility": status, "dye_nominal_uM": row["nominal_dye_uM"]})
    return pd.DataFrame(rows)
