"""Checks for schedule missingness, independent batch ranks, and citrate pairing."""

import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biological_sensitivity import _citrate_endpoint_effects, _decision, _effect_results, _scenario
from acidophile_data.peak_response import observed_peak


class BiologicalSensitivityTests(unittest.TestCase):
    def test_thinning_preserves_missing_earlier_read(self):
        curve = pd.DataFrame(dict(
            days_mean=[0., 1., 2., 3., 4.], ordinal=list(range(5)),
            biotic_change=[0., np.nan, .7, .8, 1.1], abiotic_change=[0.] * 5,
        ))
        alternate = _scenario(curve, "every second acquisition")
        self.assertEqual(alternate.ordinal.tolist(), [0, 1, 2, 4])
        self.assertEqual(_decision(alternate), "qualifying")
        self.assertTrue(np.isnan(observed_peak(alternate)["auc_to_max"]))
        truncated = _scenario(curve, "remove final acquisition")
        self.assertEqual(_decision(truncated), "unevaluable")

    def test_ranks_and_direction_use_batch_means(self):
        units = []
        # Two discordant replicate units for concentration A still contribute
        # only their common batch mean when ranking the three concentrations.
        for batch, values in [("early", [[80., 120.], [110.], [125.]]),
                              ("late", [[125.], [110.], [100.]])]:
            for condition, measurements in zip(("a", "b", "c"), values):
                for i, value in enumerate(measurements):
                    units.append(dict(
                        format="plate", species="AF", stratum="salt", condition_id=condition,
                        series_id="series", batch_id=batch, unit_id=f"{batch}-{condition}-{i}",
                        response_percent=value, delay_d=np.nan, reference_ids="control",
                        final_read_maximum=False,
                    ))
        summary, batches, conditions, pairs = _effect_results(pd.DataFrame(units))
        auc = summary.loc[summary.metric.eq("relative AUC")].iloc[0]
        self.assertEqual(auc.median_spearman_rho, -1.)
        self.assertEqual(auc.direction_agreement_n, 1)
        self.assertEqual(len(batches), 6)
        self.assertEqual(pairs.loc[pairs.metric.eq("relative AUC"), "common_condition_n"].item(), 3)
        self.assertTrue(conditions.loc[conditions.metric.eq("peak-time delay"), "between_batch_sd"].isna().all())

    def test_plate_only_citrate_dose_requires_its_own_terminal_reference(self):
        experiment = "AF_AFM1_gallium_sulfate_OD0p001_citrate_ratios"
        condition_rows, plate_rows = [], []
        for column, dose in [(2, 0.), (4, 3.)]:
            metadata = dict(species="AF", experiment_id=experiment, group_id="g", column=column)
            condition_rows.append(metadata | {"component_final_mM_json": json.dumps({"citric_acid": dose})})
            for ordinal, time in [(0, 0.), (1, 1.), (2, 2.)]:
                # The terminal reference is missing. Its finite preceding event
                # cannot replace it in either the net response or abiotic signal.
                value = np.nan if column == 2 and ordinal == 2 else float(ordinal)
                plate_rows.append(metadata | dict(batch_id="b", plate_id="p", signal="raw_A380",
                    ordinal=ordinal, days_mean=time, net_change=value, abiotic_change=value))
        flask_rows = []
        for dose in [0., 1.]:
            flask_rows.append(dict(
                kind="biotic", organism="AF", medium="AFM1", incubation_temperature_C=30,
                batch="b", nominal_citrate_mM=dose, nominal_Ga_mM=0., nominal_Cu_mM=0.,
                nominal_Al_mM=0., nominal_V_mM=0., nominal_NaCl_mM=0., nominal_In_mM=0.,
                flask_ID=f"f-{dose}", inoculum_culture_id="starter", sampling_event=1,
                days_mean=7., net_change=dose,
            ))
        results, pairs = _citrate_endpoint_effects({
            "af_iron_condition_coverage": pd.DataFrame(condition_rows),
            "af_iron_plate_trajectories": pd.DataFrame(plate_rows),
            "validation_flask_points": pd.DataFrame(flask_rows),
        })
        plate = pairs[pairs.readout.str.startswith("plate")]
        self.assertEqual(set(plate.citric_acid_mM), {3.})
        self.assertTrue(plate.paired_effect.isna().all())
        self.assertEqual(set(plate.pair_status), {"nonfinite_selected_endpoint"})
        self.assertEqual(results.loc[results.readout.str.startswith("flask") &
            results.Ga_mM.eq(0) & results.citric_acid_mM.eq(1), "mean_effect"].item(), 1.)


if __name__ == "__main__":
    unittest.main()
