"""Instrument ceiling checks for biological and optical microplate analyses."""

from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from acidophile_data.biological import _well_values
from acidophile_data.plate_comparisons import _safe_channel
from acidophile_data.plate_reader import PlateData
from acidophile_data.publication_calibration import _channel
from acidophile_data.publication_ratio_support import classify_channel_ratios


class AbsorbanceLimitTests(unittest.TestCase):
    def setUp(self):
        values = np.full((3, 8, 12), np.nan)
        values[:, 0, :3] = [[3.99, 4.00, 4.01], [1., 1., 1.], [.1, .1, .1]]
        self.plate = PlateData(
            path=Path("boundary.txt"), mode="absorbance", wavelengths_nm=(545, 485, 750),
            values=values, exported_mask=np.isfinite(values),
            parse_error_mask=np.zeros_like(values, dtype=bool),
            censored_high_mask=np.zeros_like(values, dtype=bool),
            temperature_c=None, saved_at=None, original_filename=None,
        )
        self.wells = ("A1", "A2", "A3")

    def test_biological_channels_and_ratios_include_unflagged_four_au(self):
        channels = pd.DataFrame([_well_values(self.plate, well) for well in self.wells])
        expected = [True, True, False]
        self.assertEqual(channels.A545_eligible.tolist(), expected)
        bare = channels[["A545", "A485", "A750"]]
        for source in (channels, bare, bare.assign(A545_eligible=True)):
            self.assertEqual(_safe_channel(source, "A545").notna().tolist(), expected)
            ratios = classify_channel_ratios(source)
            self.assertEqual(ratios.raw_ok.tolist(), expected)
            self.assertEqual(ratios.corrected_ok.tolist(), expected)
        flagged = channels.assign(A545_eligible=False)
        self.assertFalse(_safe_channel(flagged, "A545").notna().any())
        self.assertFalse(classify_channel_ratios(flagged).raw_ok.any())

    def test_optical_calibration_includes_unflagged_four_au(self):
        self.assertEqual(
            [_channel(self.plate, well, 545) for well in self.wells],
            [(3.99, "eligible"), (4.00, "eligible"), (4.01, "numeric_ceiling")],
        )
