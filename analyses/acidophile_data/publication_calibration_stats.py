"""Calibration stats."""
from __future__ import annotations
import numpy as np


Array = np.ndarray


def four_parameter_logistic(x: Array, bottom: float, span: float, midpoint: float, slope: float) -> Array:
    exponent = np.clip((np.asarray(x) - midpoint) * slope, -300.0, 300.0)
    return bottom + span / (1.0 + np.power(10.0, exponent))
