"""Calibration fourpl."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Sequence
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit
from scipy.stats import f as f_distribution, t as t_distribution


PARAMETERS = ("bottom", "span", "midpoint", "slope")


TECHNICAL_SCOPE = "conditional_on_prepared_standards_and_observed_technical_units"


def _model(x: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    bottom, span, midpoint, slope = parameters
    return bottom + span * expit(-np.log(10.0) * (x - midpoint) * slope)


def _jacobian(x: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    _, span, midpoint, slope = parameters
    q = expit(-np.log(10.0) * (x - midpoint) * slope)
    d = q * (1.0 - q) * np.log(10.0)
    return np.column_stack((np.ones(len(x)), q, span * slope * d,
                            -span * (x - midpoint) * d))


def ferric_response(x: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    """Blank plus stretched-exponential response, exactly defined at zero.

    The ``slope`` field reports the fourth parameter, a dimensionless shape coefficient.
    """
    x = np.asarray(x, dtype=float)
    bottom, span, midpoint, exponent = parameters
    q = np.zeros_like(x)
    positive = x > 0
    power = np.exp(np.clip(exponent * (np.log(x[positive]) - np.log(midpoint)), -745, 700))
    q[positive] = -np.expm1(-np.log(2.) * power)
    q[x < 0] = np.nan
    return bottom + span * q


def ferric_response_jacobian(x: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    _, span, midpoint, exponent = parameters
    positive = x > 0
    log_ratio = np.zeros_like(x)
    log_ratio[positive] = np.log(x[positive]) - np.log(midpoint)
    power = np.zeros_like(x)
    power[positive] = np.log(2.) * np.exp(np.clip(exponent * log_ratio[positive], -745, 700))
    q = -np.expm1(-power)
    d = power * np.exp(-power)
    return np.column_stack((np.ones(len(x)), q, -span * exponent / midpoint * d,
                            span * log_ratio * d))


@dataclass
class FourPLInference:
    x: np.ndarray
    y: np.ndarray
    observation_indices: np.ndarray
    estimates: np.ndarray
    covariance: np.ndarray
    ci_low: np.ndarray
    ci_high: np.ndarray
    free_mask: np.ndarray
    parameter_identified: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    influence: np.ndarray
    normalized_syx_influence: np.ndarray
    sse: float
    df_residual: int
    r_squared: float
    sy_x: float
    normalized_sy_x: float
    converged: bool
    identified: bool
    flags: tuple[str, ...]
    diagnostics: dict[str, Any]

    def predict(self, x: Sequence[float]) -> np.ndarray:
        return _model(np.asarray(x, dtype=float), self.estimates)

    def as_dict(self) -> dict[str, Any]:
        out = {
            "n": len(self.y), "n_levels": len(np.unique(self.x)),
            "n_free_parameters": int(self.free_mask.sum()),
            "df_residual": self.df_residual, "sse": self.sse,
            "fit_rmse": float(np.sqrt(self.sse / len(self.y))) if len(self.y) else np.nan,
            "r_squared": self.r_squared, "sy_x": self.sy_x,
            "normalized_sy_x": self.normalized_sy_x,
            "converged": self.converged, "fit_identified": self.identified,
            "identification_flags": ";".join(self.flags),
            "covariance_method": "HC3_linearized_t_interval",
            "uncertainty_scope": TECHNICAL_SCOPE,
            **self.diagnostics,
        }
        for i, name in enumerate(PARAMETERS):
            variance = self.covariance[i, i]
            out.update({name: float(self.estimates[i]),
                        f"{name}_se": float(np.sqrt(max(0., variance)))
                        if self.free_mask[i] and np.isfinite(variance) else np.nan,
                        f"{name}_ci_low": float(self.ci_low[i]),
                        f"{name}_ci_high": float(self.ci_high[i]),
                        f"{name}_identified": bool(self.parameter_identified[i]),
                        f"{name}_fixed": not bool(self.free_mask[i])})
        span = abs(self.estimates[1])
        out["bottom_se_relative_to_span"] = out["bottom_se"] / span if span > 0 else np.nan
        out["span_se_relative_to_span"] = out["span_se"] / span if span > 0 else np.nan
        out["top"] = float(self.estimates[0] + self.estimates[1])
        v = np.asarray([1., 1., 0., 0.])
        top_se = np.sqrt(max(0., float(v @ self.covariance @ v))) if np.isfinite(self.covariance).all() else np.nan
        critical = t_distribution.ppf((1 + self.diagnostics.get("confidence_level", .95)) / 2,
                                      self.df_residual) if self.df_residual > 0 else np.nan
        out.update(top_se=top_se, top_ci_low=out["top"] - critical * top_se,
                   top_ci_high=out["top"] + critical * top_se)
        return out


class FerricResponseInference(FourPLInference):
    def predict(self, x: Sequence[float]) -> np.ndarray:
        return ferric_response(np.asarray(x, dtype=float), self.estimates)


def fit_ferric_response(x: Sequence[float], y: Sequence[float], *,
                      confidence: float = .95, max_nfev: int = 4000) -> FerricResponseInference:
    """Fit all nonnegative concentrations, including blank, with HC3 diagnostics.

    Positive span, half-response concentration and exponent define a saturating
    increasing response. The blank offset is estimated, allowing raw or
    background-corrected observations. No response-based range selection occurs.
    """
    xx, yy = np.asarray(x, float), np.asarray(y, float)
    keep = np.isfinite(xx) & np.isfinite(yy) & (xx >= 0)
    indices = np.flatnonzero(keep)
    xa, ya = xx[keep], yy[keep]
    n, p = len(xa), 4
    estimates = np.full(4, np.nan)
    covariance = np.full((4, 4), np.nan)
    low, high = np.full(4, np.nan), np.full(4, np.nan)
    influence = np.full((n, 4), np.nan)
    identified = np.zeros(4, bool)
    flags, converged = [], False
    diagnostics = dict(model_family='blank_plus_stretched_exponential',
        model_equation='bottom + span*(1-exp(-log(2)*(x/midpoint)**shape))',
        slope_parameter_meaning='dimensionless_response_shape',
        zero_concentration_n=int(np.sum(xa == 0)),
        excluded_negative_concentration_n=int(np.sum(np.isfinite(xx) & (xx < 0))),
        confidence_level=confidence)
    if len(np.unique(xa)) >= 5 and n > p and np.ptp(ya) > 0 and np.any(xa > 0):
        xscale = float(np.median(xa[xa > 0]))
        yscale = float(np.ptp(ya))
        ybase = float(np.mean(ya[xa == 0])) if np.any(xa == 0) else float(np.min(ya))
        z, v = xa / xscale, (ya - ybase) / yscale
        candidates = []
        for half in (np.median(z[z > 0]), np.quantile(z[z > 0], .75)):
            for exponent in (.8, 1.5):
                initial = [0., 1., half, exponent]
                result = least_squares(lambda a: ferric_response(z, a) - v, initial,
                    jac=lambda a: ferric_response_jacobian(z, a),
                    bounds=([-np.inf, 1e-10, 1e-10, 1e-6], [np.inf, np.inf, np.inf, np.inf]),
                    x_scale='jac', max_nfev=max_nfev)
                if np.isfinite(result.fun).all():
                    candidates.append(result)
        if candidates:
            result = min(candidates, key=lambda a: np.sum(a.fun ** 2))
            a = result.x
            scale_map = np.array([yscale, yscale, xscale, 1.])
            estimates = a * scale_map
            estimates[0] += ybase
            converged = bool(result.success)
            jac = ferric_response_jacobian(z, a)
            u, singular, vt = np.linalg.svd(jac, full_matrices=False)
            rank = int(np.linalg.matrix_rank(jac))
            condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else np.inf
            diagnostics.update(jacobian_rank=rank, scaled_jacobian_condition=condition)
            if not converged:
                flags.append('optimizer_did_not_converge')
            if rank == p and condition <= 1e8:
                inverse = (vt.T / singular) @ u.T
                leverage = np.sum(u ** 2, axis=1)
                residual = v - ferric_response(z, a)
                influence = ((inverse * (residual / np.maximum(1 - leverage, 1e-6))).T
                             * scale_map)
                covariance = influence.T @ influence
                halfwidth = t_distribution.ppf((1 + confidence) / 2, n - p) * np.sqrt(np.maximum(0, np.diag(covariance)))
                low, high = estimates - halfwidth, estimates + halfwidth
                identified[:] = converged
                weak = [halfwidth[0] > 2 * yscale,
                        halfwidth[1] > 2 * yscale or low[1] <= 0,
                        halfwidth[2] > np.ptp(xa) or not xa.min() <= estimates[2] <= xa.max(),
                        low[3] <= 0 or halfwidth[3] > estimates[3]]
                for i, value in enumerate(weak):
                    if value:
                        identified[i] = False
                        flags.append(f'{PARAMETERS[i]}_weakly_identified')
            else:
                flags.append('singular_or_ill_conditioned_jacobian')
    else:
        flags.append('insufficient_finite_levels_or_response_range')
    fitted = ferric_response(xa, estimates) if np.isfinite(estimates).all() else np.full(n, np.nan)
    residuals = ya - fitted
    sse, df = float(np.sum(residuals ** 2)), n - p
    syx = float(np.sqrt(sse / df)) if df > 0 else np.nan
    total = float(np.sum((ya - np.mean(ya)) ** 2)) if n else np.nan
    r2 = 1 - sse / total if total > 0 else np.nan
    span = estimates[1]
    normalized = syx / span if span > 0 else np.nan
    normalized_influence = np.full(n, np.nan)
    diagnostics.update(_lack_of_fit(xa, ya, sse, p))
    diagnostics['fit_rmse'] = float(np.sqrt(sse / n)) if n else np.nan
    if np.any(xa == 0):
        diagnostics['blank_mean_residual'] = float(np.mean(residuals[xa == 0]))
    if np.any(xa > 0):
        lowest = xa[xa > 0].min()
        diagnostics['lowest_positive_reference'] = float(lowest)
        diagnostics['lowest_positive_mean_residual'] = float(np.mean(residuals[xa == lowest]))
    return FerricResponseInference(xa, ya, indices, estimates, covariance, low, high,
        np.ones(4, bool), identified, fitted, residuals, influence, normalized_influence,
        sse, df, float(r2), syx, normalized, converged, bool(identified.all()), tuple(flags), diagnostics)


def _lack_of_fit(x: np.ndarray, y: np.ndarray, sse: float, p: int) -> dict[str, Any]:
    levels = np.unique(x)
    pure_error = float(sum(np.sum((y[x == level] - np.mean(y[x == level])) ** 2)
                           for level in levels))
    pure_df, lack_df = len(y) - len(levels), len(levels) - p
    stat = p_value = np.nan
    if pure_df > 0 and lack_df > 0 and pure_error > 0:
        stat = max(0., sse - pure_error) / lack_df / (pure_error / pure_df)
        p_value = float(f_distribution.sf(stat, lack_df, pure_df))
    return {"pure_error_sse": pure_error, "pure_error_df": pure_df,
            "lack_of_fit_df": lack_df, "lack_of_fit_f": stat,
            "lack_of_fit_p": p_value,
            "lack_of_fit_method": "approximate_nonlinear_F_requires_homoscedastic_independent_technical_errors"}


def fit_fourpl(x: Sequence[float], y: Sequence[float], *, expected: str | None = None,
               fixed_bottom: float | None = None, initial: Sequence[float] | None = None,
               confidence: float = .95, max_nfev: int = 3000) -> FourPLInference:
    """Unbounded deterministic 4PL fit with scaled optimization and analytic J.

    ``fixed_bottom`` constrains the lower asymptote to the specified value.
    Fixed parameters have no estimated SE/CI. Parameters are ordered
    bottom, span, midpoint, slope. Increasing curves have negative slope.
    Nonfinite observations are recorded through ``observation_indices``.
    """
    if expected not in (None, "increasing", "decreasing"):
        raise ValueError("expected must be increasing, decreasing or None")
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie between zero and one")
    xa, ya = np.asarray(x, float), np.asarray(y, float)
    if xa.ndim != 1 or ya.ndim != 1 or len(xa) != len(ya):
        raise ValueError("x and y must be matching one-dimensional arrays")
    keep = np.isfinite(xa) & np.isfinite(ya)
    indices = np.flatnonzero(keep)
    xa, ya = xa[keep], ya[keep]
    free = np.asarray([fixed_bottom is None, True, True, True])
    p, n = int(free.sum()), len(ya)
    flags: list[str] = []
    estim = np.full(4, np.nan)
    if fixed_bottom is not None:
        estim[0] = fixed_bottom
    covariance = np.full((4, 4), np.nan)
    low = high = np.full(4, np.nan)
    identified = np.zeros(4, bool)
    influence = np.full((n, 4), np.nan)
    normalized_influence = np.full(n, np.nan)
    xmin, xmax = (float(xa.min()), float(xa.max())) if n else (np.nan, np.nan)
    xrange = xmax - xmin
    yrange = float(np.ptp(ya)) if n else np.nan
    sufficient = n > p and len(np.unique(xa)) >= p and xrange > 0 and yrange > 0
    converged = False
    diagnostics: dict[str, Any] = {}
    if sufficient:
        center, ycenter = (xmin + xmax) / 2, float(np.min(ya))
        z, v = (xa - center) / xrange, (ya - ycenter) / yrange
        direction = -1. if expected == "increasing" else 1.
        if expected is None:
            direction = -1. if np.cov(xa, ya)[0, 1] > 0 else 1.
        fixed_scaled = (fixed_bottom - ycenter) / yrange if fixed_bottom is not None else None
        starts = []
        if initial is not None:
            initial_array = np.asarray(initial, float)
            if initial_array.shape != (4,):
                raise ValueError("initial must contain bottom, span, midpoint, slope")
            starts.append(np.asarray([(initial_array[0] - ycenter) / yrange,
                                      initial_array[1] / yrange,
                                      (initial_array[2] - center) / xrange,
                                      initial_array[3] * xrange]))
        # Four fixed starts avoid random optimizer behavior and cover broad transitions.
        starts += [np.asarray([0., 1., m, direction * s])
                   for m, s in ((0., 4.), (-.25, 8.), (.25, 8.), (0., 1.))]
        best = None
        for start in starts:
            if fixed_scaled is not None:
                start[0] = fixed_scaled

            def expand(q: np.ndarray) -> np.ndarray:
                full = start.copy()
                full[free] = q
                return full

            try:
                result = least_squares(lambda q: _model(z, expand(q)) - v, start[free],
                                       jac=lambda q: _jacobian(z, expand(q))[:, free],
                                       method="trf", max_nfev=max_nfev,
                                       ftol=1e-10, xtol=1e-10, gtol=1e-10)
                if np.isfinite(result.cost) and (best is None or result.cost < best[0]):
                    best = (result.cost, result, expand(result.x))
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                continue
        if best is not None:
            _, result, scaled = best
            # Equivalent positive-span orientation; never changes the response.
            if scaled[1] < 0 and fixed_bottom is None:
                scaled = np.asarray([scaled[0] + scaled[1], -scaled[1], scaled[2], -scaled[3]])
            estim = np.asarray([ycenter + yrange * scaled[0], yrange * scaled[1],
                                center + xrange * scaled[2], scaled[3] / xrange])
            if fixed_bottom is not None and scaled[1] <= 0:
                flags.append("fixed_bottom_is_not_lower_asymptote_nonpositive_span")
            converged = bool(result.success)
            j_scaled = _jacobian(z, scaled)[:, free]
            singular = np.linalg.svd(j_scaled, compute_uv=False)
            condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else np.inf
            diagnostics["scaled_jacobian_condition"] = condition
            rank = int(np.sum(singular > singular[0] * max(j_scaled.shape) * np.finfo(float).eps))
            diagnostics["jacobian_rank"] = rank
            if not converged:
                flags.append("optimizer_did_not_converge")
            if rank < p or condition > 1e8:
                flags.append("singular_or_ill_conditioned_jacobian")
            else:
                # SVD inverse avoids squaring the Jacobian condition number.
                u, s, vt = np.linalg.svd(j_scaled, full_matrices=False)
                inverse_j = (vt.T / s) @ u.T
                leverage = np.sum(u ** 2, axis=1)
                residual = v - _model(z, scaled)
                scale_map = np.asarray([yrange, yrange, xrange, 1 / xrange])
                inf_free = (inverse_j * (residual / np.maximum(1 - leverage, 1e-6))).T
                influence = np.zeros((n, 4))
                influence[:, free] = inf_free * scale_map[free]
                covariance = influence.T @ influence
                se = np.sqrt(np.maximum(0., np.diag(covariance)))
                critical = t_distribution.ppf((1 + confidence) / 2, n - p)
                low, high = estim - critical * se, estim + critical * se
                low[~free] = np.nan
                high[~free] = np.nan
                identified = free.copy() if converged else np.zeros(4, bool)
                # These are conservative reporting gates, not optimizer constraints.
                half_width = critical * se
                gates = [half_width[0] > 2 * yrange,
                         half_width[1] > 2 * yrange or low[1] <= 0,
                         half_width[2] > xrange or not xmin <= estim[2] <= xmax,
                         half_width[3] > abs(estim[3]) or low[3] <= 0 <= high[3]]
                for i, weak in enumerate(gates):
                    if free[i] and weak:
                        identified[i] = False
                        flags.append(f"{PARAMETERS[i]}_weakly_identified")
                if len(np.unique(xa)) <= p:
                    flags.append("no_level_degrees_of_freedom_for_lack_of_fit")
                derivative_sign = -estim[1] * estim[3]
                sign_matches = expected is None or (derivative_sign > 0 if expected == "increasing" else derivative_sign < 0)
                if not sign_matches:
                    flags.append("fitted_direction_opposes_expected_direction")
                diagnostics["sampled_transition_10_to_90"] = bool(
                    np.min(expit(-np.log(10.) * (xa - estim[2]) * estim[3])) <= .1
                    and np.max(expit(-np.log(10.) * (xa - estim[2]) * estim[3])) >= .9)
        else:
            flags.append("optimization_failed")
    else:
        flags.append("insufficient_finite_levels_or_response_range")
    fitted = _model(xa, estim) if np.isfinite(estim).all() else np.full(n, np.nan)
    residuals = ya - fitted
    sse = float(np.sum(residuals ** 2))
    df = n - p
    syx = float(np.sqrt(sse / df)) if df > 0 else np.nan
    total = float(np.sum((ya - np.mean(ya)) ** 2)) if n else np.nan
    r2 = 1 - sse / total if total > 0 else np.nan
    span = abs(estim[1])
    normalized = syx / span if span > 0 else np.nan
    if np.isfinite(influence).all() and span > 0 and syx > 0 and df > 0:
        syx_inf = (residuals ** 2 - sse / n) / (2 * syx * df)
        normalized_influence = syx_inf / span - syx * np.sign(estim[1]) / span ** 2 * influence[:, 1]
    diagnostics.update(_lack_of_fit(xa, ya, sse, p))
    diagnostics["confidence_level"] = confidence
    diagnostics["identification_gate_policy"] = (
        "full-rank scaled Jacobian condition<=1e8; bottom/span CI halfwidth<=2 observed response ranges; "
        "positive-span CI; midpoint inside sampled x with CI halfwidth<=sampled x range; "
        "slope CI excludes zero; reporting gates only, never optimization bounds")
    return FourPLInference(xa, ya, indices, estim, covariance, low, high, free,
                           identified, fitted, residuals, influence, normalized_influence,
                           sse, df, float(r2), syx, normalized, converged,
                           bool(np.all(identified[free])), tuple(flags), diagnostics)
