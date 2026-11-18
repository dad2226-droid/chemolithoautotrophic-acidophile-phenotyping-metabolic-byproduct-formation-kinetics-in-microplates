# Supplementary Table S5. Cerimetry against measured total iron

The ferrous-standard slopes are below one in 0.001 and 1 M sulfuric acid. The ferric-standard series includes a constant cerimetry response. Tests compare the slope with zero and one and the intercept with zero.

**Figure panels** 1g · [All Supplementary tables](README.md)

**Column definitions**

- `series_ID` — Chemical standard series, either ferric standards in 1 M H₂SO₄ or ferrous standards in 0.001 or 1 M H₂SO₄.
- `n` — Number of reference levels used in the ordinary least-squares regression, each represented by a total Fe concentration and a mean of three Fe²⁺ titrations.
- `x_min_M` — Lowest measured total Fe concentration included in the regression, in M.
- `x_max_M` — Highest measured total Fe concentration included in the regression, in M.
- `slope` — Fitted increase in Fe²⁺ measured by cerimetry per unit increase in measured total Fe (M/M).
- `slope_ci_low` — Lower 95% Student t confidence limit for the regression slope (M/M), treating measured total Fe as fixed.
- `slope_ci_high` — Upper 95% Student t confidence limit for the regression slope (M/M), treating measured total Fe as fixed.
- `intercept_M` — Fitted Fe²⁺ concentration at zero total Fe, in M.
- `intercept_ci_low_M` — Lower 95% Student t confidence limit for the fitted intercept, in M.
- `intercept_ci_high_M` — Upper 95% Student t confidence limit for the fitted intercept, in M.
- `sy_x_M` — Residual standard deviation of measured Fe²⁺ about the fitted line, in M.
- `residual_df` — Residual degrees of freedom, equal to the number of fitted reference levels minus two for the estimated intercept and slope.
- `slope_vs_zero_p` — Unadjusted two-sided p-value for the hypothesis that the slope equals zero.
- `slope_zero_p_holm` — P-value for a zero slope after Holm adjustment across all estimable slope and intercept tests in the three series.
- `slope_vs_one_p` — Unadjusted two-sided p-value for the hypothesis that the slope equals one (a 1 M/M response).
- `slope_one_p_holm` — P-value for a unit slope after Holm adjustment across all estimable slope and intercept tests in the three series.
- `intercept_vs_zero_p` — Unadjusted two-sided p-value for the hypothesis that the intercept equals 0 M Fe²⁺.
- `intercept_zero_p_holm` — P-value for a zero intercept after Holm adjustment across all estimable slope and intercept tests in the three series.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| series&#95;<wbr>ID | n | x&#95;<wbr>min&#95;<wbr>M | x&#95;<wbr>max&#95;<wbr>M | slope | slope&#95;<wbr>ci&#95;<wbr>low | slope&#95;<wbr>ci&#95;<wbr>high | intercept&#95;<wbr>M | intercept&#95;<wbr>ci&#95;<wbr>low&#95;<wbr>M | intercept&#95;<wbr>ci&#95;<wbr>high&#95;<wbr>M | sy&#95;<wbr>x&#95;<wbr>M | residual&#95;<wbr>df | slope&#95;<wbr>vs&#95;<wbr>zero&#95;<wbr>p | slope&#95;<wbr>zero&#95;<wbr>p&#95;<wbr>holm | slope&#95;<wbr>vs&#95;<wbr>one&#95;<wbr>p | slope&#95;<wbr>one&#95;<wbr>p&#95;<wbr>holm | intercept&#95;<wbr>vs&#95;<wbr>zero&#95;<wbr>p | intercept&#95;<wbr>zero&#95;<wbr>p&#95;<wbr>holm |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ferric | 12 | 0.0 | 0.0976 | 0.0 |  |  | 0.0006 |  |  | 0.0 | 10 |  |  |  |  |  |  |
| ferrous&#95;<wbr>0.001M | 12 | 0.0 | 0.09305706 | 0.9174446328393027 | 0.8733187137003241 | 0.9615705519782812 | 0.001511557679955488 | 0.0001490451921759518 | 0.002874070167735024 | 0.0018341867441568627 | 10 | 5.291313132232571e-13 | 3.1747878793395426e-12 | 0.0019221376286089264 | 0.0076885505144357056 | 0.032998761902939454 | 0.06599752380587891 |
| ferrous&#95;<wbr>1M | 12 | 0.0 | 0.09847928 | 0.9249817469901457 | 0.8802016017329427 | 0.9697618922473487 | 0.001131046464387251 | -0.00033892551217802776 | 0.0026010184409525298 | 0.001979514576775998 | 10 | 5.647111226784471e-13 | 3.1747878793395426e-12 | 0.003892984542090639 | 0.011678953626271917 | 0.11722327694549835 | 0.11722327694549835 |
