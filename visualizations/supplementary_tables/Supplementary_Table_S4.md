# Supplementary Table S4. Acid-specific ferric calibrations

The fitted A380 response approaches a plateau. The acid-specific model gives baseline and span in AU, half-response concentration in mM total Fe and a dimensionless shape parameter.

**Figure panels** 1f · [All Supplementary tables](README.md)

**Column definitions**

- `acid_M` — Nominal H₂SO₄ concentration in the dye-free ferric calibration series, in M.
- `plate_n` — Number of distinct physical microplates contributing at least one eligible standard mean at this acid concentration.
- `plate_mean_n` — Number of plate-by-standard mean A380 values paired with a measured, nonnegative total Fe concentration for the fit.
- `level_n` — Number of distinct measured total Fe concentrations represented among the fitted standard means.
- `fit_rmse` — Root mean squared difference between observed and fitted A380 values, in absorbance units (AU).
- `r_squared` — Fraction of variation in the observed standard means explained by the fitted ferric response.
- `identified` — True when all four fitted parameters pass the numerical identification checks, including parameter precision and coverage of the half-response concentration.
- `converged` — True when the numerical optimizer reported a solution for the stretched-exponential ferric response model.
- `identification_flags` — Diagnostic reasons for an unavailable or weakly determined fit, with a blank indicating that no flag was raised.
- `ci_df` — Degrees of freedom for 95% Student t intervals using HC3 uncertainty clustered by physical plate, equal to the contributing plate count minus one.
- `bottom` — Fitted A380 at zero total Fe, in AU.
- `bottom_ci_low` — Lower 95% plate-clustered confidence limit for the zero-Fe baseline (AU).
- `bottom_ci_high` — Upper 95% plate-clustered confidence limit for the zero-Fe baseline (AU).
- `span` — Fitted increase in A380 from the zero-Fe baseline to the upper plateau, in AU.
- `span_ci_low` — Lower 95% plate-clustered confidence limit for the increase from baseline to plateau (AU).
- `span_ci_high` — Upper 95% plate-clustered confidence limit for the increase from baseline to plateau (AU).
- `midpoint` — Measured total Fe concentration (mM) at which the model reaches half the fitted increase above baseline.
- `midpoint_ci_low` — Lower 95% plate-clustered confidence limit for the half-response total Fe concentration (mM).
- `midpoint_ci_high` — Upper 95% plate-clustered confidence limit for the half-response total Fe concentration (mM).
- `slope` — Dimensionless exponent controlling the shape of the stretched-exponential ferric response curve.
- `slope_ci_low` — Lower 95% plate-clustered confidence limit for the dimensionless ferric response exponent.
- `slope_ci_high` — Upper 95% plate-clustered confidence limit for the dimensionless ferric response exponent.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| acid&#95;<wbr>M | plate&#95;<wbr>n | plate&#95;<wbr>mean&#95;<wbr>n | level&#95;<wbr>n | fit&#95;<wbr>rmse | r&#95;<wbr>squared | identified | converged | identification&#95;<wbr>flags | ci&#95;<wbr>df | bottom | bottom&#95;<wbr>ci&#95;<wbr>low | bottom&#95;<wbr>ci&#95;<wbr>high | span | span&#95;<wbr>ci&#95;<wbr>low | span&#95;<wbr>ci&#95;<wbr>high | midpoint | midpoint&#95;<wbr>ci&#95;<wbr>low | midpoint&#95;<wbr>ci&#95;<wbr>high | slope | slope&#95;<wbr>ci&#95;<wbr>low | slope&#95;<wbr>ci&#95;<wbr>high |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.001 | 3 | 36 | 30 | 0.055089962384382844 | 0.9982387234531114 | True | True |  | 2 | 0.02004460902504711 | -0.006109889300484308 | 0.04619910735057853 | 3.1106978102284693 | 3.0023897017743915 | 3.219005918682547 | 32.59797379874302 | 28.037280201085053 | 37.158667396400986 | 1.076313836947537 | 0.9828508900369325 | 1.1697767838581414 |
| 0.01 | 3 | 36 | 31 | 0.07857364268814852 | 0.9964450106290237 | True | True |  | 2 | 0.033728693132323496 | 0.023346187002210222 | 0.04411119926243677 | 3.1024364113191623 | 3.0288807903564523 | 3.175992032281872 | 36.80656930183721 | 18.768679152487124 | 54.84445945118729 | 1.1532858138637734 | 1.0628894597590213 | 1.2436821679685255 |
| 0.1 | 3 | 36 | 23 | 0.08617225954076646 | 0.9955867808882631 | True | True |  | 2 | 0.03840056464278879 | 0.016878325858599295 | 0.05992280342697829 | 3.048946946940279 | 3.014094137767077 | 3.083799756113481 | 37.60205913057858 | 18.267353546096793 | 56.93676471506037 | 1.167951666883524 | 1.1237452164336976 | 1.2121581173333502 |
| 1.0 | 3 | 36 | 34 | 0.08852966244148087 | 0.9953177546131359 | True | True |  | 2 | 0.03977920264547603 | 0.03465902333190044 | 0.044899381959051625 | 3.053985874057111 | 2.949537113706142 | 3.15843463440808 | 40.94801186322335 | 22.76699551908585 | 59.12902820736085 | 1.171700221488205 | 1.03897621896307 | 1.3044242240133401 |
