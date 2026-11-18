# Supplementary Table S1. pH indicator calibration

Four-parameter logistic fits describe the measured and scatter-corrected ratio response where replicate groups are complete. Bottom and span are in ratio units, midpoint is pH and slope is per pH unit.

**Figure panels** 1b–d · [All Supplementary tables](README.md)

**Column definitions**

- `dye_uM` — Nominal thymol blue concentration in the prepared standard, in µM.
- `readout` — R is A545/A485, and R&#42; is (A545 − A750)/(A485 − A750), using absorbance at the indicated wavelengths in nm.
- `plate_n` — Number of distinct physical microplates contributing at least one eligible standard mean to the fit.
- `plate_mean_n` — Number of four-well standard means on individual microplates with both a finite ratio and a measured pH.
- `level_n` — Number of distinct measured pH values among the standard means used to fit the curve.
- `fit_rmse` — Root mean squared difference between the observed standard means and the fitted ratios, in ratio units.
- `r_squared` — Fraction of the variation in the observed standard means explained by the fitted curve.
- `identified` — True when all four fitted parameters pass the numerical identification checks, including parameter precision and coverage of the midpoint by the measured pH values.
- `converged` — True when the numerical optimizer reported a solution for the four-parameter logistic fit.
- `identification_flags` — Diagnostic reasons for an unavailable or weakly determined fit, such as too few finite levels or a poorly determined parameter. A blank means no diagnostic flag was raised.
- `ci_df` — Degrees of freedom for 95% Student t intervals using HC3 uncertainty clustered by physical plate, equal to the contributing plate count minus one (−1 indicates no contributing plates).
- `bottom` — Fitted lower asymptote of the indicator ratio.
- `bottom_ci_low` — Lower 95% plate-clustered confidence limit for the lower asymptote (ratio units).
- `bottom_ci_high` — Upper 95% plate-clustered confidence limit for the lower asymptote (ratio units).
- `span` — Fitted upper asymptote minus the lower asymptote, in ratio units.
- `span_ci_low` — Lower 95% plate-clustered confidence limit for the fitted ratio span (ratio units).
- `span_ci_high` — Upper 95% plate-clustered confidence limit for the fitted ratio span (ratio units).
- `midpoint` — Measured pH at which the fitted ratio is halfway between the two asymptotes.
- `midpoint_ci_low` — Lower 95% plate-clustered confidence limit for the half-response pH.
- `midpoint_ci_high` — Upper 95% plate-clustered confidence limit for the half-response pH.
- `slope` — Steepness coefficient of the base-10 logistic curve, per pH unit, with positive values indicating a ratio that falls as pH increases.
- `slope_ci_low` — Lower 95% plate-clustered confidence limit for the base-10 logistic steepness (per pH unit).
- `slope_ci_high` — Upper 95% plate-clustered confidence limit for the base-10 logistic steepness (per pH unit).
- `unavailable_plate_mean_n` — Number of scheduled plate-by-standard means without an eligible four-well ratio mean.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| dye&#95;<wbr>uM | readout | plate&#95;<wbr>n | plate&#95;<wbr>mean&#95;<wbr>n | level&#95;<wbr>n | fit&#95;<wbr>rmse | r&#95;<wbr>squared | identified | converged | identification&#95;<wbr>flags | ci&#95;<wbr>df | bottom | bottom&#95;<wbr>ci&#95;<wbr>low | bottom&#95;<wbr>ci&#95;<wbr>high | span | span&#95;<wbr>ci&#95;<wbr>low | span&#95;<wbr>ci&#95;<wbr>high | midpoint | midpoint&#95;<wbr>ci&#95;<wbr>low | midpoint&#95;<wbr>ci&#95;<wbr>high | slope | slope&#95;<wbr>ci&#95;<wbr>low | slope&#95;<wbr>ci&#95;<wbr>high | unavailable&#95;<wbr>plate&#95;<wbr>mean&#95;<wbr>n |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | R&#42; | 0 | 0 | 0 |  |  | False | False | insufficient&#95;<wbr>finite&#95;<wbr>levels&#95;<wbr>or&#95;<wbr>response&#95;<wbr>range | -1 |  |  |  |  |  |  |  |  |  |  |  |  | 72 |
| 1 | R | 3 | 72 | 23 | 0.03250798643155862 | 0.977957904232591 | True | True |  | 2 | 0.8794228564319742 | 0.8306170115212755 | 0.9282287013426728 | 0.6097165039995909 | 0.3088664178634017 | 0.91056659013578 | 1.931405753540558 | 1.6662274498598086 | 2.1965840572213073 | 0.8417573567485513 | 0.690961999498722 | 0.9925527139983807 | 0 |
| 10 | R&#42; | 3 | 72 | 23 | 0.15802655623610792 | 0.995713178448507 | True | True |  | 2 | -0.060711442644031874 | -0.153219802621195 | 0.03179691733313125 | 8.200006089777194 | 7.160799010749973 | 9.239213168804413 | 1.4329583627009461 | 1.249646839119225 | 1.6162698862826672 | 0.7817440511088445 | 0.7583382092086051 | 0.8051498930090839 | 0 |
| 10 | R | 3 | 72 | 23 | 0.06958926035773849 | 0.9917151864727864 | True | True |  | 2 | 0.5309224467051882 | 0.4847129712455341 | 0.5771319221648422 | 2.4494699500970034 | 1.9691240920245707 | 2.929815808169436 | 1.6375884816163462 | 1.3923401990168678 | 1.8828367642158246 | 0.7277508431893078 | 0.5944287983571521 | 0.8610728880214635 | 0 |
