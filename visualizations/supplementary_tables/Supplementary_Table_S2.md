# Supplementary Table S2. Scatter correction and pH-indicator precision

At 10 µM thymol blue, correction increases response span and within-standard SD. The span-to-SD contrast has Holm-adjusted p = 0.0505. The estimates and their intervals use log₂(corrected/measured). The geometric ratios give the corresponding folds.

**Figure panels** 1c–d · [All Supplementary tables](README.md)

**Column definitions**

- `metric` — Within-plate quantity selected from the range of standard means, the pooled within-standard well SD, or the range divided by that SD.
- `raw_mean` — Arithmetic mean of the plate-specific metric calculated with the measured absorbance ratio A545/A485.
- `corrected_mean` — Arithmetic mean of the plate-specific metric after A750 is subtracted from both A545 and A485 before forming the ratio.
- `plate_n` — Number of physical plates with a finite log₂ ratio of the corrected metric to the measured metric.
- `estimate` — Mean across plates of log₂(corrected metric/measured metric), where zero denotes equal metrics.
- `sd` — Sample standard deviation among the plate-specific log₂(corrected/measured) values.
- `ci_low` — Lower limit of the 95% Student t confidence interval for the mean plate log₂ ratio.
- `ci_high` — Upper limit of the 95% Student t confidence interval for the mean plate log₂ ratio.
- `p_value` — Unadjusted two-sided one-sample t-test p-value for a mean plate log₂ ratio of zero.
- `p_holm` — P-value after Holm adjustment across the tests of span, pooled within-standard SD and their quotient.
- `multiplicity_family_n` — Number of estimable metric tests included in the Holm adjustment.
- `status` — Indicates whether the mean log₂ ratio and its interval could be estimated from the contributing physical plates.
- `geometric_ratio` — Geometric mean of the plate-specific corrected/measured ratios, equal to 2 raised to the mean log₂ ratio.
- `ratio_ci_low` — Lower 95% confidence limit for the corrected/measured fold change, obtained by exponentiating the lower log₂ limit with base 2.
- `ratio_ci_high` — Upper 95% confidence limit for the corrected/measured fold change, obtained by exponentiating the upper log₂ limit with base 2.

| metric | raw&#95;<wbr>mean | corrected&#95;<wbr>mean | plate&#95;<wbr>n | estimate | sd | ci&#95;<wbr>low | ci&#95;<wbr>high | p&#95;<wbr>value | p&#95;<wbr>holm | multiplicity&#95;<wbr>family&#95;<wbr>n | status | geometric&#95;<wbr>ratio | ratio&#95;<wbr>ci&#95;<wbr>low | ratio&#95;<wbr>ci&#95;<wbr>high |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| technical SD | 0.03282522919972126 | 0.16346484899832625 | 3 | 2.3242919980737815 | 0.22064038281946113 | 1.7761909023769284 | 2.8723930937706346 | 0.002990308250131149 | 0.007257920715513632 | 3 | estimated | 5.008199379785926 | 3.4252063433012223 | 7.322788326823538 |
| span/technical SD | 66.87051813826257 | 44.2026629035848 | 3 | -0.5955747910973491 | 0.24098261628761505 | -1.194208796093673 | 0.0030592138989747886 | 0.05047663812861834 | 0.05047663812861834 | 3 | estimated | 0.6617807394072253 | 0.4370260564098669 | 1.0021227353081137 |
| response span | 2.174278226213571 | 7.199943978310476 | 3 | 1.728717206976432 | 0.14754319965464174 | 1.3621995806020282 | 2.095234833350836 | 0.0024193069051712105 | 0.007257920715513632 | 3 | estimated | 3.314329888653536 | 2.5707682871647655 | 4.272957102227597 |
