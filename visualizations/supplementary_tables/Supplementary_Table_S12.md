# Supplementary Table S12. Paired scatter-correction effects in growth media

Across seven estimable comparisons, Holm-adjusted p ≥ 0.891 for changes in span divided by within-standard SD. Each comparison uses the same eligible wells for both readouts. The effect is log₂(corrected/measured).

**Figure panels** 4d · [All Supplementary tables](README.md)

**Column definitions**

- `medium` — Growth-medium code for the standards. St-BS and St-YS denote Starkey medium with brown and yellow elemental sulfur, respectively.
- `readout` — Identifies the measured signal paired with its scatter-corrected counterpart. The pair is A380 and A380 − A750 for iron, or A545/A485 and (A545 − A750)/(A485 − A750) for pH.
- `scheduled_plate_n` — Number of physical plates assigned to this medium and measured-versus-corrected readout pair.
- `unavailable_plate_n` — Number of scheduled plates without a finite log₂ ratio of corrected to measured Q. Both plate Q values must be finite and positive.
- `common_well_n` — Total physical wells eligible for both readouts across scheduled plates, including plates where a Q contrast could not be estimated.
- `raw_Q_mean` — Arithmetic mean of available plate Q values for the measured signal. For each plate, Q is the span of common-well standard means divided by pooled within-standard technical SD. Single-well means may contribute to span, while only replicated levels supply variance.
- `corrected_Q_mean` — Arithmetic mean of available plate Q values for the scatter-corrected signal, calculated on the same eligible wells and using the same span and pooled-SD rules as the measured signal.
- `plate_n` — Number of physical microplates with a finite log₂ Q contrast.
- `estimate` — Mean across plates of log₂(corrected Q / measured Q). A positive value indicates greater response span relative to technical scatter after correction.
- `sd` — Sample standard deviation of the log₂ Q contrast across physical microplates, in log₂-fold units.
- `ci_low` — Lower limit of the two-sided 95% Student’s t confidence interval for the mean log₂ Q contrast, in log₂-fold units.
- `ci_high` — Upper limit of the two-sided 95% Student’s t confidence interval for the mean log₂ Q contrast, in log₂-fold units.
- `p_value` — Unadjusted two-sided p-value for a one-sample t test of whether the mean log₂ Q contrast is zero.
- `p_holm` — P-value after Holm adjustment across the seven estimable medium–readout contrasts.
- `multiplicity_family_n` — Number of contrasts with finite p-values in the Holm adjustment family (the seven estimable medium–readout contrasts).
- `status` — Whether the plate-level mean and its interval could be estimated. A single contributing plate permits a mean only.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| medium | readout | scheduled&#95;<wbr>plate&#95;<wbr>n | unavailable&#95;<wbr>plate&#95;<wbr>n | common&#95;<wbr>well&#95;<wbr>n | raw&#95;<wbr>Q&#95;<wbr>mean | corrected&#95;<wbr>Q&#95;<wbr>mean | plate&#95;<wbr>n | estimate | sd | ci&#95;<wbr>low | ci&#95;<wbr>high | p&#95;<wbr>value | p&#95;<wbr>holm | multiplicity&#95;<wbr>family&#95;<wbr>n | status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| AFM1 | A380 | 3 | 0 | 72 | 54.36846913366426 | 49.15739717421824 | 3 | -0.21588972616885962 | 0.3747304269739694 | -1.1467717115552003 | 0.714992259217481 | 0.4234707264527847 | 1.0 | 7 | estimated |
| AFM3 | A380 | 3 | 0 | 72 | 112.52269229915186 | 117.62944633093453 | 3 | 0.04034410053668503 | 0.07252915179360722 | -0.13982830063507876 | 0.2205165017084488 | 0.43697796904019465 | 1.0 | 7 | estimated |
| AFM4 | R | 3 | 2 | 10 | 4.588455197652633 | 3.034427654410635 | 1 | -0.596584096121474 |  |  |  |  |  | 7 | single&#95;<wbr>plate&#95;<wbr>no&#95;<wbr>interval |
| AFM5 | R | 3 | 0 | 72 | 16.044411447960165 | 8.94919894925514 | 3 | -0.6886484049050408 | 0.5400082364207596 | -2.030103229653638 | 0.6528064198435561 | 0.1578285026931752 | 0.9469710161590512 | 7 | estimated |
| F2S | A380 | 3 | 0 | 72 | 12.075919603507984 | 20.542324183030594 | 3 | 0.631955662003572 | 0.43293687468731024 | -0.44351915521450214 | 1.7074304792216461 | 0.12725641641924498 | 0.8907949149347149 | 7 | estimated |
| F2S | R | 3 | 0 | 72 | 3.7113187165136168 | 3.7249119890998728 | 3 | -0.04231358557757417 | 0.2562089444512869 | -0.6787718865767615 | 0.5941447154216132 | 0.8017450980198793 | 1.0 | 7 | estimated |
| St-BS | R | 3 | 0 | 72 | 12.496908717558483 | 12.409115799129822 | 3 | -0.015365046915682226 | 0.20408539715851134 | -0.5223412784146838 | 0.49161118458331937 | 0.9081817210563584 | 1.0 | 7 | estimated |
| St-YS | R | 3 | 0 | 64 | 4.173664046224405 | 2.7324199803715428 | 3 | -0.7002582647699581 | 0.570728838017628 | -2.118027294472993 | 0.7175107649330769 | 0.1674886970856234 | 0.9469710161590512 | 7 | estimated |
