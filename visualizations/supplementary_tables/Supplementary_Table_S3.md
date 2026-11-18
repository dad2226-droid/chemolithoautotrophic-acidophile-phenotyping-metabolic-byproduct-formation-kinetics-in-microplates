# Supplementary Table S3. pH prediction in acid standards and ferric mixtures

Scatter correction reduces withheld-acid prediction RMSE. Both readouts have higher error in ferric mixtures over the common pH support. The paired SD comparison is a ratio. Paired bias and RMSE contrasts are corrected minus measured.

**Figure panels** 1b–d · [All Supplementary tables](README.md)

**Column definitions**

- `acid_M` — Nominal H₂SO₄ concentration (M) of the ferric mixture that defines the comparison, including the matching pH range used for the held-out acid standards.
- `matrix` — Prediction target, either acid standards on the microplate withheld from calibration or ferric mixtures evaluated with the acid calibration.
- `readout_or_comparison` — Measured A545/A485 (R), its A750-corrected counterpart (R&#42;), or the paired comparison of these two readouts for the same plate and standard.
- `metric` — Bias is predicted minus reference pH, RMSE is root mean squared prediction error, and technical SD is well-to-well ratio scatter converted to pH by the local calibration slope.
- `unit` — pH units for errors, SDs and paired differences, or a dimensionless fold for the paired corrected/measured SD ratio.
- `attempted_n` — All physical plate–standard observations considered before pH-range filtering for a single readout, or observations within the shared pH support of both readouts for a paired comparison.
- `eligible_n` — Physical plate–standard observations on shared pH support with a finite prediction for the selected readout, or finite predictions for both readouts in paired comparisons.
- `support_reference_pH_min` — Lowest measured reference pH retained within the common bins for the row’s readout or readout pair, before requiring a finite inverse prediction.
- `support_reference_pH_max` — Highest measured reference pH retained within the common bins for the row’s readout or readout pair, before requiring a finite inverse prediction.
- `support_shared_bin_n` — Number of occupied 0.25-pH-unit reference bins shared by the acid standards and ferric mixtures for this readout or pair.
- `estimate` — Metric computed with equal weight per shared pH bin, with paired rows giving corrected minus measured bias or RMSE, or corrected divided by measured technical SD.
- `ci_low` — Lower 95% bootstrap confidence limit (2.5th percentile) from resampling whole wells within each plate and standard and refitting the calibration.
- `ci_high` — Upper 95% bootstrap confidence limit (97.5th percentile) from resampling whole wells within each plate and standard and refitting the calibration.
- `bootstrap_requested_n` — Number of whole-well bootstrap resamples requested for this metric.
- `bootstrap_success_n` — Number of resamples yielding a finite value of this metric, with limits reported only when at least 20 and at least 80% of requested resamples succeed.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| acid&#95;<wbr>M | matrix | readout&#95;<wbr>or&#95;<wbr>comparison | metric | unit | attempted&#95;<wbr>n | eligible&#95;<wbr>n | support&#95;<wbr>reference&#95;<wbr>pH&#95;<wbr>min | support&#95;<wbr>reference&#95;<wbr>pH&#95;<wbr>max | support&#95;<wbr>shared&#95;<wbr>bin&#95;<wbr>n | estimate | ci&#95;<wbr>low | ci&#95;<wbr>high | bootstrap&#95;<wbr>requested&#95;<wbr>n | bootstrap&#95;<wbr>success&#95;<wbr>n |
| ---: | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.001 | held-out acid | R&#42; | bias | pH | 72 | 24 | 1.05 | 2.65 | 7 | -0.0018661666665977755 | -0.005735919845185099 | 0.0022350192723012627 | 199 | 199 |
| 0.001 | ferric transfer | R&#42; | bias | pH | 36 | 34 | 1.04 | 2.7 | 7 | 0.5870272363623654 | 0.5688841588979695 | 0.604172985009148 | 199 | 199 |
| 0.001 | held-out acid | R | bias | pH | 72 | 24 | 1.05 | 2.65 | 7 | -0.006614691850537335 | -0.0100475655230268 | -0.003221656013552747 | 199 | 199 |
| 0.001 | ferric transfer | R | bias | pH | 36 | 34 | 1.04 | 2.7 | 7 | 0.6550174982647274 | 0.6443744785031428 | 0.6657515889874701 | 199 | 199 |
| 0.01 | held-out acid | R&#42; | bias | pH | 72 | 13 | 0.97 | 1.77 | 5 | -0.018312467745949988 | -0.027190281982893826 | -0.009942109671555628 | 199 | 199 |
| 0.01 | ferric transfer | R&#42; | bias | pH | 36 | 36 | 0.93 | 1.79 | 5 | 1.0071772632869185 | 0.9874587138809139 | 1.0250937266606446 | 199 | 199 |
| 0.01 | held-out acid | R | bias | pH | 72 | 13 | 0.97 | 1.77 | 5 | 0.019937032846134527 | 0.014227345795948821 | 0.026141237827515377 | 199 | 199 |
| 0.01 | ferric transfer | R | bias | pH | 36 | 36 | 0.93 | 1.79 | 5 | 1.0862365225797743 | 1.0727385252496628 | 1.1081903977216985 | 199 | 199 |
| 0.1 | held-out acid | R&#42; | bias | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | R&#42; | bias | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | held-out acid | R | bias | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | R | bias | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | R&#42; | bias | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | R&#42; | bias | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | R | bias | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | R | bias | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.001 | held-out acid | R&#42; | RMSE | pH | 72 | 24 | 1.05 | 2.65 | 7 | 0.04475134971875612 | 0.03752121939215306 | 0.06447850313538377 | 199 | 199 |
| 0.001 | ferric transfer | R&#42; | RMSE | pH | 36 | 34 | 1.04 | 2.7 | 7 | 0.8622775092211178 | 0.8449952464669352 | 0.8797968591143618 | 199 | 199 |
| 0.001 | held-out acid | R | RMSE | pH | 72 | 24 | 1.05 | 2.65 | 7 | 0.11097368196912925 | 0.10089374474151085 | 0.12501920460020552 | 199 | 199 |
| 0.001 | ferric transfer | R | RMSE | pH | 36 | 34 | 1.04 | 2.7 | 7 | 0.949861511555476 | 0.9380366933240848 | 0.9611013010378149 | 199 | 199 |
| 0.01 | held-out acid | R&#42; | RMSE | pH | 72 | 13 | 0.97 | 1.77 | 5 | 0.054152154034662905 | 0.04042126731926034 | 0.0763758674992835 | 199 | 199 |
| 0.01 | ferric transfer | R&#42; | RMSE | pH | 36 | 36 | 0.93 | 1.79 | 5 | 1.1845949296478722 | 1.1643979757921632 | 1.2092017078065882 | 199 | 199 |
| 0.01 | held-out acid | R | RMSE | pH | 72 | 13 | 0.97 | 1.77 | 5 | 0.1449947556458059 | 0.13299834072527958 | 0.1636359011504822 | 199 | 199 |
| 0.01 | ferric transfer | R | RMSE | pH | 36 | 36 | 0.93 | 1.79 | 5 | 1.2780059627118625 | 1.2602166031356385 | 1.309535888386919 | 199 | 199 |
| 0.1 | held-out acid | R&#42; | RMSE | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | R&#42; | RMSE | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | held-out acid | R | RMSE | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | R | RMSE | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | R&#42; | RMSE | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | R&#42; | RMSE | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | R | RMSE | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | R | RMSE | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.001 | held-out acid | R&#42; | technical SD | pH | 72 | 24 | 1.05 | 2.65 | 7 | 0.05517758387738698 | 0.039358723332780496 | 0.05450229826092802 | 199 | 199 |
| 0.001 | ferric transfer | R&#42; | technical SD | pH | 36 | 34 | 1.04 | 2.7 | 7 | 0.03245382437337436 | 0.020600505426132766 | 0.032901007260922516 | 199 | 199 |
| 0.001 | held-out acid | R | technical SD | pH | 72 | 24 | 1.05 | 2.65 | 7 | 0.03988922829561593 | 0.02751072203715908 | 0.03871276073884293 | 199 | 199 |
| 0.001 | ferric transfer | R | technical SD | pH | 36 | 34 | 1.04 | 2.7 | 7 | 0.038800847003423185 | 0.027966191754954378 | 0.039457148189765985 | 199 | 199 |
| 0.01 | held-out acid | R&#42; | technical SD | pH | 72 | 13 | 0.97 | 1.77 | 5 | 0.052395104076973814 | 0.03265718847606604 | 0.05296284682620577 | 199 | 199 |
| 0.01 | ferric transfer | R&#42; | technical SD | pH | 36 | 36 | 0.93 | 1.79 | 5 | 0.07262328764633123 | 0.01963296740861028 | 0.08534198530940107 | 199 | 199 |
| 0.01 | held-out acid | R | technical SD | pH | 72 | 13 | 0.97 | 1.77 | 5 | 0.0428729467710334 | 0.029202540417121007 | 0.04210485005875358 | 199 | 199 |
| 0.01 | ferric transfer | R | technical SD | pH | 36 | 36 | 0.93 | 1.79 | 5 | 0.0650779657221574 | 0.026420038899021208 | 0.0818230906813631 | 199 | 199 |
| 0.1 | held-out acid | R&#42; | technical SD | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | R&#42; | technical SD | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | held-out acid | R | technical SD | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | R | technical SD | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | R&#42; | technical SD | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | R&#42; | technical SD | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | R | technical SD | pH | 72 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | R | technical SD | pH | 36 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.001 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 24 | 24 | 1.05 | 2.65 | 7 | 0.00474852518393956 | 0.0008300928657183912 | 0.008528766518111991 | 199 | 199 |
| 0.001 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 34 | 34 | 1.04 | 2.7 | 7 | -0.067990261902362 | -0.08522221649517271 | -0.052235566113874185 | 199 | 199 |
| 0.01 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 13 | 13 | 0.97 | 1.77 | 5 | -0.03824950059208452 | -0.04864351243718466 | -0.028922767207119556 | 199 | 199 |
| 0.01 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 36 | 36 | 0.93 | 1.79 | 5 | -0.07905925929285582 | -0.09871841616690852 | -0.06342626148544366 | 199 | 199 |
| 0.1 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired bias difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.001 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 24 | 24 | 1.05 | 2.65 | 7 | -0.06622233225037313 | -0.08331969592277516 | -0.04001389714892341 | 199 | 199 |
| 0.001 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 34 | 34 | 1.04 | 2.7 | 7 | -0.08758400233435826 | -0.10429508039799391 | -0.07324038127875515 | 199 | 199 |
| 0.01 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 13 | 13 | 0.97 | 1.77 | 5 | -0.090842601611143 | -0.11714969075629082 | -0.05703609386853521 | 199 | 199 |
| 0.01 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 36 | 36 | 0.93 | 1.79 | 5 | -0.09341103306399035 | -0.11481070201531479 | -0.07397606855055384 | 199 | 199 |
| 0.1 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired RMSE difference | pH | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.001 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 24 | 24 | 1.05 | 2.65 | 7 | 1.3832702770901018 | 1.1522884479244477 | 1.7125421352340728 | 199 | 199 |
| 0.001 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 34 | 34 | 1.04 | 2.7 | 7 | 0.8364205134622742 | 0.6909940936523243 | 0.9570758197838988 | 199 | 199 |
| 0.01 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 13 | 13 | 0.97 | 1.77 | 5 | 1.2221017686699795 | 0.9163051274408044 | 1.5694272424248137 | 199 | 199 |
| 0.01 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 36 | 36 | 0.93 | 1.79 | 5 | 1.1159428055324851 | 0.5611561096523325 | 2.3854339209598407 | 199 | 199 |
| 0.1 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 0.1 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | held-out acid | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
| 1.0 | ferric transfer | corrected&#95;<wbr>minus&#95;<wbr>raw | paired SD ratio | fold | 0 | 0 |  |  | 0 |  |  |  | 199 | 0 |
