# Supplementary Table S24. Inverse ferric concentration predictions

Near saturation, inverse concentration predictions show larger errors or cannot be evaluated. Rows report the attempted and successful predictions for each withheld reference window.

**Figure panels** 1f · [All Supplementary tables](README.md)

**Column definitions**

- `acid_M` — Nominal sulfuric acid concentration (M) in the ferric standards used for this acid-specific calibration.
- `validation` — Identifies whether cross-validation withheld one physical microplate at a time or withheld one of the ten interior dilution levels across all microplates at once.
- `nominal_reference_window_mM` — Nominal total-iron concentration window (mM) defined by dilution rank. “all” includes all attempted levels in the indicated scheme.
- `reference_min_mM` — Lowest measured total-iron reference concentration (mM) among the attempted predictions in this row.
- `reference_max_mM` — Highest measured total-iron reference concentration (mM) among the attempted predictions in this row.
- `plate_n` — Number of distinct physical microplates represented among the withheld observations in this row.
- `attempted_n` — Number of withheld plate-by-standard mean absorbances for which a total-iron concentration prediction was attempted.
- `predicted_n` — Number of attempts yielding a finite concentration prediction within the acid-specific training range.
- `failed_n` — Number of attempts for which an eligible concentration prediction was unavailable.
- `outside_reference_support_n` — Number of attempts whose measured total-iron reference lay outside the concentration range used to fit the model for that acid.
- `response_or_inverse_failure_n` — Number of attempts with a reference inside the training range but no valid inverse prediction, including missing responses, fit failures, and responses outside the fitted response range.
- `converged_unidentified_fit_n` — Number of withheld observations whose training fit converged but failed the model’s parameter-identification checks. An observation may still yield a prediction.
- `bias_mM` — Mean signed concentration error (predicted minus measured total Fe, mM) over successful predictions.
- `rmse_mM` — Root mean squared concentration error (mM) over successful predictions.
- `median_absolute_error_mM` — Median absolute difference between predicted and measured total-iron concentration (mM) over successful predictions.
- `median_absolute_relative_error_percent` — Median absolute concentration error divided by the measured concentration and expressed as a percentage, using successful predictions with a positive reference.
- `relative_error_n` — Number of successful predictions with a positive measured reference, the denominator for the relative-error summary.
- `plate_rmse_min_mM` — Smallest plate-specific root mean squared concentration error (mM) among plates with at least one successful prediction.
- `plate_rmse_max_mM` — Largest plate-specific root mean squared concentration error (mM) among plates with at least one successful prediction.

| acid&#95;<wbr>M | validation | nominal&#95;<wbr>reference&#95;<wbr>window&#95;<wbr>mM | reference&#95;<wbr>min&#95;<wbr>mM | reference&#95;<wbr>max&#95;<wbr>mM | plate&#95;<wbr>n | attempted&#95;<wbr>n | predicted&#95;<wbr>n | failed&#95;<wbr>n | outside&#95;<wbr>reference&#95;<wbr>support&#95;<wbr>n | response&#95;<wbr>or&#95;<wbr>inverse&#95;<wbr>failure&#95;<wbr>n | converged&#95;<wbr>unidentified&#95;<wbr>fit&#95;<wbr>n | bias&#95;<wbr>mM | rmse&#95;<wbr>mM | median&#95;<wbr>absolute&#95;<wbr>error&#95;<wbr>mM | median&#95;<wbr>absolute&#95;<wbr>relative&#95;<wbr>error&#95;<wbr>percent | relative&#95;<wbr>error&#95;<wbr>n | plate&#95;<wbr>rmse&#95;<wbr>min&#95;<wbr>mM | plate&#95;<wbr>rmse&#95;<wbr>max&#95;<wbr>mM |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.001 | held-out interior level | all | 0.296 | 600.0 | 3 | 30 | 30 | 0 | 0 | 0 | 0 | -53.087685829087064 | 139.42059437113596 | 1.2327766741641608 | 17.277776785332346 | 30 | 133.15161825266384 | 144.70635670423727 |
| 0.001 | held-out interior level | 0–1 | 0.296 | 1.0 | 3 | 6 | 6 | 0 | 0 | 0 | 0 | 0.09956842338520376 | 0.1714186094843627 | 0.16557655044686947 | 35.25048624354019 | 6 | 0.13000408670062222 | 0.20907101318249496 |
| 0.001 | held-out interior level | &gt;1–10 | 2.96 | 10.0 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | -0.05948633200628143 | 1.2631964666673505 | 0.9796177460187545 | 16.85915858699271 | 9 | 0.6416113002212959 | 1.5015275441660572 |
| 0.001 | held-out interior level | &gt;10–100 | 29.6 | 100.0 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.6393858051279536 | 5.641789007688664 | 3.157357125970208 | 4.886017480398846 | 9 | 3.2098102941520112 | 8.079834325572982 |
| 0.001 | held-out interior level | &gt;100–1000 | 296.0 | 600.0 | 3 | 6 | 6 | 0 | 0 | 0 | 0 | -266.407846778503 | 311.6734559067927 | 274.61565238997497 | 57.488498337903195 | 6 | 297.5661553383775 | 323.5440999480076 |
| 0.01 | held-out interior level | all | 0.296 | 599.0 | 3 | 30 | 30 | 0 | 0 | 0 | 0 | -51.274265109425095 | 135.68282577917984 | 2.0178595776739736 | 16.135719470257875 | 30 | 131.61766371257497 | 139.65138994989047 |
| 0.01 | held-out interior level | 0–1 | 0.296 | 0.998 | 3 | 6 | 6 | 0 | 0 | 0 | 0 | 0.014153016587793119 | 0.14186902331021295 | 0.13094839909894807 | 22.78873790090668 | 6 | 0.1311154802394081 | 0.15050913002438734 |
| 0.01 | held-out interior level | &gt;1–10 | 2.96 | 9.979999999999999 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.038061531581419526 | 0.9651838178978936 | 0.49849319113333346 | 11.721692119275017 | 9 | 0.43726410363228735 | 1.370972862185842 |
| 0.01 | held-out interior level | &gt;10–100 | 29.6 | 99.8 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 2.1303791767211364 | 9.496645272522192 | 8.691276171936295 | 14.706891000353659 | 9 | 7.631947838364508 | 12.164813983133195 |
| 0.01 | held-out interior level | &gt;100–1000 | 296.0 | 599.0 | 3 | 6 | 6 | 0 | 0 | 0 | 0 | -259.6381396261672 | 303.17066011315285 | 260.13314987723254 | 53.40217137127227 | 6 | 294.1527503235668 | 311.91260140211483 |
| 0.1 | held-out interior level | all | 0.299 | 600.0 | 3 | 30 | 28 | 2 | 0 | 2 | 0 | -47.51590423190136 | 125.48011690531389 | 1.4024152658601778 | 15.016829995475621 | 28 | 118.64275070025214 | 137.1472053436899 |
| 0.1 | held-out interior level | 0–1 | 0.299 | 0.9990000000000001 | 3 | 6 | 5 | 1 | 0 | 1 | 0 | -0.02987298368332216 | 0.1313832003738914 | 0.15451326844326832 | 16.96948768291417 | 5 | 0.032733019939934715 | 0.16210433542566918 |
| 0.1 | held-out interior level | &gt;1–10 | 2.99 | 9.99 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.28360090045185066 | 0.8733638614450085 | 0.7517932932075171 | 13.170204131193834 | 9 | 0.8442246054313451 | 0.9177627650314698 |
| 0.1 | held-out interior level | &gt;10–100 | 29.9 | 99.9 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.8435812298790913 | 9.318828533408189 | 7.297357427034399 | 11.283923472373457 | 9 | 5.351701696720628 | 11.4786039941737 |
| 0.1 | held-out interior level | &gt;100–1000 | 299.0 | 600.0 | 3 | 6 | 5 | 1 | 0 | 1 | 0 | -268.08811854956 | 296.6744870747204 | 355.8044516432632 | 59.39974150972674 | 5 | 268.8821500054577 | 355.8044516432632 |
| 1.0 | held-out interior level | all | 0.299 | 612.0 | 3 | 30 | 29 | 1 | 0 | 1 | 0 | -52.45915333410919 | 132.72690223099065 | 1.0401100749808467 | 13.273439027884805 | 29 | 111.89199043151203 | 144.90953781147803 |
| 1.0 | held-out interior level | 0–1 | 0.299 | 1.02 | 3 | 6 | 5 | 1 | 0 | 1 | 0 | -0.15433283260101657 | 0.21323204701731308 | 0.12634733824024247 | 12.67275208026504 | 5 | 0.08993670936938511 | 0.2987548997762361 |
| 1.0 | held-out interior level | &gt;1–10 | 2.99 | 10.200000000000001 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.3706721292698869 | 0.7075038079509965 | 0.7105389639933435 | 10.432397943639385 | 9 | 0.521113511659687 | 0.8448725607537358 |
| 1.0 | held-out interior level | &gt;10–100 | 29.9 | 102.0 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 1.9105569517881835 | 11.329043905135872 | 7.801656365388624 | 13.273439027884805 | 9 | 5.012771590783576 | 15.403889989138634 |
| 1.0 | held-out interior level | &gt;100–1000 | 299.0 | 612.0 | 3 | 6 | 6 | 0 | 0 | 0 | 0 | -256.84580737594735 | 291.4665873311208 | 241.9874736481974 | 51.978798247151694 | 6 | 249.48364489552037 | 313.3053214850426 |
| 0.001 | held-out plate | all | 0.0 | 1000.0 | 3 | 36 | 33 | 3 | 1 | 2 | 0 | -67.56589211714079 | 186.76280035839912 | 0.902319600084823 | 18.068382008789975 | 30 | 114.29767254858872 | 263.83864496064393 |
| 0.001 | held-out plate | 0–1 | 0.0 | 1.0 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.07499660391690172 | 0.17372709247495366 | 0.10178477613816889 | 28.113681884185795 | 6 | 0.04095919042785725 | 0.24392804484794037 |
| 0.001 | held-out plate | &gt;1–10 | 2.96 | 10.0 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | -0.07365663048571985 | 1.1589006231729078 | 0.8893848908468485 | 12.463080633477865 | 9 | 0.6785278688071327 | 1.336768975169563 |
| 0.001 | held-out plate | &gt;10–100 | 29.6 | 100.0 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.793780477502518 | 8.15958136323497 | 3.577197019964217 | 5.3376665746383365 | 9 | 4.940225193593831 | 12.207790329293896 |
| 0.001 | held-out plate | &gt;100–1000 | 296.0 | 1000.0 | 3 | 9 | 6 | 3 | 1 | 2 | 0 | -372.80508732067483 | 437.881218811382 | 391.47821020230396 | 65.60130105900346 | 6 | 304.1393636903018 | 527.6524140203733 |
| 0.01 | held-out plate | all | 0.0 | 998.0 | 3 | 36 | 31 | 5 | 1 | 4 | 0 | -70.3831912118621 | 187.97337037059296 | 2.7669546703143784 | 21.887299635880808 | 31 | 120.30784198964793 | 267.7152171681691 |
| 0.01 | held-out plate | 0–1 | 0.0 | 0.998 | 3 | 9 | 6 | 3 | 0 | 3 | 0 | 0.014682622929627484 | 0.15682851586280155 | 0.13057519541393675 | 23.27902546330674 | 6 | 0.11056796715211192 | 0.1887419460921036 |
| 0.01 | held-out plate | &gt;1–10 | 2.96 | 9.979999999999999 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.13103711088090333 | 1.3960461757338491 | 0.6626431712108642 | 18.53399469704924 | 9 | 0.6361090332475016 | 1.8149882117220923 |
| 0.01 | held-out plate | &gt;10–100 | 29.6 | 99.8 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 1.6746043233077375 | 12.304684769101225 | 8.198854530060927 | 19.791656402081863 | 9 | 5.621080060685811 | 16.32135265095458 |
| 0.01 | held-out plate | &gt;100–1000 | 296.0 | 998.0 | 3 | 9 | 7 | 2 | 1 | 1 | 0 | -314.0311137447144 | 395.32505459649593 | 372.85840112776566 | 62.35090319862302 | 7 | 268.57134610447883 | 512.3736216135604 |
| 0.1 | held-out plate | all | 0.0 | 999.0 | 3 | 36 | 27 | 9 | 1 | 8 | 0 | -34.38129572330714 | 106.38146271338645 | 2.069958599648765 | 19.679468384978517 | 27 | 4.85942329249197 | 140.0579245650392 |
| 0.1 | held-out plate | 0–1 | 0.0 | 0.9990000000000001 | 3 | 9 | 5 | 4 | 0 | 4 | 0 | -0.028665102541473142 | 0.17646763001269014 | 0.09354641415790443 | 9.373388192174792 | 5 | 0.03942227184519875 | 0.20876105380334173 |
| 0.1 | held-out plate | &gt;1–10 | 2.99 | 9.99 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.30064420683871595 | 1.236517898080718 | 1.0607423648951926 | 18.339555314284677 | 9 | 1.024803081305266 | 1.3638462123894703 |
| 0.1 | held-out plate | &gt;10–100 | 29.9 | 99.9 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 1.6808894798756977 | 14.097788142457237 | 11.468692923549895 | 19.679468384978517 | 9 | 7.867534687287706 | 17.225973796441824 |
| 0.1 | held-out plate | &gt;100–1000 | 299.0 | 999.0 | 3 | 9 | 4 | 5 | 1 | 4 | 0 | -236.49636554925382 | 275.57066982751945 | 238.29638435526925 | 52.44500279698603 | 4 | 232.90601522247584 | 312.46307975753444 |
| 1.0 | held-out plate | all | 0.0 | 1020.0 | 3 | 36 | 28 | 8 | 1 | 7 | 0 | -39.497194377258914 | 120.30461726371122 | 1.2484657228131781 | 18.282714499405195 | 28 | 13.26896550206942 | 149.0901337380691 |
| 1.0 | held-out plate | 0–1 | 0.0 | 1.02 | 3 | 9 | 6 | 3 | 0 | 3 | 0 | -0.15295799837832863 | 0.2062186128013295 | 0.08439763770437869 | 15.862581684662336 | 6 | 0.05718601396732489 | 0.301593859491149 |
| 1.0 | held-out plate | &gt;1–10 | 2.99 | 10.200000000000001 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 0.34119106933375654 | 0.9386679899338455 | 0.8089422862804208 | 13.482371438007013 | 9 | 0.7599842709904707 | 1.115099825378574 |
| 1.0 | held-out plate | &gt;10–100 | 29.9 | 102.0 | 3 | 9 | 9 | 0 | 0 | 0 | 0 | 2.3878274209252672 | 16.333357231053437 | 7.743715297606187 | 22.113302273195508 | 9 | 5.515194672620979 | 21.639367524789957 |
| 1.0 | held-out plate | &gt;100–1000 | 299.0 | 1020.0 | 3 | 9 | 4 | 5 | 1 | 4 | 0 | -282.3912152463277 | 317.3485599928393 | 285.0794889842415 | 60.37134904313401 | 4 | 301.22198732096143 | 332.69435133988674 |
