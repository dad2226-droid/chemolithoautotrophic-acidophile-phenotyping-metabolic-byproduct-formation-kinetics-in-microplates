# Supplementary Table S10. Edge–interior differences in storage drift

None of the 12 endpoint and time-averaged edge–interior contrasts reaches significance after Holm adjustment. Effects are percentages of the initial observed response span.

**Figure panels** 2d, 3f · [All Supplementary tables](README.md)

**Column definitions**

- `assay` — Identifies the thymol-blue absorbance ratio assay for pH or the ferric absorbance assay for iron.
- `dye` — Nominal concentration of added thymol blue, in micromoles per liter (µM).
- `readout` — The optical signal used in the calculation. R = A545/A485 and R&#42; = (A545 − A750)/(A485 − A750). A380 is absorbance at 380 nm and A&#42;380 = A380 − A750.
- `metric` — Final-acquisition contrast or time-averaged contrast. The latter is the trapezoidal integral of the plate contrast divided by the observed duration, using actual acquisition times.
- `mean_endpoint_d` — Mean elapsed time, in days, from the first to the last recorded storage acquisition across the physical plates.
- `maximum_absolute_observed_percent` — Largest absolute concentration-matched plate contrast at any recorded acquisition on any plate, expressed as a percentage of that plate’s initial response span.
- `plate_n` — Number of physical microplates with a finite edge-minus-interior drift contrast.
- `estimate` — Mean across plates of the edge-well change from baseline minus the interior-well change from baseline, averaged over matched standards and scaled to percent of each plate’s initial response span. Positive values indicate greater signed drift at the edge.
- `sd` — Sample standard deviation of the edge-minus-interior drift contrast across physical microplates, as a percentage of initial response span.
- `ci_low` — Lower limit of the two-sided 95% Student’s t confidence interval for the mean edge-minus-interior drift contrast, as a percentage of initial response span.
- `ci_high` — Upper limit of the two-sided 95% Student’s t confidence interval for the mean edge-minus-interior drift contrast, as a percentage of initial response span.
- `p_value` — Unadjusted two-sided p-value for a one-sample t test of whether the mean edge-minus-interior drift contrast is zero.
- `p_holm` — P-value after Holm adjustment across the 12 endpoint and time-average contrasts.
- `multiplicity_family_n` — Number of contrasts with finite p-values in the Holm adjustment family (the 12 endpoint and time-average contrasts).
- `status` — Whether the plate-level mean and its interval could be estimated. A single contributing plate permits a mean only.

| assay | dye | readout | metric | mean&#95;<wbr>endpoint&#95;<wbr>d | maximum&#95;<wbr>absolute&#95;<wbr>observed&#95;<wbr>percent | plate&#95;<wbr>n | estimate | sd | ci&#95;<wbr>low | ci&#95;<wbr>high | p&#95;<wbr>value | p&#95;<wbr>holm | multiplicity&#95;<wbr>family&#95;<wbr>n | status |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| iron | 0 µM | A380 | final acquisition | 15.868838734567902 | 0.7226251316831941 | 3 | -0.4897437979590924 | 0.17117818356472167 | -0.914973979181138 | -0.06451361673704675 | 0.03839290487831179 | 0.42232195366142966 | 12 | estimated |
| iron | 0 µM | A380 | time-average over observed window | 15.868838734567902 | 0.7226251316831941 | 3 | -0.2795771985642442 | 0.19044425680748847 | -0.7526669588859909 | 0.19351256175750248 | 0.1260769368960057 | 1.0 | 12 | estimated |
| iron | 0 µM | A&#42;380 | final acquisition | 15.868838734567902 | 1.1521763396135531 | 3 | -0.6936038804470605 | 0.43432227879931756 | -1.7725202322657863 | 0.38531247137166524 | 0.10962463257117483 | 1.0 | 12 | estimated |
| iron | 0 µM | A&#42;380 | time-average over observed window | 15.868838734567902 | 1.1521763396135531 | 3 | -0.4677321366590766 | 0.3020217703118264 | -1.217995806060281 | 0.28253153274212767 | 0.11541322207525581 | 1.0 | 12 | estimated |
| iron | 10 µM | A380 | final acquisition | 15.868838734567902 | 1.7075531292197663 | 3 | -0.9141096638871661 | 0.6194812347396212 | -2.452986360825518 | 0.6247670330511857 | 0.12501730081511123 | 1.0 | 12 | estimated |
| iron | 10 µM | A380 | time-average over observed window | 15.868838734567902 | 1.7075531292197663 | 3 | -0.6133602749526338 | 0.6424644813711838 | -2.209330521586909 | 0.9826099716816413 | 0.24002918826393405 | 1.0 | 12 | estimated |
| iron | 10 µM | A&#42;380 | final acquisition | 15.868838734567902 | 0.8455159517297404 | 3 | -0.6067611567688107 | 0.15030825110491697 | -0.9801475517257732 | -0.23337476181184824 | 0.019848462301958133 | 0.23818154762349758 | 12 | estimated |
| iron | 10 µM | A&#42;380 | time-average over observed window | 15.868838734567902 | 0.8455159517297404 | 3 | -0.3673535946078788 | 0.2861045437496691 | -1.0780766812395541 | 0.34336949202379646 | 0.1561647126803957 | 1.0 | 12 | estimated |
| pH | 10 µM | R&#42; | final acquisition | 15.122006172839505 | 371.7030442545471 | 3 | 99.2523122967686 | 170.0736451593973 | -323.2340434185342 | 521.7386680120715 | 0.41851608166866944 | 1.0 | 12 | estimated |
| pH | 10 µM | R&#42; | time-average over observed window | 15.122006172839505 | 371.7030442545471 | 3 | -16.098784931396434 | 33.72370418849816 | -99.87311028595728 | 67.6755404231644 | 0.49527464334856375 | 1.0 | 12 | estimated |
| pH | 10 µM | R | final acquisition | 15.122006172839505 | 4.41261385605109 | 3 | 1.3297627575003805 | 2.608368528016268 | -5.149783869087529 | 7.80930938408829 | 0.47037765967744416 | 1.0 | 12 | estimated |
| pH | 10 µM | R | time-average over observed window | 15.122006172839505 | 4.41261385605109 | 3 | 0.9661833538238147 | 2.625171460963812 | -5.555104072167101 | 7.48747077981473 | 0.5890575392713775 | 1.0 | 12 | estimated |
