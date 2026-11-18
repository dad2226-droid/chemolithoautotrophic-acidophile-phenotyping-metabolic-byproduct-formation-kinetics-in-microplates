# Supplementary Table S32. Direction of sulfur indicator changes after scatter correction

The sign of the net sulfur indicator change sometimes reverses after scatter correction. Counts identify the common observations and the nonzero comparisons used to assess sign changes.

**Figure panels** 8a–f · [All Supplementary tables](README.md)

**Column definitions**

- `experiment_id` — Identifier linking the sulfur-culture experiment to the microplate records and experimental metadata.
- `cohort` — Short label for the organism, medium, treatment series, and starting inoculum OD600 encoded by the experiment identifier.
- `plate_n` — Number of physical microplates summarized for the experiment.
- `postinitial_condition_acquisition_n` — Total number of postinitial condition measurements across plates, counting each treatment condition on each plate at each acquisition once.
- `common_valid_n` — Number of postinitial condition measurements valid for both ratios on the same three biotic wells, with eligible baselines and complete matched abiotic controls.
- `nonzero_common_n` — Number of jointly valid comparisons in which the mean net change was nonzero for both measured and scatter-corrected ratios. Net change is biotic change from baseline minus the matched abiotic change.
- `opposite_sign_n` — Total number of jointly valid, nonzero comparisons in which the two within-plate mean net changes had opposite signs.
- `opposite_sign_percent` — Mean of the within-plate percentages of jointly valid, nonzero comparisons with opposite signs. Each plate with an eligible denominator receives equal weight.
- `opposite_sign_percent_min` — Lowest observed opposite-sign percentage among physical plates with a jointly valid, nonzero comparison.
- `opposite_sign_percent_max` — Highest observed opposite-sign percentage among physical plates with a jointly valid, nonzero comparison.
- `raw_comparison_yield_percent` — Mean percentage of all postinitial condition measurements on each plate with an eligible net change for the measured ratio A545/A485, weighting physical plates equally.
- `corrected_comparison_yield_percent` — Mean percentage of all postinitial condition measurements on each plate with an eligible net change for the corrected ratio (A545 − A750)/(A485 − A750), weighting physical plates equally.

| experiment&#95;<wbr>id | cohort | plate&#95;<wbr>n | postinitial&#95;<wbr>condition&#95;<wbr>acquisition&#95;<wbr>n | common&#95;<wbr>valid&#95;<wbr>n | nonzero&#95;<wbr>common&#95;<wbr>n | opposite&#95;<wbr>sign&#95;<wbr>n | opposite&#95;<wbr>sign&#95;<wbr>percent | opposite&#95;<wbr>sign&#95;<wbr>percent&#95;<wbr>min | opposite&#95;<wbr>sign&#95;<wbr>percent&#95;<wbr>max | raw&#95;<wbr>comparison&#95;<wbr>yield&#95;<wbr>percent | corrected&#95;<wbr>comparison&#95;<wbr>yield&#95;<wbr>percent |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AF&#95;<wbr>AFM1&#95;<wbr>BS&#95;<wbr>potassium&#95;<wbr>tetrathionate&#95;<wbr>OD0p001 | AF AFM1 BS tetrathionate OD 10⁻³ | 3 | 370 | 370 | 370 | 82 | 22.564102564102566 | 13.846153846153847 | 30.0 | 100.0 | 100.0 |
| AF&#95;<wbr>AFM4&#95;<wbr>aluminum&#95;<wbr>sulfate&#95;<wbr>OD0p1 | AF AFM4 Al OD 10⁻¹ | 3 | 300 | 235 | 235 | 94 | 38.74951606658924 | 23.809523809523807 | 52.4390243902439 | 100.0 | 78.33333333333333 |
| AF&#95;<wbr>AFM4&#95;<wbr>cobalt&#95;<wbr>ii&#95;<wbr>sulfate&#95;<wbr>OD0p1 | AF AFM4 Co OD 10⁻¹ | 3 | 300 | 300 | 300 | 151 | 50.333333333333336 | 43.0 | 63.0 | 100.0 | 100.0 |
| AF&#95;<wbr>AFM4&#95;<wbr>copper&#95;<wbr>ii&#95;<wbr>sulfate&#95;<wbr>OD0p1 | AF AFM4 Cu OD 10⁻¹ | 3 | 300 | 178 | 178 | 58 | 32.60828625235405 | 28.333333333333332 | 35.59322033898305 | 100.0 | 59.333333333333336 |
| AF&#95;<wbr>AFM4&#95;<wbr>gallium&#95;<wbr>sulfate&#95;<wbr>OD0p1 | AF AFM4 Ga OD 10⁻¹ | 3 | 300 | 150 | 150 | 74 | 49.333333333333336 | 46.0 | 52.0 | 100.0 | 50.0 |
| AF&#95;<wbr>AFM4&#95;<wbr>potassium&#95;<wbr>sulfate&#95;<wbr>OD0p001 | AF AFM4 K OD 10⁻³ | 3 | 300 | 205 | 205 | 106 | 51.39282578703123 | 29.850746268656714 | 72.85714285714285 | 100.0 | 68.33333333333333 |
| AF&#95;<wbr>AFM4&#95;<wbr>potassium&#95;<wbr>sulfate&#95;<wbr>OD0p01 | AF AFM4 K OD 10⁻² | 3 | 300 | 195 | 195 | 65 | 33.766982014021714 | 27.397260273972602 | 38.983050847457626 | 100.0 | 65.0 |
| AF&#95;<wbr>AFM4&#95;<wbr>potassium&#95;<wbr>sulfate&#95;<wbr>OD0p1 | AF AFM4 K OD 10⁻¹ | 3 | 300 | 186 | 186 | 107 | 57.361111111111114 | 52.22222222222223 | 62.5 | 100.0 | 62.0 |
| AF&#95;<wbr>AFM4&#95;<wbr>zinc&#95;<wbr>sulfate&#95;<wbr>OD0p1 | AF AFM4 Zn OD 10⁻¹ | 3 | 300 | 228 | 228 | 66 | 29.628681177976954 | 18.181818181818183 | 50.70422535211267 | 100.0 | 76.0 |
| AF&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>potassium&#95;<wbr>tetrathionate&#95;<wbr>OD0p001 | AF SM4 BS tetrathionate OD 10⁻³ | 3 | 390 | 378 | 378 | 148 | 39.15343915343916 | 32.53968253968254 | 49.2063492063492 | 100.0 | 96.92307692307692 |
| AF&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>sodium&#95;<wbr>thiosulfate&#95;<wbr>OD0p001 | AF SM4 BS thiosulfate OD 10⁻³ | 3 | 390 | 390 | 390 | 136 | 34.87179487179487 | 33.07692307692307 | 36.15384615384615 | 100.0 | 100.0 |
| AF&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>sodium&#95;<wbr>thiosulfate&#95;<wbr>OD0p001&#95;<wbr>and&#95;<wbr>potassium&#95;<wbr>tetrathionate | AF SM4 BS thiosulfate OD 10⁻³ +tetrathionate | 3 | 390 | 228 | 228 | 81 | 35.526315789473685 | 31.57894736842105 | 38.15789473684211 | 100.0 | 58.461538461538474 |
| AF&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>sodium&#95;<wbr>thiosulfate&#95;<wbr>OD0p001&#95;<wbr>threefold&#95;<wbr>dye | AF SM4 BS thiosulfate OD 10⁻³ 3× dye | 3 | 390 | 390 | 390 | 129 | 33.07692307692307 | 19.230769230769234 | 53.07692307692308 | 100.0 | 100.0 |
| AF&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>sulfuric&#95;<wbr>acid&#95;<wbr>OD0p001&#95;<wbr>tetrathionate | AF SM4 BS acid OD 10⁻³ tetrathionate | 3 | 420 | 420 | 420 | 195 | 46.42857142857144 | 35.0 | 55.00000000000001 | 100.0 | 100.0 |
| AF&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>sulfuric&#95;<wbr>acid&#95;<wbr>OD0p001&#95;<wbr>thiosulfate | AF SM4 BS acid OD 10⁻³ thiosulfate | 3 | 420 | 420 | 420 | 156 | 37.142857142857146 | 29.28571428571429 | 47.85714285714286 | 100.0 | 100.0 |
| AT&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>sulfuric&#95;<wbr>acid&#95;<wbr>OD0p001&#95;<wbr>tetrathionate | AT SM4 BS acid OD 10⁻³ tetrathionate | 3 | 420 | 420 | 420 | 29 | 6.904761904761905 | 4.285714285714286 | 9.285714285714286 | 100.0 | 100.0 |
| AT&#95;<wbr>SM4&#95;<wbr>BS&#95;<wbr>sulfuric&#95;<wbr>acid&#95;<wbr>OD0p001&#95;<wbr>thiosulfate | AT SM4 BS acid OD 10⁻³ thiosulfate | 3 | 420 | 420 | 420 | 42 | 10.0 | 8.571428571428571 | 10.714285714285714 | 100.0 | 100.0 |
