# Supplementary Table S8. Lid effects on indicator response

With the lid, measured and corrected pH ratios are lower. The corresponding ferric indicator contrasts have Holm-adjusted p = 0.140 and 0.193. Effects are lid-present minus lid-absent, expressed as a percentage of the paired reference span.

**Figure panels** 2a, 3a · [All Supplementary tables](README.md)

**Column definitions**

- `assay` — Identifies the thymol-blue absorbance ratio assay for pH or the ferric absorbance assay for iron.
- `dye` — Nominal concentration of added thymol blue, in micromoles per liter (µM).
- `readout` — The optical signal used in the calculation. R = A545/A485 and R&#42; = (A545 − A750)/(A485 − A750). A380 is absorbance at 380 nm and A&#42;380 = A380 − A750.
- `paired_well_n` — Total well pairs with eligible measurements of the same standard with and without the lid, summed across plates.
- `paired_level_n_min` — Smallest number of matched standard levels contributing to a single physical plate’s lid contrast.
- `paired_level_n_max` — Largest number of matched standard levels contributing to a single physical plate’s lid contrast.
- `plate_n` — Number of physical microplates with a finite normalized lid contrast.
- `estimate` — Mean across plates of 100 × (lid-present minus lid-absent response) / lid-absent response span. Within a plate, paired well differences are averaged within each standard, then equally across standards.
- `sd` — Sample standard deviation of the normalized lid contrast across physical microplates, as a percentage of the lid-absent response span.
- `ci_low` — Lower limit of the two-sided 95% Student’s t confidence interval for the mean normalized lid contrast, as a percentage of the lid-absent response span.
- `ci_high` — Upper limit of the two-sided 95% Student’s t confidence interval for the mean normalized lid contrast, as a percentage of the lid-absent response span.
- `p_value` — Unadjusted two-sided p-value for a one-sample t test of whether the mean normalized lid contrast is zero.
- `p_holm` — P-value after Holm adjustment across the four lid contrasts.
- `multiplicity_family_n` — Number of contrasts with finite p-values in the Holm adjustment family (the four lid contrasts).
- `status` — Whether the plate-level mean and its interval could be estimated. A single contributing plate permits a mean only.

| assay | dye | readout | paired&#95;<wbr>well&#95;<wbr>n | paired&#95;<wbr>level&#95;<wbr>n&#95;<wbr>min | paired&#95;<wbr>level&#95;<wbr>n&#95;<wbr>max | plate&#95;<wbr>n | estimate | sd | ci&#95;<wbr>low | ci&#95;<wbr>high | p&#95;<wbr>value | p&#95;<wbr>holm | multiplicity&#95;<wbr>family&#95;<wbr>n | status |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| iron | 0 µM | A380 | 130 | 11 | 12 | 3 | 2.3995944301804424 | 1.1614015739319252 | -0.4854870181100419 | 5.284675878470926 | 0.0699867537714602 | 0.1399735075429204 | 4 | estimated |
| iron | 0 µM | A&#42;380 | 130 | 11 | 12 | 3 | 1.2138875158733373 | 1.0874494930787462 | -1.4874867795073554 | 3.9152618112540303 | 0.19287095374401228 | 0.19287095374401228 | 4 | estimated |
| pH | 10 µM | R&#42; | 143 | 12 | 12 | 3 | -16.223470351097927 | 2.9173102914000513 | -23.470470862842138 | -8.976469839353715 | 0.010607263985828352 | 0.03182179195748505 | 4 | estimated |
| pH | 10 µM | R | 144 | 12 | 12 | 3 | -16.76176615283154 | 1.013661958545608 | -19.279842051021387 | -14.243690254641697 | 0.0012168364930186606 | 0.0048673459720746425 | 4 | estimated |
