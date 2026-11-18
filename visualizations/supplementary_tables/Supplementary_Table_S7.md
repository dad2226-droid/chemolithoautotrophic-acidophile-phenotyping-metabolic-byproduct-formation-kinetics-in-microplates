# Supplementary Table S7. Dye and scatter correction in the ferric standards

With-versus-without-dye comparisons give Holm-adjusted p = 1.00 for both A380 and scatter-corrected A380. The four planned contrasts use matched standard positions within physical plates and report effects in AU.

**Figure panels** 1h · [All Supplementary tables](README.md)

**Column definitions**

- `contrast` — Direction of the subtraction at matched ferric-standard positions, either A750-corrected minus measured A380 or absorbance with dye minus absorbance without dye.
- `context` — Dye condition held fixed for corrected-minus-measured contrasts, or readout held fixed for dyed-minus-undyed contrasts.
- `plate_n` — Number of physical microplates with a finite paired absorbance difference.
- `estimate` — Mean paired difference in absorbance units (AU), with standard-position differences averaged equally within each plate and plates weighted equally.
- `sd` — Sample standard deviation of the paired absorbance difference across physical microplates, in AU.
- `ci_low` — Lower limit of the two-sided 95% Student’s t confidence interval for the mean paired absorbance difference, in AU.
- `ci_high` — Upper limit of the two-sided 95% Student’s t confidence interval for the mean paired absorbance difference, in AU.
- `p_value` — Unadjusted two-sided p-value for a one-sample t test of whether the mean paired absorbance difference is zero.
- `p_holm` — P-value after Holm adjustment across the four ferric-standard contrasts.
- `multiplicity_family_n` — Number of contrasts with finite p-values in the Holm adjustment family (the four ferric-standard contrasts).
- `status` — Whether the plate-level mean and its interval could be estimated. A single contributing plate permits a mean only.

| contrast | context | plate&#95;<wbr>n | estimate | sd | ci&#95;<wbr>low | ci&#95;<wbr>high | p&#95;<wbr>value | p&#95;<wbr>holm | multiplicity&#95;<wbr>family&#95;<wbr>n | status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| A&#42;380 − A380 | without dye | 3 | -0.03457569444444447 | 0.008011732210281333 | -0.05447794056444916 | -0.014673448324439786 | 0.0174308150624613 | 0.0697232602498452 | 4 | estimated |
| A&#42;380 − A380 | with dye | 3 | -0.04268125000000001 | 0.014428024638365379 | -0.07852245011022638 | -0.006840049889773651 | 0.0360438364670359 | 0.10813150940110769 | 4 | estimated |
| dyed − undyed | A380 | 3 | 0.004008333333333334 | 0.017958907207574665 | -0.040604065322827704 | 0.048620731989494376 | 0.736317617177547 | 1.0 | 4 | estimated |
| dyed − undyed | A&#42;380 | 3 | -0.004097222222222208 | 0.024262590189358616 | -0.06436883749635153 | 0.05617439305190711 | 0.7974638459770749 | 1.0 | 4 | estimated |
