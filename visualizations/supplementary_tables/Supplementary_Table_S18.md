# Supplementary Table S18. L. ferriphilum plate response and time-matched A. ferrooxidans benchmark

L. ferriphilum A380 remains close to the matched abiotic response, while the independent A. ferrooxidans benchmark shows a larger change at the matched duration.

**Figure panels** 7g · [All Supplementary tables](README.md)

**Column definitions**

- `LF_plate_context` — Treatment of the Leptospirillum ferriphilum (LF) microplates. HH identifies medium without an added compound. Ga and In identify the gallium and indium sulfate treatments.
- `LF_Ga_mM` — Ga concentration in the LF microplate treatment, in mM on a Ga cation basis.
- `LF_In_mM` — In concentration in the LF microplate treatment, in mM on an In cation basis.
- `LF_batch_n` — Number of distinct LF starting batches represented by eligible microplate endpoints.
- `LF_plate_n` — Number of LF microplates with both a biotic and a matched abiotic change in absorbance at 380 nm (A380) at the final scheduled acquisition.
- `LF_inoculum_culture_n` — Number of distinct inoculum cultures represented on the eligible LF microplates.
- `LF_biotic_mean_AU` — Mean endpoint-minus-initial A380 of the inoculated LF wells, averaged within each eligible microplate and then across microplates, in absorbance units (AU).
- `LF_abiotic_mean_AU` — Mean endpoint-minus-initial A380 of the matched uninoculated wells, averaged within each eligible LF microplate and then across microplates, in AU.
- `LF_net_change_AU` — Mean across LF microplates of biotic minus matched abiotic change in A380, in AU.
- `LF_net_between_microplate_sd_AU` — Sample standard deviation of the individual LF microplates’ biotic-minus-abiotic changes in A380, in AU.
- `LF_endpoint_min_d` — Earliest final acquisition among the eligible LF microplates in this treatment, in days from the initial acquisition.
- `LF_endpoint_max_d` — Latest final acquisition among the eligible LF microplates in this treatment, in days from the initial acquisition.
- `AF_reference_day_d` — Time in days at which the Acidithiobacillus ferrooxidans (AF) reference was evaluated, set to the latest LF acquisition across treatments.
- `AF_reference_mean_AU` — Mean biotic change in A380 for no-additive AF controls, in AU. Each plate is evaluated at the reference day using an exact reading or interpolation between adjacent readings, then plates are averaged within batch and batch means receive equal weight.
- `AF_reference_batch_sd_AU` — Sample standard deviation of the AF batch means at the reference day, in AU.
- `AF_batch_n` — Number of AF starting batches contributing at least one eligible microplate to the reference.
- `AF_plate_n` — Total number of AF control microplates with an eligible A380 estimate at the reference day.
- `LF_biotic_minus_AF_reference_AU` — LF mean biotic change at its endpoint minus the time-matched AF reference mean, in AU, using biotic changes before abiotic subtraction.
- `status` — Indicates that the comparison is descriptive because the LF treatment and AF reference differ in species, cohort, and medium.

| LF&#95;<wbr>plate&#95;<wbr>context | LF&#95;<wbr>Ga&#95;<wbr>mM | LF&#95;<wbr>In&#95;<wbr>mM | LF&#95;<wbr>batch&#95;<wbr>n | LF&#95;<wbr>plate&#95;<wbr>n | LF&#95;<wbr>inoculum&#95;<wbr>culture&#95;<wbr>n | LF&#95;<wbr>biotic&#95;<wbr>mean&#95;<wbr>AU | LF&#95;<wbr>abiotic&#95;<wbr>mean&#95;<wbr>AU | LF&#95;<wbr>net&#95;<wbr>change&#95;<wbr>AU | LF&#95;<wbr>net&#95;<wbr>between&#95;<wbr>microplate&#95;<wbr>sd&#95;<wbr>AU | LF&#95;<wbr>endpoint&#95;<wbr>min&#95;<wbr>d | LF&#95;<wbr>endpoint&#95;<wbr>max&#95;<wbr>d | AF&#95;<wbr>reference&#95;<wbr>day&#95;<wbr>d | AF&#95;<wbr>reference&#95;<wbr>mean&#95;<wbr>AU | AF&#95;<wbr>reference&#95;<wbr>batch&#95;<wbr>sd&#95;<wbr>AU | AF&#95;<wbr>batch&#95;<wbr>n | AF&#95;<wbr>plate&#95;<wbr>n | LF&#95;<wbr>biotic&#95;<wbr>minus&#95;<wbr>AF&#95;<wbr>reference&#95;<wbr>AU | status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| HH | 0.0 | 0.0 | 1 | 33 | 99 | 0.3333141414141414 | 0.33334141414141416 | -2.7272727272725952e-05 | 0.026774138193602075 | 5.308263888888889 | 5.365844907407407 | 5.365844907407407 | 2.950851180910824 | 0.696777748311818 | 8 | 342 | -2.617537039496683 | descriptive&#95;<wbr>species&#95;<wbr>cohort&#95;<wbr>and&#95;<wbr>medium&#95;<wbr>confounding |
| In | 0.0 | 6.22 | 1 | 3 | 9 | 0.3360222222222222 | 0.32215555555555553 | 0.013866666666666694 | 0.007943620780922948 | 5.330289351851852 | 5.350983796296297 | 5.365844907407407 | 2.950851180910824 | 0.696777748311818 | 8 | 342 | -2.614828958688602 | descriptive&#95;<wbr>species&#95;<wbr>cohort&#95;<wbr>and&#95;<wbr>medium&#95;<wbr>confounding |
| Ga | 142.8 | 0.0 | 1 | 3 | 9 | 0.2670111111111111 | 0.26812222222222226 | -0.0011111111111111553 | 0.004018061076541125 | 5.327743055555556 | 5.34 | 5.365844907407407 | 2.950851180910824 | 0.696777748311818 | 8 | 342 | -2.683840069799713 | descriptive&#95;<wbr>species&#95;<wbr>cohort&#95;<wbr>and&#95;<wbr>medium&#95;<wbr>confounding |
