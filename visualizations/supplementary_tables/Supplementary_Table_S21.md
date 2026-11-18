# Supplementary Table S21. Scatter correction and sulfur endpoint separation

Scatter correction does not consistently increase endpoint separation scaled by abiotic SD. Signed and absolute score changes use the paired plate–formulation observations.

**Figure panels** 8a–f · [All Supplementary tables](README.md)

**Column definitions**

- `species` — Inoculated organism, identified as AF for Acidithiobacillus ferrooxidans or AT for Acidithiobacillus thiooxidans.
- `medium` — Medium formulation identifier used to group microplate treatments for the species.
- `candidate_plate_formulation_pairs` — Number of physical microplate–formulation combinations considered for pairing measured and scatter-corrected endpoint ratios. A plate can contain several formulations.
- `paired_n` — Number of microplate–formulation combinations with finite signed and magnitude scores for both readouts at the final scheduled acquisition, using the same eligible wells.
- `unavailable_pair_n` — Number of candidate microplate–formulation combinations that lack a complete pair of endpoint scores.
- `batch_n` — Number of distinct experimental batches represented by the complete pairs.
- `inoculum_culture_n` — Number of distinct inoculum cultures represented by the complete pairs.
- `median_signed_snr_change` — Median corrected-minus-measured change in dimensionless signed technical signal-to-noise score across paired microplate formulations. Each well’s net biotic-minus-abiotic ratio change is divided by the square root of the sum of the initial and endpoint abiotic variances, then signed scores are averaged within each plate formulation.
- `signed_snr_q25` — 25th percentile of the paired corrected-minus-measured differences in signed signal-to-noise score.
- `signed_snr_q75` — 75th percentile of the paired corrected-minus-measured differences in signed signal-to-noise score.
- `median_magnitude_change` — Median corrected-minus-measured change in score magnitude across paired microplate formulations. For each readout, absolute well-level scores are averaged within the plate formulation before taking the difference.
- `magnitude_q25` — 25th percentile of the paired corrected-minus-measured differences in mean absolute signal-to-noise score.
- `magnitude_q75` — 75th percentile of the paired corrected-minus-measured differences in mean absolute signal-to-noise score.
- `magnitude_increased_n` — Number of complete pairs for which correction strictly increased the mean absolute score.
- `magnitude_increased_percent` — Percentage of all complete pairs for which correction strictly increased the mean absolute score.
- `direction_agreement_percent` — Percentage of complete pairs in which the measured and corrected signed plate scores have the same sign, counting two zero scores as agreement.
- `status` — Identifies descriptive paired comparisons across formulations, some of which share the same physical microplate.

| species | medium | candidate&#95;<wbr>plate&#95;<wbr>formulation&#95;<wbr>pairs | paired&#95;<wbr>n | unavailable&#95;<wbr>pair&#95;<wbr>n | batch&#95;<wbr>n | inoculum&#95;<wbr>culture&#95;<wbr>n | median&#95;<wbr>signed&#95;<wbr>snr&#95;<wbr>change | signed&#95;<wbr>snr&#95;<wbr>q25 | signed&#95;<wbr>snr&#95;<wbr>q75 | median&#95;<wbr>magnitude&#95;<wbr>change | magnitude&#95;<wbr>q25 | magnitude&#95;<wbr>q75 | magnitude&#95;<wbr>increased&#95;<wbr>n | magnitude&#95;<wbr>increased&#95;<wbr>percent | direction&#95;<wbr>agreement&#95;<wbr>percent | status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| AF | AFM1&#95;<wbr>BS | 30 | 30 | 0 | 1 | 9 | 0.07658231475250879 | -0.3407204982012812 | 0.46468670390655586 | -0.14686687476123075 | -0.48674959236295956 | 0.03734573182918024 | 8 | 26.666666666666668 | 80.0 | descriptive&#95;<wbr>paired&#95;<wbr>indicator&#95;<wbr>effect&#95;<wbr>shared&#95;<wbr>microplates&#95;<wbr>across&#95;<wbr>conditions |
| AF | AFM4 | 240 | 183 | 57 | 1 | 69 | -1.7720468365496571 | -7.776170495649572 | -0.411713666792615 | -0.9199713307740844 | -4.759856532399518 | 0.3957181375154343 | 60 | 32.78688524590164 | 67.75956284153006 | descriptive&#95;<wbr>paired&#95;<wbr>indicator&#95;<wbr>effect&#95;<wbr>shared&#95;<wbr>microplates&#95;<wbr>across&#95;<wbr>conditions |
| AF | SM4&#95;<wbr>BS | 180 | 159 | 21 | 2 | 54 | 0.029275660299810824 | -0.3838209142354618 | 0.7079880478513716 | 0.12146717184911193 | -0.3039363194509439 | 0.5729687795731367 | 91 | 57.23270440251572 | 66.0377358490566 | descriptive&#95;<wbr>paired&#95;<wbr>indicator&#95;<wbr>effect&#95;<wbr>shared&#95;<wbr>microplates&#95;<wbr>across&#95;<wbr>conditions |
| AT | SM4&#95;<wbr>BS | 60 | 60 | 0 | 1 | 18 | -0.8532566844090996 | -3.5158002591938198 | 0.22886985074923713 | 1.111648543723847 | -0.1815944070769199 | 3.3859271498712133 | 42 | 70.0 | 91.66666666666666 | descriptive&#95;<wbr>paired&#95;<wbr>indicator&#95;<wbr>effect&#95;<wbr>shared&#95;<wbr>microplates&#95;<wbr>across&#95;<wbr>conditions |
