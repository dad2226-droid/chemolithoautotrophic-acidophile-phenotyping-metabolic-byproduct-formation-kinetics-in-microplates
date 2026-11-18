# Supplementary Table S11. Medium-specific indicator response and precision

Response span, within-standard SD, their ratio and adjacent-standard reversals vary by medium. AFM4 has too few corrected standards with replicates to estimate a span in the primary comparison.

**Figure panels** 4a–d · [All Supplementary tables](README.md)

**Column definitions**

- `medium` — Growth-medium code for the standards. St-BS and St-YS denote Starkey medium with brown and yellow elemental sulfur, respectively.
- `readout` — The optical signal used in the calculation. R = A545/A485 and R&#42; = (A545 − A750)/(A485 − A750). A380 is absorbance at 380 nm and A&#42;380 = A380 − A750.
- `plate_n` — Number of distinct physical microplates scheduled for the indicated medium and readout.
- `eligible_well_n` — Total wells with an eligible value for this readout across the physical plates, including wells at levels with only one eligible measurement.
- `scheduled_well_n` — Total wells assigned to the indicated medium and readout across the physical plates before eligibility filtering.
- `span_plate_n` — Number of plates with at least two standard levels that each have two or more eligible well responses, permitting a response span.
- `span_mean` — Arithmetic mean of the plate response spans. Each span is the highest minus the lowest replicated standard mean, in absorbance units for iron or ratio units for pH.
- `span_ci_low` — Approximate lower 95% confidence limit for the arithmetic mean plate span, using a t interval on log(mean) with a delta-method standard error, then transforming back to the readout units.
- `span_ci_high` — Approximate upper 95% confidence limit for the arithmetic mean plate span, calculated by the same log-scale method and expressed in the readout units.
- `sd_plate_n` — Number of plates with an estimable pooled within-standard technical standard deviation for this readout.
- `technical_sd_mean` — Arithmetic mean of the plate-level pooled technical standard deviations, in the readout units. Within each plate, well variances are pooled across replicated standards using their degrees of freedom.
- `sd_ci_low` — Approximate lower 95% confidence limit for the arithmetic mean pooled technical standard deviation, using the log-scale delta-method t interval in the readout units.
- `sd_ci_high` — Approximate upper 95% confidence limit for the arithmetic mean pooled technical standard deviation, using the log-scale delta-method t interval in the readout units.
- `span_sd_mean` — Mean of the individual plates’ response span divided by pooled technical standard deviation, a dimensionless signal-to-scatter ratio.
- `median_rank_rho` — Median plate-level Spearman rank correlation between reference value and standard mean response. Each correlation requires at least three replicated levels and variation in response.
- `plates_with_adjacent_reversal` — Number of plates with at least one increase in mean pH-indicator response as reference pH rises, or a decrease in mean iron response as the iron reference rises, between adjacent scheduled standards with eligible means.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| medium | readout | plate&#95;<wbr>n | eligible&#95;<wbr>well&#95;<wbr>n | scheduled&#95;<wbr>well&#95;<wbr>n | span&#95;<wbr>plate&#95;<wbr>n | span&#95;<wbr>mean | span&#95;<wbr>ci&#95;<wbr>low | span&#95;<wbr>ci&#95;<wbr>high | sd&#95;<wbr>plate&#95;<wbr>n | technical&#95;<wbr>sd&#95;<wbr>mean | sd&#95;<wbr>ci&#95;<wbr>low | sd&#95;<wbr>ci&#95;<wbr>high | span&#95;<wbr>sd&#95;<wbr>mean | median&#95;<wbr>rank&#95;<wbr>rho | plates&#95;<wbr>with&#95;<wbr>adjacent&#95;<wbr>reversal |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AFM1 | A380 | 3 | 72 | 72 | 3 | 2.702866666666667 | 2.6545946101810887 | 2.752016518740472 | 3 | 0.051678667185269776 | 0.029260147963207 | 0.09127379141773705 | 54.36846913366426 | 1.0 | 0 |
| AFM1 | A&#42;380 | 3 | 72 | 72 | 3 | 2.6641 | 2.4791291384732568 | 2.8628717640626293 | 3 | 0.06213491483296523 | 0.01971406184530518 | 0.19583724914707332 | 49.15739717421824 | 1.0 | 0 |
| AFM3 | A380 | 3 | 72 | 72 | 3 | 2.7220222222222223 | 2.650390164528458 | 2.7955902785316304 | 3 | 0.028207838395483822 | 0.008943384731094572 | 0.08896879323320227 | 112.52269229915184 | 1.0 | 0 |
| AFM3 | A&#42;380 | 3 | 72 | 72 | 3 | 2.7184333333333335 | 2.6472566534714295 | 2.7915237376350346 | 3 | 0.027775945851381722 | 0.008141216175195332 | 0.09476510036541055 | 117.62944633093444 | 1.0 | 0 |
| AFM4 | R&#42; | 3 | 10 | 72 | 0 |  |  |  | 3 | 0.15499827583487813 | 0.047130635672906936 | 0.5097420216972681 |  |  | 0 |
| AFM4 | R | 3 | 72 | 72 | 3 | 0.5541365130832491 | 0.48492203780956394 | 0.6332301920513081 | 3 | 0.017651886485436506 | 0.006224003610935455 | 0.05006248645924303 | 35.110530085116125 | -0.6666666666666669 | 3 |
| AFM5 | R&#42; | 3 | 72 | 72 | 3 | 6.9376849758649 | 6.7176319132899875 | 7.164946434340873 | 3 | 0.8176353619250053 | 0.4314402196290953 | 1.5495254143087556 | 8.949198949255138 | -0.9761904761904763 | 2 |
| AFM5 | R | 3 | 72 | 72 | 3 | 1.7747389999663745 | 1.6897026708898824 | 1.8640548850780103 | 3 | 0.14287259967107993 | 0.034597772847924015 | 0.5899969291808762 | 16.044411447960165 | -1.0 | 1 |
| F2S | A380 | 3 | 72 | 72 | 3 | 1.6530777777777779 | 1.1396576343311857 | 2.3977956686846538 | 3 | 0.15108375233430213 | 0.04981514440454848 | 0.4582200953597666 | 12.075919603507984 | 0.9285714285714287 | 3 |
| F2S | A&#42;380 | 3 | 72 | 72 | 3 | 1.6233333333333333 | 1.4186999944848575 | 1.8574829924264429 | 3 | 0.10124647384676967 | 0.029088763209290613 | 0.3523989106257572 | 20.542324183030594 | 0.9761904761904763 | 2 |
| F2S | R&#42; | 3 | 72 | 72 | 3 | 0.8320231123795665 | 0.5545977501897074 | 1.2482244280597665 | 3 | 0.23975824321446151 | 0.10943345529721343 | 0.52528739984644 | 3.7249119890998728 | -0.880952380952381 | 3 |
| F2S | R | 3 | 72 | 72 | 3 | 0.04916305181080416 | 0.036479962221689824 | 0.06625570631525324 | 3 | 0.013355548507338574 | 0.01102969420121267 | 0.01617186049566655 | 3.711318716513618 | -0.9047619047619048 | 3 |
| St-BS | R&#42; | 3 | 72 | 72 | 3 | 0.8130219544515859 | 0.7480137500756249 | 0.883679876678 | 3 | 0.06683228846068269 | 0.0437557617908075 | 0.10207923706702014 | 12.409115799129822 | -1.0 | 1 |
| St-BS | R | 3 | 72 | 72 | 3 | 0.06556964392309901 | 0.05420054166253874 | 0.07932352836934016 | 3 | 0.005282404925699317 | 0.004259066300352903 | 0.0065516241897291706 | 12.496908717558497 | -0.9761904761904763 | 2 |
| St-YS | R&#42; | 3 | 64 | 72 | 3 | 11.948721336884264 | 2.657940352689127 | 53.71525415988591 | 3 | 6.158902578086376 | 0.5578791556877084 | 67.99336483471838 | 2.6718797208201437 | -0.8285714285714287 | 3 |
| St-YS | R | 3 | 72 | 72 | 3 | 1.022370031266055 | 0.5190718110016612 | 2.013672210043418 | 3 | 0.27562911581040556 | 0.1764628394468256 | 0.4305235579376404 | 3.673709944064889 | -0.9285714285714287 | 3 |
