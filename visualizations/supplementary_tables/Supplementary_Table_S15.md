# Supplementary Table S15. Gallium response across turbidity and fluorescence readouts

The same-batch AUC and timing responses differ between plate A600, corrected A600 and fluorescence-derived apparent planktonic OD600. The AUC effect is percentage points relative to the untreated response, with zero denoting no effect. Timing is in days.

**Figure panels** 6a–d · [All Supplementary tables](README.md)

**Column definitions**

- `comparison_quantity` — Whether the effect is based on area under the curve (AUC) through the first observed maximum or on time to that maximum, each relative to the same-batch no-Ga control.
- `format` — Whether the measurements were made in microplates or shake flasks.
- `readout` — Identifies measured plate absorbance at 600 nm, same-well scatter-corrected absorbance (A600 − A750), or fluorescence-derived apparent flask OD600.
- `Ga_mM` — Added gallium(III) sulfate concentration expressed as mM Ga³⁺, using the concentration calculated from the measured stock or the nominal concentration, as applicable.
- `unit` — The effect and its uncertainty are expressed in percentage points for control-relative AUC or days for the timing difference.
- `candidate_batch_n` — Number of starting batches represented at that dose and readout before requiring a finite estimate for the specified outcome.
- `batch_n` — Number of starting batches with a finite estimate included in the mean effect.
- `physical_unit_n` — Number of treated microplates or shake flasks with a finite control-relative estimate for the specified outcome.
- `effect` — Mean of the starting-batch treatment effects. For AUC the vessel effect is 100 × (treated AUC/mean control AUC − 1). For timing it is treated time to maximum minus mean control time to maximum, in days.
- `between_batch_sd` — Sample standard deviation of the starting-batch mean effects when at least two batches contribute, in the stated unit.
- `within_batch_unit_sd` — Sample standard deviation across treated vessel effects when exactly one batch contributes and at least two vessel estimates are available, in the stated unit.
- `ci_low` — Lower 95% Student’s t confidence limit for the mean effect across independent starting batches, calculated when at least three batches contribute.
- `ci_high` — Upper 95% Student’s t confidence limit for the same mean effect, in the stated unit.
- `final_read_maximum_n` — Number of treated vessel trajectories whose first observed maximum occurred at the final recorded acquisition, counted across candidate batches.
- `status` — Indicates whether the batch count supports a confidence interval or only descriptive reporting.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| comparison&#95;<wbr>quantity | format | readout | Ga&#95;<wbr>mM | unit | candidate&#95;<wbr>batch&#95;<wbr>n | batch&#95;<wbr>n | physical&#95;<wbr>unit&#95;<wbr>n | effect | between&#95;<wbr>batch&#95;<wbr>sd | within&#95;<wbr>batch&#95;<wbr>unit&#95;<wbr>sd | ci&#95;<wbr>low | ci&#95;<wbr>high | final&#95;<wbr>read&#95;<wbr>maximum&#95;<wbr>n | status |
| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | --- | ---: | --- | --- | ---: | --- |
| AUC relative to same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 0.001 | percentage points | 1 | 1 | 3 | -0.49630523855995534 |  | 19.235692486510082 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 0.01 | percentage points | 1 | 1 | 3 | 5.97951134719419 |  | 28.971175326888897 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 0.1 | percentage points | 1 | 1 | 3 | 2.3344973999699192 |  | 45.12561894516989 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 1.0 | percentage points | 3 | 1 | 3 | 15.952273708685581 |  | 44.82085126343332 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 10.0 | percentage points | 1 | 1 | 3 | 8.955759716946844 |  | 45.787521572052746 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 100.0 | percentage points | 3 | 1 | 2 | -23.051849449867376 |  | 14.635281756873825 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 1.9680000000000002 | percentage points | 1 | 1 | 3 | 16.798852684423352 |  | 29.715901702133515 |  |  | 1 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 6.04 | percentage points | 1 | 1 | 3 | 18.12764700151827 |  | 34.73790048098651 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 12.02 | percentage points | 1 | 1 | 3 | 16.848362821623212 |  | 59.27421658934422 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 19.68 | percentage points | 1 | 1 | 3 | 15.353171404894965 |  | 34.59226498831858 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 60.4 | percentage points | 1 | 1 | 3 | 149.93929527924158 |  | 188.72858365712858 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 120.12 | percentage points | 1 | 1 | 3 | 8.256696681185716 |  | 82.05034296242303 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 196.8 | percentage points | 1 | 1 | 3 | 56.69616253586895 |  | 116.08390244688245 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 304.0 | percentage points | 1 | 1 | 3 | -32.074481865965325 |  | 48.01219255990265 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 606.0 | percentage points | 1 | 1 | 3 | -99.90922497878059 |  | 0.15722694881017132 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 1.9680000000000002 | percentage points | 1 | 1 | 3 | 2.5532111807951168 |  | 3.013369940390256 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 6.04 | percentage points | 1 | 1 | 3 | 5.914224712785213 |  | 3.498093601076641 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 12.02 | percentage points | 1 | 1 | 3 | 24.551242836042817 |  | 28.251776455619694 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 19.68 | percentage points | 1 | 1 | 3 | 17.706942224375254 |  | 25.19425388216312 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 60.4 | percentage points | 1 | 1 | 3 | 22.430719927569456 |  | 86.54157789499223 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 120.12 | percentage points | 1 | 1 | 3 | -87.76071090203507 |  | 8.99385370042734 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 196.8 | percentage points | 1 | 1 | 3 | -95.58726773353554 |  | 3.378902311558602 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 304.0 | percentage points | 1 | 1 | 3 | -98.61453630829774 |  | 1.0458226795163634 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| AUC relative to same-batch control | plate | raw&#95;<wbr>A600 | 606.0 | percentage points | 1 | 1 | 3 | -100.01550221245246 |  | 0.01555177311966108 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 0.001 | d | 1 | 1 | 3 | 0.0 |  | 0.0 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 0.01 | d | 1 | 1 | 3 | 0.0 |  | 0.0 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 0.1 | d | 1 | 1 | 3 | 0.0 |  | 0.0 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 1.0 | d | 3 | 1 | 3 | -0.1638666666666667 |  | 0.28382539233362003 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 10.0 | d | 1 | 1 | 3 | -0.1638666666666667 |  | 0.28382539233362003 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | flask | flask&#95;<wbr>OD600&#95;<wbr>apparent | 100.0 | d | 3 | 1 | 2 | 0.0 |  | 0.0 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 1.9680000000000002 | d | 1 | 1 | 3 | 0.05952932098765432 |  | 0.10310780849069359 |  |  | 1 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 6.04 | d | 1 | 1 | 3 | -0.19120756172839495 |  | 0.6479361364078011 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 12.02 | d | 1 | 1 | 3 | -0.4876658950617281 |  | 1.3875945404129217 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 19.68 | d | 1 | 1 | 3 | -0.02417052469135787 |  | 0.7925803021703517 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 60.4 | d | 1 | 1 | 3 | 0.03610725308641971 |  | 2.448504668055446 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 120.12 | d | 1 | 1 | 3 | 0.7860455246913581 |  | 2.485187606223867 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 196.8 | d | 1 | 1 | 3 | 2.71125 |  | 2.348050893704533 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 304.0 | d | 1 | 1 | 3 | 2.71125 |  | 2.348050893704533 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | A600&#95;<wbr>minus&#95;<wbr>A750 | 606.0 | d | 1 | 1 | 1 | -2.5372222222222223 |  |  |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 1.9680000000000002 | d | 1 | 1 | 3 | 0.0 |  | 0.0 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 6.04 | d | 1 | 1 | 3 | 0.0 |  | 0.0 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 12.02 | d | 1 | 1 | 3 | 0.23584104938271602 |  | 0.4084886800412248 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 19.68 | d | 1 | 1 | 3 | 0.23584104938271602 |  | 0.4084886800412248 |  |  | 2 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 60.4 | d | 1 | 1 | 3 | 0.9303742283950616 |  | 1.6114554336329374 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 120.12 | d | 1 | 1 | 3 | 0.9303742283950616 |  | 1.6114554336329374 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 196.8 | d | 1 | 1 | 3 | 0.9303742283950616 |  | 1.6114554336329374 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 304.0 | d | 1 | 1 | 3 | 0.9303742283950616 |  | 1.6114554336329374 |  |  | 3 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
| time to maximum minus same-batch control | plate | raw&#95;<wbr>A600 | 606.0 | d | 1 | 1 | 2 | -4.256481481481481 |  | 1.9560996063073963 |  |  | 0 | descriptive&#95;<wbr>fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>batches |
