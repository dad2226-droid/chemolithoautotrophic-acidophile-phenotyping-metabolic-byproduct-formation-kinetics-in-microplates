# Supplementary Table S14. Batch repeatability of integrated response and timing

The batch distributions quantify baseline consistency and show how AUC and time to first observed maximum vary among independent batches. Native AUC units are specific to each readout.

**Figure panels** 5a–d · [All Supplementary tables](README.md)

**Column definitions**

- `metric` — Outcome averaged within each starting batch, either the signed area under the response curve (AUC) to its first observed maximum or the elapsed time to that maximum.
- `cohort` — Cultivation format and added Ga concentration in mM on a cation basis. Zero identifies no-additive AFM1 controls.
- `native_unit` — Unit of the untransformed outcome, expressed in absorbance-unit days for plate AUC, mM days of apparent ferric iron for flask AUC, or days for elapsed time.
- `batch_n` — Number of starting batches with a finite mean outcome for the cohort.
- `vessel_n` — Number of distinct physical microplates or shake flasks recorded in the cohort’s batches.
- `native_mean` — Arithmetic mean of the batch means in the native unit, weighting starting batches equally.
- `native_sd` — Sample standard deviation of the batch means in the native unit.
- `cv_percent` — Between-batch coefficient of variation, calculated as 100 times the sample standard deviation divided by the absolute value of the mean.
- `log_outcome_sd` — Sample standard deviation of the natural logarithms of positive batch means.
- `centered_deviation_iqr_pp` — Difference between the 75th and 25th percentiles of the batch deviations from the cohort mean, in percentage points. Each deviation is 100 × (batch mean/cohort mean − 1).
- `final_read_maximum_vessels` — Number of physical vessels whose first observed maximum occurred at the final acquisition in at least one contributing trajectory.
- `status` — Identifies the row as a descriptive summary of dispersion across starting batches.

| metric | cohort | native&#95;<wbr>unit | batch&#95;<wbr>n | vessel&#95;<wbr>n | native&#95;<wbr>mean | native&#95;<wbr>sd | cv&#95;<wbr>percent | log&#95;<wbr>outcome&#95;<wbr>sd | centered&#95;<wbr>deviation&#95;<wbr>iqr&#95;<wbr>pp | final&#95;<wbr>read&#95;<wbr>maximum&#95;<wbr>vessels | status |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | plate 0 | AU d | 8 | 345 | 5.981685834115897 | 2.269590108231389 | 37.94231544704384 | 0.4757233476771573 | 48.73421895245983 | 32 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | flask 0 | mM d | 5 | 15 | 105.28904566666665 | 80.8944032538133 | 76.83078780096074 | 0.7270758499643435 | 90.6290561654725 | 5 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | flask 1 | mM d | 3 | 9 | 111.39657722222222 | 97.48665673002427 | 87.51315270266392 | 0.8421372652407748 | 80.38693952690217 | 3 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | flask 100 | mM d | 3 | 9 | 64.59938722222222 | 28.116982029934235 | 43.5251528520103 | 0.45167514490605604 | 43.32532325276528 | 4 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | plate 0 | d | 8 | 345 | 4.6515940835466 | 0.7551133312736825 | 16.23343132937231 | 0.15928227436260883 | 16.147535237819707 | 32 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | flask 0 | d | 5 | 15 | 4.443606666666667 | 0.9348850337519227 | 21.03887908812637 | 0.2130284613901886 | 36.708169489948844 | 5 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | flask 1 | d | 3 | 9 | 4.220155555555555 | 1.1607278158950933 | 27.504384627885862 | 0.25909867368637085 | 24.02162637501514 | 3 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | flask 100 | d | 3 | 9 | 3.8881666666666668 | 0.40196945295769765 | 10.338277327558771 | 0.10222695157480821 | 10.180033434780722 | 4 | descriptive&#95;<wbr>between&#95;<wbr>starting&#95;<wbr>batch&#95;<wbr>dispersion |
