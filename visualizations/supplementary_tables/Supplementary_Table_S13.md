# Supplementary Table S13. Between-batch variance in plates and flasks

Both plate–flask comparisons give Holm-adjusted p = 0.617. All four estimable paired flask comparisons give p = 1.00 after adjustment. The independent-format and paired-flask tests have separate families of two and six planned comparisons, respectively.

**Figure panels** 5c–d · [All Supplementary tables](README.md)

**Column definitions**

- `metric` — Outcome averaged within each starting batch, either the signed area under the response curve (AUC) to its first observed maximum or the elapsed time to that maximum.
- `comparison` — AFM1 cohorts in numerator/denominator order. The number following plate or flask is the added Ga concentration in mM on a cation basis.
- `method` — The two-sided test compares between-batch variances of log-transformed outcomes. It uses an independent-samples F test for plates versus flasks and a paired Pitman–Morgan test for two flask doses.
- `log_outcome_variance_ratio` — Variance of the natural-log batch means in the numerator cohort divided by that in the denominator cohort. Values above 1 indicate greater dispersion in the numerator cohort.
- `ci_low` — Lower unadjusted 95% confidence limit for the variance ratio, from the F distribution for independent batches or inversion of the Pitman–Morgan test for paired batches.
- `ci_high` — Upper unadjusted 95% confidence limit for the same variance ratio.
- `p_value` — Unadjusted two-sided p-value for the hypothesis that the two variances of log batch means are equal.
- `p_holm` — P-value adjusted by Holm’s method within the independent-format or paired-flask comparison family.
- `status` — Whether the variance comparison could be estimated, with the reason given when too few matched batches were available.
- `numerator_batch_n` — Number of eligible starting batches in the numerator cohort, restricted to shared batches for paired flask comparisons.
- `denominator_batch_n` — Number of eligible starting batches in the denominator cohort, restricted to shared batches for paired flask comparisons.
- `numerator_df` — Numerator degrees of freedom for the independent F test, equal to the numerator batch count minus one.
- `denominator_df` — Denominator degrees of freedom for the independent F test, equal to the denominator batch count minus one.
- `family_planned_n` — Number of planned tests in the Holm adjustment, including unavailable comparisons. The adjustment covers two plate–flask tests or six paired flask-dose tests.
- `matched_batch_n` — Number of starting batches with eligible outcomes under both flask doses in a paired comparison.
- `statistic` — Signed t statistic for the paired Pitman–Morgan test. Positive values indicate greater log-outcome variance in the numerator cohort.
- `paired_test_df` — Degrees of freedom for the Pitman–Morgan test, equal to the matched batch count minus two when at least three pairs are available.

Blank numeric cells indicate quantities that were unavailable or inapplicable for that row.

| metric | comparison | method | log&#95;<wbr>outcome&#95;<wbr>variance&#95;<wbr>ratio | ci&#95;<wbr>low | ci&#95;<wbr>high | p&#95;<wbr>value | p&#95;<wbr>holm | status | numerator&#95;<wbr>batch&#95;<wbr>n | denominator&#95;<wbr>batch&#95;<wbr>n | numerator&#95;<wbr>df | denominator&#95;<wbr>df | family&#95;<wbr>planned&#95;<wbr>n | matched&#95;<wbr>batch&#95;<wbr>n | statistic | paired&#95;<wbr>test&#95;<wbr>df |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | AFM1 plate 0 / AFM1 flask 0 | two-sided F test of natural-log native batch variance | 0.4281042047396953 | 0.04717848249292058 | 2.364245860298256 | 0.3086172201706692 | 0.6172344403413385 | estimated | 8 | 5 | 7 | 4 | 2 |  |  |  |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | AFM1 plate 0 / AFM1 flask 0 | two-sided F test of natural-log native batch variance | 0.5590615641104878 | 0.0616104114905597 | 3.0874702326359347 | 0.4710321618161222 | 0.6172344403413385 | estimated | 8 | 5 | 7 | 4 | 2 |  |  |  |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | AFM1 flask 0 / AFM1 flask 1 | two-sided paired Pitman-Morgan test of natural-log native batch variance | 1.0301417537455513 | 0.056847909233108146 | 18.667213044871723 | 0.9404463568312986 | 1.0 | estimated | 3 | 3 |  |  | 6 | 3 | 0.09382047726809854 | 1.0 |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | AFM1 flask 0 / AFM1 flask 100 | two-sided paired Pitman-Morgan test of natural-log native batch variance | 1.8020540726948162 | 0.0062068481408210485 | 523.1961226115006 | 0.7323501720136986 | 1.0 | estimated | 3 | 3 |  |  | 6 | 3 | 0.4470804397685313 | 1.0 |
| AUC&#95;<wbr>to&#95;<wbr>own&#95;<wbr>max | AFM1 flask 1 / AFM1 flask 100 | two-sided paired Pitman-Morgan test of natural-log native batch variance |  |  |  |  |  | fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>matched&#95;<wbr>starting&#95;<wbr>batches | 1 | 1 |  |  | 6 | 1 |  |  |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | AFM1 flask 0 / AFM1 flask 1 | two-sided paired Pitman-Morgan test of natural-log native batch variance | 1.0463262459508842 | 0.40166443812755637 | 2.7256548228897364 | 0.6660294418181669 | 1.0 | estimated | 3 | 3 |  |  | 6 | 3 | 0.5786856419480162 | 1.0 |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | AFM1 flask 0 / AFM1 flask 100 | two-sided paired Pitman-Morgan test of natural-log native batch variance | 3.5108971865830854 | 0.06672564735895155 | 184.73255101515898 | 0.2520116250663929 | 1.0 | estimated | 3 | 3 |  |  | 6 | 3 | 2.3928000178798046 | 1.0 |
| time&#95;<wbr>to&#95;<wbr>observed&#95;<wbr>max | AFM1 flask 1 / AFM1 flask 100 | two-sided paired Pitman-Morgan test of natural-log native batch variance |  |  |  |  |  | fewer&#95;<wbr>than&#95;<wbr>three&#95;<wbr>matched&#95;<wbr>starting&#95;<wbr>batches | 1 | 1 |  |  | 6 | 1 |  |  |
