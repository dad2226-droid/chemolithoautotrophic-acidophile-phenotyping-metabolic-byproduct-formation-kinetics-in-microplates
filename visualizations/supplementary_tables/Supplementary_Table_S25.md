# Supplementary Table S25. Withheld prediction by shared and acid-specific ferric models

Acid-specific fits do not consistently reduce prediction error at withheld microplates or concentrations. RMSE and bias describe absorbance in AU.

**Figure panels** 1f · [All Supplementary tables](README.md)

**Column definitions**

- `acid_M` — Nominal sulfuric acid concentration (M) of the withheld ferric standards being evaluated.
- `validation` — Identifies whether cross-validation withheld one physical microplate at a time or withheld one interior dilution level across all microplates at once.
- `plate_n` — Number of distinct physical microplates represented among the withheld observations at this acid concentration.
- `attempted_n` — Number of distinct withheld plate-by-standard mean absorbances submitted to both models.
- `paired_prediction_n` — Number of withheld observations with finite absorbance predictions from both models within the acid-specific training concentration range. Both error summaries use these same observations.
- `shared_rmse_AU` — Root mean squared error of predicted 380 nm absorbance (AU) for a single calibration curve fitted across all acid concentrations.
- `shared_bias_AU` — Mean signed error in 380 nm absorbance (predicted minus measured, AU) for the shared calibration curve.
- `acid_specific_rmse_AU` — Root mean squared error of predicted 380 nm absorbance (AU) for the curve fitted only to training standards at the indicated acid concentration.
- `acid_specific_bias_AU` — Mean signed error in 380 nm absorbance (predicted minus measured, AU) for the acid-specific curve.
- `rmse_change_AU` — Acid-specific RMSE minus shared-model RMSE (AU). Negative values indicate lower error for the acid-specific model.
- `squared_error_reduction_percent` — Percentage reduction in mean squared error relative to the shared model, calculated as 100 × [1 − (acid-specific RMSE / shared RMSE)²]. Negative values indicate larger acid-specific error.

| acid&#95;<wbr>M | validation | plate&#95;<wbr>n | attempted&#95;<wbr>n | paired&#95;<wbr>prediction&#95;<wbr>n | shared&#95;<wbr>rmse&#95;<wbr>AU | shared&#95;<wbr>bias&#95;<wbr>AU | acid&#95;<wbr>specific&#95;<wbr>rmse&#95;<wbr>AU | acid&#95;<wbr>specific&#95;<wbr>bias&#95;<wbr>AU | rmse&#95;<wbr>change&#95;<wbr>AU | squared&#95;<wbr>error&#95;<wbr>reduction&#95;<wbr>percent |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.001 | held-out interior level | 3 | 30 | 30 | 0.0922504497551459 | -0.04625984964670857 | 0.061990633450916995 | 0.008948482591651926 | -0.030259816304228908 | 54.84403122525829 |
| 0.01 | held-out interior level | 3 | 30 | 30 | 0.08601135111849392 | -0.005724114738419225 | 0.08733314637971283 | 0.004733109814037893 | 0.0013217952612189082 | -3.0971533636369175 |
| 0.1 | held-out interior level | 3 | 30 | 30 | 0.09819146283718196 | 0.01942543575435276 | 0.09543874752959597 | 0.002445988661978314 | -0.0027527153075859873 | 5.52824084026895 |
| 1.0 | held-out interior level | 3 | 30 | 30 | 0.1181507485650668 | 0.055989168494514255 | 0.09646110797062868 | 0.008007475701193879 | -0.02168964059443812 | 33.34518353063013 |
| 0.001 | held-out plate | 3 | 36 | 35 | 0.09324158154208553 | -0.048684886167580224 | 0.06592366271969227 | 0.002404769009682485 | -0.02731791882239326 | 50.012272357624575 |
| 0.01 | held-out plate | 3 | 36 | 35 | 0.08690547762719701 | -0.010044755624939384 | 0.11375885313651873 | 0.0020049444018891225 | 0.026853375509321722 | -71.34684356959258 |
| 0.1 | held-out plate | 3 | 36 | 35 | 0.09742093299283465 | 0.01761314332090761 | 0.12892912503989692 | 0.0010907423974420873 | 0.031508192047062275 | -75.14490247536479 |
| 1.0 | held-out plate | 3 | 36 | 35 | 0.11878522151219342 | 0.046992969715479366 | 0.13015170674378645 | 0.0005031928137933311 | 0.011366485231593032 | -20.053523648205474 |
