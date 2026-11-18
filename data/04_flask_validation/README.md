# Flask validation

Shared materials, strain and recipe definitions are linked from [00_materials_and_strains](../00_materials_and_strains/README.md). See [reading shared YAML references](../00_materials_and_strains/README.md#reading-shared-yaml-references) for the `$ref` convention.

Preparation reagent identities are linked through `preparation_reagent_id` to the shared [chemical reagent table](../00_materials_and_strains/chemicals_materials_equipment.md) and [chemical_reagents.yaml](../00_materials_and_strains/chemical_reagents.yaml).

This folder (`data/04_flask_validation/`) contains *Acidithiobacillus ferrooxidans* flask experiments and abiotic controls in AFM1 medium, including gallium, copper, sodium chloride, citric acid, vanadyl and aluminium treatments. It comprises 33 conditions, 99 physical flasks and 1,614 flask-by-sampling-event records.

## Files

| File | Contents |
|---|---|
| [flask_measurements.csv](flask_measurements.csv) | 1 row per flask and sampling event, containing reported measurements, sampling information and measurement-status fields. |
| [flask_validation.yaml](flask_validation.yaml) | Shared methods, AFM1 formulation, condition definitions, individual flask and culture assignments and the CSV data dictionary. |
| `README.md` | Instructions for interpreting and using the data. |

The CSV and YAML can be opened in ordinary spreadsheet and text applications. The YAML's `columns` section defines every CSV column; its `conditions` section identifies every experimental condition and physical flask. Follow the linked shared definitions when reading the YAML. Keep this experimental folder together with `00_materials_and_strains`.

## Experimental design

Each condition has 3 flask replicates. Sampling-event counts apply to each flask and include the initial event.

| Condition | Start date | Inoculation | Addition to AFM1 | Flask identifiers | Events per flask |
|---|---|---|---|---|---:|
| FV01 | 2025-07-11 | AF | None | DD25-4-1 to DD25-4-3 | 11 |
| FV02 | 2025-07-11 | AF | 0.001 mM Ga(III) | DD25-4-4 to DD25-4-6 | 11 |
| FV03 | 2025-07-11 | AF | 0.01 mM Ga(III) | DD25-4-7 to DD25-4-9 | 11 |
| FV04 | 2025-07-11 | AF | 0.1 mM Ga(III) | DD25-4-10 to DD25-4-12 | 11 |
| FV05 | 2025-07-11 | AF | 1 mM Ga(III) | DD25-4-13 to DD25-4-15 | 11 |
| FV06 | 2025-07-11 | AF | 10 mM Ga(III) | DD25-4-16 to DD25-4-18 | 11 |
| FV07 | 2025-07-11 | AF | 100 mM Ga(III) | DD25-4-19 to DD25-4-21 | 11 |
| FV08 | 2025-12-17 | Abiotic | None | DD25-4-22 to DD25-4-24 | 11 |
| FV09 | 2025-12-17 | Abiotic | 1 mM Ga(III) | DD25-4-25 to DD25-4-27 | 11 |
| FV10 | 2025-12-17 | Abiotic | 100 mM Ga(III) | DD25-4-28 to DD25-4-30 | 11 |
| FV11 | 2025-07-16 | AF | None | DD25-4b-1 to DD25-4b-3 | 15 |
| FV12 | 2025-07-16 | AF | 100 mM Ga(III) | DD25-4b-4 to DD25-4b-6 | 15 |
| FV13 | 2025-07-16 | AF | 100 mM Cu(II) | DD25-4b-7 to DD25-4b-9 | 15 |
| FV14 | 2025-07-16 | AF | 50 mM NaCl | DD25-4b-10 to DD25-4b-12 | 11 |
| FV15 | 2025-07-22 | AF | None | DD25-4c-1 to DD25-4c-3 | 17 |
| FV16 | 2025-07-22 | AF | 50 mM Citric acid | DD25-4c-4 to DD25-4c-6 | 24 |
| FV17 | 2025-07-22 | AF | 100 mM Ga(III) | DD25-4c-7 to DD25-4c-9 | 17 |
| FV18 | 2025-07-22 | AF | 100 mM Vanadyl(IV) | DD25-4c-13 to DD25-4c-15 | 17 |
| FV19 | 2025-07-22 | AF | 100 mM Al(III) | DD25-4c-16 to DD25-4c-18 | 17 |
| FV20 | 2025-08-20 | Abiotic | 1 mM Ga(III) | DD25-4e-1 to DD25-4e-3 | 20 |
| FV21 | 2025-08-20 | AF | None | DD25-4e-4 to DD25-4e-6 | 20 |
| FV22 | 2025-08-20 | AF | 1 mM Citric acid | DD25-4e-7 to DD25-4e-9 | 20 |
| FV23 | 2025-08-20 | AF | 6 mM Citric acid | DD25-4e-10 to DD25-4e-12 | 20 |
| FV24 | 2025-08-20 | AF | 1 mM Ga(III) | DD25-4e-13 to DD25-4e-15 | 20 |
| FV25 | 2025-08-20 | AF | 1 mM Ga(III) + 1 mM Citric acid | DD25-4e-16 to DD25-4e-18 | 20 |
| FV26 | 2025-08-20 | AF | 1 mM Ga(III) + 6 mM Citric acid | DD25-4e-19 to DD25-4e-21 | 20 |
| FV27 | 2025-08-29 | Abiotic | None | DD25-4g-1 to DD25-4g-3 | 20 |
| FV28 | 2025-08-29 | AF | None | DD25-4g-4 to DD25-4g-6 | 20 |
| FV29 | 2025-08-29 | AF | 12 mM Citric acid | DD25-4g-7 to DD25-4g-9 | 20 |
| FV30 | 2025-08-29 | AF | 18 mM Citric acid | DD25-4g-10 to DD25-4g-12 | 20 |
| FV31 | 2025-08-29 | AF | 1 mM Ga(III) | DD25-4g-13 to DD25-4g-15 | 20 |
| FV32 | 2025-08-29 | AF | 1 mM Ga(III) + 12 mM Citric acid | DD25-4g-16 to DD25-4g-18 | 20 |
| FV33 | 2025-08-29 | AF | 1 mM Ga(III) + 18 mM Citric acid | DD25-4g-19 to DD25-4g-21 | 20 |

All flasks initially contained 100 mL in the [documented culture flasks](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment) and were incubated in the New Brunswick I26 at 30 °C and 150 rpm. A 1 mL sample was removed at each sampling event. AF is identified in the [shared strain table](../00_materials_and_strains/strains.md). Biotic flask tests started at OD600 0.001. Approximately 20–50 µL of inoculum was added, with a consistent volume across tests started in the same batch. Abiotic controls received the same volume of sterile medium.

Each condition's `additives` list gives its additions and nominal concentrations in the prepared flask culture. An empty list means no addition to AFM1. Gallium concentrations describe Ga(III), while citric acid concentrations describe the amount of citric acid. Stocks were prepared by mass and final solution volume in the [documented stock-preparation vessel](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment). Measured soluble elemental concentrations are separate CSV fields.

The citric acid/Ga experiments started on 20 and 29 August 2025 and were followed through 28 August and 7 September, respectively. They include citric acid alone, Ga alone and combined conditions. The label “citrate” denotes citric acid; no citrate salt was used under that label.

## Media formulations

The shared [AFM1 formulation](../00_materials_and_strains/media_formulations.md#afm1-afm1) provides the complete recipe, preparation target and parent-formulation references. The `medium_definition` field in [flask_validation.yaml](flask_validation.yaml) references the flask-context definition in the shared folder.

Recipe amounts describe prepared medium; condition-specific additions are listed above. Preparation pH is a target, while initial chemistry values in the CSV are assigned from measured starting-medium references.

## Replication and time

The physical experimental unit is a flask. Sampling events repeatedly observe that same flask. Within each biotic condition, replicates 1, 2 and 3 were inoculated with different cultures. Across conditions started on the same day, corresponding replicate numbers used the same culture. The 84 inoculated flasks represent 15 source cultures across 5 start dates, together with 15 uninoculated abiotic flasks.

`inoculum_culture_id` identifies the starter culture for each biotic flask. For example, `AF_20250820_1` identifies the culture used for replicate 1 of the biotic tests started on 20 August 2025.

`sample_date` is a calendar date in M/D/YY format, with years in the 2000s; for example, `8/20/25` means 20 August 2025. Condition dates in the YAML use YYYY-MM-DD. `elapsed_days_mean` and `elapsed_days_sd` describe elapsed sampling times across the full experimental series. Sampling event 1 is the initial event at elapsed time 0.

## Reported measurements

The CSV contains measurements, initial-reference assignments and sampling information. At sampling event 1, pH, ORP, Fe(II), total Fe, Ga, Cu, V and Al values are assigned from measured starting-medium references. The YAML defines these assignments under `initial_chemistry_assignments` and in the corresponding column descriptions.

- `measured_pH` and `measured_ORP_mV` contain pH and oxidation-reduction potential. ORP is reported in mV versus an Ag/AgCl reference electrode containing 3 M KCl electrolyte.
- `cerimetric_FeII_equivalent_mM` contains the cerimetric result expressed as an Fe(II)-equivalent concentration. Vanadyl interferes with the Fe(II) titration; FV18 values have `cerimetry_status=vanadyl_interference`.
- `soluble_total_Fe_mM`, `soluble_Ga_mM`, `soluble_Cu_mM`, `soluble_V_mM` and `soluble_Al_mM` contain ICP-OES elemental concentrations in mM, accounting for analytical sample dilution.
- `measured_planktonic_fluorescence_RFU` and `measured_nonplanktonic_fluorescence_RFU` contain SYBR Green I fluorescence after TE-buffer blank subtraction, measured with the [assigned plate reader](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment) at 497 nm excitation and 520 nm emission using the cited method.<sup>1</sup>

The nonplanktonic fraction was recovered and measured using the cited fluorescence method.<sup>1</sup>

There are 1,545 numerical planktonic RFU observations and 588 numerical nonplanktonic RFU observations. Each is a reported result for a flask and sampling event.

## Reading the files

The CSV uses UTF-8, 1 header row, commas and a period as the decimal separator.

Use `flask_ID` and `sampling_event` as the record key. Join `condition_ID` to the YAML condition definition; `inoculum_culture_id` identifies the starter culture. Column and identifier names are case-sensitive.

## Contributors

CRediT roles for this experimental collection are listed below. See the [dataset information](../../README.md#dataset-information) and [dataset-wide contributor credits](../README.md#contributors).

| CRediT role<sup>2</sup> | Contributors |
|---|---|
| Methodology | [DAD](https://orcid.org/0000-0003-1908-7196) |
| Investigation | [DAD](https://orcid.org/0000-0003-1908-7196), [NRL](https://orcid.org/0000-0002-3824-9122) |
| Validation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [WW](https://orcid.org/0009-0002-5905-2318) |

[DAD](https://orcid.org/0000-0003-1908-7196) designed and prepared the flask experiments examining gallium, copper, sodium chloride, aluminum, vanadyl and citric acid treatments, including the citric acid-to-gallium series and abiotic controls. [DAD](https://orcid.org/0000-0003-1908-7196) prepared and measured the starting-medium references and collected the time-course pH, ORP, Fe(II) and elemental measurements. [DAD](https://orcid.org/0000-0003-1908-7196) and [NRL](https://orcid.org/0000-0002-3824-9122) measured fluorescence in the planktonic and nonplanktonic fractions.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [WW](https://orcid.org/0009-0002-5905-2318) contributed to assay validation using the flask reference measurements.

## References

1. [Dailey, D. A.](https://orcid.org/0000-0003-1908-7196); Hernández-Pagán, E.; Bailey, S.; Raffaele, N. E.; [Bavaresco, S. T.](https://orcid.org/0009-0007-4652-5733); Piatt-Price, A.; Carneiro, J. S. A.; Austin, R. N.; Doherty, C. J.; [Banta, S.](https://orcid.org/0000-0001-7885-0150) Biophytometallurgy: Biomining Metals from Plant Resources. *bioRxiv* (preprint), **2026**. DOI: [10.64898/2026.08.27.747594](https://doi.org/10.64898/2026.08.27.747594).

2. Elsevier. [CRediT Author Statement](https://www.elsevier.com/researcher/author/policies-and-guidelines/credit-author-statement) (accessed 2026-09-12).
