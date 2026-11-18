# Test-tube validation

Shared materials, strain and recipe definitions are linked from [00_materials_and_strains](../00_materials_and_strains/README.md). See [reading shared YAML references](../00_materials_and_strains/README.md#reading-shared-yaml-references) for the `$ref` convention.

Preparation reagent identities are linked through `preparation_reagent_id` to the shared [chemical reagent table](../00_materials_and_strains/chemicals_materials_equipment.md) and [chemical_reagents.yaml](../00_materials_and_strains/chemical_reagents.yaml).

This folder contains 96 tubes in 32 conditions with 3 replicates each, spanning 5 acidophile species and abiotic controls. Experiments began on 4 August 2026. HH iron cultures were followed through 11 August; elemental-sulfur and tetrathionate cultures were followed through 25 August.

This section is stored under `data/05_test_tube_validation/`.

## Files

| File | Contents |
|---|---|
| [test_tube_measurements.csv](test_tube_measurements.csv) | Individual tube measurements and sampling times, shared initial-reference assignments, and endpoint planktonic fluorescence. |
| [test_tube_validation.yaml](test_tube_validation.yaml) | Organism and culture identities, condition assignments, medium recipes, methods, data dictionary and the qualitative LF starter follow-up observation. |
| `README.md` | Guide to using these files. |

Keep this experimental folder together with `00_materials_and_strains`. The CSV and YAML, including the referenced shared definitions, describe this section. The CSV has 1 header row and 22 columns; the YAML's `columns` section defines those exact column names in the same order.

## Experimental design

| Timing block | timing_block_ID | Tube IDs | Conditions | Tubes | Sampling events after the start |
|---|---|---|---:|---:|---:|
| HH iron medium | `HH` | GT-S2-1–36 | 12 | 36 | 9 |
| SM4 with brown elemental sulfur | `brown_sulfur` | GT-S2-37–66 | 10 | 30 | 9 |
| SM4 with 5 mM potassium tetrathionate (TB-04) | `tetrathionate` | GT-S2-67–96 | 10 | 30 | 16 |

The strain codes AF, AT, LF, AFur and AC are defined in the [shared strain table](../00_materials_and_strains/strains.md). `Ab` denotes an abiotic control. The YAML assigns every physical tube to its organism, medium, treatment and incubation temperature.

Cultures contained 6 mL prepared medium in the [documented capped culture tubes](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment), with 20 µL washed, normalized inoculum added at an initial target OD600 of 0.003. Abiotic tubes received matching cell-free basal salts. Tubes were shaken at 150 rpm at their specified 30 or 37 °C. Evaporation was compensated with sterile ultrapure Milli-Q water before measuring pH or ORP directly in the tubes.

3 independent starter cultures were used per organism across treatments. `inoculum_culture_ID` identifies the starter culture for each tube, and `tube_ID` connects each tube to its sampling events. The test-tube starters were freshly thawed.

## Media formulations

The shared [media formulation tables](../00_materials_and_strains/media_formulations.md) provide complete recipes, preparation targets and parent-formulation references for `HH`, `SM4_BS_tube`, `SM4_S0` and `SM4_tetrathionate`. The local `SM4_BS` identifier selects `SM4_BS_tube`; the microplate basal formulation is a separate entry. The YAML `media` field references these tube-context definitions.

Recipe concentrations apply to the 6 mL prepared medium before inoculum addition. For TB-04 at 30 and 37 °C, this mixture contained 5 mM potassium tetrathionate. The substrate media inherit the tube basal-salt recipe. The 2× medium and treatment solutions were mixed in equal volumes; the treatment solution was prepared at pH 1.8. A 200 mM Ga stock gives 100 mM in the prepared mixture; a 20 mM In stock gives 10 mM. These treatment additions are separate from the recipes. Recorded initial mixed-medium pH values remain distinct from stock-preparation targets.

## Reading the files

The record key is `tube_ID` plus `sampling_event`. Use `condition_ID` to find the condition in the YAML. `medium_ID`, `treatment_ID` and `timing_block_ID` identify the corresponding definitions. Column names are case-sensitive.

There are 1,074 time-course rows and 96 shared initial-reference rows, giving 1,170 rows in total. The measured time course includes 749 numerical pH readings and 323 numerical ORP readings. These counts exclude the assigned initial-reference values.

`record_role` distinguishes 3 types of row:

- `shared_initial_reference`: event 0, containing the initial values assigned to each tube series.
- `measurement`: an intermediate tube-reading event.
- `endpoint`: the final time-course event for the tube’s block.

The YAML defines 21 initial references and their assignments to the 96 event-0 rows through `initial_reference_ID`.

Use `time_average_d` as the sampling time for each `timing_block_ID` and `sampling_event`. `time_sd_d` is the sample SD across tube times in that block/event, using n − 1. The 3 blocks have separate sampling schedules.

`tube_datetime` contains ISO 8601 local dates and clock times; `datetime_basis` identifies recorded, interpolated or assigned-start times. For HH event 6, times for GT-S2-2–35 are linearly interpolated between GT-S2-1 at 06:23:10 and GT-S2-36 at 07:34:15 on 8 August 2026. The block/event summary uses 2 recorded times and 34 interpolated times. Other poststart times are recorded. Event 0 is assigned time 0, with `time_average_d` and `time_sd_d` both 0.

## Endpoint fluorescence

`planktonic_background_corrected_rfu` contains 96 reported endpoint SYBR Green measurements: 36 for HH iron cultures, 30 for brown elemental-sulfur cultures and 30 for tetrathionate cultures. These rows have `fluorescence_status=reported`.

SYBR Green I dsDNA fluorescence was measured using the cited method.<sup>1</sup> Fluorescence values include TE-buffer blank subtraction.

## LF starter follow-up

The YAML contains the recorded qualitative follow-up of the original LF microplate starters. 3 replicate HH tubes were incubated at 37 °C and 150 rpm for 14 days, beginning on the day the LF microplate experiments ended. They did not become more red than abiotic controls. This record is separate from the freshly thawed starter cultures used for the numbered test-tube series.

## Contributors

CRediT roles for this experimental collection are listed below. See the [dataset information](../../README.md#dataset-information) and [dataset-wide contributor credits](../README.md#contributors).

| CRediT role<sup>2</sup> | Contributors |
|---|---|
| Methodology | [DAD](https://orcid.org/0000-0003-1908-7196) |
| Investigation | [DAD](https://orcid.org/0000-0003-1908-7196), [STB](https://orcid.org/0009-0007-4652-5733), [SVRA](https://orcid.org/0009-0002-9372-3226), [LNM](https://orcid.org/0009-0000-5291-5200) |
| Validation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [WW](https://orcid.org/0009-0002-5905-2318) |

[DAD](https://orcid.org/0000-0003-1908-7196) designed the test-tube validation experiments and designed, prepared and recorded the 14-day qualitative follow-up of the original LF microplate starters. [DAD](https://orcid.org/0000-0003-1908-7196), [STB](https://orcid.org/0009-0007-4652-5733) and [SVRA](https://orcid.org/0009-0002-9372-3226) prepared the iron, elemental-sulfur and potassium tetrathionate experiments and their abiotic controls.

[STB](https://orcid.org/0009-0007-4652-5733), [SVRA](https://orcid.org/0009-0002-9372-3226) and [DAD](https://orcid.org/0000-0003-1908-7196) measured the iron-medium experiments; [STB](https://orcid.org/0009-0007-4652-5733), [DAD](https://orcid.org/0000-0003-1908-7196) and [LNM](https://orcid.org/0009-0000-5291-5200) measured the elemental-sulfur experiments; and [LNM](https://orcid.org/0009-0000-5291-5200), [STB](https://orcid.org/0009-0007-4652-5733) and [DAD](https://orcid.org/0000-0003-1908-7196) measured the potassium tetrathionate experiments. [STB](https://orcid.org/0009-0007-4652-5733), [SVRA](https://orcid.org/0009-0002-9372-3226), [DAD](https://orcid.org/0000-0003-1908-7196) and [LNM](https://orcid.org/0009-0000-5291-5200) contributed to preparation of the pH/ORP and endpoint fluorescence assays and to pH/ORP measurements, while [STB](https://orcid.org/0009-0007-4652-5733) and [LNM](https://orcid.org/0009-0000-5291-5200) measured endpoint SYBR Green fluorescence.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [WW](https://orcid.org/0009-0002-5905-2318) contributed to assay validation using the test-tube reference measurements and LF starter follow-up.

## References

1. [Dailey, D. A.](https://orcid.org/0000-0003-1908-7196); Hernández-Pagán, E.; Bailey, S.; Raffaele, N. E.; [Bavaresco, S. T.](https://orcid.org/0009-0007-4652-5733); Piatt-Price, A.; Carneiro, J. S. A.; Austin, R. N.; Doherty, C. J.; [Banta, S.](https://orcid.org/0000-0001-7885-0150) Biophytometallurgy: Biomining Metals from Plant Resources. *bioRxiv* (preprint), **2026**. DOI: [10.64898/2026.08.27.747594](https://doi.org/10.64898/2026.08.27.747594).

2. Elsevier. [CRediT Author Statement](https://www.elsevier.com/researcher/author/policies-and-guidelines/credit-author-statement) (accessed 2026-09-12).
