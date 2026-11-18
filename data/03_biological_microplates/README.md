# Biological microplate experiments

Shared materials, strain and recipe definitions are linked from [00_materials_and_strains](../00_materials_and_strains/README.md). See [reading shared YAML references](../00_materials_and_strains/README.md#reading-shared-yaml-references) for the `$ref` convention.

This collection contains plate-reader absorbance measurements and experimental metadata for AF, AT and LF, as identified in the [shared strain table](../00_materials_and_strains/strains.md). It covers iron- and sulfur-based culture experiments, including additive concentration series, copper treatments, and matched abiotic controls. The optical indicators track two metabolic byproducts in cultures. Absorbance at 380 nm indicates apparent Fe³⁺ generation, and the thymol blue pH ratio tracks acidification associated with H⁺ release. The organism-specific experiment records define the included conditions.

This section is stored under `data/03_biological_microplates/`. Reader exports are organized in the species and condition subfolders below.

## Files

Read [biological_definitions.yaml](biological_definitions.yaml) first. It references the shared organism identities and medium formulations and records preparation methods, units, well layout, incubation conditions and measurement definitions. Then read the experiment YAML in the relevant species directory:

| File or folder | Contents |
|---|---|
| [biological_definitions.yaml](biological_definitions.yaml) | Shared preparations, well layout, incubation conditions, units, media and readout definitions. |
| [AF experiments](01_Acidithiobacillus_ferrooxidans_ATCC_23270/AF_experiments.yaml) and [AF reader exports](01_Acidithiobacillus_ferrooxidans_ATCC_23270/) | *Acidithiobacillus ferrooxidans* ATCC 23270 experiments and their microplate measurements. |
| [AT experiments](02_Acidithiobacillus_thiooxidans_ATCC_19377/AT_experiments.yaml) and [AT reader exports](02_Acidithiobacillus_thiooxidans_ATCC_19377/) | *Acidithiobacillus thiooxidans* ATCC 19377 experiments and their microplate measurements. |
| [LF experiments](03_Leptospirillum_ferriphilum_DSM_14647/LF_experiments.yaml) and [LF reader exports](03_Leptospirillum_ferriphilum_DSM_14647/) | *Leptospirillum ferriphilum* DSM 14647 experiments and their microplate measurements. |
| [acquisitions.csv](chronology/acquisitions.csv) | Export paths, physical plates, acquisition times, chronological order and checksums. |
| [condition_timecourse_index.csv](chronology/condition_timecourse_index.csv) | Experiment, physical-plate and column assignments with condition values and units. |
| [Chronology instructions](chronology/README.md) and [load_timecourses.py](chronology/load_timecourses.py) | Instructions and a Python reader for joining conditions to acquisitions. |

The YAML files are plain-text documents. Named fields describe their values, indentation groups related information, and lists contain repeated records. The species experiment documents list `acquisition_ids`, which resolve through the acquisition registry. Each species document's `definitions_id` connects it to `biological_definitions.yaml`. Media, preparation and readout identifiers select the corresponding shared definitions.

Each treatment compound's `reagent_id` selects its entry in the shared [chemical reagent table](../00_materials_and_strains/chemicals_materials_equipment.md) and [chemical_reagents.yaml](../00_materials_and_strains/chemical_reagents.yaml). The table gives hydrate forms and supplier/bottle information for the treatment compounds. Nominal-versus-measured concentration assignments are documented in the shared definitions; N in the inventory marks nominal additive stocks. Citric acid is the citrate component in the combination experiments.

## Reading the files

The `.xls` files contain UTF-16LE tab-delimited text with a byte-order mark. Import them as text using UTF-16 encoding and a tab delimiter; their `.xls` extension does not indicate a binary Excel workbook. In a spreadsheet program, use its text-import option and select these settings if the format is not recognized automatically.

Instrument headers and wavelength labels identify each measurement block. Exports contain 380 and 750 nm channels, with additional wavelength channels in some files.

Biological export filenames use the embedded filename followed by `__saved_YYYYMMDDTHHMMSS.xls`. The timestamp is the recorded local `Date Last Saved` value. Exports were saved immediately after acquisition, so this field supplies the acquisition time.

The [acquisition registry](chronology/acquisitions.csv) gives each export's path, acquisition ID, physical plate identity, embedded filename, timestamp, chronological index and SHA-256 checksum. Paths are relative to the directory containing `data/`. Species YAML files list each physical plate's `acquisition_ids` in chronological order. Resolve those IDs through the registry to identify each physical plate and export.

The [condition timecourse index](chronology/condition_timecourse_index.csv) maps each experiment, plate and column to its conditions and units. Join it to the acquisition registry using `physical_plate_id` to obtain the readings in each time course. The [chronology instructions and loader](chronology/README.md) provide a working example.

`native_read_label` records the embedded label; `read_index_zero_based` is the chronological position within a physical plate, starting at 0. Acquisition timestamps determine elapsed time.

## Experimental design

The tables list all 88 experiments. A dosing series is a plate group with 1 set of column conditions. Microplate counts are unique physical plates within each experiment; a physical plate can contribute to several experiments through different columns. The linked species metadata give every dose, assigned column and acquisition.

### Acidithiobacillus ferrooxidans ATCC 23270

| Comparison | Medium | Initial OD600 | Readout; final well volume | Dosing series | Microplates |
|---|---|---:|---|---:|---:|
| Sodium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Cobalt(II) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 4 | 12 |
| Indium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Iron(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Potassium bromide concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 4 | 12 |
| Titanyl sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Zinc sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Vanadyl sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Ammonium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Gallium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 4 | 12 |
| Potassium iodide concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Potassium chloride concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Citric acid concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Potassium chloride concentration series (nominal pH 1.4) | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Sodium gluconate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Dysprosium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Yttrium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Europium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Neodymium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Magnesium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Terbium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Gadolinium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Cerium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Samarium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Lanthanum(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Zirconium(IV) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Erbium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Potassium chloride concentration series (nominal pH 1.2) | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Potassium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Sodium chloride concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Manganese(II) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Monopotassium phosphate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Iron(III) chloride concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Potassium chloride concentration series (nominal pH 1.0) | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Praseodymium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Ytterbium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Sucrose concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Cerium(IV) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Potassium nitrate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Cadmium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Nickel(II) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Chromium(III) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Potassium fluoride concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Potassium thiocyanate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Potassium dichromate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Potassium perrhenate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Potassium tetrathionate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Sodium selenite concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Sodium sulfite concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Sodium formate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Iron(II) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Sodium molybdate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Sodium tungstate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Sodium chromate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Lithium sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Citric acid and gallium combinations | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 2 | 6 |
| Silver sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Aluminum sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Potassium tetrathionate concentration series | AFM1_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Aluminum sulfate concentration series | AFM4 | 0.1 | pH indicator; 210 µL | 1 | 3 |
| Zinc sulfate concentration series | AFM4 | 0.1 | pH indicator; 210 µL | 1 | 3 |
| Cobalt(II) sulfate concentration series | AFM4 | 0.1 | pH indicator; 210 µL | 1 | 3 |
| Potassium sulfate concentration series | AFM4 | 0.1 | pH indicator; 210 µL | 1 | 3 |
| Potassium sulfate concentration series | AFM4 | 0.01 | pH indicator; 210 µL | 1 | 3 |
| Potassium sulfate concentration series | AFM4 | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Gallium sulfate concentration series | AFM4 | 0.1 | pH indicator; 210 µL | 1 | 3 |
| Sulfuric acid series with 20 mM sodium thiosulfate | SM4_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Sulfuric acid series with 20 mM potassium tetrathionate | SM4_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Sodium thiosulfate concentration series | SM4_BS | 0.001 | pH indicator; 230 µL | 1 | 3 |
| Sodium thiosulfate and potassium tetrathionate mixtures | SM4_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Potassium tetrathionate concentration series | SM4_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Sodium thiosulfate concentration series | SM4_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Potassium iodate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Copper(II) sulfate concentration series | AFM1 | 0.001 | Fe(III) absorbance; 200 µL | 3 | 9 |
| Copper(II) sulfate concentration series | AFM4 | 0.1 | pH indicator; 210 µL | 1 | 3 |

### Acidithiobacillus thiooxidans ATCC 19377

| Comparison | Medium | Initial OD600 | Readout; final well volume | Dosing series | Microplates |
|---|---|---:|---|---:|---:|
| Sulfuric acid series with 20 mM potassium tetrathionate | SM4_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |
| Sulfuric acid series with 20 mM sodium thiosulfate | SM4_BS | 0.001 | pH indicator; 210 µL | 1 | 3 |

### Leptospirillum ferriphilum DSM 14647

| Comparison | Medium | Initial OD600 | Readout; final well volume | Dosing series | Microplates |
|---|---|---:|---|---:|---:|
| Potassium sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Zinc sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Iron(II) sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Iron(III) sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Gallium sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Indium sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Aluminum sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Cobalt(II) sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Vanadyl sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Cerium(III) sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |
| Cerium(IV) sulfate concentration series | HH | 0.001 | Fe(III) absorbance; 200 µL | 1 | 3 |

An `experiment_id` identifies a comparison. Its `plate_groups` list the physical `plate_ids` and the `column_conditions` for those plates. Each column condition specifies the added compound or compounds and their concentrations; `fixed_components`, when present, apply throughout the experiment in addition to the medium formulation.

The `plates` section identifies each physical plate, its `batch_id`, and its ordered `acquisition_ids`. The registry's `plate_id` matches the physical plate key in the species YAML; `physical_plate_id` is a unique registry identifier for that same plate. Repeated measurements retain the same well identities.

The earliest recorded reading of each physical plate is time 0. Elapsed time is calculated from acquisition timestamps. The timestamp field and its interpretation are specified in `biological_definitions.yaml`. Recorded wavelength labels identify the measurements, and reader temperatures remain distinct from incubation temperatures.

## Well assignments and replication

For each assigned column, rows E, F and G contain biotic wells and rows B, C and D contain condition-matched abiotic wells. Within each microplate, E, F and G were inoculated from 3 separately grown cultures. Other replicate microplates used independent sets of 3 cultures. Independent microplates were the biological replicates. Culture numbers 1, 2 and 3 apply within each microplate. Abiotic wells match the plate and column condition.

A physical plate can support several experiments through different column assignments. A shared control well remains 1 physical well wherever it is referenced. `blank_condition_columns` identify columns without an assigned experimental condition; 0-dose controls are explicitly represented by numerical 0 in `column_conditions`. Peripheral water-filled wells are defined separately in the shared layout.

## Concentrations and preparation volumes

Additive concentrations and initial OD600 values use the 200 µL preparation basis. Each plate group states its concentration unit. The shared definitions distinguish nominal concentrations from concentrations determined using elemental measurements and specify the chemical basis, including any explicit component-level exception. Added iron is additional to the iron already present in the medium.

Iron experiments contain no indicator dye and retain a 200 µL well volume. Standard pH preparations add 10 µL of dye stock, giving 210 µL; 3-fold-dye preparations add 30 µL, giving 230 µL. The shared preparation definitions give the volume factors for converting concentrations and OD600 from the 200 µL preparation basis to the final-well volume. Medium formulations, additive-stock pH, and experiment-specific pH values are separate quantities defined in the YAML files.

AFM1 additive stocks use pH 1.80 unless the experiment specifies an exception. This value describes the additive stock before mixing.

Potassium fluoride, potassium thiocyanate, potassium dichromate, sodium chromate, sodium selenite, sodium thiosulfate, potassium tetrathionate, ammonium sulfate, citric acid, sucrose and sulfuric acid have nominal concentrations. These assignments apply to additive screens and sulfur-substrate experiments; see `concentration_conventions.nominal_reagents` in the shared definitions. Titanyl sulfate and zirconium(IV) sulfate stocks were measured by ICP-OES in 2% (w/v) HCl, and their concentrations use the compound basis.

The shared [ICP-OES method and wavelengths](../00_materials_and_strains/experimental_definitions.yaml) specify the analytical matrix and processing order: subtract the analytical blank, apply the analytical dilution factor, then convert to elemental molarity and, where applicable, compound molarity. Concentrations include the analytical dilution correction. The additive stocks passed the stock-only and stock–AFM1 solubility and stability checks defined in `stock_preparation.solubility_screen`.

All nominal stock concentrations, including citric acid and sucrose, were calculated from the mass added and final prepared volume in the [documented stock-preparation vessel](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment), followed by the stated mixing and dilution volumes.

The initial pH for the identified SM4 experiments is 3.95 ± 0.03 (mean ± SD, 3 measurements). Initial pH was measured within the ranges stated in the methods. The normalized inoculum suspension OD600 of 5.0 ± 0.1 also denotes mean ± SD from 3 measurements.

## Media formulations

The shared [media formulation tables](../00_materials_and_strains/media_formulations.md) provide the complete recipes and parent-formulation references. This collection uses `AFM1`, `AFM1_BS`, `AFM4`, `SM4_BS_microplate` and `HH`. The local `SM4_BS` identifier selects `SM4_BS_microplate`; the test-tube formulation is a separate entry.

Recipe amounts describe nominal 1× medium on the 200 µL preparation basis before dye addition. Medium is dispensed at 2× strength for equal-volume mixing. Substrate and treatment additions remain in the species experiment YAMLs. The `media` field in `biological_definitions.yaml` references the corresponding shared definitions while retaining these local IDs.

The experiment metadata give the initial-mixture pH values.

## Contributors

CRediT roles for this experimental collection are listed below. See the [dataset information](../../README.md#dataset-information) and [dataset-wide contributor credits](../README.md#contributors).

| CRediT role<sup>1</sup> | Contributors |
|---|---|
| Conceptualization | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) |
| Methodology | [DAD](https://orcid.org/0000-0003-1908-7196), [LNM](https://orcid.org/0009-0000-5291-5200), [STB](https://orcid.org/0009-0007-4652-5733) |
| Investigation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [NRL](https://orcid.org/0000-0002-3824-9122), [LNM](https://orcid.org/0009-0000-5291-5200) |
| Validation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [WW](https://orcid.org/0009-0002-5905-2318) |

[DAD](https://orcid.org/0000-0003-1908-7196) conceived the thymol blue pH indicator approach, and [MAH](https://orcid.org/0009-0009-0165-2098) conceived the Fe(III) absorbance indicator approach applied in the biological microplate experiments.

[DAD](https://orcid.org/0000-0003-1908-7196) designed the AFM1, AFM4 and LF experiments; [DAD](https://orcid.org/0000-0003-1908-7196), [LNM](https://orcid.org/0009-0000-5291-5200) and [STB](https://orcid.org/0009-0007-4652-5733) designed the AF/AT soluble-sulfur series in SM4 basal salts.

[DAD](https://orcid.org/0000-0003-1908-7196) and [NRL](https://orcid.org/0000-0002-3824-9122) prepared and measured the AFM1 additive screen, with [MAH](https://orcid.org/0009-0009-0165-2098) contributing to a subset of those series. [LNM](https://orcid.org/0009-0000-5291-5200) and [DAD](https://orcid.org/0000-0003-1908-7196) prepared and measured the AFM4 series, the basal-AFM1 potassium tetrathionate series and the AF/AT soluble-sulfur experiments in SM4 basal salts. [DAD](https://orcid.org/0000-0003-1908-7196) and [MAH](https://orcid.org/0009-0009-0165-2098) prepared the LF experiments, and [MAH](https://orcid.org/0009-0009-0165-2098) and [LNM](https://orcid.org/0009-0000-5291-5200) collected their microplate measurements. [DAD](https://orcid.org/0000-0003-1908-7196) and [NRL](https://orcid.org/0000-0002-3824-9122) prepared the biological additive stocks and inocula; [DAD](https://orcid.org/0000-0003-1908-7196) performed the stock elemental measurements and initial inoculum pH and OD600 measurements.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [WW](https://orcid.org/0009-0002-5905-2318) contributed to validation of the indicator approaches using the biological microplate data.

## References

1. Elsevier. [CRediT Author Statement](https://www.elsevier.com/researcher/author/policies-and-guidelines/credit-author-statement) (accessed 2026-09-12).
