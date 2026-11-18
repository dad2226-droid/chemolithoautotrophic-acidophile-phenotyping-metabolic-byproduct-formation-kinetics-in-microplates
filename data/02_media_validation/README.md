# Media validation

Shared materials, strain and recipe definitions are linked from [00_materials_and_strains](../00_materials_and_strains/README.md). See [reading shared YAML references](../00_materials_and_strains/README.md#reading-shared-yaml-references) for the `$ref` convention.

Abiotic medium tests characterize how the optical indicators of the metabolic byproducts Fe³⁺ and H⁺ respond to added ferric iron and to pH in each medium.

This collection contains endpoint absorbance measurements testing ratiometric pH and iron(III) indicators in acidophile media. It comprises 12 abiotic plates, 936 measured wells and 5,616 individual absorbance values. Each well was measured at 6 wavelengths. The F2S plates support both indicator tests.

See the [dataset overview](../README.md) for the other experimental collections. This folder provides the underlying measurements and their annotations.

## Files and experimental design

| Source exports | Media and source subfolder | Measurements | Measured region per plate |
|---|---|---|---|
| `Plate1-4-1-2026.xls`, `Plate2-4-1-2026.xls`, `Plate3-4-1-2026.xls` | [F2S](03_ratiometric_pH_and_iron_III_indicators/01_F2S_medium/) | Iron(III) additions and pH series with/without dye | Rows A–H, columns 1–9: 72 wells |
| `Plate4-4-1-2026.xls`, `Plate5-4-1-2026.xls`, `Plate6-4-1-2026.xls` | [AFM1 and AFM3](02_iron_III_indicator/01_AFM1_and_AFM3_media/) | Iron(III) additions without dye | Rows A–H, columns 1–6: 48 wells |
| `Plate7-4-1-2026.xls`, `Plate8-4-1-2026.xls`, `Plate9-4-1-2026.xls` | [AFM4 and AFM5](01_ratiometric_pH_indicator/01_AFM4_and_AFM5_media/) | pH series with/without dye | Rows A–H, columns 1–12: 96 wells |
| `Plate10-4-1-2026.xls`, `Plate11-4-1-2026.xls`, `Plate12-4-1-2026.xls` | [Starkey with brown or colloidal sulfur](01_ratiometric_pH_indicator/02_Starkey_with_brown_and_colloidal_sulfur_media/) | pH series with/without dye | Rows A–H, columns 1–12: 96 wells |

The supporting files are:

- [media_validation.yaml](media_validation.yaml): experimental design, preparation volumes, known medium recipes, complete well assignments, acquisition information, exact source filenames, checksums, and reference-measurement associations.
- [reference_measurements.csv](reference_measurements.csv): measured pH references and the measured iron concentration of the iron(III) stock. Its column definitions and measurement methods are included in the YAML.
- `README.md`: this guide.

## Reading the files

The `.xls` files contain UTF-16LE tab-delimited text, with a byte-order mark. They are not binary Excel workbooks. Import them as tab-delimited text using UTF-16 encoding when software does not recognize the contents automatically.

Every file contains 1 raw endpoint absorbance block. The wavelength order in the export is 380, 430, 485, 545, 428, 750 nm. Absorbance is reported in AU. Preserve the recorded numerical values and instrument tokens.

Physical plates are identified by the exact source filename and the YAML `physical_plate_id`. The `Original Filename` field records the export filename. Acquisition times are recorded under `Date Last Saved` and are dated 1 April 2026. Recorded reader temperatures range from 23.8 to 24.2 °C.

The [assigned plate reader](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment) is identified through the YAML reader settings. Measurements followed 30 minutes at approximately 23 °C.

The measured plate regions are listed above.

## Well assignments and replication

Rows A–H correspond to 8 dilution levels, in descending stock fraction: 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 0. Each 3-column block contains technical replicate wells at each row's level.

| Plates | Columns 1–3 | Columns 4–6 | Columns 7–9 | Columns 10–12 |
|---|---|---|---|---|
| 1–3 | F2S iron(III) series; no dye | F2S pH; dye | F2S pH; no dye | Unmeasured |
| 4–6 | AFM1 iron(III) series; no dye | AFM3 iron(III) series; no dye | Unmeasured | Unmeasured |
| 7–9 | AFM4 pH; no dye | AFM4 pH; dye | AFM5 pH; no dye | AFM5 pH; dye |
| 10–12 | Brown-sulfur Starkey pH; no dye | Brown-sulfur Starkey pH; dye | Colloidal-sulfur Starkey pH; no dye | Colloidal-sulfur Starkey pH; dye |

The 4 experimental families used the same calibration stocks, with 3 physical plate replicates per family and 3 technical wells per level on each plate. Physical plates are the replicate unit. The YAML links the 4 families to their stock set. The wells were abiotic.

The YAML's `groups` section encodes these assignments explicitly. Its `experiments` section assigns groups to physical plates, and `acquisitions` connects each physical plate to its export. The `matched_group_id` associates dye and no-dye groups within the same medium; comparisons should retain physical plate and dilution-level pairing.

## Preparation and concentration basis

Each well received 100 µL of 2× medium plus 100 µL of an acid or ferric standard. Dye-containing pH wells received an additional 10 µL of 200 µM thymol blue. Thus, undyed wells contain 200 µL, and dyed wells contain 210 µL. All ferric-series wells are undyed. Plates were non-UV-transparent and were measured without lids.

The nominal dye label is 10 µM. The shared medium formulations describe the standard 1× medium before the additional dye solution.

For the acid series, nominal 1 M sulfuric acid was serially diluted with water. Row H contains water in place of an acid standard. Acid contributed by this series is separate from acid already present in the medium.

For the iron(III) series, iron(III) sulfate stock contained 0.2094 M total iron measured by ICP-OES and was diluted in 1 M sulfuric acid. Row H contains the acid solvent without iron(III) addition. The medium recipes specify their iron content.

## Media formulations

The complete recipes, preparation targets and parent-formulation references are in the shared [media formulation tables](../00_materials_and_strains/media_formulations.md). This collection uses `AFM1`, `AFM3`, `F2S`, `AFM4`, `AFM5`, `Starkey_brown` and `Starkey_colloidal`; its YAML `media` field references their context-specific stock and standard-medium definitions.

All calibration media were dispensed as 2× stocks. The listed 1× recipes describe medium before the extra dye solution. Stock-preparation quantities remain available through the shared definitions. The 2 Starkey sulfur conditions use distinct preparation identities; the colloidal condition is also labelled “Yellow Starkey.” Sulfur-containing stocks were aggressively vortexed before dispensing. Experimental acid/ferric additions and their reference measurements remain defined in this folder.

## Using the reference table

`reference_measurements.csv` is UTF-8 comma-separated text with 1 header row and 49 reference records:

- 8 final-mixture pH values associated with F2S.
- 8 final-mixture pH values for AFM4 and 8 separately measured values for AFM5.
- 8 final-mixture pH values for brown-sulfur Starkey and 8 separately measured values for colloidal-sulfur Starkey.
- 8 measured pH values for the acid/DI-water solutions before medium mixing.
- 1 measured total-iron concentration for the ferric stock.

The CSV contains measured reference values; medium recipes and stock-preparation details are recorded in the YAML.

The CSV columns are `reference_ID`, `reference_set_ID`, `level_ID`, `plate_row`, `quantity`, `value`, `unit`, `measurement_basis`, and `method_ID`. `reference_ID` identifies each reference value; `reference_set_ID` identifies its series, and `method_ID` selects its method definition in the YAML. Reference sets list 8 identifiers in plate-row order A–H. pH was measured separately for each medium before dye addition. Each medium’s pH references are associated with its dye conditions and physical plates.

The pH-meter and ICP-OES method definitions are referenced through the YAML from the shared folder. For the Fe-stock reference, a 1 mL sample was centrifuged at 3,000 × g for 10 seconds; 750 µL of supernatant was transferred to a new tube and centrifuged at 17,000 × g for 1 minute. The clarified supernatant was diluted for ICP-OES analysis. The reported elemental concentration describes the starting ferric stock.

## Preservation and reuse

Preserve source exports exactly as recorded. Do not automatically remove edge wells; they contain assigned conditions and technical replicates.

File sizes and SHA-256 checksums are included with the acquisition entries. Match each `source_filename` to its acquisition entry and use the physical-plate, well and reference identifiers to connect the measurements.

## Contributors

CRediT roles for this experimental collection are listed below. See the [dataset information](../../README.md#dataset-information) and [dataset-wide contributor credits](../README.md#contributors).

| CRediT role<sup>1</sup> | Contributors |
|---|---|
| Conceptualization | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) |
| Methodology | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) |
| Investigation | [MAH](https://orcid.org/0009-0009-0165-2098), [DAD](https://orcid.org/0000-0003-1908-7196) |
| Validation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [WW](https://orcid.org/0009-0002-5905-2318) |

[DAD](https://orcid.org/0000-0003-1908-7196) conceived the thymol blue pH indicator approach, and [MAH](https://orcid.org/0009-0009-0165-2098) conceived the Fe(III) absorbance indicator approach evaluated in these media.

[DAD](https://orcid.org/0000-0003-1908-7196) and [MAH](https://orcid.org/0009-0009-0165-2098) designed the indicator-validation experiments in AFM1, AFM3, F2S, AFM4, AFM5 and the 2 Starkey sulfur media. [MAH](https://orcid.org/0009-0009-0165-2098) and [DAD](https://orcid.org/0000-0003-1908-7196) prepared the media-validation experiments and reference solutions, and [MAH](https://orcid.org/0009-0009-0165-2098) collected the plate-reader measurements. [MAH](https://orcid.org/0009-0009-0165-2098) and [DAD](https://orcid.org/0000-0003-1908-7196) measured reference pH and Fe(II), and [DAD](https://orcid.org/0000-0003-1908-7196) performed total-Fe reference measurements by ICP-OES.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [WW](https://orcid.org/0009-0002-5905-2318) contributed to validation of the indicator approaches using the medium-validation measurements.

## References

1. Elsevier. [CRediT Author Statement](https://www.elsevier.com/researcher/author/policies-and-guidelines/credit-author-statement) (accessed 2026-09-12).
