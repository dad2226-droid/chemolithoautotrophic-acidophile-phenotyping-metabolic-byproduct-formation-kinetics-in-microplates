# Indicator calibration

Shared materials, strain and recipe definitions are linked from [00_materials_and_strains](../00_materials_and_strains/README.md). See [reading shared YAML references](../00_materials_and_strains/README.md#reading-shared-yaml-references) for the `$ref` convention.

This section contains abiotic absorbance measurements used to evaluate the thymol-blue pH indicator and the iron(III) absorbance indicator. These chemical standards calibrate the optical indicators used to follow the metabolic byproducts H⁺ and Fe³⁺ in cultures. The measurements include calibration concentration series, spectra, repeated storage measurements, plate-material and temperature comparisons, gas-exchange membrane measurements, and paired readings with and without a lid. Iron(II) preparations provide specificity controls.

See the [dataset overview](../README.md) for the other experimental collections. This folder provides the underlying measurements and their annotations.

## Files

| Subfolder | Experimental documentation | Source data |
|---|---|---|
| [01_ratiometric_pH_indicator](01_ratiometric_pH_indicator/) | [ratiometric_pH_indicator.yaml](01_ratiometric_pH_indicator/ratiometric_pH_indicator.yaml) | 96 reader exports covering initial calibrations, spectra, plate configurations and pH storage measurements. |
| [02_iron_III_indicator](02_iron_III_indicator/) | [iron_III_indicator.yaml](02_iron_III_indicator/iron_III_indicator.yaml) | 201 reader exports covering iron concentration series in different acid backgrounds, spectra, plate configurations, the final iron(III) storage series and initial iron(II) measurements. |
| [03_lid_comparison](03_lid_comparison/) | [lid_comparison.yaml](03_lid_comparison/lid_comparison.yaml) | [acid-iron.xlsx](03_lid_comparison/acid-iron.xlsx), containing 12 measured sheets: paired lid and no-lid readings of 3 pH plates and 3 iron(III) plates. |

The shared [reference_measurements.csv](reference_measurements.csv) contains the measured total-Fe and Fe(II) reference values for the final iron(III) and iron(II) preparations. Reference pH values and the configuration-series Fe concentrations are included directly in the corresponding YAML files. The YAML files contain the preparation details, reference methods, well assignments and acquisition associations needed to interpret these data.

There are 298 optical source files: 297 reader exports and 1 workbook. The YAML files describe 309 acquisitions because each of the workbook's 12 measured sheets is a separate acquisition.

## Reading the files

The `.xls` files are UTF-16LE tab-delimited reader exports containing instrument headers and wavelength-specific measurement blocks. `acid-iron.xlsx` is an Excel workbook with sheet names, instrument settings, channel times and well coordinates.

The pH/iron and paired-lid reader assignments are documented in the shared [materials and equipment table](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment) and referenced by the YAML `readouts`.

Endpoint acquisitions measure absorbance at 380, 428, 430, 485, 545 and 750 nm. Channels are identified by their recorded wavelength. Spectra span 200–750 nm in 10 nm steps; the pH scans ending at 700 nm are identified individually in the metadata.

Open a YAML file in a text editor to read its labelled sections. `experiments` describes preparation and replication; `reference_sets` provides calibration values and well assignments; `physical_plates` identifies the prepared plates; and `acquisitions` associates each export or workbook sheet with an experiment, physical plate, reference set and readout. Identifiers connect these sections. Match each exact `source_filename`, including capitalization, within its corresponding subfolder; workbook acquisitions also specify `source_sheet`.

The plate uses rows A–H and columns 1–12. A `values_by_column` list follows its accompanying `plate_columns` list in order. In the initial pH calibrations, rows A–D and E–H are different pH ranges and both contain dye. In the pH standard spectra and storage series, the final iron calibrations, and the iron configuration experiments, rows A–D are undyed and rows E–H contain dye. Each 4-row group supplies 4 technical wells per concentration. Use the explicit reference set for every acquisition.

## Experimental design

Each acquisition is assigned to its physical plate through `physical_plate_id`. Well coordinates and wavelengths identify measurements within each plate.

The configuration experiments used 3 independently prepared standard series, prepared on separate days. Corresponding UV-transparent and non-UV-transparent plates used the same standards within a preparation run. The 20 µM and 200 µM pH dye-stock experiments used separate physical plates. The initial pH calibrations used 3 independently prepared 24-standard series and 3 UV-transparent plates for each dye-stock concentration. The `acid-final` endpoint and spectral exports come from 3 non-UV-transparent plates containing 3 independently prepared 12-standard series, with 1 plate per preparation and 4 technical wells per condition per plate. Spectra and endpoint readings were acquired from these plates. The final iron calibration series used 1 preparation per series distributed among 3 technical plates. These replication structures are recorded separately in the YAML files.

For the configuration sequence, the same plates were read after 30 minutes at 23 °C, 2 hours at 4 °C, 2 hours at 30 °C, 2 hours at 75 °C, and overnight at 23 °C with the [documented sealing film](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment) applied before that incubation. The [documented gas-exchange membranes](../00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment) were then attached for an additional reading. These are prescribed incubation conditions; recorded reader temperatures and timestamps are separate fields.

The pH temperature and plate-type comparison used 3 UV-transparent and 3 non-UV-transparent plates, with 4 technical wells per condition per plate. The nominal 10 µM initial pH calibration provides the initial UV-transparent measurements. The same 6 plates were used for the membrane comparison after overnight incubation at 23 °C with lids and Parafilm. Paired measurements were acquired immediately before and after membrane attachment.

The pH lid comparison used 3 independently prepared 12-standard series distributed among 3 separate UV-transparent plates, with 4 technical wells per condition per plate. Each plate was measured without and with a non-UV-transparent low-evaporation lid on the same day. The pH storage experiment used 3 independently prepared 12-standard series distributed among 3 separate non-UV-transparent plates, with 4 technical wells per condition per plate. These plates were stored on the benchtop for approximately 15 days with low-evaporation lids, which were removed for measurements. The pH lid-comparison and storage experiments used the same 12 nominal sulfuric-acid levels.

The final iron(III) endpoint series and its corresponding spectra were measured on the same physical plates. The iron(III) lid comparison used 3 separate UV-transparent plates. Each plate was measured without and with a lid; these plates were distinct from those used for the final iron(III) calibration and storage series. The lid experiment used 1 standard preparation, with reference concentrations assigned in `reference_sets.lid_iron_standards`. During the final iron(III) storage series, plates were stored at 23 °C with low-evaporation lids obtained from sacrificial donor non-UV-transparent plates. Those lids were removed for measurements. Measurement lids were present in `acid-iron.xlsx` sheets ending in `-lid`; `-nolid` sheets and the other indicator-calibration acquisitions were measured without lids.

The final iron(II) controls include initial endpoint readings. The 1 M sulfuric-acid series also includes spectra.

Final-series filename suffixes are sequential reading indices. Date Last Saved records the acquisition time. Elapsed time is calculated from the recorded timestamps and the time origin specified for each experiment.

## Reference measurements and concentration basis

Undyed wells contained 200 µL of standard. Dyed wells contained that standard plus 10 µL of dye solution, giving 210 µL. The usual 200 µM dye stock produces a calculated final concentration of approximately 9.52 µM; the 20 µM stock produces approximately 0.952 µM. The YAML files distinguish these calculated values from the nominal 10 and 1 µM labels.

Iron reference concentrations and stated sulfuric-acid concentrations describe the 200 µL standard. Iron values are expressed on an elemental-iron basis. Iron calibration values represent measured total iron unless identified as iron(II) titration measurements. Configuration-series iron arrays are in M. Reference solutions were prepared by serial dilution, and the reference concentrations were measured after dilution. For the final reference series, the CSV contains 1 total-Fe value and 3 Fe(II) titration values per reference sample.

`reference_measurements.csv` is UTF-8 comma-separated text with a byte-order mark and 1 header row. Read it with a BOM-aware UTF-8 decoder so that its first column is recognized as `measurement_ID`. It has 144 concentration rows for 36 reference samples: 36 total-iron values and 108 iron(II) values. Its 3 `series_ID` values are `ferric`, `ferrous_1M` and `ferrous_0.001M`. Select measurements using `series_ID`, `plate_column` and `analyte`; `measurement_ID` identifies an individual row. `reference_replicate` identifies iron(II) reference replicates.

`measured_concentration_M` is the concentration before dye addition, accounting for analytical sample dilutions. Undyed wells contain 200 µL at this concentration. Dyed wells contain 210 µL, giving a final-well concentration of 200/210 times the reference concentration. The CSV’s `optical_final_volume_uL` value of 210 applies to dyed wells.

Final pH reference values were measured after dye addition.

## Data associations

The YAML files document source filenames, header labels and acquisition associations. Each acquisition identifies its source file.

## Contributors

CRediT roles for this experimental collection are listed below. See the [dataset information](../../README.md#dataset-information) and [dataset-wide contributor credits](../README.md#contributors).

| CRediT role<sup>1</sup> | Contributors |
|---|---|
| Conceptualization | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) |
| Methodology | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) |
| Investigation | [MAH](https://orcid.org/0009-0009-0165-2098), [SVRA](https://orcid.org/0009-0002-9372-3226), [NRL](https://orcid.org/0000-0002-3824-9122), [DAD](https://orcid.org/0000-0003-1908-7196) |
| Validation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [WW](https://orcid.org/0009-0002-5905-2318) |

[DAD](https://orcid.org/0000-0003-1908-7196) conceived the thymol blue pH indicator approach evaluated in these calibrations, and [MAH](https://orcid.org/0009-0009-0165-2098) conceived the Fe(III) absorbance indicator approach.

[DAD](https://orcid.org/0000-0003-1908-7196) designed the pH calibrations, iron configuration tests, ferrous specificity controls and lid comparison; [MAH](https://orcid.org/0009-0009-0165-2098) designed the final ferric-indicator calibration, spectra and storage series. [MAH](https://orcid.org/0009-0009-0165-2098) and [DAD](https://orcid.org/0000-0003-1908-7196) prepared the experiments and reference solutions, and [MAH](https://orcid.org/0009-0009-0165-2098) collected optical measurements across the calibration sets.

[MAH](https://orcid.org/0009-0009-0165-2098), [SVRA](https://orcid.org/0009-0002-9372-3226), [NRL](https://orcid.org/0000-0002-3824-9122) and [DAD](https://orcid.org/0000-0003-1908-7196) measured the final ferric series; [MAH](https://orcid.org/0009-0009-0165-2098), [SVRA](https://orcid.org/0009-0002-9372-3226) and [NRL](https://orcid.org/0000-0002-3824-9122) measured the ferrous specificity controls, and [MAH](https://orcid.org/0009-0009-0165-2098) and [NRL](https://orcid.org/0000-0002-3824-9122) performed the lid comparison. [MAH](https://orcid.org/0009-0009-0165-2098) and [DAD](https://orcid.org/0000-0003-1908-7196) measured reference pH and Fe(II), and [DAD](https://orcid.org/0000-0003-1908-7196) performed total-Fe measurements by ICP-OES.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [WW](https://orcid.org/0009-0002-5905-2318) contributed to validation of the indicator approaches using the calibration measurements.

## References

1. Elsevier. [CRediT Author Statement](https://www.elsevier.com/researcher/author/policies-and-guidelines/credit-author-statement) (accessed 2026-09-12).
