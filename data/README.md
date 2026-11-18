# Chemolithoautotrophic acidophile phenotyping via a miniaturized kinetic assay for metabolic byproduct formation

The two metabolic byproducts monitored in cultures are Fe³⁺ from Fe²⁺ oxidation and H⁺ associated with sulfur oxidation. Their optical indicators are absorbance at 380 nm and the thymol blue pH ratio, respectively. Abiotic standards characterize how each indicator responds to its chemical target.

This collection provides the underlying measurements and experimental annotations for testing the tolerance of some iron- and/or sulfur-oxidizing acidophiles to some more (non-)metals<sup>1</sup> using microplate-based absorbance measurements. It includes indicator calibration, medium validation, and flask and test-tube culture experiments. The documentation identifies preparations, conditions, experimental units, reference measurements and sampling times so that the measurement records can be interpreted together with their experimental context.

See the [dataset information](../README.md#dataset-information) for citation order, contacts, licensing, and funding.

## Organization

The `data/` directory contains a shared materials-and-strains folder and five experimental folders. Each family has a README and YAML documentation alongside its data. Folder number prefixes control display order; experimental identifiers define the relationships among measurements.

The shared [materials and strains catalog](00_materials_and_strains/README.md) provides tables and machine-readable data for chemicals, materials, equipment, strains and medium formulations. Reagent identities are distinct from experiment-specific concentrations and concentration-determination methods. Experiment-specific assignments are documented with each experimental family.

| Folder within `data/` | Contents | Documentation to read first |
|---|---|---|
| [00_materials_and_strains](00_materials_and_strains/) | Shared chemicals, materials, equipment, strains and media formulations, with tables and YAML records. | [README.md](00_materials_and_strains/README.md). |
| [01_indicator_calibration](01_indicator_calibration/) | Thymol-blue pH and iron absorbance calibrations, spectra, configuration and storage measurements, and paired lid/no-lid readings. | The family [README.md](01_indicator_calibration/README.md), followed by the pH, iron or lid-comparison YAML. |
| [02_media_validation](02_media_validation/) | Abiotic measurements in iron- and sulfur-containing media, with pH and iron reference values. | [README.md](02_media_validation/README.md) and [media_validation.yaml](02_media_validation/media_validation.yaml). |
| [03_biological_microplates](03_biological_microplates/) | AF, AT, and LF culture measurements, additive concentration series and corresponding abiotic wells. | [README.md](03_biological_microplates/README.md), [biological_definitions.yaml](03_biological_microplates/biological_definitions.yaml), then the relevant organism's experiment YAML. |
| [04_flask_validation](04_flask_validation/) | AF flask cultures and abiotic controls, including citric acid–gallium combinations, chemistry measurements and fluorescence. | [README.md](04_flask_validation/README.md) and [flask_validation.yaml](04_flask_validation/flask_validation.yaml). |
| [05_test_tube_validation](05_test_tube_validation/) | Acidophile and abiotic tube cultures, pH/ORP time courses, endpoint fluorescence and documented starter-culture observations. | [README.md](05_test_tube_validation/README.md) and [test_tube_validation.yaml](05_test_tube_validation/test_tube_validation.yaml). |

Organism codes and strain designations are defined in the shared [strain table](00_materials_and_strains/strains.md). Experimental metadata select the applicable strain and preserve starter-culture identities.

## Reading the collection

The collection contains 7,416 `.xls` reader exports, 1 `.xlsx` reader workbook and 4 `.csv` measurement tables, accompanied by 11 Markdown documents, 15 YAML files, 2 chronology CSV tables and 1 time-course reader. CSV files can be opened in spreadsheet applications; Markdown and YAML files can be read in a text editor. Follow each family’s README for the reader-export format and data dictionary.

Start with the README for the experimental family of interest. Its YAML documents contain named fields and lists that can be read in a text editor. They define the applicable methods, units, concentrations, well assignments, replication and identifiers. Shared definitions use the [documented `$ref` convention](00_materials_and_strains/README.md#reading-shared-yaml-references); retain `00_materials_and_strains` with any experimental family. The CSV column definitions describe the corresponding data tables.

Reader exports retain their supplied formats. The family READMEs explain how their measurement blocks and wavelength labels are organized. The indicator-calibration and media-validation families each have a `reference_measurements.csv`; these are separate reference collections and should be used with their own family's YAML. `flask_measurements.csv` and `test_tube_measurements.csv` contain 1 row per flask or tube and sampling event, as defined in their documentation. Test-tube rows also explicitly identify shared initial-reference assignments.

Use identifiers to connect measurements to their definitions. Physical plate, flask and tube identifiers connect each experimental unit to its time-course observations. Starter-culture identities and reference assignments are defined in the family metadata.

Biological acquisition times and condition assignments are provided in [chronology](03_biological_microplates/chronology/README.md). The acquisition registry connects conditions and physical plates to their reader exports in chronological order.

## Values and units

CSV tables use a period as the decimal separator. Units and column definitions are given in the corresponding YAML files.

Use the stated unit and concentration basis for each experiment. Stock, preparation and final-well concentrations are defined separately. Elemental-analysis concentrations account for analytical dilution. Initial flask chemistry values are assigned from measured starting-medium references; the test-tube initial-reference rows use the assignments in their YAML.

All nominal stock concentrations were calculated from the mass added and the final prepared volume in the [documented stock-preparation vessel](00_materials_and_strains/chemicals_materials_equipment.md#materials-and-equipment). Subsequent mixing and dilution volumes determine the stated preparation or final-culture concentration; nominal concentrations are distinct from elemental measurements.

All ORP values are reported in mV versus Ag/AgCl with 3 M KCl reference electrolyte.

All deposited fluorescence values have the TE Buffer background-blank subtracted. The fluorescence methods and bibliography are provided in the flask and test-tube documentation.

Sampling times and time origins are defined in each experimental family. Recorded and interpolated times are identified in the metadata.

## Contributors

| Initials | Name | Affiliation | ORCID |
|---|---|---|---|
| [DAD](https://orcid.org/0000-0003-1908-7196) | [Daniel A. Dailey](https://orcid.org/0000-0003-1908-7196) | Department of Chemical Engineering, Columbia University | [0000-0003-1908-7196](https://orcid.org/0000-0003-1908-7196) |
| [MAH](https://orcid.org/0009-0009-0165-2098) | [Mallori A. Herishko](https://orcid.org/0009-0009-0165-2098) | Department of Chemical Engineering, Columbia University | [0009-0009-0165-2098](https://orcid.org/0009-0009-0165-2098) |
| [WW](https://orcid.org/0009-0002-5905-2318) | [Wei Wang](https://orcid.org/0009-0002-5905-2318) | Department of Chemical Engineering, Columbia University | [0009-0002-5905-2318](https://orcid.org/0009-0002-5905-2318) |
| [LNM](https://orcid.org/0009-0000-5291-5200) | [Lauren N. Murray](https://orcid.org/0009-0000-5291-5200) | Department of Chemical Engineering, Columbia University | [0009-0000-5291-5200](https://orcid.org/0009-0000-5291-5200) |
| [STB](https://orcid.org/0009-0007-4652-5733) | [Sofia T. Bavaresco](https://orcid.org/0009-0007-4652-5733) | Department of Chemical Engineering, Columbia University | [0009-0007-4652-5733](https://orcid.org/0009-0007-4652-5733) |
| [NRL](https://orcid.org/0000-0002-3824-9122) | [Natalie R. Ling](https://orcid.org/0000-0002-3824-9122) | Department of Chemical Engineering, Columbia University | [0000-0002-3824-9122](https://orcid.org/0000-0002-3824-9122) |
| [SVRA](https://orcid.org/0009-0002-9372-3226) | [Sebastian Valentino Rojas Abregu](https://orcid.org/0009-0002-9372-3226) | Department of Biomedical Engineering, Columbia University | [0009-0002-9372-3226](https://orcid.org/0009-0002-9372-3226) |
| [SB](https://orcid.org/0000-0001-7885-0150) | [Scott Banta](https://orcid.org/0000-0001-7885-0150) | Department of Chemical Engineering, Columbia University | [0000-0001-7885-0150](https://orcid.org/0000-0001-7885-0150) |

The following roles summarize contributions to the dataset using the CRediT taxonomy.<sup>2</sup> Experimental design is represented by Methodology; experimental preparation and measurements are represented by Investigation. The family READMEs specify contributions to each experimental collection.

| CRediT role | Contributors |
|---|---|
| Conceptualization | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) |
| Methodology | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [LNM](https://orcid.org/0009-0000-5291-5200), [STB](https://orcid.org/0009-0007-4652-5733) |
| Investigation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [LNM](https://orcid.org/0009-0000-5291-5200), [STB](https://orcid.org/0009-0007-4652-5733), [NRL](https://orcid.org/0000-0002-3824-9122), [SVRA](https://orcid.org/0009-0002-9372-3226) |
| Validation | [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098), [WW](https://orcid.org/0009-0002-5905-2318) |
| Data curation | [DAD](https://orcid.org/0000-0003-1908-7196), [WW](https://orcid.org/0009-0002-5905-2318) |
| Software | [WW](https://orcid.org/0009-0002-5905-2318), [DAD](https://orcid.org/0000-0003-1908-7196) |
| Resources | [SB](https://orcid.org/0000-0001-7885-0150) |
| Supervision | [SB](https://orcid.org/0000-0001-7885-0150) |
| Project administration | [SB](https://orcid.org/0000-0001-7885-0150) |
| Funding acquisition | [SB](https://orcid.org/0000-0001-7885-0150) |

[DAD](https://orcid.org/0000-0003-1908-7196) conceived the use of thymol blue as the pH indicator in this assay, and [MAH](https://orcid.org/0009-0009-0165-2098) conceived the Fe(III) absorbance indicator approach.

[DAD](https://orcid.org/0000-0003-1908-7196) contributed experimental design across all 5 families; [DAD](https://orcid.org/0000-0003-1908-7196) and [MAH](https://orcid.org/0009-0009-0165-2098) co-designed media validation, and [MAH](https://orcid.org/0009-0009-0165-2098) designed the final ferric-indicator series. [DAD](https://orcid.org/0000-0003-1908-7196), [LNM](https://orcid.org/0009-0000-5291-5200) and [STB](https://orcid.org/0009-0007-4652-5733) co-designed the AF/AT soluble-sulfur series in SM4 basal salts.

[MAH](https://orcid.org/0009-0009-0165-2098) and [DAD](https://orcid.org/0000-0003-1908-7196) prepared calibration and media-validation experiments; [MAH](https://orcid.org/0009-0009-0165-2098) collected their optical data, with [SVRA](https://orcid.org/0009-0002-9372-3226) and [NRL](https://orcid.org/0000-0002-3824-9122) contributing to the final iron studies and [NRL](https://orcid.org/0000-0002-3824-9122) to the lid comparison. [DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [NRL](https://orcid.org/0000-0002-3824-9122) contributed to AFM1 screening; [LNM](https://orcid.org/0009-0000-5291-5200) and [DAD](https://orcid.org/0000-0003-1908-7196) prepared and measured the AFM4 and soluble-sulfur series, while [MAH](https://orcid.org/0009-0009-0165-2098) and [LNM](https://orcid.org/0009-0000-5291-5200) measured the LF series prepared by [DAD](https://orcid.org/0000-0003-1908-7196) and [MAH](https://orcid.org/0009-0009-0165-2098). [DAD](https://orcid.org/0000-0003-1908-7196) performed the flask chemical measurements, and [DAD](https://orcid.org/0000-0003-1908-7196) and [NRL](https://orcid.org/0000-0002-3824-9122) measured flask fluorescence. [DAD](https://orcid.org/0000-0003-1908-7196), [STB](https://orcid.org/0009-0007-4652-5733) and [SVRA](https://orcid.org/0009-0002-9372-3226) prepared the tube cultures, with [STB](https://orcid.org/0009-0007-4652-5733), [SVRA](https://orcid.org/0009-0002-9372-3226), [DAD](https://orcid.org/0000-0003-1908-7196) and [LNM](https://orcid.org/0009-0000-5291-5200) contributing measurements and [STB](https://orcid.org/0009-0007-4652-5733) and [LNM](https://orcid.org/0009-0000-5291-5200) measuring endpoint fluorescence; [DAD](https://orcid.org/0000-0003-1908-7196) carried out the original LF starter follow-up.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [WW](https://orcid.org/0009-0002-5905-2318) contributed to assay validation using the calibration, medium-validation, biological microplate, flask and test-tube datasets.

[DAD](https://orcid.org/0000-0003-1908-7196) and [WW](https://orcid.org/0009-0002-5905-2318) curated the data and experimental metadata, and [WW](https://orcid.org/0009-0002-5905-2318) and [DAD](https://orcid.org/0000-0003-1908-7196) developed the supporting software. [SB](https://orcid.org/0000-0001-7885-0150) provided resources and supervision, administered the project and acquired funding.

## References

1. Tuovinen, O. H.; Niemelä, S. I.; Gyllenberg, H. G. Tolerance of *Thiobacillus ferrooxidans* to Some Metals. *Antonie van Leeuwenhoek* 1971, *37*, 489–496. [DOI: 10.1007/BF02218519](https://doi.org/10.1007/BF02218519).

2. Elsevier. [CRediT Author Statement](https://www.elsevier.com/researcher/author/policies-and-guidelines/credit-author-statement) (accessed 2026-09-12).

