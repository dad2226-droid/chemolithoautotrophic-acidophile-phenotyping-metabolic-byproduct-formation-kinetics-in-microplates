# Chemolithoautotrophic acidophile phenotyping via a miniaturized kinetic assay for metabolic byproduct formation

This repository contains the experimental data, analysis code, and figures and tables for the manuscript. The 96-well microplate assay tracks two metabolic byproducts of substrate oxidation. Absorbance at 380 nm indicates apparent ferric iron (Fe³⁺) generation during Fe²⁺ oxidation. The thymol blue absorbance ratio indicates acidification associated with H⁺ release during sulfur oxidation. Chemical standards establish how the indicators respond, and measurements from shake flasks and culture tubes provide independent comparisons for the microplate cultures.

## Dataset information

| Field | Information |
|---|---|
| Dataset title | Chemolithoautotrophic acidophile phenotyping via a miniaturized kinetic assay for metabolic byproduct formation |
| Dataset creators, in citation order | [Daniel A. Dailey](https://orcid.org/0000-0003-1908-7196); [Wei Wang](https://orcid.org/0009-0002-5905-2318); [Scott Banta](https://orcid.org/0000-0001-7885-0150) |
| Corresponding contacts | [Scott Banta](https://orcid.org/0000-0001-7885-0150); [Daniel A. Dailey](https://orcid.org/0000-0003-1908-7196) |
| Contact emails, in the same order | sbanta@columbia.edu; dad2226@columbia.edu |
| Data and documentation license | [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) (accessed 2026-09-12). SPDX identifier `CC-BY-4.0`. |
| Funders | U.S. National Science Foundation (NSF) |
| Grants and awards | U.S. National Science Foundation: DGE-2437839 |

## Repository contents

| Folder | Contents |
|---|---|
| [data](data/README.md) | Instrument exports, chemical reference measurements, culture measurements, and experimental records. Start here to look up a strain, medium, treatment, sampling time, or experimental unit. |
| [analyses](analyses/README.md) | Python workflows for indicator calibration, measurement precision, culture time courses, and comparisons across plates, flasks, and tubes. The guide maps each analysis to its manuscript figure panels. |
| [visualizations](visualizations/README.md) | Captioned figures, Markdown tables, plotted-value CSVs, and code to generate them from the analysis results and experimental records. |

## Manuscript figures and tables

- [Figures 1–8](visualizations/figures/README.md) show the artwork with its full panel caption. Each figure has a single-page PDF and a CSV of the plotted values.
- [Main-text Tables 1–6](visualizations/main_tables/README.md) describe the experimental design. Select an experiment ID within a table to find the measurement files and metadata for each component.
- [Supplementary Tables S1–S32](visualizations/supplementary_tables/README.md) give the numerical results supporting the Technical Validation claims and identify the associated figure panels.

## Reproduce the results

From the repository root, use a Python 3.12 environment and run the commands below. Choose a new directory for the analysis results. If you already have a complete run, pass its directory to the final command.

```sh
python -m pip install -r visualizations/requirements.txt
python analyses/run.py --output analyses/results/manuscript
python visualizations/run.py --results analyses/results/manuscript
```

The analysis excludes individual absorbance readings above 4 AU and instrument-flagged over-range readings before calculating indicators. An unflagged reading of exactly 4 AU is eligible. Each analysis section writes a `result_guide.csv` that maps figure panels to the numerical outputs.

## CRediT contribution statement

The contributors below are credited for the dataset. [Daniel A. Dailey](https://orcid.org/0000-0003-1908-7196) and [Mallori A. Herishko](https://orcid.org/0009-0009-0165-2098) contributed equally to the manuscript.

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

For this dataset, Methodology covers experimental design, and Investigation covers experimental preparation and measurements. The role names follow the [CRediT taxonomy](https://credit.niso.org/).

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

In the experimental contributions below, AF denotes *Acidithiobacillus ferrooxidans*, AT denotes *Acidithiobacillus thiooxidans*, and LF denotes *Leptospirillum ferriphilum*.

[DAD](https://orcid.org/0000-0003-1908-7196) conceived the use of thymol blue as the pH indicator in this assay, and [MAH](https://orcid.org/0009-0009-0165-2098) conceived the Fe(III) absorbance indicator approach.

[DAD](https://orcid.org/0000-0003-1908-7196) contributed experimental design across all 5 families. [DAD](https://orcid.org/0000-0003-1908-7196) and [MAH](https://orcid.org/0009-0009-0165-2098) co-designed media validation, and [MAH](https://orcid.org/0009-0009-0165-2098) designed the final ferric-indicator series. [DAD](https://orcid.org/0000-0003-1908-7196), [LNM](https://orcid.org/0009-0000-5291-5200) and [STB](https://orcid.org/0009-0007-4652-5733) co-designed the AF/AT soluble-sulfur series in SM4 basal salts.

[MAH](https://orcid.org/0009-0009-0165-2098) and [DAD](https://orcid.org/0000-0003-1908-7196) prepared calibration and media-validation experiments. [MAH](https://orcid.org/0009-0009-0165-2098) collected their optical data, with [SVRA](https://orcid.org/0009-0002-9372-3226) and [NRL](https://orcid.org/0000-0002-3824-9122) contributing to the final iron studies and [NRL](https://orcid.org/0000-0002-3824-9122) to the lid comparison.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [NRL](https://orcid.org/0000-0002-3824-9122) contributed to AFM1 screening. [LNM](https://orcid.org/0009-0000-5291-5200) and [DAD](https://orcid.org/0000-0003-1908-7196) prepared and measured the AFM4 and soluble-sulfur series, while [MAH](https://orcid.org/0009-0009-0165-2098) and [LNM](https://orcid.org/0009-0000-5291-5200) measured the LF series prepared by [DAD](https://orcid.org/0000-0003-1908-7196) and [MAH](https://orcid.org/0009-0009-0165-2098).

[DAD](https://orcid.org/0000-0003-1908-7196) performed the flask chemical measurements, and [DAD](https://orcid.org/0000-0003-1908-7196) and [NRL](https://orcid.org/0000-0002-3824-9122) measured flask fluorescence. [DAD](https://orcid.org/0000-0003-1908-7196), [STB](https://orcid.org/0009-0007-4652-5733) and [SVRA](https://orcid.org/0009-0002-9372-3226) prepared the tube cultures, with [STB](https://orcid.org/0009-0007-4652-5733), [SVRA](https://orcid.org/0009-0002-9372-3226), [DAD](https://orcid.org/0000-0003-1908-7196) and [LNM](https://orcid.org/0009-0000-5291-5200) contributing measurements and [STB](https://orcid.org/0009-0007-4652-5733) and [LNM](https://orcid.org/0009-0000-5291-5200) measuring endpoint fluorescence. [DAD](https://orcid.org/0000-0003-1908-7196) carried out the LF starter follow-up.

[DAD](https://orcid.org/0000-0003-1908-7196), [MAH](https://orcid.org/0009-0009-0165-2098) and [WW](https://orcid.org/0009-0002-5905-2318) contributed to assay validation using the calibration, medium-validation, biological microplate, flask and test-tube datasets.

[DAD](https://orcid.org/0000-0003-1908-7196) and [WW](https://orcid.org/0009-0002-5905-2318) curated the data and experimental metadata, and [WW](https://orcid.org/0009-0002-5905-2318) and [DAD](https://orcid.org/0000-0003-1908-7196) developed the supporting software. [SB](https://orcid.org/0000-0001-7885-0150) provided resources and supervision, administered the project and acquired funding.

## Funding

[Daniel A. Dailey](https://orcid.org/0000-0003-1908-7196) was partially supported by the National Science Foundation Graduate Research Fellowship Program (Grant No. DGE-2437839), the Tau Beta Pi Fellowship, and the Bourbon County Educational Fund.

Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views or policies of the National Science Foundation or the U.S. Government.

## License
The authors license the original data, analysis and plotting scripts, figures, tables, and documentation in this deposit under the Creative Commons Attribution 4.0 International License (CC BY 4.0).