# Manuscript analyses

Run the analyses in this folder to calculate the numerical results for Figures 1–8 and the accompanying technical validation from the files in [`../data/`](../data/). The workflow reads the instrument exports, chemical reference measurements, culture records and experimental metadata directly. The culture analyses use absorbance at 380 nm and the thymol blue pH ratio as indicators of the metabolic byproducts Fe³⁺ and H⁺, respectively. Its outputs are CSV files that can be used for statistical review and plotting.

## Run from the repository root

Use Python 3.12. The commands below create a virtual environment, install the versioned dependencies and run both sections.

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r analyses/requirements.txt
python analyses/run.py
```

On Windows PowerShell, activate the environment with `.venv\Scripts\Activate.ps1`. To run one section or choose a destination, use either of these commands. Each destination must be a new directory outside `data/`.

```sh
python analyses/run.py --section optical --output analyses/results/optical-run
python analyses/run.py --section biological --output analyses/results/biological-run
```

## Contents

| File or directory | Use |
|---|---|
| [`run.py`](run.py) | Run the selected analyses, write the result tables and record input and code checksums. |
| [`optics_workflow.py`](optics_workflow.py) | Prepare optical standards and references, fit calibration models, validate predictions and evaluate handling, storage and medium effects. |
| [`biology_workflow.py`](biology_workflow.py) | Prepare microplate, shake-flask and culture-tube data and calculate biological screening, repeatability and cross-format comparisons. |
| [`biological_sensitivity.py`](biological_sensitivity.py) | Calculate treatment-effect agreement and rank correlations across batches, acquisition-schedule sensitivity and matched citrate endpoint contrasts. The biological workflow calls these functions. |
| [`acidophile_data/`](acidophile_data/) | Read the source data and implement the indicator calculations, statistical models, eligibility rules and plate-, vessel- and batch-level summaries used by the workflows. |
| [`metadata/`](metadata/) | Define the fluorescence conversion used for shake flasks and the baseline conditions used for the soluble-sulfur comparisons. |
| [`requirements.txt`](requirements.txt) | Install the Python packages and versions used for the analysis. |
| [`tests/`](tests/) | Check absorbance cutoffs, missing observations, batch-level treatment ranks and endpoint matching on small examples. Run with `python -m unittest discover -s analyses/tests`. |

## Analyses by manuscript figure

| Section | Figure | Measurements and results |
|---|---|---|
| Optical | 1a–d | Measured and scatter-corrected thymol-blue ratios, four-parameter logistic pH calibration, correction effects on span and precision, withheld-preparation pH prediction and transfer to ferric mixtures. |
| Optical | 1e–h | Ferric–ferrous spectral contrasts, 380/428 nm comparisons, acid-specific Weibull calibration, withheld-plate and withheld-concentration prediction, adjacent-standard resolution and cerimetry versus total-iron regression. |
| Optical | 2–3 | Lid, plate material, temperature and membrane contrasts, storage changes, complete-standard sensitivity, edge–interior differences and availability of corrected ratios. |
| Optical | 4 | Response span, pooled within-standard SD, span-to-SD ratio, reversals, paired scatter-correction effects and dye-background contrasts at 485, 545 and 750 nm. |
| Biological | 5 | Apparent oxidation screen, signed area under the curve (AUC) to the first observed maximum, time to that maximum, between-batch variability, variance tests, treatment-effect repeatability and sensitivity to the acquisition schedule. |
| Biological | 6 and 7a–f | Measured and scatter-corrected absorbance, fluorescence-derived apparent planktonic cell density, matched citrate endpoints, final net 380/428 nm contrasts, dose matching, and plate–flask comparisons of response, timing and trajectory shape. |
| Biological | 7g–i and 8 | *L. ferriphilum* endpoints with time-matched *A. ferrooxidans* benchmarks, fresh-starter tube contrasts, soluble-sulfur indicator trajectories, endpoint changes relative to abiotic variation and paired ratio eligibility. |

## Read the results

The default output is a new dated directory under `analyses/results/`. It contains `optical/` and/or `biological/` for the selected sections. In each section, begin with **`result_guide.csv`** to find tables by analysis or figure panel. **`table_index.csv`** records the filename, row and column counts, and SHA-256 for each CSV. The run also writes `source_files.csv` with input-file hashes and `run.json` with software versions, code and specification hashes, selected sections and execution times. The destination becomes available when the selected calculations finish and the input integrity check passes.

Technical wells contribute to within-standard precision or are averaged within biological microplates. Biological comparisons retain the identities of independently inoculated microplates, flasks, tubes and batches. For the primary iron comparison, microplate biotic change in absorbance and shake-flask apparent ferric generation are evaluated within their own formats. Net comparisons use the matched abiotic condition where specified. Missing or ineligible measurements can yield blank statistics. Use the accompanying status and sample-size fields when interpreting those rows. Vanadyl-interfered titrations are ineligible for quantitative iron comparisons. Fluorescence-derived values are apparent OD-equivalents based on the *A. ferrooxidans* calibration.

For control-relative AUC, the workflow requires all assigned control areas to be estimable and their mean to be positive. The separate iron-screen criterion requires a biotic increase of at least 1 AU and a matched abiotic increase below 1 AU at the same acquisition. Citrate endpoint comparisons use zero citrate at the same gallium concentration, including doses measured only in microplates.

The absorbance cutoff applies to each measured channel before calculating ratios, differences or trajectory summaries. Unflagged numeric readings of exactly 4 AU are eligible. Readings above 4 AU and instrument-flagged over-range readings are excluded.
