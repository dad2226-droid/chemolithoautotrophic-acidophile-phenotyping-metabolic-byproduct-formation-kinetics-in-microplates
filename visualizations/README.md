# Figures and tables

| Manuscript item | Contents |
|---|---|
| [Figures 1–8](figures/README.md) | Each figure appears with its full panel caption, a downloadable PDF and the plotted values. |
| [Main-text Tables 1–6](main_tables/README.md) | Each experimental overview appears with its caption and notes. Select an experiment ID in the table to reach the underlying records for each component. |
| [Supplementary Tables S1–S32](supplementary_tables/README.md) | Each table presents the numerical results alongside the Technical Validation claim they support and identifies the associated figure panels. |

Each table is a GitHub Flavored Markdown file containing its caption and pipe-delimited data. The six experimental overview tables include component identifiers and relative links into the bundled data folder. Each figure’s plotted-values CSV identifies the panel, series and type of observation or summary, including uncertainty bounds and inset values where applicable.

## Reproduce the files

From the project root, use Python 3.12 and run the commands below. The example uses a new analysis output directory. For an existing complete run, pass its path to `--results` in the final command. The figure renderer uses Arial when installed and selects an available sans-serif font otherwise.

```sh
python -m pip install -r visualizations/requirements.txt
python analyses/run.py --output analyses/results/manuscript
python visualizations/run.py --results analyses/results/manuscript
```

To generate a separate copy, add `--output path/to/output` to the last command. The [analysis guide](../analyses/README.md) describes the measurements, statistical methods, eligibility rules and full numerical outputs.
