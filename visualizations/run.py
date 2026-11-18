"""Build manuscript figures, experimental overview tables and Supplementary tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import tempfile
from collections.abc import Mapping
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib import _mathtext, font_manager, pyplot as plt
from matplotlib.mathtext import MathTextParser
from matplotlib.text import _get_text_metrics_with_cache_impl
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "analyses"))

COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#332288", "#000000"]


class ArialScriptOffsets(_mathtext.FontConstantsBase):
    """Set subscript and superscript baselines for the small chemical labels."""

    sup1 = .385
    sub1 = .14
    sub2 = .245


def configure_style():
    font_path = font_manager.findfont(
        font_manager.FontProperties(family=["Arial", "Liberation Sans", "DejaVu Sans"])
    )
    family = font_manager.FontProperties(fname=font_path).get_name()
    if family == "Arial" and _mathtext._font_constant_mapping.get(family) is not ArialScriptOffsets:
        _mathtext._font_constant_mapping[family] = ArialScriptOffsets
        MathTextParser._parse_cached.cache_clear()
        _get_text_metrics_with_cache_impl.cache_clear()
    plt.rcParams.update({
        "font.family": family, "font.size": 6, "axes.titlesize": 7, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "legend.title_fontsize": 6, "legend.frameon": False,
        "axes.linewidth": .5, "lines.linewidth": .8,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.width": .5, "ytick.major.width": .5, "xtick.minor.width": .5,
        "ytick.minor.width": .5, "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.minor.size": 1.6, "ytick.minor.size": 1.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.fontset": "custom", "mathtext.rm": family, "mathtext.it": family + ":italic",
        "mathtext.bf": family + ":bold", "mathtext.cal": family, "mathtext.sf": family, "mathtext.tt": family,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "axes.prop_cycle": plt.cycler(color=COLORS),
    })


class AnalysisTables(Mapping):
    """Read each named CSV on demand and verify it against its analysis index."""

    def __init__(self, results):
        self.entries = {}
        self.frames = {}
        for section in ("optical", "biological"):
            directory = results / section
            with (directory / "table_index.csv").open(encoding="utf-8", newline="") as stream:
                for row in csv.DictReader(stream):
                    key = row["table"]
                    if key == "result_guide":
                        continue
                    if key in self.entries or row["file"] != key + ".csv" or Path(key).name != key:
                        raise ValueError(f"Invalid or duplicate analysis index entry: {key}")
                    self.entries[key] = (directory / row["file"], row["sha256"])

    def __getitem__(self, key):
        if key not in self.frames:
            path, expected = self.entries[key]
            with path.open("rb") as stream:
                if hashlib.file_digest(stream, "sha256").hexdigest() != expected:
                    raise ValueError(f"Analysis table checksum mismatch: {path}")
                stream.seek(0)
                self.frames[key] = pd.read_csv(stream, low_memory=False, float_precision="round_trip")
        return self.frames[key]

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)


def render(results, output):
    from build_main_tables import build as write_main_tables
    from captioned_figure import artwork_bounds, save_captioned_pdf
    from plotted_values import write_plotted_values
    from supplementary import write_tables
    configure_style()
    tables = AnalysisTables(results)
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output) as tmp:
        bundle = Path(tmp)
        (bundle / "README.md").write_bytes((HERE / "README.md").read_bytes())
        figures = bundle / "figures"
        figures.mkdir()
        (figures / "README.md").write_bytes((HERE / "figures" / "README.md").read_bytes())
        for number in range(1, 9):
            first = (number - 1) // 2 * 2 + 1
            module = importlib.import_module(f"figure_{first}_{first + 1}")
            figure = getattr(module, f"figure_{number}")(tables)
            try:
                bounds = artwork_bounds(figure)
                stem = figures / f"Figure_{number}"
                caption = HERE / "figures" / f"Figure_{number}.md"
                save_captioned_pdf(figure, caption, stem.with_suffix(".pdf"), bounds)
                figure.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches=bounds, pad_inches=0)
                write_plotted_values(figure, number, stem.with_suffix(".csv"))
                stem.with_suffix(".md").write_bytes(caption.read_bytes())
            finally:
                plt.close(figure)
        write_tables(results, bundle / "supplementary_tables")
        write_main_tables(HERE.parent / "data", bundle / "main_tables",
                          link_directory=output / "main_tables")
        for file in sorted(p for p in bundle.rglob("*") if p.is_file()):
            target = output / file.relative_to(bundle)
            target.parent.mkdir(exist_ok=True, parents=True)
            file.replace(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True,
                        help="Directory containing both the optical and biological analysis outputs.")
    parser.add_argument("--output", type=Path, default=HERE,
                        help="Directory to receive figures/, main_tables/ and supplementary_tables/ (default: visualizations/).")
    args = parser.parse_args()
    results = args.results.resolve()
    output = args.output.resolve()
    if output.is_relative_to(HERE.parent / "data") or output.is_relative_to(results):
        parser.error("Choose an output directory outside data/ and the analysis results.")
    render(results, output)


if __name__ == "__main__":
    main()
