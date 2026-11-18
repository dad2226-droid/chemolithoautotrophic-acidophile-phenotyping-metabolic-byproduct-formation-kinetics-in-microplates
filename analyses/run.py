"""Calculate manuscript analyses from the neighboring data directory."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import tempfile
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

from acidophile_data.provenance import assert_sources_unchanged, source_manifest, write_tables


ANALYSES = Path(__file__).resolve().parent
PACKAGES = ("numpy", "pandas", "scipy", "PyYAML", "openpyxl")


def code_checksums() -> dict[str, str]:
    paths = [*ANALYSES.glob("*.py"), *ANALYSES.joinpath("acidophile_data").glob("*.py"),
             *ANALYSES.joinpath("metadata").glob("*"), ANALYSES / "requirements.txt"]
    return {str(path.relative_to(ANALYSES)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths) if path.is_file()}


def calculate(section: str, root: Path):
    if section == "optical":
        from optics_workflow import run
    else:
        from biology_workflow import run
    return run(root)


def run_analysis(root: Path, output: Path, sections: list[str]) -> None:
    data = (root / "data").resolve()
    output = output.resolve()
    if not data.is_dir():
        raise ValueError(f"Missing data directory: {data}")
    if output.is_relative_to(data) or output.exists():
        raise ValueError("Choose a new output directory outside data/: " + str(output))
    output.parent.mkdir(parents=True, exist_ok=True)
    before = source_manifest(root)
    started = datetime.now(timezone.utc)
    with tempfile.TemporaryDirectory(prefix=".analysis-", dir=output.parent) as tmp:
        bundle = Path(tmp) / "results"
        bundle.mkdir()
        for section in sections:
            print(f"Calculating {section} analyses", flush=True)
            tables = calculate(section, root)
            guide = tables["result_guide"]
            if (guide.empty or guide.table.duplicated().any()
                    or set(guide.table) != set(tables) - {"result_guide"}):
                raise ValueError(f"Incomplete or duplicate entries in the {section} result guide")
            index = write_tables(tables, bundle / section)
            print(f"Wrote {len(index)} {section} CSV files", flush=True)
            del tables
        assert_sources_unchanged(before, source_manifest(root))
        before.to_csv(bundle / "source_files.csv", index=False)
        record = {
            "started_utc": started.isoformat(), "completed_utc": datetime.now(timezone.utc).isoformat(),
            "sections": sections, "python": platform.python_version(),
            "dependencies": {name: version(name) for name in PACKAGES},
            "code_and_specification_sha256": code_checksums(),
            "data_manifest_sha256": hashlib.sha256((bundle / "source_files.csv").read_bytes()).hexdigest(),
        }
        (bundle / "run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        bundle.rename(output)
    print(f"Results: {output}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--section", choices=("all", "optical", "biological"), default="all",
                        help="Select which analyses to calculate (default: all)")
    parser.add_argument("--output", type=Path,
                        help="New results directory (default: a dated folder under analyses/results/)")
    args = parser.parse_args()
    sections = ["optical", "biological"] if args.section == "all" else [args.section]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run_analysis(ANALYSES.parent, args.output or ANALYSES / "results" / stamp, sections)


if __name__ == "__main__":
    main()
