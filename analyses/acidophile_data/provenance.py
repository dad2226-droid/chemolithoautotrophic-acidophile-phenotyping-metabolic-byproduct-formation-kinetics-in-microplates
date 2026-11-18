"""Provenance."""
from __future__ import annotations
import hashlib
from pathlib import Path
import pandas as pd


def source_manifest(root: Path) -> pd.DataFrame:
    """Hash every deposited file, including its experimental documentation."""
    root = Path(root).resolve()
    data = root / "data"
    if not data.is_dir():
        raise ValueError(f"Missing canonical data directory: {data}")
    rows = []
    for path in sorted(data.rglob("*")):
        if not path.is_file() or path.name == ".DS_Store":
            continue
        if path.is_symlink() or not path.resolve().is_relative_to(data.resolve()):
            raise ValueError(f"Input symlinks are not supported: {path}")
        payload = path.read_bytes()
        rows.append({"source": path.relative_to(root).as_posix(),
                     "size_bytes": len(payload),
                     "sha256": hashlib.sha256(payload).hexdigest()})
    if not rows:
        raise ValueError("The deposited data directory contains no files")
    return pd.DataFrame(rows)


def assert_sources_unchanged(before: pd.DataFrame, after: pd.DataFrame) -> None:
    if not before.equals(after):
        raise ValueError("Deposited source files changed during the analysis")


def write_tables(tables: dict[str, pd.DataFrame], destination: Path) -> pd.DataFrame:
    destination.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, frame in sorted(tables.items()):
        if not name or Path(name).name != name:
            raise ValueError(f"Unsafe output table name: {name}")
        path = destination / f"{name}.csv"
        temporary = path.with_name(path.name + ".partial")
        frame.to_csv(temporary, index=False, lineterminator="\n")
        temporary.replace(path)
        rows.append({"table": name, "file": path.name, "rows": len(frame),
                     "columns": len(frame.columns),
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    index = pd.DataFrame(rows)
    temporary = destination / "table_index.csv.partial"
    index.to_csv(temporary, index=False)
    temporary.replace(destination / "table_index.csv")
    return index
