"""Metadata."""
from __future__ import annotations
from pathlib import Path
from typing import Any
import yaml


def load_metadata(path: str | Path) -> Any:
    """Resolve ``{$ref: relative.yaml#/path/to/value}`` without changing source files.

    Reference objects have exactly one key. Paths are relative to the file that
    contains the reference; fragments use JSON Pointer escaping. References stay
    within the enclosing data deposit (or the input directory for standalone
    metadata). Each call reads a fresh snapshot and returns independent objects.
    """
    source = Path(path).resolve()
    boundary = next((p for p in source.parents if p.name == "data"), source.parent)
    documents: dict[Path, Any] = {}

    def document(filename: Path) -> Any:
        if filename not in documents:
            documents[filename] = yaml.safe_load(filename.read_text(encoding="utf-8-sig"))
        return documents[filename]

    def select(filename: Path, pointer: str) -> Any:
        value = document(filename)
        if not pointer:
            return value
        if not pointer.startswith("/"):
            raise ValueError(f"Invalid metadata JSON pointer: {filename}#{pointer}")
        try:
            for segment in pointer[1:].split("/"):
                key = segment.replace("~1", "/").replace("~0", "~")
                if isinstance(value, list):
                    if not key.isdigit() or (key != "0" and key.startswith("0")):
                        raise ValueError("Invalid array index")
                    value = value[int(key)]
                else:
                    value = value[key]
        except (KeyError, IndexError, TypeError, ValueError) as error:
            raise ValueError(f"Unresolved metadata reference: {filename}#{pointer}") from error
        return value

    def resolve(value: Any, filename: Path, active: tuple) -> Any:
        if isinstance(value, dict):
            if "$ref" in value:
                reference = value["$ref"]
                if len(value) != 1 or not isinstance(reference, str):
                    raise ValueError(f"A metadata reference must be a sole string $ref: {filename}")
                relative, _, pointer = reference.partition("#")
                if ":" in relative or Path(relative).is_absolute():
                    raise ValueError(f"Metadata reference must use a local relative path: {reference}")
                target = (filename.parent / relative).resolve() if relative else filename
                if not target.is_relative_to(boundary):
                    raise ValueError(f"Metadata reference escapes the data directory: {reference}")
                identity = (target, pointer)
                if identity in active:
                    raise ValueError(f"Cyclic metadata reference: {target}#{pointer}")
                return resolve(select(target, pointer), target, (*active, identity))
            return {key: resolve(child, filename, active) for key, child in value.items()}
        if isinstance(value, list):
            return [resolve(child, filename, active) for child in value]
        return value

    return resolve(document(source), source, ((source, ""),))
