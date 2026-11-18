"""Plate reader."""
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
import numpy as np


ROWS = "ABCDEFGH"


N_COLUMNS = 12


@dataclass(frozen=True)
class PlateData:
    """One plate-reader export."""

    path: Path
    mode: str
    wavelengths_nm: tuple[float, ...]
    values: np.ndarray  # wavelength x row x column
    exported_mask: np.ndarray
    parse_error_mask: np.ndarray
    censored_high_mask: np.ndarray
    temperature_c: float | None
    saved_at: datetime | None
    original_filename: str | None

    def wavelength_index(self, wavelength_nm: float, tolerance_nm: float = 0.01) -> int:
        distances = np.abs(np.asarray(self.wavelengths_nm) - wavelength_nm)
        index = int(np.argmin(distances))
        if distances[index] > tolerance_nm:
            raise KeyError(
                f"{self.path}: wavelength {wavelength_nm:g} nm is absent; "
                f"available={self.wavelengths_nm}"
            )
        return index

    def value(self, well: str, wavelength_nm: float) -> float:
        wavelength_index, row, column = self._indices(well, wavelength_nm)
        return float(self.values[wavelength_index, row, column])

    def _indices(self, well: str, wavelength_nm: float) -> tuple[int, int, int]:
        match = re.fullmatch(r"([A-Ha-h])(1[0-2]|[1-9])", well.strip())
        if not match:
            raise ValueError(f"Invalid 96-well identifier: {well!r}")
        row = ROWS.index(match.group(1).upper())
        column = int(match.group(2)) - 1
        return self.wavelength_index(wavelength_nm), row, column

    def is_parse_error(self, well: str, wavelength_nm: float) -> bool:
        wavelength_index, row, column = self._indices(well, wavelength_nm)
        return bool(self.parse_error_mask[wavelength_index, row, column])

    def is_censored_high(self, well: str, wavelength_nm: float) -> bool:
        wavelength_index, row, column = self._indices(well, wavelength_nm)
        return bool(self.censored_high_mask[wavelength_index, row, column])


def _decode(path: Path) -> str:
    raw = path.read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        return raw.decode("utf-16")
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        # A few reader versions omit the UTF-16 byte-order mark.
        return raw.decode("utf-16-le")


def _cell(value: str) -> tuple[float, bool, bool]:
    """Return numeric value, parse-error flag, and high-censor flag."""

    text = value.strip()
    upper = text.upper()
    if not text or upper in {"NAN", "NA", "N/A"}:
        return np.nan, False, False
    if upper in {"OVER", "OVERFLW", "***"}:
        return np.nan, False, True
    if upper == "#ERR":
        return np.nan, True, False
    try:
        return float(text.replace(",", "")), False, False
    except ValueError:
        return np.nan, True, False


def _number(value: str) -> float:
    return _cell(value)[0]


def _footer_metadata(text: str) -> tuple[str | None, datetime | None]:
    match = re.search(
        r"Original Filename:\s*(.*?);\s*Date Last Saved:\s*(.*?)\s*$",
        text,
        flags=re.MULTILINE,
    )
    if not match:
        return None, None
    original = match.group(1).strip()
    timestamp = match.group(2).strip()
    parsed = None
    for fmt in (
        "%m/%d/%Y %I:%M:%S %p",
        "%m/%d/%Y %I:%M %p",
        "%m/%d/%Y %H:%M:%S",
    ):
        try:
            parsed = datetime.strptime(timestamp, fmt)
            break
        except ValueError:
            pass
    return original, parsed


def _parse_endpoint(path: Path, text: str, lines: list[str], header: list[str]) -> PlateData:
    wavelengths = tuple(float(x) for x in header[15].split())
    first_column = int(header[16])
    column_count = int(header[17])
    first_row = int(header[19])
    row_count = int(header[20])
    values = np.full((len(wavelengths), len(ROWS), N_COLUMNS), np.nan, dtype=float)
    exported_mask = np.zeros(values.shape, dtype=bool)
    parse_error_mask = np.zeros(values.shape, dtype=bool)
    censored_high_mask = np.zeros(values.shape, dtype=bool)

    candidate_lines: list[list[str]] = []
    for line in lines[3:]:
        if line.startswith("~End"):
            break
        cells = line.split("\t")
        # A row of OVER/#ERR/NA is still an exported row. Requiring a
        # finite numeric cell would drop it and shift every later well.
        if any(cell.strip() for cell in cells[2:]):
            candidate_lines.append(cells)
    if len(candidate_lines) != row_count:
        raise ValueError(
            f"{path}: expected {row_count} populated plate rows, found {len(candidate_lines)}"
        )

    temperature = _number(candidate_lines[0][1])
    for local_row, cells in enumerate(candidate_lines[:row_count]):
        plate_row = first_row - 1 + local_row
        for wavelength_index in range(len(wavelengths)):
            group_start = 2 + wavelength_index * (N_COLUMNS + 1)
            for plate_column in range(first_column - 1, first_column - 1 + column_count):
                cell_index = group_start + plate_column
                if cell_index < len(cells):
                    value, parse_error, censored_high = _cell(cells[cell_index])
                    values[wavelength_index, plate_row, plate_column] = value
                    exported_mask[wavelength_index, plate_row, plate_column] = True
                    parse_error_mask[wavelength_index, plate_row, plate_column] = parse_error
                    censored_high_mask[wavelength_index, plate_row, plate_column] = censored_high

    original, saved_at = _footer_metadata(text)
    return PlateData(
        path=path,
        mode="Endpoint",
        wavelengths_nm=wavelengths,
        values=values,
        exported_mask=exported_mask,
        parse_error_mask=parse_error_mask,
        censored_high_mask=censored_high_mask,
        temperature_c=float(temperature) if np.isfinite(temperature) else None,
        saved_at=saved_at,
        original_filename=original,
    )


def _parse_spectrum(path: Path, text: str, lines: list[str], header: list[str]) -> PlateData:
    first_column = int(header[16])
    column_count = int(header[17])
    first_row = int(header[19])
    row_count = int(header[20])
    wavelengths: list[float] = []
    blocks: list[list[list[str]]] = []
    index = 3
    while index < len(lines):
        if lines[index].startswith("~End"):
            break
        cells = lines[index].split("\t")
        wavelength = _number(cells[0] if cells else "")
        if not np.isfinite(wavelength):
            index += 1
            continue
        block = [cells]
        for offset in range(1, row_count):
            if index + offset >= len(lines):
                break
            block.append(lines[index + offset].split("\t"))
        if len(block) != row_count:
            raise ValueError(f"{path}: incomplete spectrum block at {wavelength:g} nm")
        wavelengths.append(float(wavelength))
        blocks.append(block)
        index += row_count

    expected_count = int(header[8])
    if len(blocks) != expected_count:
        raise ValueError(
            f"{path}: expected {expected_count} spectrum blocks, found {len(blocks)}"
        )
    values = np.full((len(wavelengths), len(ROWS), N_COLUMNS), np.nan, dtype=float)
    exported_mask = np.zeros(values.shape, dtype=bool)
    parse_error_mask = np.zeros(values.shape, dtype=bool)
    censored_high_mask = np.zeros(values.shape, dtype=bool)
    temperature = _number(blocks[0][0][1]) if blocks else np.nan
    for wavelength_index, block in enumerate(blocks):
        for local_row, cells in enumerate(block):
            plate_row = first_row - 1 + local_row
            for plate_column in range(first_column - 1, first_column - 1 + column_count):
                cell_index = 2 + plate_column
                if cell_index < len(cells):
                    value, parse_error, censored_high = _cell(cells[cell_index])
                    values[wavelength_index, plate_row, plate_column] = value
                    exported_mask[wavelength_index, plate_row, plate_column] = True
                    parse_error_mask[wavelength_index, plate_row, plate_column] = parse_error
                    censored_high_mask[wavelength_index, plate_row, plate_column] = censored_high

    original, saved_at = _footer_metadata(text)
    return PlateData(
        path=path,
        mode="Spectrum",
        wavelengths_nm=tuple(wavelengths),
        values=values,
        exported_mask=exported_mask,
        parse_error_mask=parse_error_mask,
        censored_high_mask=censored_high_mask,
        temperature_c=float(temperature) if np.isfinite(temperature) else None,
        saved_at=saved_at,
        original_filename=original,
    )


def read_plate(path: str | Path) -> PlateData:
    """Read one plate-reader export.

    Raises a descriptive exception for malformed exports.  Missing or overflow
    cells are retained as ``NaN`` rather than silently removed.
    """

    path = Path(path)
    text = _decode(path)
    lines = text.splitlines()
    if len(lines) < 3 or not lines[0].startswith("##BLOCKS="):
        raise ValueError(f"{path}: not a recognized plate-reader export")
    header = lines[1].split("\t")
    if len(header) < 21 or header[3] != "PlateFormat":
        raise ValueError(f"{path}: malformed PlateFormat header")
    mode = header[4]
    if mode == "Endpoint":
        return _parse_endpoint(path, text, lines, header)
    if mode == "Spectrum":
        return _parse_spectrum(path, text, lines, header)
    raise ValueError(f"{path}: unsupported reader mode {mode!r}")
