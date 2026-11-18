"""Load chronological biological condition records using the acquisition registry.

Uses only the Python standard library. Source exports are checked against their
SHA256 values and are never modified.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import hashlib
from pathlib import Path
import sys


CHRONOLOGY_DIRECTORY = Path('data/03_biological_microplates/chronology')
ACQUISITION_FIELDS = {
    'acquisition_id', 'physical_plate_id', 'species_id', 'plate_id', 'batch_id',
    'native_original_filename', 'native_read_label', 'date_last_saved_text',
    'timestamp_local_iso', 'timestamp_timezone', 'read_index_zero_based',
    'elapsed_seconds', 'elapsed_hours', 'elapsed_days',
    'plate_first_timestamp_local_iso', 'path', 'file_sha256', 'metadata_path',
}
CONDITION_FIELDS = {
    'timecourse_id', 'experiment_id', 'physical_plate_id', 'species_id',
    'plate_id', 'column', 'number_of_acquisitions',
}


def default_root() -> Path:
    """Return the deposit root containing this script."""
    return Path(__file__).resolve().parents[3]


def read_table(path: Path, required: set[str]) -> list[dict]:
    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        unresolved_fields = required - set(reader.fieldnames or [])
        if unresolved_fields:
            raise ValueError(f'{path.name}: Required columns: {sorted(unresolved_fields)}')
        return list(reader)


def export_path(root: Path, relative_path: str) -> Path:
    relative = Path(relative_path)
    if relative.is_absolute():
        raise ValueError(f'Expected a deposit-relative export path: {relative_path}')
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root):
        raise ValueError(f'Export path leaves the deposit: {relative_path}')
    return resolved


def load_timecourses(experiment: str, root: str | Path | None = None,
                     plate: str | None = None,
                     column: int | None = None) -> list[dict]:
    """Return ordered condition/acquisition rows with source hashes checked.

    ``plate`` is an exact ``plate_id`` from the condition index or organism YAML.
    The join uses ``physical_plate_id``; native filename labels are not parsed.
    Returned ``path`` values are relative to the deposit root.
    """
    root = Path(root).resolve() if root is not None else default_root()
    directory = root / CHRONOLOGY_DIRECTORY
    acquisitions = read_table(directory / 'acquisitions.csv', ACQUISITION_FIELDS)
    conditions = read_table(directory / 'condition_timecourse_index.csv', CONDITION_FIELDS)
    selected = [row for row in conditions
                if row['experiment_id'] == experiment
                and (plate is None or row['plate_id'] == plate)
                and (column is None or int(row['column']) == column)]
    if not selected:
        raise ValueError('No condition records match the experiment, plate and column.')

    selected_ids = {row['physical_plate_id'] for row in selected}
    by_plate = collections.defaultdict(list)
    acquisition_ids = set()
    for row in acquisitions:
        if row['physical_plate_id'] not in selected_ids:
            continue
        if row['acquisition_id'] in acquisition_ids:
            raise ValueError(f'Acquisition IDs must be unique: {row["acquisition_id"]}')
        acquisition_ids.add(row['acquisition_id'])
        row['read_index_zero_based'] = int(row['read_index_zero_based'])
        row['elapsed_seconds'] = int(row['elapsed_seconds'])
        for field in ['elapsed_hours', 'elapsed_days']:
            row[field] = float(row[field])
        by_plate[row['physical_plate_id']].append(row)

    # Verify every physical export, even when several column views share it.
    for physical_id in selected_ids:
        reads = by_plate[physical_id]
        if not reads:
            raise ValueError(f'No acquisitions for physical plate {physical_id}')
        reads.sort(key=lambda row: row['read_index_zero_based'])
        if [r['read_index_zero_based'] for r in reads] != list(range(len(reads))):
            raise ValueError(f'Nonconsecutive chronological indices: {physical_id}')
        timestamps = [datetime.datetime.fromisoformat(r['timestamp_local_iso']) for r in reads]
        if any(b <= a for a, b in zip(timestamps, timestamps[1:])):
            raise ValueError(f'Timestamps do not strictly increase: {physical_id}')
        for row, timestamp in zip(reads, timestamps):
            elapsed = int((timestamp - timestamps[0]).total_seconds())
            if row['elapsed_seconds'] != elapsed:
                raise ValueError(f'Elapsed time mismatch: {row["acquisition_id"]}')
            path = export_path(root, row['path'])
            if hashlib.sha256(path.read_bytes()).hexdigest() != row['file_sha256']:
                raise ValueError(f'Export checksum mismatch: {row["path"]}')

    result = []
    joins = set()
    for condition in selected:
        reads = by_plate[condition['physical_plate_id']]
        if int(condition['number_of_acquisitions']) != len(reads):
            raise ValueError(f'Acquisition count mismatch: {condition["timecourse_id"]}')
        for acquisition in reads:
            for field in ['physical_plate_id', 'plate_id', 'species_id']:
                if condition[field] != acquisition[field]:
                    raise ValueError(f'Identity mismatch in {condition["timecourse_id"]}: {field}')
            key = (condition['timecourse_id'], acquisition['acquisition_id'])
            if key in joins:
                raise ValueError(f'Condition/acquisition keys must be unique: {key}')
            joins.add(key)
            result.append({**condition, **acquisition,
                           'column': int(condition['column']),
                           'number_of_acquisitions': len(reads)})
    result.sort(key=lambda row: (row['experiment_id'], row['plate_id'],
                                 row['column'], row['read_index_zero_based']))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=None,
                        help='Deposit root; defaults to the root containing this script.')
    parser.add_argument('--experiment', required=True,
                        help='Exact experiment_id from the organism YAML or condition index.')
    parser.add_argument('--plate',
                        help='Optional exact plate_id from the YAML or condition index.')
    parser.add_argument('--column', type=int, choices=range(1, 13),
                        help='Optional physical plate column.')
    parser.add_argument('--output', type=Path,
                        help='Output CSV path; omitted means write CSV to standard output.')
    args = parser.parse_args()
    try:
        rows = load_timecourses(args.experiment, args.root, args.plate, args.column)
        if args.output:
            with args.output.open('w', encoding='utf-8', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
