"""Biological chronology."""
from __future__ import annotations
from collections import defaultdict
import csv
from datetime import datetime
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


DIRECTORY = Path('data/03_biological_microplates/chronology')


ACQUISITION_FIELDS = {
    'acquisition_id', 'physical_plate_id', 'species_id', 'plate_id', 'batch_id',
    'native_original_filename', 'native_read_label', 'timestamp_local_iso',
    'date_last_saved_text', 'timestamp_timezone',
    'read_index_zero_based', 'elapsed_seconds', 'elapsed_hours', 'elapsed_days',
    'plate_first_timestamp_local_iso', 'path', 'file_sha256', 'metadata_path',
}


CONDITION_FIELDS = {
    'timecourse_id', 'experiment_id', 'group_id', 'physical_plate_id', 'species_id',
    'plate_id', 'column', 'column_conditions_json', 'concentration_unit',
    'volume_basis_uL', 'all_column_component_values_zero', 'fixed_components_json',
    'medium_id', 'initial_OD600_200_uL', 'initial_pH_json', 'nominal_initial_pH_json',
    'preparation_id', 'readout_id', 'biotic_rows_json', 'abiotic_rows_json',
    'number_of_acquisitions', 'first_timestamp_local_iso',
    'last_timestamp_local_iso', 'period_days',
}


def _read_table(path: Path, required: set[str]) -> list[dict]:
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'{path.name}: required chronology columns missing: {sorted(missing)}')
        rows = list(reader)
    if not rows:
        raise ValueError(f'{path.name}: empty chronology table')
    return rows


def _unique(rows: list[dict], field: str) -> None:
    values = [row[field] for row in rows]
    if any(not value for value in values) or len(values) != len(set(values)):
        raise ValueError(f'Chronology {field} values must be nonempty and unique')


def _path(repository: Path, relative: str) -> Path:
    candidate = Path(relative)
    if candidate.is_absolute() or '..' in candidate.parts:
        raise ValueError(f'Expected a safe deposit-relative path: {relative}')
    result = (repository / candidate).resolve()
    if not result.is_relative_to((repository / 'data/03_biological_microplates').resolve()):
        raise ValueError(f'Chronology path leaves biological data: {relative}')
    return result


@lru_cache(maxsize=8192)
def _hash_file(path: str, size: int, mtime_ns: int) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _check_hash(path: Path, expected: str) -> None:
    state = path.stat()
    actual = _hash_file(str(path), state.st_size, state.st_mtime_ns)
    if actual != expected:
        raise ValueError(f'Export checksum mismatch: {path}')


def _json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(',', ':'), default=str)


def _check_value(row: dict, field: str, expected, *, numeric=False, json_value=False) -> None:
    value = json.loads(row[field]) if json_value else float(row[field]) if numeric else row[field]
    matches = (np.isfinite(value) and np.isclose(value, expected, rtol=1e-12, atol=1e-12)
               if numeric else value == expected)
    if not matches:
        identity = row.get('timecourse_id', row.get('acquisition_id'))
        raise ValueError(f'Chronology/YAML mismatch for {identity}: {field}')


def load_chronology(repository: Path, definitions: dict, documents: list,
                    read_export) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return verified acquisitions and condition assignments for these species."""
    repository = repository.resolve()
    directory = repository / DIRECTORY
    acquisitions = _read_table(directory / 'acquisitions.csv', ACQUISITION_FIELDS)
    conditions = _read_table(directory / 'condition_timecourse_index.csv', CONDITION_FIELDS)
    for field in ('acquisition_id', 'path'):
        _unique(acquisitions, field)
    _unique(conditions, 'timecourse_id')
    known_species = set(definitions['species'])
    if ({r['species_id'] for r in acquisitions + conditions} - known_species):
        raise ValueError('Chronology contains an unknown species')
    selected = {document['species_id'] for _, document in documents}
    acquisitions = [row for row in acquisitions if row['species_id'] in selected]
    conditions = [row for row in conditions if row['species_id'] in selected]
    by_plate = defaultdict(list)
    physical_ids = {}
    for row in acquisitions:
        key = (row['species_id'], row['plate_id'])
        physical_id = row['physical_plate_id']
        if not physical_id or physical_id in physical_ids and physical_ids[physical_id] != key:
            raise ValueError('Chronology physical_plate_id must identify exactly one microplate')
        physical_ids[physical_id] = key
        by_plate[key].append(row)

    expected_conditions = {}
    expected_plates = set()
    records = []
    for document_path, document in documents:
        species = document['species_id']
        registry_reference = document.get('acquisition_registry')
        if not registry_reference or (document_path.parent / registry_reference).resolve() != directory / 'acquisitions.csv':
            raise ValueError(f'{document_path}: unresolved acquisition_registry')
        memberships = defaultdict(list)
        for experiment in document['experiments']:
            for group in experiment['plate_groups']:
                for plate_id in group['plate_ids']:
                    memberships[plate_id].append({'experiment_id': experiment['experiment_id'],
                                                 'group_id': group['group_id']})
                    for column, condition in group['column_conditions'].items():
                        key = (species, experiment['experiment_id'], group['group_id'], plate_id, int(column))
                        if key in expected_conditions:
                            raise ValueError(f'Duplicate YAML condition assignment: {key}')
                        expected_conditions[key] = (experiment, group, condition)
        for plate_id, plate in document['plates'].items():
            key = (species, plate_id)
            expected_plates.add(key)
            if not plate.get('batch_id') or plate_id not in memberships:
                raise ValueError(f'{plate_id}: missing batch or experiment membership')
            reads = sorted(by_plate.get(key, []), key=lambda row: int(row['read_index_zero_based']))
            if not reads:
                raise ValueError(f'{plate_id}: no registered acquisitions')
            if [int(r['read_index_zero_based']) for r in reads] != list(range(len(reads))):
                raise ValueError(f'{plate_id}: nonconsecutive chronological indices')
            if len({r['physical_plate_id'] for r in reads}) != 1:
                raise ValueError(f'{plate_id}: multiple physical_plate_id values')
            if plate.get('acquisition_ids') != [r['acquisition_id'] for r in reads]:
                raise ValueError(f'{plate_id}: acquisition IDs/order differ between YAML and chronology')
            timestamps = [datetime.fromisoformat(r['timestamp_local_iso']) for r in reads]
            if any(b <= a for a, b in zip(timestamps, timestamps[1:])):
                raise ValueError(f'{plate_id}: timestamps do not strictly increase')
            first = timestamps[0]
            cultures = {row: f"{reads[0]['physical_plate_id']}_{row}"
                        for row in definitions['plate_layout']['default_biotic_rows']}
            for row, timestamp in zip(reads, timestamps):
                _check_value(row, 'batch_id', plate['batch_id'])
                _check_value(row, 'metadata_path', document_path.relative_to(repository).as_posix())
                _check_value(row, 'plate_first_timestamp_local_iso', first.isoformat())
                elapsed = (timestamp - first).total_seconds()
                for field, expected in [('elapsed_seconds', elapsed), ('elapsed_hours', elapsed / 3600),
                                        ('elapsed_days', elapsed / 86400)]:
                    _check_value(row, field, expected, numeric=True)
                path = _path(repository, row['path'])
                if not path.is_file():
                    raise FileNotFoundError(f'Missing registered acquisition: {row["path"]}')
                _check_hash(path, row['file_sha256'])
                if row['acquisition_id'] != 'ACQ_' + row['file_sha256'][:24]:
                    raise ValueError(f'Acquisition ID/checksum mismatch: {row["acquisition_id"]}')
                reading = read_export(path)
                if reading.saved_at != timestamp:
                    raise ValueError(f'Export acquisition timestamp differs from chronology: {row["path"]}')
                if reading.original_filename != row['native_original_filename']:
                    raise ValueError(f'Export native label differs from chronology: {row["path"]}')
                records.append({
                    'species': species, 'species_name': definitions['species'][species]['scientific_name'],
                    'plate_id': plate_id, 'physical_plate_id': row['physical_plate_id'],
                    'biological_replicate_id': row['physical_plate_id'],
                    'biological_replicate_unit': definitions['culture_preparation']['biological_replicate_unit'],
                    'culture_id_scope': definitions['plate_layout']['biotic_culture_number_scope'],
                    'batch_id': plate['batch_id'], 'acquisition_id': row['acquisition_id'],
                    'acquisition_label': row['native_original_filename'],
                    'native_original_filename': row['native_original_filename'],
                    'native_read_label': row['native_read_label'],
                    'native_ordinal': int(row['native_read_label']) if row['native_read_label'].isdigit() else np.nan,
                    'ordinal': int(row['read_index_zero_based']),
                    'read_index_zero_based': int(row['read_index_zero_based']),
                    'source': row['path'], 'file_sha256': row['file_sha256'],
                    'source_metadata': row['metadata_path'],
                    'timestamp': timestamp, 'time_zero': first,
                    'timestamp_timezone': row['timestamp_timezone'],
                    'date_last_saved_text': row['date_last_saved_text'],
                    'days': elapsed / 86400, 'hours': elapsed / 3600,
                    'reader_temperature_C': reading.temperature_c,
                    'incubation_temperature_C': definitions['species'][species]['incubation_temperature_C'],
                    'culture_ids_json': _json(cultures),
                    'memberships_json': _json(memberships[plate_id]),
                    'blank_condition_columns_json': _json(plate.get('blank_condition_columns', [])),
                    'wavelengths_nm_json': _json(reading.wavelengths_nm),
                    'plate_has_ordinal_zero': True, 'baseline_ordinal': 0,
                    'plate_has_native_ordinal_zero': any(r['native_read_label'] == '0' for r in reads),
                })
        actual_paths = {p.resolve() for p in document_path.parent.rglob('*.xls')}
        expected_paths = {_path(repository, r['path']) for r in acquisitions if r['species_id'] == species}
        if actual_paths != expected_paths:
            raise ValueError(f'{species}: biological export inventory differs from acquisition registry')
    if set(by_plate) != expected_plates:
        raise ValueError('Chronology/YAML microplate inventory differs')

    indexed_keys = set()
    for row in conditions:
        key = (row['species_id'], row['experiment_id'], row['group_id'], row['plate_id'], int(row['column']))
        if key in indexed_keys or key not in expected_conditions:
            raise ValueError(f'Unexpected or duplicate condition-index assignment: {key}')
        indexed_keys.add(key)
        experiment, group, condition = expected_conditions[key]
        reads = sorted(by_plate[(row['species_id'], row['plate_id'])], key=lambda r: int(r['read_index_zero_based']))
        for field, expected in [('physical_plate_id', reads[0]['physical_plate_id']),
                                ('concentration_unit', group['concentration_unit']),
                                ('medium_id', experiment['medium_id']),
                                ('preparation_id', experiment['preparation_id']),
                                ('readout_id', experiment['readout_id']),
                                ('first_timestamp_local_iso', reads[0]['timestamp_local_iso']),
                                ('last_timestamp_local_iso', reads[-1]['timestamp_local_iso'])]:
            _check_value(row, field, expected)
        for field, expected in [('column_conditions_json', condition),
                                ('fixed_components_json', experiment.get('fixed_components', [])),
                                ('initial_pH_json', experiment.get('initial_pH')),
                                ('nominal_initial_pH_json', experiment.get('nominal_initial_pH')),
                                ('biotic_rows_json', definitions['plate_layout']['default_biotic_rows']),
                                ('abiotic_rows_json', definitions['plate_layout']['default_abiotic_rows'])]:
            _check_value(row, field, expected, json_value=True)
        for field, expected in [('volume_basis_uL', group['volume_basis_uL']),
                                ('initial_OD600_200_uL', experiment['initial_OD600_200_uL']),
                                ('number_of_acquisitions', len(reads)),
                                ('period_days', float(reads[-1]['elapsed_days']))]:
            _check_value(row, field, expected, numeric=True)
        if row['all_column_component_values_zero'].lower() != str(all(float(v) == 0 for v in condition.values())).lower():
            raise ValueError(f'Zero-addition flag differs from YAML: {row["timecourse_id"]}')
    if indexed_keys != set(expected_conditions):
        raise ValueError('Condition index omits YAML experiment/microplate/column assignments')
    registry = pd.DataFrame(records).sort_values(['species', 'plate_id', 'ordinal']).reset_index(drop=True)
    index = pd.DataFrame(conditions)
    for field in ('column', 'number_of_acquisitions'):
        index[field] = index[field].astype(int)
    return registry, index.sort_values(['species_id', 'experiment_id', 'group_id', 'plate_id', 'column']).reset_index(drop=True)
