# Biological acquisition chronology

These files connect biological conditions to the chronological measurements of each physical microplate. They cover the deposited AF, AT and LF reader exports.

## Files

| File | Contents |
|---|---|
| `acquisitions.csv` | 1 row per deposited reader export: acquisition and physical plate identities, native label, timestamp, chronological index, elapsed time, file path and SHA256. |
| `condition_timecourse_index.csv` | 1 row per experiment, physical plate and assigned column, with condition values, units, fixed components, well rows and measurement context. |
| `load_timecourses.py` | Python standard-library loader that joins the 2 tables and verifies the selected source files. |

Paths in the acquisition registry are relative to the directory containing `data/`. The organism YAML files define the experimental groups and their column assignments. Values and units in the condition index retain their documented preparation basis.

## Identities and ordering

`acquisition_id` is a file key formed from `ACQ_` and the first 24 hexadecimal characters of the export's SHA-256 checksum. The full checksum is stored in `file_sha256`. The registry records the acquisition time and chronological position separately.

`physical_plate_id` identifies a physical plate independently of its acquisition labels. `plate_id` identifies that plate in the organism YAML and condition index. Join on `physical_plate_id` to associate acquisitions with physical plates.

`native_original_filename` and `native_read_label` give the labels recorded by the instrument. `read_index_zero_based` orders the measurements of each physical plate from 0.

The recorded `Date Last Saved` values are acquisition timestamps because exports were saved immediately after acquisition. Times use the recorded local clock. Elapsed seconds, hours and days are measured from the first deposited acquisition of each physical plate. Filenames combine the native label and timestamp; ordering comes from the registry.

An export appears 1 time in `acquisitions.csv`. The same acquisition can support several column timecourses or experiment views. Shared control columns remain associated with each relevant comparison. The `fixed_components_json` and `medium_id` fields specify components present throughout each experiment.

The loader returns all acquisitions for the selected experiment, plate and column.

## Load a timecourse

Python 3.10 or later is required. No additional packages are needed. Run the following command from the directory containing `data/`.

```bash
python data/03_biological_microplates/chronology/load_timecourses.py \
  --experiment AF_AFM1_potassium_iodide_OD0p001 \
  --plate HT3-21-A --column 2 --output timecourse.csv
```

Omit `--plate` or `--column` to include all corresponding plate or column timecourses for the experiment. `--experiment` accepts an exact `experiment_id`; `--plate` accepts an exact YAML `plate_id`. Omit `--output` to write CSV to standard output. Use `--root` to specify the directory containing `data/` if the loader is stored elsewhere.

The output contains condition fields joined to acquisition fields. Rows are grouped by experiment, plate and column, then ordered by chronological index. Each selected export's path and SHA256 are checked 1 time. The loader also checks acquisition counts, identities, chronological indices, timestamps and elapsed seconds.

The function can be imported when this directory is on the Python import path:

```python
from pathlib import Path
from load_timecourses import load_timecourses

root = Path.cwd()  # Run from the directory containing data/.
rows = load_timecourses(
    experiment='AF_AFM1_potassium_iodide_OD0p001',
    root=root,
    plate='HT3-21-A',
    column=2,
)
raw_export = (root / rows[0]['path']).read_text(encoding='utf-16')
```

The loader reads and verifies exports.
