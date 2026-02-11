# ANN Rapid Predictor — Run Scripts

This directory contains the main prediction script (`dti_ann_LUT.py`) that
uses a trained ANN to predict neural fiber activation thresholds.

## Usage

```bash
python run/dti_ann_LUT.py <electrode_file> <tract_file> <model_path> <output_json> <centering> <tract_type> <X> <Y> <Z> <mode>
```

### Arguments

| # | Argument | Description |
|---|----------|-------------|
| 1 | `electrode_file` | Path to the electrode voltage data file (FEM export). |
| 2 | `tract_file` | Path to the tractography file (one fiber per line). |
| 3 | `model_path` | Path to a trained ANN model directory (e.g. `models/ann_19_reg`). |
| 4 | `output_json` | Output path for the results JSON file. |
| 5 | `centering` | Centering strategy: `ssd` (second spatial derivative) or `ec` (extracellular potential). |
| 6 | `tract_type` | Type of tract data: `dti` (real DTI tracts) or `artificial`. |
| 7–9 | `X Y Z` | **Electrode center coordinates** — where the electrode should be placed in the fiber coordinate space. The FEM bounding box is re-centered here. Alternatively, pass a single string `anisotropic` or `isotropic` to use the legacy default positions. |
| 10 | `mode` | Prediction mode: `reg` (regression) or `class` (classification). Must match the model type. |

### Electrode Center Placement

When you provide explicit X, Y, Z coordinates, the script:
1. Reads the FEM bounding box from the electrode file.
2. Re-centers the FEM grid so the middle of the bounding box is at the given coordinates.
3. Applies an internal offset (`fem_offset`) so that fiber positions are correctly translated back to the original FEM grid for voltage lookups.

This lets you use the same electrode FEM file with any tractography file, regardless of the coordinate system the tracts were generated in.

### Examples

**Custom electrode placement** (recommended for arbitrary tract files):
```bash
python run/dti_ann_LUT.py \
    electrodes/medtronic_3387/monopolar/3387_anisotropic_monopolar_01-23.txt \
    example_tracks/whole_brain.txt \
    models/ann_19_reg \
    run/results.json \
    ssd dti 38 50 30 reg
```
This places the electrode center at coordinates (38, 50, 30) in the fiber's
coordinate space.

**Legacy conductivity mode** (for files already aligned with the FEM grid):
```bash
python run/dti_ann_LUT.py \
    electrodes/medtronic_3387/monopolar/3387_anisotropic_monopolar_01-23.txt \
    example_tracks/L_DRTT_voxel.txt \
    models/ann_17_class \
    run/results.json \
    ssd dti anisotropic class
```

### Output

The script writes a JSON file containing:
- `valid_inds` — indices of fibers that passed all filtering criteria.
- `problem_inds` — indices of fibers that were skipped.
- Per-pulse-width threshold predictions keyed by pulse width in seconds (e.g. `"0.06"`).

## Supporting Script — process_DTI.py

`process_DTI.py` handles tract processing, spline interpolation, and filtering.
It is called internally by `dti_ann_LUT.py` and does not need to be run
directly. It supports the `custom_center` parameter for electrode placement and
uses list comprehensions for ragged-array operations.
