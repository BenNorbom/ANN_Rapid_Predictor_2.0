# Visualization Scripts

This directory contains three visualization scripts for tractography activation results.
All produce interactive 3D HTML plots rendered with Plotly.

---

## plot_tracts_fast.py

Fast 3D visualization designed to handle 100,000+ fibers efficiently by merging
all fiber geometry into just two WebGL traces (activated and inactive).

### Usage

```bash
python graphing/plot_tracts_fast.py --tract <TRACT_FILE> --results <RESULTS_JSON> --output <OUTPUT_DIR> [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--tract FILE` | Path to the tract text file (one fiber per line: `x1 y1 z1 x2 y2 z2 ...`). |
| `--results FILE` | Path to the JSON results file produced by `dti_ann_LUT.py`. |
| `--output DIR` | Directory where HTML plots will be saved. |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--activation_threshold VOLTS` | `3.0` | Voltage (V) below which a fiber is considered **activated**. |
| `--electrode_center X Y Z` | `0 0 0` | X Y Z position of the electrode for display. |
| `--electrode_config CONFIG` | none | Contact configuration string (e.g. `01-23`, `+012-3`, `-0+1-2+3`). `-` = cathode (red), `+` = anode (blue), unmarked = inactive (grey). Each sign applies only to the immediately following digit. |
| `--downsample N` | `1` | Keep every Nth point per fiber to speed up rendering. |
| `--show_axes` | off | Show X/Y/Z axes in the 3D scene. |
| `--all_fibers` | off | Plot ALL fibers at full length, ignoring the `valid_inds` filter. Fibers without a threshold are drawn as inactive. |

### Electrode Rendering

The script renders a detailed 3D mesh electrode model (Medtronic 3387 geometry)
with four contacts, a rounded tip, insulating gaps, and a shaft. Contacts are
coloured according to `--electrode_config`:
- **Contact 0** = tip (distal end)
- **Contact 3** = shaft (proximal end)

### Pulse Width Auto-Detection

The script automatically reads whatever pulse widths are present in the results
JSON file. There is no need to specify them — if the simulation ran one pulse
width, one plot is generated; if it ran all 17, you get 17 plots.

### Examples

```bash
# Basic usage
python graphing/plot_tracts_fast.py \
    --tract example_tracks/whole_brain.txt \
    --results run/results.json \
    --output output_viz/

# Custom settings with electrode config
python graphing/plot_tracts_fast.py \
    --tract example_tracks/whole_brain.txt \
    --results run/results.json \
    --output output_viz/ \
    --activation_threshold 2.5 \
    --electrode_center 38 50 30 \
    --electrode_config 01-23 \
    --downsample 3

# Show all fibers regardless of FEM bounds
python graphing/plot_tracts_fast.py \
    --tract example_tracks/whole_brain.txt \
    --results run/results.json \
    --output output_viz/ \
    --all_fibers --show_axes
```

---

## plot_tracts.py

Original visualization script with per-fiber rendering and detailed 3D electrode
mesh. Better suited for smaller tract files or when you need single pulse-width
interactive plots.

### Usage

```bash
python graphing/plot_tracts.py --tract <TRACT_FILE> --results <RESULTS_JSON> --output <OUTPUT_DIR> [options]
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--tract FILE` | required | Path to the fiber tract file. |
| `--results FILE` | required | Path to the results JSON file. |
| `--output DIR` | required | Directory to save output files. |
| `--voltage VOLTS` | `3.0` | Activation threshold voltage. |
| `--cond TYPE` | `anisotropic` | Conductivity type: `anisotropic` or `isotropic`. |
| `--show_axes` | off | Show XYZ axes. |
| `--filter_indices FILE` | none | Path to a file listing the indices of fibers that were simulated. |
| `--all_fibers` | off | Plot ALL fibers, ignoring the filter_indices file. |
| `--interactive_pw INDEX` | none | Render a single pulse width (index) interactively. |
| `--electrode_config CONFIG` | none | Electrode contact configuration string. |

---

## plot_tracts_bundles.py

Bundle-aware visualization that colours each tract bundle differently.
Requires a manifest JSON (produced by `move_whole_brain_tracts.py`) to identify
which fibers belong to which bundle.

### Usage

```bash
python graphing/plot_tracts_bundles.py \
    --tract <MERGED_TRACT_FILE> --results <RESULTS_JSON> \
    --manifest <MANIFEST_JSON> --output <OUTPUT_DIR> [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--tract FILE` | Path to the merged tract file (whole-brain + bundles). |
| `--results FILE` | Path to the JSON results file from `dti_ann_LUT.py`. |
| `--manifest FILE` | Path to the manifest JSON from `move_whole_brain_tracts.py`. |
| `--output DIR` | Directory where HTML plots will be saved. |

### Optional Arguments

Same as `plot_tracts_fast.py`: `--activation_threshold`, `--electrode_center`,
`--electrode_config`, `--downsample`, `--show_axes`, `--all_fibers`.

### Bundle Colour Scheme

| Bundle | Activated | Inactive |
|--------|-----------|----------|
| Whole-brain | Red | Black (transparent) |
| DRTT | Lime green | Dark green |
| ML | Cyan | Teal |
| PRT | Magenta | Purple |

### Example

```bash
python graphing/plot_tracts_bundles.py \
    --tract merged_tracts.txt \
    --results run/results.json \
    --manifest merged_tracts_manifest.json \
    --output output_viz_whole_brain/ \
    --electrode_center 38 50 30 \
    --electrode_config 01-23 \
    --all_fibers --show_axes
```

---

## Viewing the Output

All scripts generate `.html` files in your output directory (e.g.
`activation_pw_00_60us.html`).

**Option A: Local server (recommended for many files)**
```bash
cd output_viz
python -m http.server
```
Then open [http://localhost:8000](http://localhost:8000) in your browser.

**Option B: Open directly**
Double-click any `.html` file to open it in your browser.

### 3D Navigation
- **Left-click + drag** — Rotate
- **Right-click + drag** — Zoom
- **Middle-click + drag** — Pan
- **Double-click** — Reset view
