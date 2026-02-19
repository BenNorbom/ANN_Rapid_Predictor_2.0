# Visualization Scripts

This directory contains two visualization scripts for tractography activation results.
Both produce interactive 3D HTML plots rendered with Plotly.

> **Windows / PowerShell Note:** PowerShell does not support `\` for line continuation. Replace `\` with a backtick `` ` `` at the end of each line, or run the command on a single line.

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
| `--electrode FILE` | none | Path to the electrode FEM export file (`.txt`, COMSOL format). When provided, an electric-field isosurface is rendered. Toggle it on/off in the plot legend. |
| `--electrode_center X Y Z` | `0 0 0` | X Y Z voxel position of the electrode center. Used to re-center the FEM grid so the field aligns with the fiber tracts. |
| `--electrode_config CONFIG` | none | Contact configuration string (e.g. `01-23`, `+012-3`, `-0+1-2+3`). `-` = cathode (red), `+` = anode (blue), unmarked = inactive (grey). |
| `--field_subsample N` | `4` | Keep every Nth point per axis when loading the electrode grid. Increase this value to reduce memory usage and rendering time for large FEM files. |
| `--downsample N` | `1` | Keep every Nth point per fiber to speed up rendering. |
| `--show_axes` | off | Show toggleable X/Y/Z axes in the 3D scene. |
| `--all_fibers` | off | Plot ALL fibers at full length, ignoring the `valid_inds` filter. Fibers without a threshold are drawn as inactive. |

### Electric Field Visualization

When `--electrode` is provided, the script loads the COMSOL FEM export file
(columns: `x y z V`), sub-samples it by `--field_subsample`, and renders a
**Plotly Volume** trace with a symmetric-log (symlog) voltage colour map
(`RdBu` scale: blue = positive/anode, red = negative/cathode).
The symlog transform compresses the large dynamic range of the field so both
strong and weak regions are visible simultaneously; values near zero are
rendered transparent. Hovering over the volume shows the actual voltage in V.
Toggle **Electric Field** on/off from the plot legend.

### Toggleable Legend Groups

Each plot has three legend groups you can click to show/hide:
- **Activated** — fibers below the activation threshold (shown by default)
- **Inactive** — fibers above the threshold (hidden by default; click to show)
- **Electric Field** — the isosurface (shown by default when `--electrode` is given)

### Pulse Width Auto-Detection

The script automatically reads whatever pulse widths are present in the results
JSON file. There is no need to specify them — if the simulation ran one pulse
width, one plot is generated; if it ran all 17, you get 17 plots.

### Examples

```bash
# Basic usage (no electric field)
python graphing/plot_tracts_fast.py \
    --tract example_tracks/L_DRTT_voxel.txt \
    --results run/results.json \
    --output output_viz/

# With electric field overlay
python graphing/plot_tracts_fast.py \
    --tract example_tracks/L_DRTT_voxel.txt \
    --results run/results.json \
    --output output_viz/ \
    --electrode "electrodes/directed/monopolar/electrode_file.txt" \
    --electrode_center 167 223 147 \
    --field_subsample 4

# Show all fibers + axes
python graphing/plot_tracts_fast.py \
    --tract example_tracks/L_DRTT_voxel.txt \
    --results run/results.json \
    --output output_viz/ \
    --all_fibers --show_axes
```

---

## plot_tracts.py

Original visualization script with per-fiber rendering. Better suited for smaller
tract files or when you need single pulse-width interactive plots.

### Usage

```bash
python graphing/plot_tracts.py --tract <TRACT_FILE> --results <RESULTS_JSON> --output <OUTPUT_DIR> [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--tract FILE` | Path to the fiber tract file. |
| `--results FILE` | Path to the results JSON file. |
| `--output DIR` | Directory to save output files. |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--voltage VOLTS` | `3.0` | Activation threshold voltage. |
| `--cond TYPE` | `anisotropic` | Conductivity type: `anisotropic` or `isotropic`. |
| `--electrode FILE` | none | Path to the electrode FEM export file (`.txt`, COMSOL format). Renders an electric-field isosurface in the 3D scene. |
| `--electrode_center X Y Z` | none | X Y Z voxel position of the electrode center for field re-centering. |
| `--electrode_config CONFIG` | none | Electrode contact configuration string (same format as above). |
| `--field_subsample N` | `4` | Subsample the electrode grid by keeping every Nth point per axis. |
| `--show_axes` | off | Show toggleable XYZ axes. |
| `--filter_indices FILE` | none | Path to a JSON or text file listing the indices of fibers that were simulated. |
| `--all_fibers` | off | Plot ALL fibers, ignoring the `--filter_indices` file. |
| `--interactive_pw INDEX` | none | Render a single pulse width by index interactively. |

### Examples

```bash
# Basic usage
python graphing/plot_tracts.py \
    --tract example_tracks/L_DRTT_voxel.txt \
    --results run/results.json \
    --output output_viz/

# With electric field overlay
python graphing/plot_tracts.py \
    --tract example_tracks/L_DRTT_voxel.txt \
    --results run/results.json \
    --output output_viz/ \
    --electrode "electrodes/directed/monopolar/electrode_file.txt" \
    --electrode_center 167 223 147 \
    --field_subsample 4

# Render a single pulse width interactively
python graphing/plot_tracts.py \
    --tract example_tracks/L_DRTT_voxel.txt \
    --results run/results.json \
    --output output_viz/ \
    --interactive_pw 0
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
