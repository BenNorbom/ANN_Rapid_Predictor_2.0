# Visualization Scripts

This directory contains three visualization scripts for tractography activation results.
All produce interactive 3D HTML plots using [Plotly](https://plotly.com/python/).

> **Windows / PowerShell Note:** PowerShell does not support `\` for line continuation. Replace `\` with a backtick `` ` `` at the end of each line, or run the command on a single line.

---

## Script Overview

| Script | Best For | Fiber Rendering |
|--------|----------|-----------------|
| `plot_tracts_fast.py` | Large datasets (100k+ fibers) | 2 merged WebGL traces |
| `plot_tracts.py` | Small datasets, single-PW interactive | Per-fiber traces |
| `plot_tracts_bundles.py` | Named bundles (DRTT, ML, PTR) with distinct colours | Bundle-aware merged traces |

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
| `--electrode FILE` | none | Path to the electrode FEM export file (`.txt`, COMSOL format). When provided, an electric-field volume is rendered. |
| `--electrode_center X Y Z` | `0 0 0` | X Y Z voxel position of the electrode center. Used to re-center the FEM grid so the field aligns with the fiber tracts. |
| `--electrode_config CONFIG` | none | Contact configuration string (e.g. `01-23`, `+012-3`, `-0+1-2+3`). `-` = cathode (red), `+` = anode (blue), unmarked = inactive (grey). |
| `--field_subsample N` | `4` | Keep every Nth point per axis when loading the electrode grid. Increase to reduce memory and rendering time. |
| `--downsample N` | `1` | Keep every Nth point per fiber when reading the tract file. |
| `--simplify N` | `3` | Keep every Nth point per fiber for rendering. Set to 1 for full resolution. |
| `--show_axes` | off | Show toggleable X/Y/Z axes in the 3D scene. |
| `--all_fibers` | off | Plot ALL fibers, ignoring the `valid_inds` filter. Fibers without a threshold are drawn as inactive. |

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
| `--electrode FILE` | none | Path to the electrode FEM export file. Renders an electric-field volume in the 3D scene. |
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

## plot_tracts_bundles.py

Bundle-aware visualization that highlights named fiber bundles (DRTT, ML, PTR)
in distinct colours. The script reads a **combined tract file** (multiple tract
files concatenated end-to-end) and a **manifest JSON** that records which line
ranges belong to each bundle.

### Preparing the Combined Tract File

Concatenate the individual tract files and create a manifest that records where
each bundle starts and ends. An example combined file and manifest are included
in `example_tracks/`.

#### PowerShell

```powershell
$drtt = "example_tracks/L_DRTT_voxel.txt"
$ml   = "example_tracks/L_ML_voxel.txt"
$ptr  = "example_tracks/L_PTR_voxel.txt"
$out  = "example_tracks/L_combined_DRTT_ML_PTR_voxel.txt"

# Count lines (= number of fibers)
$drttCount = (Get-Content $drtt | Measure-Object -Line).Lines
$mlCount   = (Get-Content $ml   | Measure-Object -Line).Lines
$ptrCount  = (Get-Content $ptr  | Measure-Object -Line).Lines

# Concatenate
Get-Content $drtt, $ml, $ptr | Set-Content $out

# Build manifest JSON
$manifest = @{
    bundles = [ordered]@{
        L_DRTT_voxel = @{ start = 0;           end = $drttCount;                        count = $drttCount }
        L_ML_voxel   = @{ start = $drttCount;   end = $drttCount + $mlCount;             count = $mlCount   }
        L_PTR_voxel  = @{ start = $drttCount + $mlCount; end = $drttCount + $mlCount + $ptrCount; count = $ptrCount }
    }
}
$manifest | ConvertTo-Json -Depth 4 | Set-Content "example_tracks/L_combined_DRTT_ML_PTR_voxel_manifest.json"
```

#### Bash

```bash
drtt="example_tracks/L_DRTT_voxel.txt"
ml="example_tracks/L_ML_voxel.txt"
ptr="example_tracks/L_PTR_voxel.txt"
out="example_tracks/L_combined_DRTT_ML_PTR_voxel.txt"

drtt_n=$(wc -l < "$drtt")
ml_n=$(wc -l < "$ml")
ptr_n=$(wc -l < "$ptr")

cat "$drtt" "$ml" "$ptr" > "$out"

ml_start=$drtt_n
ptr_start=$((ml_start + ml_n))

cat > "example_tracks/L_combined_DRTT_ML_PTR_voxel_manifest.json" <<EOF
{
  "bundles": {
    "L_DRTT_voxel": { "start": 0,          "end": $drtt_n,                      "count": $drtt_n },
    "L_ML_voxel":   { "start": $ml_start,   "end": $((ml_start + ml_n)),         "count": $ml_n   },
    "L_PTR_voxel":  { "start": $ptr_start,  "end": $((ptr_start + ptr_n)),       "count": $ptr_n  }
  }
}
EOF
```

#### Manifest Format

```json
{
  "bundles": {
    "L_DRTT_voxel": { "start": 0,   "end": 100, "count": 100 },
    "L_ML_voxel":   { "start": 100, "end": 200, "count": 100 },
    "L_PTR_voxel":  { "start": 200, "end": 300, "count": 100 }
  }
}
```

- **`start`** — Zero-based index of the first fiber in this bundle.
- **`end`** — One-past-the-last index (Python-style half-open range).
- **`count`** — Number of fibers in the bundle.

> Bundle names must match the keys in `BUNDLE_STYLES` inside the script
> (`L_DRTT_voxel`, `L_ML_voxel`, `L_PTR_voxel`). Unrecognised names fall back
> to an orange colour scheme.

### Usage

```bash
python graphing/plot_tracts_bundles.py \
    --tract <COMBINED_TRACT_FILE> \
    --results <RESULTS_JSON> \
    --manifest <MANIFEST_JSON> \
    --output <OUTPUT_DIR> \
    [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--tract FILE` | Path to the combined tract file. |
| `--results FILE` | Path to the JSON results file from `dti_ann_LUT.py`. |
| `--manifest FILE` | Path to the manifest JSON (see above). |
| `--output DIR` | Directory where HTML plots will be saved. |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--activation_threshold VOLTS` | `3.0` | Voltage (V) below which a fiber is considered **activated**. |
| `--electrode FILE` | none | Path to the electrode FEM export file (`.txt`, COMSOL format). Renders an electric-field volume. |
| `--electrode_center X Y Z` | `0 0 0` | X Y Z voxel position of the electrode center for re-centering the FEM grid. |
| `--electrode_config CONFIG` | none | Contact configuration string (e.g. `01-23`, `+012-3`, `-0+1-2+3`). |
| `--field_subsample N` | `4` | Keep every Nth point per axis when loading the electrode grid. |
| `--downsample N` | `1` | Keep every Nth point per fiber when reading the tract file. |
| `--wb_simplify N` | `3` | Keep every Nth point of non-bundle fibers for rendering. Bundle fibers are always rendered at full resolution. |
| `--show_axes` | off | Show toggleable X/Y/Z axes in the 3D scene. |
| `--all_fibers` | off | Plot ALL fibers, ignoring the `valid_inds` filter from the results file. |

### Bundle Colour Scheme

| Bundle | Activated | Inactive | Opacity (active / inactive) |
|--------|-----------|----------|-----------------------------|
| DRTT (`L_DRTT_voxel`) | Lime `#00FF00` | Dark Green `#006400` | 1.0 / 0.35 |
| ML (`L_ML_voxel`) | Cyan `#00FFFF` | Teal `#008080` | 1.0 / 0.35 |
| PTR (`L_PTR_voxel`) | Magenta `#FF00FF` | Purple `#800080` | 1.0 / 0.35 |

Unrecognised bundle names use an orange fallback (`#FFA500` / `#8B4500`).

### Examples

```bash
# Using the included example combined file
python graphing/plot_tracts_bundles.py \
    --tract example_tracks/L_combined_DRTT_ML_PTR_voxel.txt \
    --results run/results.json \
    --manifest example_tracks/L_combined_DRTT_ML_PTR_voxel_manifest.json \
    --output output_viz_bundles/ \
    --all_fibers --show_axes

# With electric field overlay
python graphing/plot_tracts_bundles.py \
    --tract example_tracks/L_combined_DRTT_ML_PTR_voxel.txt \
    --results run/results.json \
    --manifest example_tracks/L_combined_DRTT_ML_PTR_voxel_manifest.json \
    --output output_viz_bundles/ \
    --electrode "electrodes/directed/monopolar/electrode_file.txt" \
    --electrode_center 167 213 147 \
    --electrode_config 01-23 \
    --all_fibers --show_axes
```

---

## Electric Field Visualization — Technical Details

All three scripts share the same electric field rendering pipeline. This section
describes the implementation in detail.

### Data Loading (`load_electrode_field`)

The electrode FEM export file is a COMSOL text file with columns `x y z V`.
The loader:

1. Parses unique x, y, z coordinates to determine the 3-D grid dimensions.
2. Reshapes the flat voltage array into a `(nx, ny, nz)` volume. The raw data
   is stored in z-outer / y-middle / x-inner order (matching COMSOL's default
   export), so it is transposed to `(x, y, z)` axis order via `np.transpose`.
3. **Sub-samples** the grid by keeping every Nth coordinate along each axis
   (`--field_subsample`). This reduces the point count by N³, which is critical
   for interactive rendering of grids that may have millions of points.
4. **Re-centers** the FEM grid if `--electrode_center` is provided. The offset
   between the original FEM midpoint and the requested center is applied to all
   axis coordinates, so the field aligns with the fiber data in MNI/voxel space.

### Symmetric-Log (symlog) Transform

Electrode fields span many orders of magnitude — voltages near the contact
surface can be > 1 V while the field boundary is < 1 µV. A linear colour map
would make everything outside the immediate contact region invisible.

The transform works in three stages:

1. **Threshold**: A linear threshold `linthresh = 1e-6 V` defines the boundary
   between linear and log regions.
2. **Log-compress**: For `|V| > linthresh`, the value is mapped to
   `sign(V) * (1 + log10(|V| / linthresh))`. Values below the threshold are
   mapped linearly: `V / linthresh`. This preserves the sign (cathode vs anode)
   while compressing the dynamic range.
3. **Normalise**: The result is divided by its absolute maximum to produce
   values in `[-1, 1]`.

### Adaptive Colour Scale

The colour scale and opacity mapping adapt to the polarity of the field:

| Field Type | Detected When | Colour Scale | Opacity Behaviour |
|------------|---------------|--------------|-------------------|
| **Negative only** (monopolar cathode) | `vmax ≤ 0` | Red gradient (dark → light) | Opaque at most negative, transparent at 0 |
| **Positive only** (monopolar anode) | `vmin ≥ 0` | Blue gradient (light → dark) | Transparent at 0, opaque at most positive |
| **Mixed** (bipolar / multipolar) | both signs present | Diverging `RdBu` | Symmetric: opaque at extremes, transparent at 0 |

This means a purely negative (cathode) field will **not** show any blue, and
vice versa.

### Opacity Mapping

A wide transparent "dead zone" (`dead = 0.30` of the normalised range) around
zero ensures that near-zero background values are fully invisible. The opacity
ramps through control points (shown for the negative-only case):

| Normalised value | Opacity | Description |
|------------------|---------|-------------|
| 0.0 (strongest)  | 1.0     | Fully opaque at the most extreme voltage |
| 0.45             | 0.20    | Mid-range: semi-transparent |
| 0.70 (dead edge) | 0.02    | Nearly invisible |
| 1.0 (zero)       | 0.0     | Fully transparent |

For the diverging (mixed polarity) case the ramp is symmetric around 0.5.

### Plotly Volume Trace

The field is rendered with `plotly.graph_objects.Volume`, which performs
marching-cubes isosurfacing on the client side in WebGL. Key parameters:

- **`surface_count=50`** — 50 concentric isosurfaces are drawn, giving a
  smooth volumetric appearance.
- **`caps`** — Disabled on all three axes to avoid opaque bounding-box faces.
- **`customdata`** — The raw (untransformed) voltage values are attached so
  the hover tooltip shows the actual voltage in V, not the symlog value.
- **`colorbar`** — Positioned at `x=0.02` (left side) to avoid overlapping
  with the fiber legend on the right.

### Toggle Button

A Plotly `updatemenus` button labelled **Toggle Electric Field** is added to
the top-left of each figure. Clicking it toggles the visibility of the Volume
trace using `restyle` with `args` / `args2`, providing a true on/off toggle
independent of the legend.

### Pulse Width Auto-Detection

All scripts automatically read whatever pulse widths are present in the results
JSON file. There is no need to specify them — if the simulation ran one pulse
width, one plot is generated; if it ran all 17, you get 17 plots.

---

## Interactive Legend Groups

Each plot has clickable legend groups you can show/hide:
- **Activated** — fibers below the activation threshold (shown by default)
- **Inactive** — fibers above the threshold (hidden by default; click to show)
- **Electric Field** — the volume rendering (shown by default when `--electrode` is given; also togglable via the button)

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
