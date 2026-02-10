# Neural Tract Visualization Setup Guide

This guide will help you set up and run the 3D neural tract visualization tool.

## Prerequisites

- Python 3.10 or higher
- pip (Python package installer)

## Step 1: Create and Activate Virtual Environment

```powershell
# Create a new virtual environment
python -m venv venvs/viz

# Activate the virtual environment
# On Windows:
.\venvs\viz\Scripts\Activate.ps1
# On Linux/Mac:
# source venvs/viz/bin/activate
```

## Step 2: Install Required Packages

```powershell
# Install required packages
python -m pip install numpy plotly kaleido
```

## Step 3: Prepare Your Data Files

You need two input files:
1. Tract file (`.xyz` or `.xyz.ini` format) containing fiber coordinates
2. Results JSON file containing threshold data

Example file structure:
```
your_project_folder/
├── processed_tract.xyz.ini    # Fiber tract data
├── results.json               # Threshold results
└── plot_tracts_plotly.py     # Visualization script
```

## Step 4: Run the Visualization Script

The arguments are now named for clarity.

```powershell
# Basic usage
python plot_tracts.py --tract path/to/tract.txt --results path/to/results.json --voltage 3.0 --cond anisotropic --output output_folder

# Handling filtered fibers (recommended):
# Pass the results.json file to --filter_indices if your simulation skipped some fibers
python plot_tracts.py --tract path/to/tract.txt --results path/to/results.json --output output_folder --filter_indices path/to/results.json

# View a specific pulse width interactively only
python plot_tracts.py --tract path/to/tract.txt --results path/to/results.json --output output_folder --interactive_pw 0
```

### Arguments

- `--tract`: **(Required)** Path to the fiber tract file (text format).
- `--results`: **(Required)** Path to the JSON file with simulation results.
- `--output`: **(Required)** Directory to save output images and HTML files.
- `--voltage`: Threshold voltage (V) to determine activation. Default is `3.0`.
- `--cond`: Conductivity type (`anisotropic` or `isotropic`). Default is `anisotropic`.
- `--filter_indices`: Path to a file containing valid fiber indices (often the same as your results file).
- `--show_axes`: Flag to show XYZ axes in the 3D scene.
- `--interactive_pw`: Index of a single pulse width to render (0-based).

## Step 5: View the Visualizations

The script generates `.html` files in your output directory (e.g., `activation_pw_00.html`, `activation_pw_01.html`).

**Option A: Use a local server (Recommended)**
This is the best way to browse all results at once. It essentially creates a mini-website folder for your files.

1.  Open your terminal/command prompt.
2.  Navigate to your output folder:
    ```powershell
    cd output_folder
    ```
3.  Start the simple Python server:
    ```powershell
    python -m http.server
    ```
    *(If that fails, try: `python -m http.server 8000 --bind 127.0.0.1`)*

4.  Open your web browser and go to: [http://localhost:8000](http://localhost:8000)
    You will see a directory listing. Click any file to view correctly sorted results (e.g., `09` appears before `10`).

**Option B: Open directly**
Simply double-click the `.html` files in your file explorer. They will open in your default web browser where you can rotate and zoom the 3D model.

## Interacting with the 3D Visualization

- Left click + drag: Rotate the view
- Right click + drag: Zoom in/out
- Middle click + drag: Pan the view
- Double click: Reset the view

## Output Files

The script generates two types of files for each pulse width:
1. Interactive HTML files (`activation_pw_X.html`)
2. Static PNG images (`activation_pw_X.png`)
3. Summary JSON file (`activation_summary.json`)

## Troubleshooting

1. If packages fail to install:
```powershell
python -m pip install --upgrade pip
python -m pip install numpy plotly kaleido --no-cache-dir
```

2. If the web browser won't open local HTML files directly:
Always use the local server method (Step 5) to view the visualizations.

3. If you see import errors:
Make sure you've activated the virtual environment and installed all required packages.

## Script Options

Optional parameters:
- `--interactive_pw INDEX`: Generate visualization for a single pulse width index
  Example: `python plot_tracts_plotly.py processed_tract.xyz.ini results.json 3.0 anisotropic output_folder --interactive_pw 0`
 - `--show_axes`: Show the XYZ axes in the Plotly 3D scene. By default the exporter hides axes for a cleaner visualization. Example:
   `python plot_tracts_plotly.py processed_tract.xyz.ini results.json 3.0 anisotropic output_folder --show_axes`

## Example File Formats

1. Tract file format (`.xyz` or `.xyz.ini`):
```
x1 y1 z1 x2 y2 z2 x3 y3 z3  # Points for fiber 1
x1 y1 z1 x2 y2 z2           # Points for fiber 2
...
```

2. Results JSON format:
```json
{
  "0.06": {                  # Pulse width in seconds (60μs)
    "0": threshold_value,    # Threshold for fiber 0
    "1": threshold_value,    # Threshold for fiber 1
    ...
  },
  "0.075": {                 # Next pulse width (75μs)
    ...
  }
}
```

open interactive images with this link...
cd output_folder; python -m http.server 8000
http://localhost:8000/index.html

## Plotly exporter — recent changes

The Plotly exporter (`plot_tracts_plotly.py`) has been updated. If you are using it, please note the following changes and options:

- **Lead coordinates (Mayavi-derived)**: the script now uses the following `leftLeadPos` by default to match the Mayavi visualizations:

  ```python
  leftLeadPos = [[167, 161], [223, 222], [143, 159]]
  ```

- **Dependencies**: ensure `plotly` and `kaleido` are installed in your environment:

  ```powershell
  python -m pip install -U plotly kaleido
  ```