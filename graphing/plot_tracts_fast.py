"""
Fast 3D visualization for large tractography datasets.

Reads ANN prediction results and the corresponding tract file, then generates
interactive Plotly HTML plots showing activated vs. inactive fibers with an
optional electric-field isosurface from the electrode FEM data.

Key design: Only TWO Plotly traces are created for fibers (one red "Activated",
one grey "Inactive"), with individual fibers merged into single coordinate arrays
separated by NaN.  This lets the browser render 100k+ fibers in real-time.
"""

import os
import sys
import json
import argparse
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import time


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def mkdirp(path):
    os.makedirs(path, exist_ok=True)


def read_tract_file(path, downsample=1):
    """Read a tract text file.  Each line = one fiber: x1 y1 z1 x2 y2 z2 ...
    Returns a list of numpy arrays, each with shape (N, 3)."""
    print(f"Reading tract file: {path} ...")
    t0 = time.time()
    fibers = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            nums = np.fromstring(line, sep=" ")
            if nums.size % 3 != 0:
                continue
            coords = nums.reshape(-1, 3)
            if downsample > 1:
                coords = coords[::downsample]
            fibers.append(coords)
    print(f"Loaded {len(fibers)} fibers in {time.time() - t0:.2f}s")
    return fibers


def load_results(path):
    """Load a results JSON written by dti_ann_LUT.py.

    Returns
    -------
    data : dict          - the raw JSON
    pulse_widths : list   - sorted list of pulse-width keys found (as floats, in seconds)
    valid_inds : list     - list of valid fiber indices (may be empty)
    """
    with open(path) as f:
        data = json.load(f)

    valid_inds = data.get("valid_inds", [])

    # Discover pulse widths present in the file.  They are stored as string
    # keys like "0.06" (seconds).  Skip non-numeric keys.
    pw_keys = []
    for k in data:
        try:
            pw_keys.append(float(k))
        except ValueError:
            pass
    pw_keys.sort()

    return data, pw_keys, [int(i) for i in valid_inds]


def get_thresholds_for_pw(data, pw_key):
    """Return a list of threshold values for a given pulse-width key."""
    pw_data = data.get(str(pw_key), data.get(pw_key, {}))
    if isinstance(pw_data, list):
        return pw_data
    if isinstance(pw_data, dict):
        indices = sorted(int(k) for k in pw_data)
        return [pw_data[str(i)] for i in indices]
    return []


# ---------------------------------------------------------------------------
# Electric field loading (subsampled for visualisation)
# ---------------------------------------------------------------------------

def load_electrode_field(path, subsample=4, electrode_center=None):
    """Load a COMSOL electrode export and return subsampled 3D field data.

    Parameters
    ----------
    path : str
        Path to the electrode .txt file (COMSOL export: x y z V).
    subsample : int
        Keep every Nth unique coordinate along each axis to reduce data volume.
    electrode_center : tuple or None
        (cx, cy, cz) — if provided, the FEM grid is re-centered so its
        midpoint sits at these coordinates (same logic as dti_ann_LUT.py).

    Returns
    -------
    X, Y, Z : 3-D ndarrays (meshgrid)
    V : 3-D ndarray of electric potential (mV, absolute value, log-scaled)
    """
    print(f"Loading electrode field: {path} ...")
    t0 = time.time()

    x_coords, y_coords, z_coords, potentials = [], [], [], []
    x_prev = y_prev = z_prev = None

    with open(path) as f:
        for line in f:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            xv = round(float(parts[0]), 3)
            yv = round(float(parts[1]), 3)
            zv = round(float(parts[2]), 3)
            vv = float(parts[3])
            potentials.append(vv)
            if xv != x_prev and xv not in x_coords:
                x_coords.append(xv)
            if yv != y_prev and yv not in y_coords:
                y_coords.append(yv)
            if zv != z_prev and zv not in z_coords:
                z_coords.append(zv)
            x_prev, y_prev, z_prev = xv, yv, zv

    nx, ny, nz = len(x_coords), len(y_coords), len(z_coords)
    print(f"  Grid: {nx} x {ny} x {nz} = {nx*ny*nz:,} points")

    # Reshape — data is stored z-outer, y-middle, x-inner (like FEM.py)
    V_full = np.array(potentials, dtype=np.float64).reshape((nz, ny, nx))
    V_full = np.transpose(V_full, (2, 1, 0))  # -> (nx, ny, nz)

    x_arr = np.array(x_coords)
    y_arr = np.array(y_coords)
    z_arr = np.array(z_coords)

    # Subsample
    xs = x_arr[::subsample]
    ys = y_arr[::subsample]
    zs = z_arr[::subsample]
    Vs = V_full[::subsample, ::subsample, ::subsample]

    # Apply centering offset (same as dti_ann_LUT.py)
    if electrode_center is not None:
        orig_cx = (x_arr[0] + x_arr[-1]) / 2.0
        orig_cy = (y_arr[0] + y_arr[-1]) / 2.0
        orig_cz = (z_arr[0] + z_arr[-1]) / 2.0
        xs = xs + (electrode_center[0] - orig_cx)
        ys = ys + (electrode_center[1] - orig_cy)
        zs = zs + (electrode_center[2] - orig_cz)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

    print(f"  Subsampled to {len(xs)} x {len(ys)} x {len(zs)} = {Vs.size:,} points "
          f"in {time.time() - t0:.1f}s")

    return X, Y, Z, Vs


# ---------------------------------------------------------------------------
# Batched trace building (the performance trick)
# ---------------------------------------------------------------------------

def build_merged_trace(fibers, indices):
    """Merge many fibers into one (x, y, z) tuple with NaN separators."""
    if not indices:
        return None, None, None
    selected = [fibers[i] for i in indices]
    total = sum(f.shape[0] for f in selected) + len(selected)
    merged = np.full((total, 3), np.nan, dtype=np.float32)
    cur = 0
    for fib in selected:
        n = fib.shape[0]
        merged[cur : cur + n, :] = fib
        cur += n + 1
    return merged[:, 0], merged[:, 1], merged[:, 2]


# ---------------------------------------------------------------------------
# Electrode config parsing
# ---------------------------------------------------------------------------

def parse_electrode_config(config_str):
    """Parse an electrode configuration string like '01-23', '+012-3', '-0+1-2+3'.

    Returns a dict mapping contact number (0-3) to polarity:
        '+' = anode (blue), '-' = cathode (red), None = inactive (grey).

    Rules:
        - A '+' or '-' sign applies to the NEXT digit only.
        - Digits without a preceding sign are inactive.
    """
    contact_polarity = {0: None, 1: None, 2: None, 3: None}
    current_sign = None
    for ch in config_str:
        if ch in '+-':
            current_sign = ch
        elif ch.isdigit():
            contact_polarity[int(ch)] = current_sign
            current_sign = None  # sign consumed — reset
    return contact_polarity


def contact_color(polarity):
    """Return a display color for a contact given its polarity."""
    if polarity == '-':
        return '#CC0000'   # red  (cathode)
    elif polarity == '+':
        return '#0066CC'   # blue (anode)
    return '#404040'       # grey (inactive)


def contact_label(index, polarity):
    if polarity == '-':
        return f'Contact {index} (cathode -)'
    elif polarity == '+':
        return f'Contact {index} (anode +)'
    return f'Contact {index} (inactive)'


# ---------------------------------------------------------------------------
# Electric field rendering
# ---------------------------------------------------------------------------

def add_electric_field(fig, X, Y, Z, V):
    """Add continuous volumetric electric-field rendering to the figure.

    A symmetric-log (symlog) transform compresses the huge dynamic range
    while preserving sign.  The result is mapped to a diverging colour
    scale (red = negative / cathode, blue = positive / anode) with
    opacity that fades smoothly to fully transparent near zero.
    """
    Vc = V.copy()
    Vc[np.isnan(Vc)] = 0.0

    # --- symmetric log transform ---
    # linthresh sets the "floor": values below this become transparent.
    # Above linthresh the scale is logarithmic to compress the huge range.
    linthresh = 1e-6  # V
    sign = np.sign(Vc)
    mag  = np.abs(Vc)
    mag_safe = np.maximum(mag, linthresh)  # avoid log10(0)
    slog = np.where(
        mag <= linthresh,
        Vc / linthresh,
        sign * (1.0 + np.log10(mag_safe / linthresh)),
    )

    # Normalise to [-1, 1]
    absmax = np.max(np.abs(slog))
    if absmax > 0:
        slog /= absmax

    # Compute actual voltage extremes for colorbar tick labels
    vmin = float(np.nanmin(Vc))
    vmax = float(np.nanmax(Vc))

    # Fraction of the normalised range that corresponds to the transparent
    # "dead zone" around zero.  Values whose |slog| is below this fraction
    # will be invisible; above it, opacity ramps up.
    dead = 0.15  # ~first-decade of the log range -> transparent

    fig.add_trace(go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=slog.flatten(),
        customdata=Vc.flatten(),
        hovertemplate=(
            "x: %{x:.1f}<br>"
            "y: %{y:.1f}<br>"
            "z: %{z:.1f}<br>"
            "V: %{customdata:.4g} V"
            "<extra>Electric Field</extra>"
        ),
        isomin=-1.0,
        isomax=1.0,
        opacity=0.5,
        surface_count=50,
        colorscale='RdBu',
        opacityscale=[
            [0.0,          1.0],   # strong negative -> opaque
            [0.5 - dead,   0.08],
            [0.5,          0.0],   # zero -> fully transparent
            [0.5 + dead,   0.08],
            [1.0,          1.0],   # strong positive -> opaque
        ],
        caps=dict(x_show=False, y_show=False, z_show=False),
        showscale=True,
        colorbar=dict(
            title="V", x=1.02, len=0.6,
            tickvals=[-1, 0, 1],
            ticktext=[f"{vmin:.3g}", "0", f"{vmax:.3g}"],
        ),
        name="Electric Field",
        visible=True,
        legendgroup="efield",
        showlegend=True,
    ))


# ---------------------------------------------------------------------------
# Scene rendering
# ---------------------------------------------------------------------------

def render(fibers, thresholds, voltage, electrode_center, title="",
           show_axes=False, electrode_config=None, field_data=None):
    """Build a Plotly figure with activated/inactive fibers and optional E-field."""
    fig = go.Figure()

    active, inactive = [], []
    for i in range(len(fibers)):
        thr = float("inf")
        if i < len(thresholds):
            try:
                thr = float(thresholds[i])
            except (TypeError, ValueError):
                pass
        (active if thr < voltage else inactive).append(i)

    print(f"  Activated: {len(active)},  Inactive: {len(inactive)}")

    ax, ay, az = build_merged_trace(fibers, active)
    if ax is not None:
        fig.add_trace(go.Scatter3d(
            x=ax, y=ay, z=az, mode="lines",
            line=dict(color="red", width=2), opacity=1.0,
            name="Activated", connectgaps=False,
            legendgroup="activated", showlegend=True,
        ))

    ix, iy, iz = build_merged_trace(fibers, inactive)
    if ix is not None:
        fig.add_trace(go.Scatter3d(
            x=ix, y=iy, z=iz, mode="lines",
            line=dict(color="black", width=1), opacity=0.15,
            name="Inactive", connectgaps=False,
            legendgroup="inactive", showlegend=True,
            visible="legendonly",  # inactive fibers hidden by default
        ))

    # Electric field (replaces lead mesh)
    if field_data is not None:
        X, Y, Z, V = field_data
        add_electric_field(fig, X, Y, Z, V)

    # Axes toggle visibility
    show = show_axes
    axis_style = dict(
        visible=True,
        showticklabels=show,
        showgrid=show,
        zeroline=show,
        title="" if not show else None,
    )
    hide = dict(visible=False, showticklabels=False, showgrid=False, zeroline=False)

    scene = dict(
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)),
        aspectmode="data",
        xaxis=axis_style if show else hide,
        yaxis=axis_style if show else hide,
        zaxis=axis_style if show else hide,
    )

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center"),
        scene=scene,
        showlegend=True,
        legend=dict(
            itemsizing="constant",
            title="Click to toggle",
        ),
        margin=dict(l=0, r=0, t=30, b=0),
    )
    return fig, len(active)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Fast 3D fiber-activation visualizer for ANN Rapid Predictor results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Plot all pulse widths found in the results file
  python plot_tracts_fast.py --tract tracts.txt --results results.json --output plots/

  # Only show fibers activated below 2.5 V
  python plot_tracts_fast.py --tract tracts.txt --results results.json --output plots/ --activation_threshold 2.5

  # Downsample points per fiber for faster rendering of huge files
  python plot_tracts_fast.py --tract tracts.txt --results results.json --output plots/ --downsample 3
""",
    )

    parser.add_argument(
        "--tract", required=True, metavar="FILE",
        help="Path to the tract text file (one fiber per line, x1 y1 z1 x2 y2 z2 ...).",
    )
    parser.add_argument(
        "--results", required=True, metavar="FILE",
        help="Path to the JSON results file produced by dti_ann_LUT.py.",
    )
    parser.add_argument(
        "--output", required=True, metavar="DIR",
        help="Directory where HTML plots will be saved.",
    )
    parser.add_argument(
        "--activation_threshold", type=float, default=3.0, metavar="VOLTS",
        help="Voltage (V) below which a fiber is considered activated.  Default: 3.0",
    )
    parser.add_argument(
        "--electrode_center", type=float, nargs=3, default=None, metavar=("X", "Y", "Z"),
        help="X Y Z position of the electrode center for display.  "
             "If omitted, defaults to (0, 0, 0).",
    )
    parser.add_argument(
        "--downsample", type=int, default=1, metavar="N",
        help="Keep every Nth point per fiber to speed up rendering.  Default: 1 (no downsampling).",
    )
    parser.add_argument(
        "--show_axes", action="store_true",
        help="Show X/Y/Z axes in the 3D scene.",
    )
    parser.add_argument(
        "--all_fibers", action="store_true",
        help="Plot ALL fibers at full length, ignoring the FEM valid_inds filter. "
             "Fibers without a threshold are drawn as inactive.",
    )
    parser.add_argument(
        "--electrode_config", type=str, default=None, metavar="CONFIG",
        help="Electrode contact configuration string, e.g. '01-23', '+012-3', '-0+1-2+3'. "
             "'-' marks cathodes (red), '+' marks anodes (blue), unmarked contacts are grey.",
    )
    parser.add_argument(
        "--electrode", type=str, default=None, metavar="FILE",
        help="Path to the electrode FEM export file (.txt, COMSOL format). "
             "When provided, an electric-field isosurface is shown instead of a "
             "lead mesh.  Toggle it on/off in the legend.",
    )
    parser.add_argument(
        "--field_subsample", type=int, default=4, metavar="N",
        help="Subsample the electrode grid by keeping every Nth point per axis. "
             "Default: 4.  Increase to reduce memory / rendering time.",
    )

    args = parser.parse_args()
    mkdirp(args.output)

    electrode_center = tuple(args.electrode_center) if args.electrode_center else (0, 0, 0)
    electrode_config = parse_electrode_config(args.electrode_config) if args.electrode_config else None

    # Load electric field if an electrode file was provided
    field_data = None
    if args.electrode:
        field_data = load_electrode_field(
            args.electrode,
            subsample=args.field_subsample,
            electrode_center=electrode_center,
        )

    # 1. Load results JSON - auto-detect pulse widths
    print(f"Loading results: {args.results}")
    data, pw_keys, valid_inds = load_results(args.results)

    if not pw_keys:
        print("ERROR: No pulse-width data found in the results file.")
        sys.exit(1)

    print(f"Found {len(pw_keys)} pulse width(s): " +
          ", ".join(f"{pw*1000:.0f} us" for pw in pw_keys))

    # 2. Read fibers
    fibers = read_tract_file(args.tract, downsample=args.downsample)

    # 3. Filter to valid indices (if present in results) or show all
    if args.all_fibers:
        print(f"--all_fibers: showing all {len(fibers)} fibers (no FEM filter).")
    elif valid_inds:
        print(f"Filtering to {len(valid_inds)} valid fibers (from results file)...")
        try:
            fibers = [fibers[i] for i in valid_inds]
        except IndexError:
            print("ERROR: valid_inds references fiber indices that don't exist in the tract file.")
            sys.exit(1)
        print(f"  {len(fibers)} fibers after filtering.")

    # 4. Generate one plot per pulse width
    summary = {}
    for idx, pw in enumerate(pw_keys):
        pw_us = pw * 1000  # convert seconds -> microseconds for display
        title = f"Pulse Width: {pw_us:.0f} us  |  Threshold: {args.activation_threshold} V"
        print(f"\nPlotting PW {pw_us:.0f} us  ({idx + 1}/{len(pw_keys)}) ...")

        thresholds = get_thresholds_for_pw(data, pw)
        fig, n_active = render(
            fibers, thresholds, args.activation_threshold,
            electrode_center, title, args.show_axes,
            electrode_config=electrode_config,
            field_data=field_data,
        )

        out_html = os.path.join(args.output, f"activation_pw_{idx:02d}_{pw_us:.0f}us.html")
        fig.write_html(out_html)
        print(f"  Saved {out_html}")

        summary[f"{pw_us:.0f}"] = n_active

    # 5. Write summary
    summary_path = os.path.join(args.output, "activation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(
            {"activation_threshold_V": args.activation_threshold, "activated_counts": summary},
            f, indent=2,
        )
    print(f"\nDone. Summary -> {summary_path}")


if __name__ == "__main__":
    main()
