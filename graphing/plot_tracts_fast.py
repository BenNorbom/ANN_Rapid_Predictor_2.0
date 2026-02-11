"""
Fast 3D visualization for large tractography datasets.

Reads ANN prediction results and the corresponding tract file, then generates
interactive Plotly HTML plots showing activated vs. inactive fibers.

Key design: Only TWO Plotly traces are created (one red "Activated", one grey
"Inactive"), with individual fibers merged into single coordinate arrays
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
# Electrode rendering (detailed 3D mesh)
# ---------------------------------------------------------------------------

def plot_electrode(fig, leftLeadPos, config=None, max_fiber_z=None):
    """Add a detailed 3D mesh electrode visualisation to the Plotly figure.

    Parameters
    ----------
    leftLeadPos : list of lists
        [[x0, x1], [y0, y1], [z0, z1]] defining the lead axis.
    config : dict or None
        Mapping {0: '+'/'-'/None, ...} from parse_electrode_config().
        If None, all contacts drawn as inactive grey.
    max_fiber_z : float or None
        If provided, cap the lead shaft so it doesn't extend above this z.
    """
    if config is None:
        config = {0: None, 1: None, 2: None, 3: None}

    lead_color = '#A6A6A6'

    scale = 3.0
    radius = (1.27 / 2) * scale
    step = np.pi / 32
    contact_height = 1.5 * scale

    x_temp = leftLeadPos[0]
    y_temp = leftLeadPos[1]
    z_temp = leftLeadPos[2]

    # Build rotation matrix from lead axis
    x_trans = np.subtract(x_temp, x_temp[0])
    y_trans = np.subtract(y_temp, y_temp[0])
    z_trans = np.subtract(z_temp, z_temp[0])

    vectorz = [x_trans[1], y_trans[1], z_trans[1]]
    uvz = vectorz / np.linalg.norm(vectorz)
    uvx = np.cross(uvz, [0, 1, 0])
    if np.linalg.norm(uvx) == 0:
        uvx = np.array([1.0, 0.0, 0.0])
    else:
        uvx = uvx / np.linalg.norm(uvx)
    uvy = np.cross(uvz, uvx)

    rotation_matrix = np.array([
        [uvx[0], uvx[1], uvx[2], 0],
        [uvy[0], uvy[1], uvy[2], 0],
        [uvz[0], uvz[1], uvz[2], 0],
        [0.0, 0.0, 0.0, 1.0],
    ])

    def create_cylinder_mesh(height, position, color, is_tip=False):
        """Create a cylinder mesh segment and add as a Mesh3d trace."""
        if is_tip:
            phi, theta = np.meshgrid(
                np.arange(np.pi / 2, np.pi + step, step),
                np.arange(0, 2 * np.pi + step, step),
            )
        else:
            height_step = np.arctan2(height, radius)
            phi, theta = np.meshgrid(
                np.arange(np.pi / 2, np.pi / 2 + height_step + height_step, height_step),
                np.arange(0, 2 * np.pi + step, step),
            )

        x = np.sin(phi) * np.cos(theta) * radius
        y = np.sin(phi) * np.sin(theta) * radius
        z = np.cos(phi) * radius

        if not is_tip:
            mult = 1.0 / np.cos(height_step)
            x[:, 1] *= mult
            y[:, 1] *= mult
            z[:, 1] *= mult

        z = np.add(z, height + position)

        # Rotate and translate into world coordinates
        x_rot = np.empty_like(x)
        y_rot = np.empty_like(y)
        z_rot = np.empty_like(z)
        for i, j in np.ndindex(x.shape):
            pt = np.array([x[i, j], y[i, j], z[i, j], 1.0])
            rp = rotation_matrix.dot(pt)
            x_rot[i, j] = rp[0] + x_temp[0]
            y_rot[i, j] = rp[1] + y_temp[0]
            z_rot[i, j] = rp[2] + z_temp[0]

        # Build triangle faces
        verts = np.column_stack((x_rot.flatten(), y_rot.flatten(), z_rot.flatten()))
        tri_i, tri_j, tri_k = [], [], []
        nrows, ncols = x_rot.shape
        for r in range(nrows - 1):
            for c in range(ncols - 1):
                base = r * ncols + c
                tri_i.extend([base, base + 1])
                tri_j.extend([base + 1, base + ncols + 1])
                tri_k.extend([base + ncols, base + ncols])

        fig.add_trace(go.Mesh3d(
            x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
            i=tri_i, j=tri_j, k=tri_k,
            color=color, opacity=1.0,
            name='Electrode', showlegend=False,
        ))

    # Rounded tip
    tip_height = contact_height - radius
    create_cylinder_mesh(tip_height, -1 * tip_height, lead_color, is_tip=True)
    create_cylinder_mesh(tip_height, -1 * tip_height, lead_color)

    # 4 contacts + 3 gaps = 7 segments (tip → shaft)
    # Contact 0 is at the tip (distal), Contact 3 is nearest the shaft (proximal)
    positions = [contact_height * i for i in range(7)]
    colors = [
        contact_color(config.get(0)),  # Contact 0 (tip / distal)
        lead_color,                    # Gap
        contact_color(config.get(1)),  # Contact 1
        lead_color,                    # Gap
        contact_color(config.get(2)),  # Contact 2
        lead_color,                    # Gap
        contact_color(config.get(3)),  # Contact 3 (shaft / proximal)
    ]
    for pos, col in zip(positions, colors):
        create_cylinder_mesh(contact_height, pos, col)

    # Shaft beyond the contacts
    shaft_pos = positions[-1] + contact_height
    default_shaft = 100.0 * scale
    shaft_length = default_shaft
    if max_fiber_z is not None:
        uvz_vec = np.array([rotation_matrix[2, 0], rotation_matrix[2, 1], rotation_matrix[2, 2]])
        uvz_z = uvz_vec[2]
        if uvz_z > 1e-6:
            max_allowed = (max_fiber_z - z_temp[0]) / uvz_z - shaft_pos
            shaft_length = max(0.0, min(default_shaft, max_allowed - 0.1 * scale))
    create_cylinder_mesh(max(0.0, shaft_length), shaft_pos, lead_color)


# ---------------------------------------------------------------------------
# Scene rendering
# ---------------------------------------------------------------------------

def render(fibers, thresholds, voltage, electrode_center, title="", show_axes=False, electrode_config=None):
    """Build a Plotly figure with activated (red) and inactive (grey) fibers."""
    fig = make_subplots(rows=1, cols=1, specs=[[{"type": "scene"}]])

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
        ))

    ix, iy, iz = build_merged_trace(fibers, inactive)
    if ix is not None:
        fig.add_trace(go.Scatter3d(
            x=ix, y=iy, z=iz, mode="lines",
            line=dict(color="black", width=1), opacity=0.15,
            name="Inactive", connectgaps=False,
        ))

    # Compute max fiber Z so we can cap the lead shaft height
    max_fiber_z = None
    if fibers:
        try:
            max_fiber_z = max(f[:, 2].max() for f in fibers)
        except Exception:
            pass

    # Build leftLeadPos from electrode center: two-point axis along Z
    cx, cy, cz = electrode_center
    leftLeadPos = [[cx, cx], [cy, cy], [cz - 5, cz + 5]]
    plot_electrode(fig, leftLeadPos, config=electrode_config, max_fiber_z=max_fiber_z)

    scene = dict(camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)), aspectmode="data")
    if not show_axes:
        hide = dict(visible=False, showticklabels=False, showgrid=False, zeroline=False)
        scene.update(xaxis=hide, yaxis=hide, zaxis=hide)

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center"),
        scene=scene, showlegend=True,
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

    args = parser.parse_args()
    mkdirp(args.output)

    electrode_center = tuple(args.electrode_center) if args.electrode_center else (0, 0, 0)
    electrode_config = parse_electrode_config(args.electrode_config) if args.electrode_config else None

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
