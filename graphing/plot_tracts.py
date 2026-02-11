import os
import sys
import json
import math
import argparse
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

pulse_widths = [60, 75, 90, 105, 120, 135, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500]


def parse_electrode_config(config_str):
    """Parse an electrode configuration string like '01-23', '+012-3', '-0+1-2+3'.

    Returns a dict mapping contact number (0-3) to polarity:
        '+' = anode (blue), '-' = cathode (red), None = inactive (grey).
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


def polarity_to_color(polarity, lead_color, inactive_color):
    """Map polarity to display color."""
    if polarity == '-':
        return '#CC0000'   # red  (cathode)
    elif polarity == '+':
        return '#0066CC'   # blue (anode)
    return inactive_color  # grey (inactive)

def mkdirp(dir):
    if not os.path.exists(dir):
        os.makedirs(dir, exist_ok=True)

def read_tract_file(path):
    fibers = []
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            nums = [float(x) for x in parts]
            if len(nums) % 3 != 0:
                raise ValueError(f"Line {i+1} in tract file does not have a 3N number of floats")
            pts = [(nums[j], nums[j+1], nums[j+2]) for j in range(0, len(nums), 3)]
            fibers.append(pts)
    return fibers

def load_valid_indices(path):
    """Load a list of indices from a JSON [0, 1, 3...] or text file (newlines)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Filter indices file not found: {path}")

    # Try JSON first
    try:
        with open(path, 'r') as f:
            data = json.load(f)
            if isinstance(data, list):
                return [int(x) for x in data]
            elif isinstance(data, dict) and 'valid_inds' in data:
                return [int(x) for x in data['valid_inds']]
    except Exception:
        pass

    # Fallback to text lines
    indices = []
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                indices.append(int(float(parts[0])))
    return indices

def load_thresholds(path):
    with open(path) as f:
        data = json.load(f)
    thresholds = []
    for pw in pulse_widths:
        key = str(pw/1000)
        row = []
        if key in data:
            # Handle new format where thresholds are in dict or list form under the PW key
            pw_data = data[key]
            # Assumes keys are string numeric indices, find max key to determine list length or iterate
            # For robustness, we'll try to reconstruct based on available keys
            if isinstance(pw_data, dict):
                # Sort keys numerically to ensure order [0, 1, 2...]
                indices = sorted([int(k) for k in pw_data.keys()])
                for idx in indices:
                    row.append(pw_data[str(idx)])
            elif isinstance(pw_data, list):
                 row = pw_data
        thresholds.append(row)
    return thresholds

def plot_electrode(fig, leftLeadPos, max_fiber_z=None, electrode_config=None):
    """Add a detailed electrode visualization to the plotly figure.
    
    Parameters
    ----------
    electrode_config : dict or None
        Mapping {0: '+'/'-'/None, ...} from parse_electrode_config().
        If None, uses the legacy hardcoded layout (contact 2 active).
    """
    # Colors (match the Mayavi visualization as closely as possible)
    lead_color = '#A6A6A6'  # ~ (0.65, 0.65, 0.65)
    inactive_contact_color = '#404040'  # ~ (0.25, 0.25, 0.25)

    # Basic parameters (same geometry as the Mayavi lead)
    # Use a scale factor to make the lead larger if desired
    scale = 3.0  # increase this to make the lead even larger
    radius = (1.27 / 2) * scale
    # finer angular resolution for higher-definition lead mesh
    step = np.pi / 32
    # contact/cylinder nominal height in the original code is 1.5; scale it
    contact_height = 1.5 * scale

    # Get original DBS lead coordinates and create rotation matrix
    x_temp = leftLeadPos[0]
    y_temp = leftLeadPos[1]
    z_temp = leftLeadPos[2]

    # Translate original lead to origin for rotation basis
    x_trans = np.subtract(x_temp, x_temp[0])
    y_trans = np.subtract(y_temp, y_temp[0])
    z_trans = np.subtract(z_temp, z_temp[0])

    vectorz = [x_trans[1], y_trans[1], z_trans[1]]
    uvz = vectorz / np.linalg.norm(vectorz)
    uvx = np.cross(uvz, [0, 1, 0])
    # In degenerate case uvx may be zero-length; guard against that
    if np.linalg.norm(uvx) == 0:
        uvx = np.array([1.0, 0.0, 0.0])
    else:
        uvx = uvx / np.linalg.norm(uvx)
    uvy = np.cross(uvz, uvx)

    rotation_matrix = np.array([
        [uvx[0], uvx[1], uvx[2], 0],
        [uvy[0], uvy[1], uvy[2], 0],
        [uvz[0], uvz[1], uvz[2], 0],
        [0.0, 0.0, 0.0, 1.0]
    ])

    def create_cylinder_mesh(height, position, color, is_tip=False):
        """Create a rounded/rolled cylinder segment (or tip) and add as a Mesh3d trace."""
        if is_tip:
            phi, theta = np.meshgrid(
                np.arange(np.pi/2, np.pi + step, step),
                np.arange(0, 2 * np.pi + step, step)
            )
        else:
            height_step = np.arctan2(height, radius)
            phi, theta = np.meshgrid(
                np.arange(np.pi/2, np.pi/2 + height_step + height_step, height_step),
                np.arange(0, 2 * np.pi + step, step)
            )

        x = np.sin(phi) * np.cos(theta) * radius
        y = np.sin(phi) * np.sin(theta) * radius
        z = np.cos(phi) * radius

        if not is_tip:
            mult = 1.0 / np.cos(height_step)
            # widen the rim as in the Mayavi code
            x[:, 1] *= mult
            y[:, 1] *= mult
            z[:, 1] *= mult

        # shift along the lead axis
        z = np.add(z, height + position)

        # Rotate and translate vertices into world coordinates
        x_rot = np.empty(shape=x.shape)
        y_rot = np.empty(shape=y.shape)
        z_rot = np.empty(shape=z.shape)

        for i, j in np.ndindex(x.shape):
            point = np.array([x[i, j], y[i, j], z[i, j], 1.0])
            rot_point = rotation_matrix.dot(point)
            x_rot[i, j] = rot_point[0] + x_temp[0]
            y_rot[i, j] = rot_point[1] + y_temp[0]
            z_rot[i, j] = rot_point[2] + z_temp[0]

        # Build triangle faces
        vertices = np.column_stack((x_rot.flatten(), y_rot.flatten(), z_rot.flatten()))
        tri_i = []
        tri_j = []
        tri_k = []
        nrows, ncols = x_rot.shape
        for r in range(nrows - 1):
            for c in range(ncols - 1):
                base = r * ncols + c
                tri_i.extend([base, base + 1])
                tri_j.extend([base + 1, base + ncols + 1])
                tri_k.extend([base + ncols, base + ncols])

        fig.add_trace(go.Mesh3d(
            x=vertices[:, 0],
            y=vertices[:, 1],
            z=vertices[:, 2],
            i=tri_i,
            j=tri_j,
            k=tri_k,
            color=color,
            opacity=1.0,
            name='Electrode',
            showlegend=False
        ))

    # Tip (rounded)
    tip_height = contact_height - radius
    create_cylinder_mesh(tip_height, -1 * tip_height, lead_color, is_tip=True)
    create_cylinder_mesh(tip_height, -1 * tip_height, lead_color)

    # Contacts and gaps (positions and colors mirror the Mayavi layout)
    positions = [contact_height * i for i in range(7)]
    
    if electrode_config is not None:
        # Build colors from electrode_config
        # positions 0,2,4,6 are contacts; tip→shaft = contact 0,1,2,3
        colors = []
        contact_order = [0, 1, 2, 3]  # tip (distal) to shaft (proximal)
        contact_slot = 0
        for slot in range(7):
            if slot % 2 == 0:
                idx = contact_order[contact_slot]
                colors.append(polarity_to_color(electrode_config.get(idx), lead_color, inactive_contact_color))
                contact_slot += 1
            else:
                colors.append(lead_color)  # gap
    else:
        # Legacy hardcoded: contact 2 (index 4 in positions) is active
        colors = [
            inactive_contact_color,  # Contact 0 (tip / distal)
            lead_color,              # Gap
            inactive_contact_color,  # Contact 1
            lead_color,              # Gap
            '#990000',               # Contact 2 (active, legacy)
            lead_color,              # Gap
            inactive_contact_color,  # Contact 3 (shaft / proximal)
        ]

    for pos, col in zip(positions, colors):
        create_cylinder_mesh(contact_height, pos, col)

    # Shaft (long rounded cylinder beyond the contacts)
    shaft_pos = positions[-1] + contact_height
    # compute shaft length so lead does not extend above max_fiber_z (if provided)
    default_shaft_length = 100.0 * scale
    shaft_length = default_shaft_length
    if max_fiber_z is not None:
        # uvz is unit vector along lead axis; component in world Z is uvz[2]
        try:
            uvz = np.array([rotation_matrix[2,0], rotation_matrix[2,1], rotation_matrix[2,2]])
            uvz_z = uvz[2]
        except Exception:
            uvz_z = 0.0

        # if the lead axis has a positive Z component, limit the shaft so its highest point
        # does not exceed max_fiber_z. Add a small margin.
        if uvz_z > 1e-6:
            # highest z contributed by shaft end approx = z_temp[0] + uvz_z*(shaft_pos + shaft_length)
            max_allowed = (max_fiber_z - z_temp[0]) / uvz_z - shaft_pos
            if max_allowed < 0:
                # no room for shaft; set to zero
                shaft_length = 0.0
            else:
                shaft_length = min(default_shaft_length, max_allowed - 0.1 * scale)

    create_cylinder_mesh(max(0.0, shaft_length), shaft_pos, lead_color)

def plot_electrode_cube(fig, bounds, color='blue', opacity=0.35):
    # Deprecated - now using detailed electrode visualization
    pass

def render_scene_plotly(fibers, thresholds_row, voltage_limit, leftLeadPos, title='', show_axes=False, electrode_config=None):
    """Render the scene using Plotly."""
    
    # Create figure with subplots
    fig = make_subplots(
        rows=1, cols=1,
        specs=[[{'type': 'scene'}]]
    )

    # Determine max fiber Z so we can cap the lead height
    max_fiber_z = None
    if fibers:
        try:
            max_fiber_z = max(p[2] for fib in fibers for p in fib)
        except Exception:
            max_fiber_z = None

    # Plot each fiber (densify for higher-definition)
    def densify_points(xs, ys, zs, factor=4):
        if factor <= 1:
            return xs, ys, zs
        new_xs = []
        new_ys = []
        new_zs = []
        for i in range(len(xs) - 1):
            x0, x1 = xs[i], xs[i+1]
            y0, y1 = ys[i], ys[i+1]
            z0, z1 = zs[i], zs[i+1]
            for t in np.linspace(0, 1, factor, endpoint=False):
                new_xs.append(x0 + (x1 - x0) * t)
                new_ys.append(y0 + (y1 - y0) * t)
                new_zs.append(z0 + (z1 - z0) * t)
        # append last point
        new_xs.append(xs[-1])
        new_ys.append(ys[-1])
        new_zs.append(zs[-1])
        return new_xs, new_ys, new_zs

    activated = []
    for i, fib in enumerate(fibers):
        thr = math.inf
        if i < len(thresholds_row) and thresholds_row[i] is not None:
            try:
                thr = float(thresholds_row[i])
            except Exception:
                thr = math.inf
                
        # Determine color based on threshold
        color = 'red' if thr < voltage_limit else 'black'
        if thr < voltage_limit:
            activated.append(i)
            
        # Extract x, y, z coordinates
        xs = [p[0] for p in fib]
        ys = [p[1] for p in fib]
        zs = [p[2] for p in fib]

        # densify points for higher-definition lines
        xs, ys, zs = densify_points(xs, ys, zs, factor=4)
        
        # Add line trace for the fiber
        fig.add_trace(
            go.Scatter3d(
                x=xs,
                y=ys,
                z=zs,
                mode='lines',
                line=dict(color=color, width=2),
                name=f'Fiber {i}',
                showlegend=False
            )
        )

    # Add detailed electrode visualization (cap height at max_fiber_z)
    plot_electrode(fig, leftLeadPos, max_fiber_z=max_fiber_z, electrode_config=electrode_config)

    # Build scene dict; respect show_axes flag
    scene_dict = dict(
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.5)
        ),
        aspectmode='data'
    )

    if not show_axes:
        # hide axes, ticks, grid and zeroline
        scene_dict.update(
            xaxis=dict(visible=False, showticklabels=False, showgrid=False, zeroline=False),
            yaxis=dict(visible=False, showticklabels=False, showgrid=False, zeroline=False),
            zaxis=dict(visible=False, showticklabels=False, showgrid=False, zeroline=False)
        )

    # Update layout
    fig.update_layout(
        title=dict(
            text=title,
            x=0.5,
            xanchor='center'
        ),
        scene=scene_dict,
        showlegend=False,
        margin=dict(l=0, r=0, t=30, b=0)
    )
    
    return fig, activated

def plot_activation_plotly(fibers, pulse_widths, thresholds, voltage_limit, out_folder, leftLeadPos, interactive_pw=None, show_axes=False, electrode_config=None):
    """Generate Plotly visualizations for fiber activation."""
    os.makedirs(out_folder, exist_ok=True)
    
    n_fibers = len(fibers)
    n_thr_rows = len(thresholds[0]) if thresholds else 0
    if n_thr_rows < n_fibers:
        print(f"Warning: thresholds provided for {n_thr_rows} fibers but tract has {n_fibers}")

    activation_summary = {}

    if interactive_pw is not None:
        # Single interactive plot
        pw_idx = int(interactive_pw)
        if pw_idx < 0 or pw_idx >= len(pulse_widths):
            raise ValueError(f"--interactive_pw {pw_idx} out of range [0, {len(pulse_widths)-1}]")
        
        pw = pulse_widths[pw_idx]
        row = thresholds[pw_idx] if pw_idx < len(thresholds) else []
        title = f"Pulse width: {pw} μs (index {pw_idx})"
        
        print(f"Rendering interactive view for PW index {pw_idx} (PW={pw})...")
        fig, activated = render_scene_plotly(
            fibers, row, voltage_limit, leftLeadPos,
            title=title, show_axes=show_axes, electrode_config=electrode_config
        )
        
        # Save HTML for interactive viewing (zero-padded index for correct sorting)
        out_html = os.path.join(out_folder, f"activation_pw_{pw_idx:02d}.html")
        fig.write_html(out_html)
        print(f"Saved interactive plot to {out_html}")
        
        # Also save static image
        out_png = os.path.join(out_folder, f"activation_pw_{pw_idx:02d}.png")
        fig.write_image(out_png)
        print(f"Saved static image to {out_png}")
        
        activation_summary[str(pw)] = activated
        
    else:
        # Batch processing
        for pw_idx, pw in enumerate(pulse_widths):
            row = thresholds[pw_idx] if pw_idx < len(thresholds) else []
            title = f"Pulse width: {pw} μs (index {pw_idx})"
            
            print(f"Rendering PW index {pw_idx} (PW={pw})...")
            fig, activated = render_scene_plotly(
                fibers, row, voltage_limit, leftLeadPos,
                title=title, show_axes=show_axes, electrode_config=electrode_config
            )
            
            # Save both interactive HTML and static PNG
            # Use :02d format to ensure files like pw_02 come before pw_10 in directory listing
            out_html = os.path.join(out_folder, f"activation_pw_{pw_idx:02d}.html")
            out_png = os.path.join(out_folder, f"activation_pw_{pw_idx:02d}.png")
            
            fig.write_html(out_html)
            fig.write_image(out_png)
            print(f"Saved {out_html} and {out_png} (activated: {len(activated)})")
            
            activation_summary[str(pw)] = activated

    # Write activation summary
    summary_path = os.path.join(out_folder, 'activation_summary.json')
    with open(summary_path, 'w') as f:
        json.dump({
            'pulse_widths': pulse_widths,
            'voltage_limit': voltage_limit,
            'activated': activation_summary
        }, f, indent=2)
    print(f"Wrote activation summary to {summary_path}")

def main():
    parser = argparse.ArgumentParser(description="Plot activation by pulse width using Plotly (interactive 3D in browser).")
    
    # Required named arguments
    parser.add_argument("--tract", required=True, dest="tract_file",
                        help="Path to the fiber tract file (.txt)")
    parser.add_argument("--results", required=True, dest="thresholds_json",
                        help="Path to the simulation results JSON file")
    parser.add_argument("--output", required=True, dest="out_folder",
                        help="Directory to save output images and HTML files")
    
    # Optional arguments with defaults
    parser.add_argument("--voltage", type=float, dest="voltage_limit", default=3.0,
                        help="Threshold voltage (V) to determine activation (default: 3.0)")
    parser.add_argument("--cond", choices=["anisotropic", "isotropic"], dest="conductivity", default="anisotropic",
                        help="Conductivity type (default: anisotropic)")
    
    # Extra features
    parser.add_argument("--show_axes", action="store_true",
                       help="Show XYZ axes in the Plotly 3D scene (default: hidden)")
    parser.add_argument("--filter_indices", type=str, default=None,
                        help="Path to a text/JSON file listing the indices of fibers that were simulated (handles mismatches).")
    parser.add_argument("--all_fibers", action="store_true",
                        help="Plot ALL fibers at full length, ignoring the filter_indices file. "
                             "Fibers without a threshold are drawn as inactive.")
    parser.add_argument("--interactive_pw", type=int, default=None,
                       help="Render a single pulse width (index) in an interactive plot")
    parser.add_argument("--electrode_config", type=str, default=None,
                        help="Electrode contact configuration, e.g. '01-23', '+012-3', '-0+1-2+3'. "
                             "'-' = cathode (red), '+' = anode (blue), unmarked = inactive (grey).")
    
    args = parser.parse_args()

    # create a subfolder so plotly outputs don't mix with other exporters
    base_out = args.out_folder
    mkdirp(base_out)

    # Use the Mayavi-derived lead coordinates provided by the user
    leftLeadPos = [[167, 161], [223, 222], [143, 159]]

    print(f"Reading {args.tract_file}...")
    fibers = read_tract_file(args.tract_file)
    print(f"Read {len(fibers)} fibers")
    
    if args.all_fibers:
        print(f"--all_fibers: showing all {len(fibers)} fibers (no FEM filter).")
    elif args.filter_indices:
        print(f"Loading filter indices from {args.filter_indices}...")
        valid_indices = load_valid_indices(args.filter_indices)
        
        # Filter and reorder fibers to match the sequential results in thresholds_json
        # Result "0" -> Fiber valid_indices[0], Result "1" -> Fiber valid_indices[1], etc.
        try:
            filtered_fibers = [fibers[idx] for idx in valid_indices]
            print(f"Filtered fibers: kept {len(filtered_fibers)} of {len(fibers)} original fibers.")
            fibers = filtered_fibers
        except IndexError as e:
            print(f"Error: filter index {e} is out of bounds for the loaded tract file.")
            sys.exit(1)

    electrode_config = parse_electrode_config(args.electrode_config) if args.electrode_config else None

    print(f"Reading {args.thresholds_json}...")
    thresholds = load_thresholds(args.thresholds_json)
    if len(thresholds) != len(pulse_widths):
        print(f"Warning: expected {len(pulse_widths)} pulse widths but got {len(thresholds)}")

    plot_activation_plotly(
        fibers, pulse_widths, thresholds, args.voltage_limit,
        base_out, leftLeadPos, args.interactive_pw, show_axes=args.show_axes,
        electrode_config=electrode_config
    )

if __name__ == "__main__":
    main()