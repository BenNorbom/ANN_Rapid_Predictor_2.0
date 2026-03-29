'''
    This script predicts the activation of the DTI-tractography based fiber trajectories
    using an ANN, across all desired pulse widths. The results are written out to a json
    LUT, to be compared with result LUTs from other predictor models.
'''

# Imports
import sys
import numpy as np
import json
import math
import time
import os

# Add parent directory to path for custom_types
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.dirname(__file__))

# Import local modules from run directory
import process_DTI
import fiber_DTI
import FEM
import ann_predict_lib

# Import custom types from parent directory
try:
    from custom_types import *
except ImportError:
    print("Warning: Could not import custom_types, using defaults")
    class Tract:
        def __init__(self, val):
            self.value = val
    class Conductivity:
        def __init__(self, val):
            self.value = val


def print_progress(current, total, prefix='', bar_length=40):
    """Print a simple progress bar to the terminal."""
    fraction = current / total if total > 0 else 1.0
    filled = int(bar_length * fraction)
    bar = '\u2588' * filled + '\u2591' * (bar_length - filled)
    percent = fraction * 100
    sys.stdout.write(f'\r{prefix} |{bar}| {percent:5.1f}% ({current}/{total})')
    sys.stdout.flush()
    if current >= total:
        sys.stdout.write('\n')
        sys.stdout.flush()


# COMMAND LINE INPUTS: 
electrode1File = sys.argv[1]
tractFile = sys.argv[2]
ANN_model = sys.argv[3]
output_json = sys.argv[4]
centering_strategy = sys.argv[5]
tractType = Tract(sys.argv[6])

# Parse conductivity or explicit coordinates
custom_center = None
try:
    # Try to parse 3 float coordinates
    cx = float(sys.argv[7])
    cy = float(sys.argv[8])
    cz = float(sys.argv[9])
    custom_center = [cx, cy, cz]
    
    # If successful, we assume the next argument is regression mode
    # Defaulting conductivity to ANISOTROPIC for model bounds purposes if not specified? 
    # Usually ANN models are trained with specific conductivity assumptions.
    # The 'anisotropic' vs 'isotropic' arg mainly affects the FEM bounds (Active Region).
    # We will default to Anisotropic bounds as they are generally wider/standard for this project.
    electrode1FileConductivity = Conductivity("anisotropic") 
    regression_mode = sys.argv[10]
    
except ValueError:
    # Fallback to original behavior: arg 7 is conductivity string
    electrode1FileConductivity = Conductivity(sys.argv[7])
    regression_mode = sys.argv[8]

ANN_REGRESSION = 1 if regression_mode.lower() == "reg" else 0
EXP_EXTRAPOLATE = True


pulse_widths = [60, 75, 90, 105, 120, 135, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500]
node_to_node = 0.5

print("Loading FEM grid...")
fem_load_start = time.time()
fem = FEM.FEMgrid(electrode1File)
grid_e1 = fem.get3dGrid()
fem_load_time = time.time() - fem_load_start
print(f"FEM grid loaded in {fem_load_time:.2f}s")
fem_bounds_original = fem.getFEMBounds()

# When custom_center is provided, the 3 coordinates specify the CENTER
# of the FEM bounding box in the fiber's coordinate space.
# Fibers are NOT moved. Instead, FEM bounds are re-centered and voltage
# lookups translate fiber coords back to the original FEM grid space.
fem_offset = [0.0, 0.0, 0.0]  # offset to subtract from fiber coords when querying FEM grid

if custom_center is not None:
    # Compute original FEM center
    orig_cx = (fem_bounds_original[0][0] + fem_bounds_original[0][1]) / 2.0
    orig_cy = (fem_bounds_original[1][0] + fem_bounds_original[1][1]) / 2.0
    orig_cz = (fem_bounds_original[2][0] + fem_bounds_original[2][1]) / 2.0
    
    # Half-widths
    hw_x = (fem_bounds_original[0][1] - fem_bounds_original[0][0]) / 2.0
    hw_y = (fem_bounds_original[1][1] - fem_bounds_original[1][0]) / 2.0
    hw_z = (fem_bounds_original[2][1] - fem_bounds_original[2][0]) / 2.0
    
    # New center is the user-specified custom_center
    new_cx, new_cy, new_cz = custom_center
    
    # Shift FEM bounds to be centered at custom_center
    fem_bounds = [
        [new_cx - hw_x, new_cx + hw_x],
        [new_cy - hw_y, new_cy + hw_y],
        [new_cz - hw_z, new_cz + hw_z],
    ]
    
    # Offset: to query the original FEM grid, subtract this from fiber coords
    fem_offset = [new_cx - orig_cx, new_cy - orig_cy, new_cz - orig_cz]
    
    print(f"Original FEM center: ({orig_cx}, {orig_cy}, {orig_cz})")
    print(f"New FEM center (custom): ({new_cx}, {new_cy}, {new_cz})")
    print(f"Shifted FEM bounds to: X={fem_bounds[0]}, Y={fem_bounds[1]}, Z={fem_bounds[2]}")
    print(f"FEM query offset: ({fem_offset[0]}, {fem_offset[1]}, {fem_offset[2]})")
else:
    fem_bounds = fem_bounds_original

print("Processing DTI tracts...")
dti_start = time.time()
test_dti = process_DTI.DTI_tracts(
    tractFile, 
    fem_bounds, 
    node_to_node,
    tractType,
    electrode1FileConductivity,
    custom_center=custom_center
)
dti_time = time.time() - dti_start
print(f"DTI tracts processed in {dti_time:.2f}s")

print("Loading ANN model...")
ann_load_start = time.time()
ANN_model = ann_predict_lib.ANN(ANN_model)
ann_load_time = time.time() - ann_load_start
print(f"ANN model loaded in {ann_load_time:.2f}s")
ann_hparam_dict = ANN_model.get_hparam_dict()
num_ecs = ann_hparam_dict["num_ecs"]
num_fsds = ann_hparam_dict["num_fsds"]
num_ssds = ann_hparam_dict["num_ssds"]

# Use command line input for regression mode instead of model hparam
ANN_REGRESSION = 1 if regression_mode.lower() == "reg" else 0

# max operation is to ensure that NUMBER_OF_SPATIAL_VALUES valid even if e.g. ssd only used
NUMBER_OF_SPATIAL_VALUES = max(ann_hparam_dict["num_ecs"], ann_hparam_dict["num_fsds"], ann_hparam_dict["num_ssds"])

xNodeComp, yNodeComp, zNodeComp = test_dti.getNodeCompPos()
orig_indices = test_dti.getOriginalFiberIndices()

input_array = []
problem_inds = []
fem_warning_count = 0
too_close_both_count = 0
too_close_extrapolated_count = 0
extrapolation_boundary_fail_count = 0
normal_boundary_fail_count = 0
feature_shape_warning_count = 0

def extrapolate_values(values, center_ind):
    num = len(values) - 1 # ex: 11-1 = 10
    num_right = num - center_ind # ex: 10-5 = 5
    num_left = center_ind # ex: 5

    missing_left = NUMBER_OF_SPATIAL_VALUES//2 + 1 - num_left # +1 for two-sided fsds & sds calc
    missing_right = NUMBER_OF_SPATIAL_VALUES//2 + 1 - num_right # +1 for two-sided fsds & sds calc

    def fit_k_from_boundary(side_vals, eps=1e-12):
        k = 5
        y = np.log(np.abs(side_vals[:k]) + eps)
        d = np.arange(k, dtype=float)       # 0 = boundary, increasing inward
        slope, _ = np.polyfit(d, y, 1)
        return slope

    # Left: use outermost-left inward
    k_left = fit_k_from_boundary(values[:])
    v_left_edge = values[0]
    A_left = abs(v_left_edge)
    s_left = -1.0 if v_left_edge == 0 else np.sign(v_left_edge)

    # Right: use outermost-right inward (reverse so boundary is index 0)
    k_right = fit_k_from_boundary(list(values[::-1]))
    v_right_edge = values[-1]
    A_right = abs(v_right_edge)
    s_right = -1.0 if v_right_edge == 0 else np.sign(v_right_edge)

    if missing_left > 0:
        # farthest new left first (largest d), then closer …
        left_add = [s_left * (A_left * math.exp(-k_left * d)) for d in range(missing_left, 0, -1)]
        values = left_add + values

    if missing_right > 0:
        right_add = [s_right * (A_right * math.exp(-k_right * d)) for d in range(1, missing_right + 1)]
        values = values + right_add

    return values, center_ind + missing_left - 1

total_fibers = len(xNodeComp)
fiber_start_time = time.time()
print(f"\nProcessing {total_fibers} fibers...")

for fib in range(total_fibers):
    print_progress(fib + 1, total_fibers, prefix='Processing fibers')
    # Get compartmental EC potentials from previously made 3d-grid
    fiberLVoltages = []
    for i in range(len(xNodeComp[fib])):
        try:
                    
            # Translate fiber coords back to original FEM grid space for voltage lookup
            query_x = xNodeComp[fib][i] - fem_offset[0]
            query_y = yNodeComp[fib][i] - fem_offset[1]
            query_z = zNodeComp[fib][i] - fem_offset[2]
            fiberLVoltages.append(float(grid_e1( [query_x, query_y, query_z] )))
                    
        except Exception as e:
            fem_warning_count += 1
            pass
                
    fiberLVoltages_nodes = fiberLVoltages
    fiber_dti = fiber_DTI.Fiber(fiberLVoltages_nodes)
    len_ecs = len(fiber_dti.ecs)
    len_fsds = len(fiber_dti.fsds)
    len_ssds = len(fiber_dti.ssds)
    center_ind = fiber_dti.getCenterInd(centering_strategy, NUMBER_OF_SPATIAL_VALUES)

    ssds = []
    ecs = []
    fsds = []

    input_sizes = [num_ecs, num_fsds, num_ssds]

    too_close = not fiber_dti.isValidCenterInd
    fiber_11 = fiber_DTI.Fiber(fiberLVoltages_nodes)
    fiber_11.getCenterInd(centering_strategy, 11)
    too_close_11 = not fiber_11.isValidCenterInd

    # center of sent in data (min_ecs or max_ssd) is too close to the edge of the fiber
    ecs_at_nodes_around_center = []
    fsds_at_nodes_around_center = []
    ssds_at_nodes_around_center = []
    if too_close:
        if too_close_11:
            problem_inds.append(fib)
            too_close_both_count += 1
            ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ecs, center_ind)
            fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_fsds, center_ind)
            ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ssds, center_ind)
        else:
            too_close_extrapolated_count += 1
            ecs_extrapolated, center_ind = extrapolate_values(fiberLVoltages_nodes, center_ind)
            fiber_dti = fiber_DTI.Fiber(ecs_extrapolated)
            #center_ind = fiber_dti.getCenterInd(centering_strategy, NUMBER_OF_SPATIAL_VALUES)
            ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ec", num_ecs, center_ind)
            fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("fsd", num_fsds, center_ind)
            ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ssd", num_ssds, center_ind)
            
            # Validate boundaries using full NUMBER_OF_SPATIAL_VALUES so CNN models
            # (which may have num_fsds=0 / num_ssds=0) detect the same problem fibers as ANNs
            ecs_check = fiber_dti.getValAroundCenterNode("ec", NUMBER_OF_SPATIAL_VALUES, center_ind)
            fsd_check = fiber_dti.getValAroundCenterNode("fsd", NUMBER_OF_SPATIAL_VALUES, center_ind)
            ssd_check = fiber_dti.getValAroundCenterNode("ssd", NUMBER_OF_SPATIAL_VALUES, center_ind)
            
            # Check if any None values were returned (boundary issue after extrapolation)
            if (None in ecs_check or None in fsd_check or None in ssd_check):
                problem_inds.append(fib)
                extrapolation_boundary_fail_count += 1
                ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ecs, center_ind)
                fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_fsds, center_ind)
                ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ssds, center_ind)
    else:
        ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ec", num_ecs, center_ind)
        fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("fsd", num_fsds, center_ind)
        ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ssd", num_ssds, center_ind)
        
        # Validate boundaries using full NUMBER_OF_SPATIAL_VALUES so CNN models
        # (which may have num_fsds=0 / num_ssds=0) detect the same problem fibers as ANNs
        ecs_check = fiber_dti.getValAroundCenterNode("ec", NUMBER_OF_SPATIAL_VALUES, center_ind)
        fsd_check = fiber_dti.getValAroundCenterNode("fsd", NUMBER_OF_SPATIAL_VALUES, center_ind)
        ssd_check = fiber_dti.getValAroundCenterNode("ssd", NUMBER_OF_SPATIAL_VALUES, center_ind)
        
        # Check if any None values were returned (boundary issue on normal fibers)
        if (None in ecs_check or None in fsd_check or None in ssd_check):
            problem_inds.append(fib)
            normal_boundary_fail_count += 1
            ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ecs, center_ind)
            fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_fsds, center_ind)
            ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ssds, center_ind)

    for pwd in pulse_widths:

        test_features = [pwd/1000]
        test_features.extend(ecs_at_nodes_around_center)
        test_features.extend(fsds_at_nodes_around_center)
        test_features.extend(ssds_at_nodes_around_center)

        if np.array(test_features).shape[0] != (1 + num_ecs + num_fsds + num_ssds):
            feature_shape_warning_count += 1

        input_array.append(test_features)

# Post-processing summary
fiber_time = time.time() - fiber_start_time
print(f"Fiber processing took {fiber_time:.2f}s")
print(f"  Fibers processed: {total_fibers}")
print(f"  Problem fibers: {len(problem_inds)}")
if fem_warning_count > 0:
    print(f"  FEM out-of-range warnings: {fem_warning_count}")
if too_close_both_count > 0:
    print(f"  Too close (both thresholds): {too_close_both_count}")
if too_close_extrapolated_count > 0:
    print(f"  Extrapolated fibers: {too_close_extrapolated_count}")
if extrapolation_boundary_fail_count > 0:
    print(f"  Extrapolation failures: {extrapolation_boundary_fail_count}")
if normal_boundary_fail_count > 0:
    print(f"  Normal boundary failures: {normal_boundary_fail_count}")
if feature_shape_warning_count > 0:
    print(f"  Feature shape mismatches: {feature_shape_warning_count}")

# Check for None values in input
none_count = sum(1 for sublist in input_array for value in sublist if value is None)
if none_count > 0:
    print(f"  WARNING: {none_count} None values found in input array!")

if len(input_array) == 0:
    print("\nERROR: No valid fibers found within range of the electrode.")
    print("This means all fibers were outside the FEM bounding box after shifting.")
    print("Check that your center coordinates place the electrode inside the fiber bundle.")
    result_dict = {"problem_inds": problem_inds, "valid_inds": [], "warning": "No fibers in range"}
    with open(output_json, 'w') as outfile:
        json.dump(result_dict, outfile, indent=4)
    print(f"Wrote empty results to {output_json}")
    sys.exit(0)

print("\nRunning ANN predictions...")
if ANN_REGRESSION == 1:
    prediction_start_time = time.time()
    ANN_prediction = ANN_model.batch_predict_threshold_reg(input_array)
    prediction_stop_time = time.time()
else:
    prediction_start_time = time.time()
    ANN_prediction = ANN_model.batch_predict_threshold(input_array, 100)
    prediction_stop_time = time.time()

print(f"ANN predictions ({len(input_array)} total) took {prediction_stop_time - prediction_start_time:.2f}s")

result_dict = {}

total_fibs = len(xNodeComp) 
good_fibs = total_fibs - len(problem_inds)
pulse_width_count = len(pulse_widths)

# Map filtered-fiber indices back to original tract-file line numbers
problem_orig = [orig_indices[i] for i in problem_inds]
problem_set = set(problem_inds)
valid_orig = [orig_indices[i] for i in range(total_fibs) if i not in problem_set]

result_dict["problem_inds"] = problem_orig
result_dict["valid_inds"] = valid_orig

for i in range(len(ANN_prediction)):
    pulse_width = pulse_widths[int(i % pulse_width_count)] / 1000
    dti_index = int(i / pulse_width_count)
    orig_fiber_idx = orig_indices[dti_index]
    
    if pulse_width not in result_dict:
        result_dict[pulse_width] = {}

    result_dict[pulse_width][orig_fiber_idx] = float(ANN_prediction[i]) 

with open(output_json, 'w') as outfile:
    json.dump(result_dict, outfile, indent=4)

######################################################################################################