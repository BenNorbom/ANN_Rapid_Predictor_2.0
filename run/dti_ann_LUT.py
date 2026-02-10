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


# COMMAND LINE INPUTS: 
electrode1File = sys.argv[1]
tractFile = sys.argv[2]
ANN_model = sys.argv[3]
output_json = sys.argv[4]
centering_strategy = sys.argv[5]
tractType = Tract(sys.argv[6])
electrode1FileConductivity = Conductivity(sys.argv[7])
regression_mode = sys.argv[8]  # "reg" or "class"
EXP_EXTRAPOLATE = True


pulse_widths = [60, 75, 90, 105, 120, 135, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500]
node_to_node = 0.5

fem = FEM.FEMgrid(electrode1File)
grid_e1 = fem.get3dGrid()
fem_bounds = fem.getFEMBounds()

test_dti = process_DTI.DTI_tracts(
    tractFile, 
    fem_bounds, 
    node_to_node,
    tractType,
    electrode1FileConductivity
)
ANN_model = ann_predict_lib.ANN(ANN_model)
ann_hparam_dict = ANN_model.get_hparam_dict()
num_ecs = ann_hparam_dict["num_ecs"]
num_fsds = ann_hparam_dict["num_fsds"]
num_ssds = ann_hparam_dict["num_ssds"]

# Use command line input for regression mode instead of model hparam
ANN_REGRESSION = 1 if regression_mode.lower() == "reg" else 0

# max operation is to ensure that NUMBER_OF_SPATIAL_VALUES valid even if e.g. ssd only used
NUMBER_OF_SPATIAL_VALUES = max(ann_hparam_dict["num_ecs"], ann_hparam_dict["num_fsds"], ann_hparam_dict["num_ssds"])

xNodeComp, yNodeComp, zNodeComp = test_dti.getNodeCompPos()

input_array = []
problem_inds = []

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

for fib in range(len(xNodeComp)):
    # Get compartmental EC potentials from previously made 3d-grid
    fiberLVoltages = []
    for i in range(len(xNodeComp[fib])):
        try:
                    
            fiberLVoltages.append(float(grid_e1( [xNodeComp[fib][i], yNodeComp[fib][i], zNodeComp[fib][i]] )))
                    
        except Exception as e:
            print("WARNING: 3d-position out of COMSOL range! X = " + str(xNodeComp[fib][i]) + ", Y = " + str(yNodeComp[fib][i]) + ", Z = " + str(zNodeComp[fib][i]))
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
    print(f"Input sizes for fiber {fib}: {input_sizes}")

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
            print("too_close_11 and too_close")
            ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ecs, center_ind)
            fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_fsds, center_ind)
            ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ssds, center_ind)
        else:
            print("too_close but not too_close_11")
            ecs_extrapolated, center_ind = extrapolate_values(fiberLVoltages_nodes, center_ind)
            fiber_dti = fiber_DTI.Fiber(ecs_extrapolated)
            #center_ind = fiber_dti.getCenterInd(centering_strategy, NUMBER_OF_SPATIAL_VALUES)
            ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ec", num_ecs, center_ind)
            fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("fsd", num_fsds, center_ind)
            ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ssd", num_ssds, center_ind)
            
            # Check if any None values were returned (boundary issue after extrapolation)
            if (None in ecs_at_nodes_around_center or None in fsds_at_nodes_around_center or None in ssds_at_nodes_around_center):
                problem_inds.append(fib)
                print("Extrapolated fiber still has boundary issues - marking as problem fiber")
                ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ecs, center_ind)
                fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_fsds, center_ind)
                ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ssds, center_ind)
    else:
        print("neither")
        ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ec", num_ecs, center_ind)
        fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("fsd", num_fsds, center_ind)
        ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("ssd", num_ssds, center_ind)
        
        # Check if any None values were returned (boundary issue on normal fibers)
        if (None in ecs_at_nodes_around_center or None in fsds_at_nodes_around_center or None in ssds_at_nodes_around_center):
            problem_inds.append(fib)
            print("Normal fiber has boundary issues - marking as problem fiber")
            ecs_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ecs, center_ind)
            fsds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_fsds, center_ind)
            ssds_at_nodes_around_center = fiber_dti.getValAroundCenterNode("err", num_ssds, center_ind)

    for pwd in pulse_widths:

        test_features = [pwd/1000]
        test_features.extend(ecs_at_nodes_around_center)
        test_features.extend(fsds_at_nodes_around_center)
        test_features.extend(ssds_at_nodes_around_center)

        if np.array(test_features).shape[0] != (1 + num_ecs + num_fsds + num_ssds):
            print(f"Fiber {fib}, Pulse width {pwd} us: Features: {np.array(test_features).shape}")

        input_array.append(test_features)

# Flatten input_array and iterate over its elements
flattened_input = [item for sublist in input_array for item in sublist]
for value in flattened_input:
    if value is None:
        print("VALUE IS NONE!")
if ANN_REGRESSION == 1:
    prediction_start_time = time.time()
    ANN_prediction = ANN_model.batch_predict_threshold_reg(input_array)
    prediction_stop_time = time.time()
else:
    prediction_start_time = time.time()
    ANN_prediction = ANN_model.batch_predict_threshold(input_array, 100)
    prediction_stop_time = time.time()

print()
print(str(len(input_array)) + " ANN predictions took " + str(prediction_stop_time - prediction_start_time) + " s")
print()

result_dict = {}

total_fibs = len(xNodeComp) 
good_fibs = total_fibs - len(problem_inds)
pulse_width_count = len(pulse_widths)

result_dict["problem_inds"] = problem_inds
result_dict["valid_inds"] = [i for i in range(total_fibs) if i not in problem_inds]

for i in range(len(ANN_prediction)):
    pulse_width = pulse_widths[int(i % pulse_width_count)] / 1000
    dti_index = int(i / pulse_width_count)
    
    if pulse_width not in result_dict:
        result_dict[pulse_width] = {}

    result_dict[pulse_width][dti_index] = float(ANN_prediction[i]) 

with open(output_json, 'w') as outfile:
    json.dump(result_dict, outfile, indent=4)

######################################################################################################