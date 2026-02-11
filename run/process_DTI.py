'''
    This file enables the processing of diffusion tensor imaging (DTI) tractography based
    fiber tracts. From an input txt containing of 3d coordinate points, fiber
    tracts are individually resampled and spline interpolated to obtain realistic
    fiber trajectories. These tracts are then rediscretized to yield the center
    positions of each axon compartment along the fiber as specified by the MRG axon
    model. The coordinates along the center axis of the DBS lead are used to remove
    tracts that intersect with the lead or encapsulating scar tissue.
'''

import numpy as np
from scipy.interpolate import splprep
from scipy.interpolate import splev
import sys
import os

# Add parent directory to path for custom_types
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

# Import custom types
from custom_types import Tract, Conductivity

# Import fiber_DTI from local directory
import fiber_DTI

class DTI_tracts:    
    
    anisotropic_leftLeadPos = [[167,161],[223,222],[143,159]]
    add_x_shift = -10
    add_y_shift = 20
    add_z_shift = 0
    
    def __init__(
        self, 
        tracts_file: str, 
        fem_bounds, 
        node_to_node: float, 
        tract: Tract, 
        conductivity: Conductivity,
        shift_fibers_to_origin = False,
        spline_smoothing = 0.2,
        node_count_requirement = 21,
        add_shift_fibers = False,
        custom_center = None
    ):
        """
        Initialize the DTI_tracts class, load fiber tracts, apply spline interpolation,
        and filter fibers based on DBS lead intersection.

        Parameters
        ----------
        tracts_file : str
            Path to the input file containing 3D coordinates of fiber tracts.
        fem_bounds : Bounds3D
            Bounds of the finite element model (FEM) in [min, max] format for each axis.
        node_to_node : float
            Distance between adjacent nodes of the fibers.
        tract : Tract
            DTI or Artificial
        conducitivity : Conductivity
            Anisotropic or Isotropic
        shift_fibers_to_origin: bool
            translate the DBS lead and surrounding pathways to the origin
        spline_smoothing : float
            smoothing constant chosen to minimize smoothing and retain small scale curvature
        node_count_requirement : int
            Set a lower limit on the length of fibers to make predictions on
        add_shift_fibers : int
            Translate the pathways independent of the DBS lead
        custom_center : list [x, y, z]
            Optional. If provided, overrides the default lead position and shifts fibers
            such that this center becomes (0,0,0).

        Attributes
        ----------
        node_to_node : float
            Distance between adjacent nodes of the fibers.
        shift_fibers_to_origin: bool
            Translate the DBS lead and surrounding pathways to the origin
        spline_smoothing : float
            Smoothing constant chosen to minimize smoothing and retain small scale curvature
        node_count_requirement : int
            Set a lower limit on the length of fibers to make predictions on
        add_shift_fibers : int
            Translate the pathways independent of the DBS lead
        take_every : int
            Sample the DTI points at some interval as they are voxel-based, causing inaccurate jaggedness.
            Can use '1' for artificial tracts, '3' for DTI tracts.
        leftLeadPos : Bounds3D
            Center of DBS
        xCompPos, yCompPos, zCompPos : List[List[float]]
            Sampled (with take_every specification) lists of nodal compartment positions for each fiber tract
        xCompPos, yCompPos, zCompPos : List[List[float]]
            Spline interpolated lists of nodal compartment positions for each fiber tract
        pre_trunc_xNodalComp, pre_trunc_yNodalComp, pre_trunc_zNodalComp : List[List[float]]
            Pre-truncation nodal compartment positions for each fiber tract
        xNodalComp, yNodalComp, zNodalComp : List[List[float]]
            Final nodal compartment positions for each fiber tract after truncation
        """
        self.node_to_node = node_to_node
        self.shift_fibers_to_origin = shift_fibers_to_origin
        self.spline_smoothing = spline_smoothing
        self.node_count_requirement = node_count_requirement
        self.add_shift_fibers = add_shift_fibers
        self.take_every = 3 if tract==Tract.DTI else 1
        
        if custom_center is not None:
             # Create a virtual vertical lead segment centered at custom_center
             # Length 10mm (arbitrary, only used for intersection checks)
             cx, cy, cz = custom_center
             self.leftLeadPos = [[cx, cx], [cy, cy], [cz-5, cz+5]]
        else:
            self.leftLeadPos = [[167,161],[223,222],[143,159]] if conductivity==Conductivity.ANISOTROPIC else [[0,0],[0,0],[0,10]]

        self.xCompPos = []
        self.yCompPos = []
        self.zCompPos = []

        #### Read in x,y,z positions from txt tract file, sample in desired interval
        with open(tracts_file) as file:
            for line in file:
                coords = [float(k) for k in line.split()]
                ## a single tract is stored on a single line in the following format:
                ## x1 y1 z1 x2 y2 z2 ... xn yn zn
                thisFiberX = coords[0::3]
                thisFiberY = coords[1::3]
                thisFiberZ = coords[2::3]

                ## Sample the positions of each tract to reduce noise/jaggedness
                self.xCompPos.append(thisFiberX[::self.take_every])
                self.yCompPos.append(thisFiberY[::self.take_every])
                self.zCompPos.append(thisFiberZ[::self.take_every])

        #### Spline interpolate the tracts:
        spline_xPosL = []
        spline_yPosL = []
        spline_zPosL = []

        num_true_pts = 1000
        u_fine = np.linspace(0,1,num_true_pts)
        for i in range(len(self.xCompPos)):

            ## Use a very precise (n = 1000) spline to determine the length of the fiber accurately
            tck, u = splprep([self.xCompPos[i],self.yCompPos[i],self.zCompPos[i]], s=self.spline_smoothing)
            x_fine, y_fine, z_fine = splev(u_fine, tck)
            tract_length = self.getLength(x_fine, y_fine, z_fine)
            rounded_tract_length = round(tract_length * 2) / 2

            ## With this length, get the positions of the node compartments based on the internodal distance
            num_step_points = int(rounded_tract_length / self.node_to_node) + 1
            u_step = np.linspace(0,1,num_step_points)
            x_step, y_step, z_step = splev(u_step, tck)

            spline_xPosL.append(x_step)
            spline_yPosL.append(y_step)
            spline_zPosL.append(z_step)

        # No fiber shifting when custom_center is used — fibers stay in their original space.
        # The FEM bounds and voltage lookups are shifted instead (handled in dti_ann_LUT.py).

        if custom_center is None and self.shift_fibers_to_origin:
            # Original shift-to-origin logic (used when conductivity is anisotropic, no custom center)
            dx = self.anisotropic_leftLeadPos[0][0]
            dy = self.anisotropic_leftLeadPos[1][0]
            dz = self.anisotropic_leftLeadPos[2][0]
            self.leftLeadPos = [[0,0],[0,0],[0,10]]
            spline_xPosL = [x - dx for x in spline_xPosL]
            spline_yPosL = [y - dy for y in spline_yPosL]
            spline_zPosL = [z - dz for z in spline_zPosL]

        if self.add_shift_fibers:
            # Also fix this for ragged arrays
            spline_xPosL = [x - self.add_x_shift for x in spline_xPosL]
            spline_yPosL = [y - self.add_y_shift for y in spline_yPosL]
            spline_zPosL = [z - self.add_z_shift for z in spline_zPosL]

        self.xCompPos = spline_xPosL
        self.yCompPos = spline_yPosL
        self.zCompPos = spline_zPosL

        #### Throw out tracts that pass through the DBS lead or surrounding 0.5 mm thick scar tissue
        dist_xPosL = []
        dist_yPosL = []
        dist_zPosL = []

        ## Set up interpolated line down center of DBS lead
        num_true_pts = 100
        tck, u = splprep([self.leftLeadPos[0],self.leftLeadPos[1],self.leftLeadPos[2]], s=.25, k=1)
        u_fine = np.linspace(0,1,num_true_pts)
        x_fine, y_fine, z_fine = splev(u_fine, tck)

        lead_tract_sep = 1.135 # 0.635 # Should be radius of lead OR scar tissue
        for i in range(len(self.xCompPos)): # For each tract in list of tracts
            remove = False
            for j in range(len(self.xCompPos[i])): # For each node in the tract
                for k in range(len(x_fine)): # Loop down center of DBS lead, make sure node is not within scar tissue
                    distance = np.sqrt((self.xCompPos[i][j] - x_fine[k])**2 + (self.yCompPos[i][j] - y_fine[k])**2 + (self.zCompPos[i][j] - z_fine[k])**2)
                    if distance <= lead_tract_sep:
                        remove = True
                        break

                if remove == True:
                    break

            if remove == False:
                dist_xPosL.append(self.xCompPos[i])
                dist_yPosL.append(self.yCompPos[i])
                dist_zPosL.append(self.zCompPos[i])

        self.xCompPos = dist_xPosL
        self.yCompPos = dist_yPosL
        self.zCompPos = dist_zPosL

        self.pre_trunc_xNodalComp = self.xCompPos
        self.pre_trunc_yNodalComp = self.yCompPos
        self.pre_trunc_zNodalComp = self.zCompPos

        #### Truncate tracts to fit within FEM bounds, ensure they meet node-count requirement
        trunc_xPosL = []
        trunc_yPosL = []
        trunc_zPosL = []

        for i in range(len(self.xCompPos)):
            tempFiberX = []
            tempFiberY = []
            tempFiberZ = []

            any_node_in_flag = False
            for j in range(len(self.xCompPos[i])): # For every nodal coordinate per tract, remove if not within FEM bounds
                if self.xCompPos[i][j] >= fem_bounds[0][0] and self.xCompPos[i][j] <= fem_bounds[0][1] and self.yCompPos[i][j] >= fem_bounds[1][0] and self.yCompPos[i][j] <= fem_bounds[1][1] and self.zCompPos[i][j] >= fem_bounds[2][0] and self.zCompPos[i][j] <= fem_bounds[2][1]:
                    tempFiberX.append(self.xCompPos[i][j])
                    tempFiberY.append(self.yCompPos[i][j])
                    tempFiberZ.append(self.zCompPos[i][j])
                    any_node_in_flag = True
                elif any_node_in_flag == False: # Allows for starting FEM bounds being outside to NOT cause no output
                    pass
                else: # Ensures no discontinuities
                    print("WARNING: Possible data loss from fibers going outisde of FEM bounds")
                    break

            """
            # WARNING: If starting from outside FEM bounds nothing will output
            j = 0
            while j < len(self.xCompPos[i]) and self.xCompPos[i][j] >= fem_bounds[0][0] and self.xCompPos[i][j] <= fem_bounds[0][1] and self.yCompPos[i][j] >= fem_bounds[1][0] and self.yCompPos[i][j] <= fem_bounds[1][1] and self.zCompPos[i][j] >= fem_bounds[2][0] and self.zCompPos[i][j] <= fem_bounds[2][1]:
                tempFiberX.append(self.xCompPos[i][j])
                tempFiberY.append(self.yCompPos[i][j])
                tempFiberZ.append(self.zCompPos[i][j])
                j = j + 1
            """
            """
            # WARNING: Can lead to discontinuities if tracts go in and out of FEM bounds
            for j in range(len(self.xCompPos[i])): # For every nodal coordinate per tract, remove if not within FEM bounds
                if self.xCompPos[i][j] >= fem_bounds[0][0] and self.xCompPos[i][j] <= fem_bounds[0][1] and self.yCompPos[i][j] >= fem_bounds[1][0] and self.yCompPos[i][j] <= fem_bounds[1][1] and self.zCompPos[i][j] >= fem_bounds[2][0] and self.zCompPos[i][j] <= fem_bounds[2][1]:
                    tempFiberX.append(self.xCompPos[i][j])
                    tempFiberY.append(self.yCompPos[i][j])
                    tempFiberZ.append(self.zCompPos[i][j])
            """


            if len(tempFiberX) >= self.node_count_requirement: # Ensure that the remaning tracts has the minimum number of nodes
                trunc_xPosL.append(tempFiberX)
                trunc_yPosL.append(tempFiberY)
                trunc_zPosL.append(tempFiberZ)

        self.xCompPos = trunc_xPosL
        self.yCompPos = trunc_yPosL
        self.zCompPos = trunc_zPosL

        self.xNodalComp = self.xCompPos
        self.yNodalComp = self.yCompPos
        self.zNodalComp = self.zCompPos

        print("Number of tracts that meet all critera: " + str(len(self.xNodalComp)))

    def getAllComps(self, lin_comp_pos):
        """Get the compartmental positions for each fiber tract.

        Parameters
        ----------
        lin_comp_pos : list
            List of positions for compartmental nodes.

        Returns
        -------
        None
        """
        self.xAllComp = []
        self.yAllComp = []
        self.zAllComp = []

        #### Linearly interpolate between nodal points to get internodal compartment coordinates
        for i in range(len(self.xNodalComp)):
            xCompTemp = []
            yCompTemp = []
            zCompTemp = []
            for j in range(0,len(self.xNodalComp[i])-1):
                for k in range(len(lin_comp_pos)):
                    xCompTemp.append(self.xNodalComp[i][j] + ((self.xNodalComp[i][j+1] - self.xNodalComp[i][j]) * lin_comp_pos[k] / self.node_to_node))
                    yCompTemp.append(self.yNodalComp[i][j] + ((self.yNodalComp[i][j+1] - self.yNodalComp[i][j]) * lin_comp_pos[k] / self.node_to_node))
                    zCompTemp.append(self.zNodalComp[i][j] + ((self.zNodalComp[i][j+1] - self.zNodalComp[i][j]) * lin_comp_pos[k] / self.node_to_node))              
                    
            #add last nodal point 
            xCompTemp.append(self.xNodalComp[i][len(self.xNodalComp[i])-1])
            yCompTemp.append(self.yNodalComp[i][len(self.yNodalComp[i])-1])
            zCompTemp.append(self.zNodalComp[i][len(self.zNodalComp[i])-1])

            self.xAllComp.append(xCompTemp)
            self.yAllComp.append(yCompTemp)
            self.zAllComp.append(zCompTemp)


    def getNodeCompPos(self):
        """Get the compartmental positions of nodes.
        
        Returns
        -------
        tuple
            xNodalComp, yNodalComp, zNodalComp
        """
        return self.xNodalComp, self.yNodalComp, self.zNodalComp

    def getPreTruncNodeCompPos(self):
        """Get the pre-truncation compartmental positions of nodes.
        
        Returns
        -------
        tuple
            pre_trunc_xNodalComp, pre_trunc_yNodalComp, pre_trunc_zNodalComp
        """
        return self.pre_trunc_xNodalComp, self.pre_trunc_yNodalComp, self.pre_trunc_zNodalComp

    def getAllCompPos(self):
        """Get the compartmental positions of all nodes.

        Returns
        -------
        tuple
            xAllComp, yAllComp, zAllComp
        """
        return self.xAllComp, self.yAllComp, self.zAllComp

    def getLength(self,xs,ys,zs):
        """Compute the length of a fiber tract.

        Parameters
        ----------
        xs : list
            List of x-coordinates of the fiber tract.
        ys : list
            List of y-coordinates of the fiber tract.
        zs : list
            List of z-coordinates of the fiber tract.

        Returns
        -------
        float
            Length of the fiber tract.
        """
        l = 0
        for i in range(len(xs) - 1):
            dl = np.sqrt((xs[i + 1] - xs[i])**2 + (ys[i + 1] - ys[i])**2 + (zs[i + 1] - zs[i])**2)
            l = l + dl
        return l

    def getTractCount(self):
        """Get the number of tracts.

        Returns
        -------
        int
            Number of tracts.
        """
        return int(len(self.xNodalComp))

    def getNodeCount(self, index):
        """Get the number of nodes in a specific tract.

        Parameters
        ----------
        index : int 
            Index of the tract.

        Returns
        -------
        int
            Number of nodes in the specified tract.
        """
        return int(len(self.xNodalComp[index]))

    def getLeadCoordinates(self):
        """Get the coordinates of the DBS lead.

        Returns
        -------
        list
            Coordinates of the DBS lead.
        """
        return self.leftLeadPos
    
    def getCompartmentalEcPotentials(self, fem3dGrid, fib_ind):
        """Get compartmental EC potentials from previously made 3d-grid for a specific fiber.

        Parameters
        ----------
        fem3dGrid : function
            Function to compute the 3D grid.
        fib_ind : int
            Index of the fiber.

        Returns
        -------
        list
            List of extracellular potentials for the specified fiber.
        """
        fiberLVoltages = []
        grid_e1 = fem3dGrid
        xNodeComp, yNodeComp, zNodeComp = self.getNodeCompPos()

        for nodal_ind in range(len(xNodeComp[fib_ind])):
            try:
                fiberLVoltages.append(float(grid_e1( [xNodeComp[fib_ind][nodal_ind], yNodeComp[fib_ind][nodal_ind], zNodeComp[fib_ind][nodal_ind]] )))
                        
            except Exception as e:
                print("WARNING: 3d-position out of COMSOL range! X = " + str(xNodeComp[fib_ind][nodal_ind]) + ", Y = " + str(yNodeComp[fib_ind][nodal_ind]) + ", Z = " + str(zNodeComp[fib_ind][nodal_ind]))
                pass
        
        return fiberLVoltages

    def getEcsAtNodes(self, ecs_all, compartDiv):
        """Get extracellular potentials at Nodes of Ranvier for a given fiber.

        This function samples the extracellular potentials at the nodes of Ranvier
        based on the specified compartment division.
        
        Parameters
        ----------
        ecs_all : list
            A list of extracellular potentials.
        compartDiv : int
            The compartment division factor.
            
        Returns
        -------
        list
            A list of extracellular potentials at Nodes of Ranvier.
        """
        return [k for k in ecs_all[::compartDiv]]