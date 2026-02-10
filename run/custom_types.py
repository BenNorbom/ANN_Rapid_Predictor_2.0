from enum import Enum
import numpy as np

class ShiftMode(Enum):
    """Enumeration to represent shifting modes."""
    ORIGIN = 'origin'
    INDEPENDENT = 'independent'

class Tract(Enum):
    """Enumeration to represent the type of tract."""
    ARTIFICIAL = 'artificial'
    DTI = 'dti'

class Conductivity(Enum):
    """Enumeration to represent conductivity mode."""
    ISOTROPIC = 'isotropic'
    ANISOTROPIC = 'anisotropic'

class Bounds3D:
    """Custom type to represent FEM bounds with validation."""
    def __init__(self, x_bound, y_bound, z_bound):
        self.x_bound = self._validate_bound(x_bound, "X")
        self.y_bound = self._validate_bound(y_bound, "Y")
        self.z_bound = self._validate_bound(z_bound, "Z")

    def _validate_bound(self, bound, axis):
        if len(bound) != 2 or bound[0] >= bound[1]:
            raise ValueError(f"Invalid {axis} bound: {bound}. Must be [min, max].")
        return bound

#FIXME: Add MAC option. Currently not included bc not sure how to load neuron dlls for mac.
class OS(Enum):
    """Enumeration to represent operating systems."""
    LINUX = 'linux'
    WINDOWS = 'windows'

class DtiTypes(Enum):
    DRTT = 'drtt'
    ML = 'ml'
    TRP = 'trp'
    CSA = 'csa'
    CST = 'cst'
    TRA = 'tra'