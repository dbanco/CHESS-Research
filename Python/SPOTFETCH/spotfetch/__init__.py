"""
spotfetch: A Python module for X-ray diffraction spot analysis.

Submodules:
- io: Handles input/output operations.
- processing: Contains functions for data processing and spot detection.
- visualization: Tools for plotting and visualizing diffraction patterns.

"""

from .spotfetch import *
from .visualization import *
from .labeler import *
from .detectors import *
from .tracking import *