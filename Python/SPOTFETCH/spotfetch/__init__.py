"""
spotfetch: A Python module for X-ray diffraction spot analysis.

Submodules:
- io: Handles input/output operations.
- processing: Contains functions for data processing and spot detection.
- visualization: Tools for plotting and visualizing diffraction patterns.

Initialization:
- Sets up directory paths for module organization.
"""

import os

# Define the base directory of the module
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Define directories for submodules
SUBMODULES = {
    "processing": os.path.join(BASE_DIR, "processing"),
    "detectors": os.path.join(BASE_DIR, "detectors"),
    "labeler": os.path.join(BASE_DIR, "labeler"),
    "visualization": os.path.join(BASE_DIR, "visualization"),
}

# Add directories to Python path (optional if they are directly used as packages)
# for submodule, path in SUBMODULES.items():
#     if os.path.exists(path):
        # __path__.append(path)

# Import key functionality from submodules (if needed)
# from .io import read_data, write_data
# from .processing import detect_spots, process_image
# from .visualization import plot_diffraction, show_image

# Module metadata
__version__ = "0.1.0"
__author__ = "Daniel Banco"
__email__ = "daniel.banco@tufts.edu"
# __all__ = ["detectors", "labeler", "visualization", "read_data", "write_data", 
#            "detect_spots", "process_image", "plot_diffraction", "show_image"]
