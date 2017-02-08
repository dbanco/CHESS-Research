try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

import numpy as np
from Cython.Build import cythonize

setup(
  name = 'Circulant Coordinate Ascent for LASSO',
  include_dirs = [np.get_include()],  
  ext_modules = cythonize("circulant_lasso_ca_np.pyx"),
)