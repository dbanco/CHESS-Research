try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

import numpy as np
from Cython.Build import cythonize

setup(
  name = 'Circulant FISTA for LASSO',
  include_dirs = [np.get_include()],  
  ext_modules = cythonize("circulant_fista.pyx"),
)