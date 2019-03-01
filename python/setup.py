# -*- coding: utf-8; mode: python -*-
"""
Script to build Python interface to the OptimPackLegacy library.

- Remove eventual old compilation and open Makefile
	cd src
	rm liboptimpack.a
	rm *.o
	gedit Makefile
- Add now "-fPIC" to compiler options and actually compile it
  	make
- Run Cython build
	py3 setup.py build_ext --inplace
- Test optimpack library in Python
	py3 test.py


@author: rfetick, emmt
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

setup(
    ext_modules = cythonize(
        [Extension(
            name="optimpacklegacy",
            sources=["optimpacklegacy.pyx"],
            libraries=["optimpacklegacy"],
            library_dirs=["."],
            include_dirs=["../src", numpy.get_include()],
            ## The following is to define a macro to avoid warnings about
            ## using deprecated NumPy API.
            #define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        )]
    )
)
