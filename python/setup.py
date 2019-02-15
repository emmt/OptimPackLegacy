#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 08 17:18:00 2019

Cython entry file to connect optimpacklegacy to Python

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


@author: rfetick
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

op_extension = Extension(
    name="optimpacklegacy",
    sources=["optimpacklegacy.pyx"],
    libraries=["optimpacklegacy"],
    library_dirs=["../src"],
    include_dirs=["../src",numpy.get_include()]
)
setup(
    name="optimpacklegacy",
    ext_modules=cythonize([op_extension])
)
