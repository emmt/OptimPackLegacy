# LOPPY

Legacy OptimPack in PYthon

This is the OptimPackLegacy adaptation in Python. Must be included in the OptimPackLegacy folder

## Getting started

* Build Cython executables for your OS
  * Open a terminal, go to `OptimPackLegacy/python/` and run
  * `python3 setup.py build_ext --inplace`
* Run or have a look at the two following test files. They both minimize the Rosenbrock's banana with bound constraints.
  * `test.py` (test if Cython imported functions work properly)
  * `testloppy.py` (test if Python LOPPY Optimizer class works properly)

## Content

| NAME                              | DESCRIPTION                                                              |
|-----------------------------------|--------------------------------------------------------------------------|
| README                            | This file                                                                |
| setup.py                          | Setup file to compile Cython code                                        |
| optimpacklegacy.pyx               | Prototypes of Optimpack Legacy functions for Cython                      |
| loppy.py                          | Python `Optimizer` class based on Cython `optimpacklegacy.pyx` functions |
| test.py                           | Test Cython imported functions                                           |
| testloppy.py                      | Test LOPPY, the Python optimizer class                                   |

## Revision history

15 Feb 2019 - RJLF - Add loppy optimizer class

14 Feb 2019 - RJLF - Add `-fPIC` option to GCC, clean package, test Cython call on one optimpacklegacy function
                     
08 Feb 2019 - RJLF - Creation


