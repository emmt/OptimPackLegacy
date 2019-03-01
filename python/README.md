# LOPPY

LOPPY (for *Legacy OptimPack in PYthon*) let you use VMLM-B algorithm from
OptimPackLegacy in Python.


## Installation

1. Build Cython executables for your OS.  Open a terminal, go to
   `OptimPackLegacy/python/` and run:
   ```sh
   python setup.py build_ext --inplace
   ```

2. Run or have a look at the two following test files. They both minimize the
   Rosenbrock's banana with bound constraint.
   - `test.py` tests if Cython imported functions work properly.
   - `testloppy.py` tests if Python LOPPY Optimizer class works properly.

Note: if the NumPy library cannot be found during Cython build, use the
following command before building Cython module:

```sh
cp -r ~/.local/lib/python3.6/site-packages/numpy/core/include/numpy ~/.local/include
```

where you may change `python3.6` with your actual Python folder.  This issue
was mentionned on [Github](https://github.com/andersbll/cudarray/issues/25).


## Content

| Name                  | Description                                                              |
|-----------------------|--------------------------------------------------------------------------|
| `README.md`           | This file                                                                |
| `setup.py`            | Setup file to compile Cython code                                        |
| `optimpacklegacy.pyx` | Prototypes of OptimPackLegacy functions for Cython                       |
| `loppy.py`            | Python `Optimizer` class based on Cython `optimpacklegacy.pyx` functions |
| `test.py`             | Test Cython imported functions                                           |
| `testloppy.py`        | Test LOPPY, the Python optimizer class                                   |

## Revision history

15 Feb 2019 - RJLF - Add loppy optimizer class

14 Feb 2019 - RJLF - Add `-fPIC` option to GCC, clean package, test Cython call on one optimpacklegacy function

08 Feb 2019 - RJLF - Creation
