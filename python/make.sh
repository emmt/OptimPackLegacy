#! /bin/sh

CFLAGS="-I../src"  \
LDFLAGS="-L."     \
    python setup.py build_ext -i
