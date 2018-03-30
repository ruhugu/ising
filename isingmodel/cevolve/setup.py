# -*- coding: utf-8 -*-

# This file compiles the Cython code into an extension module ".so" that can be 
# loaded from pure Python.
# In order to compile, call: python setup.py build_ext --inplace
#       build_ext tells distutils that a external module will be built
#       --inplace builds the module in the current path
# This two are both passed as arguments to the script

from __future__ import (print_function, division, absolute_import)

# distutils is a Python package that provides support for building and 
# installing modules (written in C or Python) into Python.
# https://docs.python.org/3/library/distutils.html
from distutils.core import setup, Extension

# Since we are using Cython, we also have to load addtional pieces specific 
# to Cython.
# TODO ¿REF?
from Cython.Distutils import build_ext

# We need numpy.get_include()
import numpy

# The setup function is meant to provide Distutils with the information
# needed to build, distribute or install the module.
# https://docs.python.org/2/distutils/setupscript.html
# https://docs.python.org/2/distutils/setupscript.html#describing-extension-modules
setup(
    # This is specific to Cython
    cmdclass={'build_ext': build_ext},
    # This is a list with the extension modules we use. Each extension
    # module is described through a Extension instance
    # https://docs.python.org/2/distutils/apiref.html#distutils.core.Extension
    ext_modules=[Extension("cevolve", 
                           # Source files
                           sources=["cevolve.pyx", "c_evolve.c"],
                           # other compile args for gcc
                           #    lgfortran: links with gfortran libraries
                           #    lm: links with math library
                           extra_compile_args=['-lgfortran', '-lm', '-O3'],
                           # other files to link to
                           extra_link_args=["dranxor2/dranxor2C.o"],
                           # Load numpy .h files
                           # https://github.com/numpy/numpy/blob/029863eae86b9df2de4b9a9843ca8f88c99130df/numpy/lib/utils.py#L22
                           include_dirs=[numpy.get_include()])],
)
