#!/usr/bin/env python

from distutils.core import setup, Extension

python = "/usr/include/python2.7/"
numpy = "/usr/local/lib/python2.7/dist-packages/numpy/core/include/numpy/"

setup (name = "sparsecoo",
    version = "0.1",
    ext_modules = [
        Extension (
            "sparsecoo",
            [ "pcoo.c" ],
            library_dirs = ["../lib" ],
            libraries = ['sparsefastmm',],
            include_dirs = [ python, numpy,"../src" ]
       )
    ]
)



