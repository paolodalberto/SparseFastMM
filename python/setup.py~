#!/usr/bin/env python

from distutils.core import setup, Extension

python = "/usr/include/python2.7/"
numpy = "/usr/local/lib/python2.7/dist-packages/numpy/core/include/numpy/"

setup (name = "anomaly",
    version = "0.4",
    ext_modules = [
        Extension (
            "anomaly",
            [ "anomalymodule.c" ],
            library_dirs = ["../lib", "../TimeSeries/Compression/bzip2-1.0.5"],
            libraries = ['cdl', 'bz2'],
            include_dirs = [ python, numpy,"../NonParametric/DistanceFunctions", "../Summation", "../Sort", "../TimeSeries/HoltWinters", "../NonParametric", "../TimeSeries/MovingAverage/","../Window", "../TimeSeries",
                             "../TimeSeries/Martingale", "../TimeSeries/Compression", "../TimeSeries/Kernels",
                             "../NonParametric/PValue","../Poset" ]
       )
    ]
)



