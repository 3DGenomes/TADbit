#!/usr/bin/env python

from distutils.core import setup, Extension

pytadbit_module = Extension(
    '_pytadbit',
    sources=['pytadbit_wrap.c', 'tadbit.c'],
)

setup(
   name        = 'pytadbit',
   version     = '1.0',
   author      = 'Guillaume Filion',
   description = 'Identify TADs in hi-C data',
   ext_modules = [pytadbit_module],
   py_modules  = ["pytadbit"],
)
