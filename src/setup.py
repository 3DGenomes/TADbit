#!/usr/bin/env python

from distutils.core import setup, Extension
from os import path

PATH = path.abspath(path.split(path.realpath(__file__))[0])

def main():
    # c module to find TADs
    pytadbit_module = Extension('pytadbit.tadbit_py',
                                sources=['tadbit_py.c'])
    # c++ module to align and calculate distances between 2 3D models
    eqv_rmsd_module = Extension('pytadbit.eqv_rms_drms',
                                sources=['eqv-tmscore/eqv_rms_drms_py.cpp'])

    setup(
        name        = 'pytadbit',
        version     = '1.0',
        author      = 'Guillaume Filion',
        description = 'Identify TADs in hi-C data',
        ext_modules = [pytadbit_module, eqv_rmsd_module],
        package_dir = {'pytadbit': PATH + '/../pytadbit'},
        packages    = ['pytadbit', 'pytadbit.parsers', 'pytadbit.boundary_aligner',
                       'pytadbit.tad_clustering', 'pytadbit.imp'],
        py_modules  = ["pytadbit"]
    )


if __name__ == '__main__':
    exit(main())
