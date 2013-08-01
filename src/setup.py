#!/usr/bin/env python

from distutils.core import setup, Extension
from os import path
from distutils.spawn import find_executable

PATH = path.abspath(path.split(path.realpath(__file__))[0])

def main():
    # check if MCL is installed
    if not find_executable('mcl'):
        print('\nWARNING: It is HIGHLY RECOMMENDED to have MCL installed ' +
              '(which do not seems to be).\nIf you are under Debian/Ubuntu' +
              ' just run "apt-get-install mcl".')
        follow = raw_input('\n  You still have the option to follow with the ' +
                           'installation. Do you want to follow? [y/N]')
        if follow.upper() != 'Y' :
            exit('\n    Wise choice :)\n')
    
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
