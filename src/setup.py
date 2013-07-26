#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.spawn import find_executable
from subprocess import Popen, PIPE
from os import path, chdir, getcwd
from re import sub
import fileinput

PATH = path.abspath(path.split(path.realpath(__file__))[0])

def compile_TMscore(fortran):
    flags= '' if fortran == 'g77' else '-ffast-math '
    print 'compiling TM_score_consistency with {}'.format(fortran)
    _, stderr = Popen(('{0} -static -O3 {1}-lm -o ' +
                       '{2}/consistency/TMscore_consistency ' +
                       '{2}/consistency/TMscore_consistency.f').format(
                          fortran, flags, PATH), stderr=PIPE,
                      shell=True).communicate()
    print stderr
    return stderr

def compile_eqvrmsd(cpp):
    cur = getcwd()
    chdir(PATH + '/eqv-tmscore/')
    inpf = fileinput.input('makefile', inplace=True)
    for line in inpf:
        line = sub(r'(CC\s*=\s*)[/A-Za-z0-9.+_-]+', '\\1' + cpp + ' ', line)
        print line.rstrip()
    inpf.close()
    
    Popen('make clean', shell=True).communicate()
    _, stderr = Popen('make', stderr=PIPE, shell=True).communicate()
    print stderr
    chdir(cur)
    return stderr


def main():
    # compile fortran
    for gfortran in ['gfortran', 'g77']:
        if not find_executable(gfortran):
            continue
        if compile_TMscore(gfortran):
            continue
        break
    else:
        dit = raw_input('FORTRAN COMPILER NOT FOUND (tried gfortran and g77), ' +
                        'enter a path to retry or "c" to continue (this ' +
                        'executable is OPTIONAL):\n')
        if dit != 'c':
            compile_TMscore(dit)


    # compile c++
    for cpp in ['g++', 'g++-4.4']:
        if not find_executable(cpp):
            continue
        compile_eqvrmsd(cpp)
        break
    else:
        dit = raw_input('C++ COMPILER NOT FOUND (tried g++ and g++-4.4), ' +
                        'enter a path to retry or "c" to continue (this ' +
                        'executable is OPTIONAL):\n')
        if dit != 'c':
            compile_TMscore(dit)


    pytadbit_module = Extension('pytadbit.tadbit_py',
                                sources=['tadbit_py.c'])

    setup(
        name        = 'pytadbit',
        version     = '1.0',
        author      = 'Guillaume Filion',
        description = 'Identify TADs in hi-C data',
        ext_modules = [pytadbit_module],
        package_dir = {'pytadbit': PATH + '/../pytadbit'},
        packages    = ['pytadbit', 'pytadbit.parsers', 'pytadbit.boundary_aligner',
                       'pytadbit.tad_clustering', 'pytadbit.imp'],
        py_modules  = ["pytadbit"],
        scripts=['consistency/TMscore_consistency',
                 'eqv-tmscore/eqv-drmsd-concat-xyz']
    )


if __name__ == '__main__':
    exit(main())
