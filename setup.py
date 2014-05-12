#!/usr/bin/env python

from distutils.core import setup, Extension
from os import path
from subprocess import Popen, PIPE
from distutils.spawn import find_executable


PATH = path.abspath(path.split(path.realpath(__file__))[0])

TAGS = [
    "Development Status :: 2 - Pre-Alpha",
    "Environment :: Console",
    "Environment :: X11 Applications",
    "Intended Audience :: Developers",
    "Intended Audience :: Other Audience",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License (GPL)",
    "Natural Language :: English",
    "Operating System :: MacOS",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Visualization",
    "Topic :: Software Development :: Libraries :: Python Modules",
    ]



def can_import(modname):
    'Test if a module can be imported '
    try:
        __import__(modname)
    except ImportError:
        return None
    else:
        return True

def ask(string, valid_values, default=-1, case_sensitive=False):
    """ Asks for a keyborad answer """
    v = None
    if not case_sensitive:
        valid_values = [value.lower() for value in valid_values]
    while v not in valid_values:
        v = raw_input("%s [%s]" % (string,','.join(valid_values)))
        if v == '' and default>=0:
            v = valid_values[default]
        if not case_sensitive:
            v = v.lower()
    return v


PATH = path.abspath(path.split(path.realpath(__file__))[0])
PYTHON_DEPENDENCIES = [
    ["numpy"     , "Numpy is required arrays, 1 dimensional interpolation and polynomial fit.", 1],
    ["scipy"     , "Required for clustering and interpolation.", 1],
    ["matplotlib", "Required fot displaying plots.", 0],
    ["IMP"       , "Required for 3D modeling.", 0]]


print "Checking dependencies..."
missing = False
for mname, msg, ex in PYTHON_DEPENDENCIES:
    if not can_import(mname):
        print "  *", mname, "cannot be found in your python installation."
        print "   ->", msg
        if ex:
            missing=True
        else:
            print ("\n  However, you can still install Tadbit and " +
                   "try to fix it afterwards.")
            if ask( "  -> Do you want to continue with the installation?",
                    ["y", "n"]) == "n":
                exit()

if missing:
    exit("Essential dependencies missing, please review and install.\n")


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
                                language = "c",
                                sources=['src/tadbit_py.c'],
                                extra_compile_args=['-std=c99'])
    # c++ module to align and calculate all distances between group of 3D models
    eqv_rmsd_module = Extension('pytadbit.eqv_rms_drms',
                                language = "c++",
                                sources=['src/3d-lib/eqv_rms_drms_py.cpp',
                                         'src/3d-lib/matrices.cc',
                                         'src/3d-lib/3dStats.cpp',
                                         'src/3d-lib/align.cpp'],
                                extra_compile_args=["-ffast-math"])
    # c++ module to align a pair of 3D models
    aligner3d_module = Extension('pytadbit.aligner3d',
                                 language = "c++",
                                 sources=['src/3d-lib/align_py.cpp',
                                          'src/3d-lib/matrices.cc',
                                          'src/3d-lib/3dStats.cpp',
                                          'src/3d-lib/align.cpp'],
                                 extra_compile_args=["-ffast-math"])
    # c++ module to align and calculate consistency of a group of 3D models
    consistency_module = Extension('pytadbit.consistency',
                                   language = "c++",
                                   sources=['src/3d-lib/consistency_py.cpp',
                                            'src/3d-lib/matrices.cc',
                                            'src/3d-lib/3dStats.cpp',
                                            'src/3d-lib/align.cpp'],
                                   extra_compile_args=["-ffast-math"])
    # c++ module to get centroid of a group of 3D models
    centroid_module = Extension('pytadbit.centroid',
                                language = "c++",
                                sources=['src/3d-lib/centroid_py.cpp',
                                         'src/3d-lib/matrices.cc',
                                         'src/3d-lib/3dStats.cpp',
                                         'src/3d-lib/align.cpp'],
                                extra_compile_args=["-ffast-math"])

    # UPDATE version number
    version_full = open(path.join(PATH, '_pytadbit', '_version.py')
                        ).readlines()[0].split('=')[1]
    version_full.strip().replace('"', '')
    version      = '.'.join(version_full.split('.')[:-1])
    revision     = version_full.split('.')[-1]
    # try to use git to check if version number matches
    git_revision = git_version = None
    try:
        git_revision, err = Popen(['git', 'describe'], stdout=PIPE,
                                  stderr=PIPE).communicate()
        git_status, err2 = Popen(['git', 'status'], stdout=PIPE,
                                stderr=PIPE).communicate()
        plus = 'nothing to commit' not in git_status
        git_version  = git_revision.split('-')[0]
        git_revision = str(int(git_revision.split('-')[1]) + plus)
    except OSError:
        git_revision = revision
        git_version  = version
    else:
        if err:
            git_revision = revision
            git_version  = version
    # update version number and write it to _version.py and README files
    revision = git_revision
    version  = git_version
    version_full = '.'.join([version, revision])
    out = open(path.join(PATH, '_pytadbit', '_version.py'), 'w')
    out.write('__version__ = "%s"' % version_full)
    out.close()
    lines = []
    for line in open(path.join(PATH, 'README.rst')):
        if line.startswith('* Current version: '):
            line = '* Current version: ' + version_full
        lines.append(line.strip())
    out = open(path.join(PATH, 'README.rst'), 'w')
    out.write('\n'.join(lines))
    out.close()

    setup(
        name         = 'TADbit',
        version      = version_full,
        author       = 'Davide Bau, Francois Serra, Guillaume Filion and Marc Marti-Renom',
        author_email = 'serra.francois@gmail.com',
        ext_modules  = [pytadbit_module, eqv_rmsd_module, centroid_module,
                        consistency_module, aligner3d_module],
        package_dir  = {'pytadbit': PATH + '/_pytadbit'},
        packages     = ['pytadbit', 'pytadbit.parsers',
                        'pytadbit.boundary_aligner', 'pytadbit.utils',
                        'pytadbit.tad_clustering', 'pytadbit.imp'],
        # py_modules   = ["pytadbit"],
        platforms = "OS Independent",
        license = "GPLv3",
        description  = 'Identification, analysis and modelling of topologically associating domains from Hi-C data',
        long_description = (open("README.rst").read() +
                            open("doc/source/install.rst").read()),
        classifiers  = TAGS,
        provides     = ["pytadbit"],
        keywords     = ["testing"],
        url          = 'https://github.com/3DGenomes/tadbit',
        download_url = 'https://github.com/3DGenomes/tadbit/tarball/master',
    )


if __name__ == '__main__':
    exit(main())
