

from pytadbit._version import __version__
from pytadbit.tadbit import tadbit, batch_tadbit
from pytadbit.chromosome import Chromosome
from pytadbit.experiment import Experiment
from pytadbit.chromosome import load_chromosome
from pytadbit.imp.structuralmodels import StructuralModels
from pytadbit.imp.structuralmodels import load_structuralmodels
from pytadbit.imp.impmodel import load_impmodel_from_cmm
from pytadbit.imp.impmodel import load_impmodel_from_xyz
from pytadbit.imp.impmodel import IMPmodel
try:
    from pytadbit.imp.impoptimizer import IMPoptimizer
except ImportError:
    from warnings import warn
    warn('IMP not found, check PYTHONPATH\n')


def get_dependencies_verison():
    versions = {}
    try:
        import IMP
        versions['IMP'] = IMP.kernel.get_module_version()
        # versions['IMP'] += '(random seed indexed 1 = %s)' % ()
    except ImportError:
        versions['IMP'] = 'Not found'
    try:
        import scipy
        versions['scipy'] = scipy.__version__
    except ImportError:
        versions['scipy'] = 'Not found'
    try:
        import numpy
        versions['numpy'] = numpy.__version__
    except ImportError:
        versions['numpy'] = 'Not found'
    try:
        import matplotlib
        versions['matplotlib'] = matplotlib.__version__
    except ImportError:
        versions['matplotlib'] = 'Not found'
    from subprocess import Popen, PIPE
    try:
        mcl, err = Popen(['mcl', '--version'], stdout=PIPE,
                         stderr=PIPE).communicate()
        versions['MCL'] = mcl.split()[1]
    except:
        versions['MCL'] = 'Not found'
    try:
        chi, err = Popen(['chimera', '--version'], stdout=PIPE,
                         stderr=PIPE).communicate()
        versions['Chimera'] = chi.strip()
    except:
        versions['Chimera'] = 'Not found'

    return '\n'.join(['%15s : %s' % (k, versions[k]) for k in
                      sorted(versions.keys())])
