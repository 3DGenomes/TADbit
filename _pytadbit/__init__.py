
from pytadbit._version import __version__
from os import environ
import locale

# make sure that we are comparing strings in the same way as the bash sort
# command used in mapping
try:
    locale.setlocale(locale.LC_ALL, '.'.join(locale.getdefaultlocale()))
except:
    environ["LANG"] = "en_US.UTF-8"
    locale.setlocale(locale.LC_ALL, '.'.join(locale.getdefaultlocale()))


def get_dependencies_version(dico=False):
    """
    Check versions of TADbit and all dependencies, as well and retieves system
    info. May be used to ensure reproductibility.
    
    :returns: string with description of versions installed
    """
    versions = {'  TADbit': __version__ + '\n\n'}
    try:
        import IMP
        try:
            versions['IMP'] = IMP.get_module_version()
            IMP.random_number_generator.seed(1)
            seed = IMP.random_number_generator()
        except AttributeError:
            versions['IMP'] = IMP.kernel.get_module_version()
            IMP.kernel.random_number_generator.seed(1)
            seed = IMP.kernel.random_number_generator()
        versions['IMP'] += ' (random seed indexed at 1 = %s)' % (seed)
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
        mcl, _ = Popen(['mcl', '--version'], stdout=PIPE,
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
    try:
        chi, err = Popen(['chimera', '--version'], stdout=PIPE,
                         stderr=PIPE).communicate()
        versions['Chimera'] = chi.strip()
    except:
        versions['Chimera'] = 'Not found'
    try:
        uname, err = Popen(['uname', '-rom'], stdout=PIPE,
                           stderr=PIPE).communicate()
        versions[' Machine'] = uname
    except:
        versions[' Machine'] = 'Not found'

    if dico:
        return versions
    else:
        return '\n'.join(['%15s : %s' % (k, versions[k]) for k in
                          sorted(versions.keys())])


from pytadbit.tadbit import tadbit, batch_tadbit
from pytadbit.chromosome import Chromosome
from pytadbit.experiment import Experiment, load_experiment_from_reads
from pytadbit.chromosome import load_chromosome
from pytadbit.imp.structuralmodels import StructuralModels
from pytadbit.imp.structuralmodels import load_structuralmodels
from pytadbit.parsers.hic_parser import load_hic_data_from_reads
from pytadbit.imp.impmodel import load_impmodel_from_cmm
from pytadbit.imp.impmodel import load_impmodel_from_xyz
from pytadbit.imp.impmodel import IMPmodel
from pytadbit.parsers.hic_parser import read_matrix
try:
    from pytadbit.imp.impoptimizer import IMPoptimizer
except ImportError:
    from warnings import warn
    warn('IMP not found, check PYTHONPATH\n')

