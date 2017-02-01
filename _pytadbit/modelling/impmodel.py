"""
25 Oct 2016


"""
from pytadbit.modelling.structuralmodel import StructuralModel
from pytadbit.utils.extraviews      import tadbit_savefig
from scipy.interpolate              import spline
from numpy                          import linspace
from warnings                       import warn
from re                             import findall, compile as compil

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')


def load_impmodel_from_cmm(f_name, rand_init=None, radius=None):
    '''
    Loads an IMPmodel object using an cmm file of the form:

    ::

        <marker_set name="1">
          <marker id="1" x="7347.50964739" y="-7743.92836303" z="-8283.39749204" r="0.00990099009901" g="0" b="0.990099009901" radius="500.0" note="1"/>
          <marker id="2" x="7647.90254377" y="-7308.1816344" z="-7387.75932893" r="0.019801980198" g="0" b="0.980198019802" radius="500.0" note="2"/>
          <link id1="1" id2="2" r="1" g="1" b="1" radius="250.0"/>
        </marker_set>

    :params f_name: path where to find the file
    :params None rand_init: IMP random initial number used to generate the model
    :param None radius: radius of each particle

    :return: IMPmodel
    '''

    if not rand_init:
        try:
            rand_init = str(int(f_name.split('.')[-2]))
        except:
            rand_init = None
    model = IMPmodel((('x', []), ('y', []), ('z', []), ('rand_init', rand_init),
                      ('index', 0), ('objfun', 0), ('radius', radius)))
    expr = compil(
        ' x="([0-9.-]+)" y="([0-9.-]+)" z="([0-9.-]+)".* radius="([0-9.]+)"')
    for xxx, yyy, zzz, radius in findall(expr, open(f_name).read()):
        model['x'].append(float(xxx))
        model['y'].append(float(yyy))
        model['z'].append(float(zzz))
    if not model['radius']:
        model['radius'] = float(radius)
    return model

def load_impmodel_from_xyz(f_name, rand_init=None, radius=None):
    """
    Loads an IMPmodel object using an xyz file of the form:

    ::

          # ID              : some identifier
          # SPECIES         : None
          # CELL TYPE       : None
          # EXPERIMENT TYPE : Hi-C
          # RESOLUTION      : 10000
          # ASSEMBLY        : None
          # CHROMOSOME      : 19
          # START           : 1
          # END             : 50
          1  19:1-10000        44.847     412.828    -162.673
          2  19:10001-20000   -55.574     396.869    -129.782

    :params f_name: path where to find the file
    :params None rand_init: IMP random initial number used to generate the model
    :param None radius: radius of each particle

    :return: IMPmodel

    """
    if not rand_init:
        try:
            rand_init = str(int(f_name.split('.')[-2]))
        except:
            rand_init = None
    model = IMPmodel((('x', []), ('y', []), ('z', []), ('rand_init', rand_init),
                      ('index', 0), ('objfun', 0), ('radius', radius)))
    expr = compil('[0-9]+\s[A-Za-z0-9_ ]+:[0-9]+-[0-9]+\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)')
    model['description'] = {}
    for line in open(f_name):
        if line.startswith('# '):
            key, val = line.strip('# ').split(':')
            model['description'][key.strip().lower()] = val.strip()
    for xxx, yyy, zzz in findall(expr, open(f_name).read()):
        model['x'].append(float(xxx))
        model['y'].append(float(yyy))
        model['z'].append(float(zzz))
    return model

def load_impmodel_from_xyz_OLD(f_name, rand_init=None, radius=None,
                               chromosome='UNKNOWN', start=0, resolution=1):
    """
    Loads an IMPmodel object using an xyz file of the form:

    ::

          p1           1      44.847     412.828    -162.673
          p2           2     -55.574     396.869    -129.782

    :params f_name: path where to find the file
    :params None rand_init: IMP random initial number used to generate the model
    :param None radius: radius of each particle

    :return: IMPmodel

    """
    if not rand_init:
        try:
            rand_init = str(int(f_name.split('.')[-2]))
        except:
            rand_init = None
    model = IMPmodel((('x', []), ('y', []), ('z', []), ('rand_init', rand_init),
                      ('objfun', None), ('radius', radius)))
    expr = compil('p[0-9]+\s+[0-9]+\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)')
    for xxx, yyy, zzz in findall(expr, open(f_name).read()):
        model['x'].append(float(xxx))
        model['y'].append(float(yyy))
        model['z'].append(float(zzz))
    model['description'] = {'chromosome':chromosome,
                            'start': start, 'resolution': resolution}

    return model

class IMPmodel(StructuralModel):
    """
    A container for the IMP modeling results. The container is a dictionary
    with the following keys:

    - log_objfun: The list of IMP objective function values
    - objfun: The final objective function value of the corresponding model
    - rand_init: Random number generator feed (needed for model reproducibility)
    - x, y, z: 3D coordinates of each particles. Each represented as a list

    """
    def __str__(self):
        try:
            return ('IMP model ranked %s (%s particles) with: \n' +
                    ' - Final objective function value: %s\n' +
                    ' - random initial value: %s\n' +
                    ' - first coordinates:\n'+
                    '        X      Y      Z\n'+
                    '  %7s%7s%7s\n'+
                    '  %7s%7s%7s\n'+
                    '  %7s%7s%7s\n') % (
                self['index'] + 1,
                len(self['x']), self['objfun'], self['rand_init'],
                int(self['x'][0]), int(self['y'][0]), int(self['z'][0]),
                int(self['x'][1]), int(self['y'][1]), int(self['z'][1]),
                int(self['x'][2]), int(self['y'][2]), int(self['z'][2]))
        except IndexError:
            return ('IMP model of %s particles with: \n' +
                    ' - Final objective function value: %s\n' +
                    ' - random initial value: %s\n' +
                    ' - first coordinates:\n'+
                    '      X    Y    Z\n'+
                    '  %5s%5s%5s\n') % (
                len(self['x']), self['objfun'], self['rand_init'],
                self['x'][0], self['y'][0], self['z'][0])


    def objective_function(self, log=False, smooth=True,
                           axe=None, savefig=None):
        """
        This function plots the objective function value per each Monte-Carlo
        step.

        :param False log: log plot
        :param True smooth: curve smoothing
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).

        """
        show = False
        if not axe:
            fig = plt.figure(figsize=(7, 7))
            axe = fig.add_subplot(111)
            show = True
            axe.patch.set_facecolor('lightgrey')
            axe.patch.set_alpha(0.4)
            axe.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
            axe.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
            axe.set_axisbelow(True)
            axe.minorticks_on() # always on, not only for log
            # remove tick marks
            axe.tick_params(axis='both', direction='out', top=False,
                            right=False, left=False, bottom=False)
            axe.tick_params(axis='both', direction='out', top=False,
                            right=False, left=False, bottom=False,
                            which='minor')
        else:
            fig = axe.get_figure()
        # text
        plt.xlabel('Iteration number')
        plt.ylabel('IMP Objective Function Value')
        plt.title('Model ' + str(self['rand_init']))
        # smooth
        nrjz = self['log_objfun'][1:]
        if smooth:
            xnew = linspace(0, len(nrjz), 10000)
            nrjz_smooth = spline(range(len(nrjz)), nrjz, xnew,
                                 order=3)
            axe.plot(xnew, nrjz_smooth, color='darkred')
        else:
            axe.plot(nrjz, color='darkred')
        # plot
        axe.plot(nrjz, color='darkred', marker='o', alpha=.5, ms=4, ls='None')
        # log
        if log:
            axe.set_yscale('log')
        if savefig:
            tadbit_savefig(savefig)
        elif show:
            plt.show()
