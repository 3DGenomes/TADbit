"""
06 Aug 2013


"""

from scipy.interpolate import spline
from numpy             import linspace
from warnings          import warn

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')


class IMPmodel(dict):
    """
    A container for the IMP modeling results. The container is a dictionary 
    with the following keys:
    
    - log_objfun: The list of IMP objective function values
    - objfun: The final objective function value of the corresponding model 
       (from log_objfun). This value will be used to rank all the generated 
       models
    - rand_init: The random number generator feed (needed for model 
       reproducibility)
    - x, y, z: The spatial 3D coordinates of each particles. Each coordinate is 
       represented as a list
    
    """
    def __repr__(self):
        try:
            return ('IMP model of %s particles with: \n' +
                    ' - Final objective function value: %s\n' +
                    ' - random initial value: %s\n' +
                    ' - first coordinates:\n'+
                    '        X      Y      Z\n'+
                    '  %7s%7s%7s\n'+
                    '  %7s%7s%7s\n'+
                    '  %7s%7s%7s\n') % (
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
            fig.savefig(savefig)
        elif show:
            plt.show()
