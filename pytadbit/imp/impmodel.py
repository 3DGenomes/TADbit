"""
06 Aug 2013


"""

from matplotlib import pyplot as plt
from scipy.interpolate import spline
from numpy import linspace


class IMPmodel(dict):
    """
    A container for IMP modelling result.

    It is a dictionnary with this keys:
    
    - log_energies: a list of enrgies found by IMP objective function
    - energy: the last energy of previous list. The one associated to the
       model, and used for model comparison.
    - rand_init: The random initial value passed to IMP. An other
       modelization with this same number should result in the same model
    - x 
    - y   }  3D coordinates of each particles (each of these is a list)
    - z 
    
    """
    def __repr__(self):
        try:
            return ('IMP model of {} particles with: \n' +
                    ' - Final energy: {}\n' +
                    ' - random initial value: {}\n' +
                    ' - first coordinates:\n'+
                    '        X      Y      Z\n'+
                    '  {:>7}{:>7}{:>7}\n'+
                    '  {:>7}{:>7}{:>7}\n'+
                    '  {:>7}{:>7}{:>7}\n').format(
                len(self['x']), self['energy'], self['rand_init'],
                int(self['x'][0]), int(self['y'][0]), int(self['z'][0]),
                int(self['x'][1]), int(self['y'][1]), int(self['z'][1]),
                int(self['x'][2]), int(self['y'][2]), int(self['z'][2]))
        except IndexError:
            return ('IMP model of {} particles with: \n' +
                    ' - Final energy: {}\n' +
                    ' - random initial value: {}\n' +
                    ' - first coordinates:\n'+
                    '      X    Y    Z\n'+
                    '  {:>5}{:>5}{:>5}\n').format(
                len(self['x']), self['energy'], self['rand_init'],
                self['x'][0], self['y'][0], self['z'][0])


    def objective_function(self, log=False, smooth=True,
                           axe=None, savefig=None):
        """
        Plots the fall in energy through the Monte-Carlo search.
        
        :param False log: to plot in log scale
        :param True smooth: to smooth the curve
        
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
        plt.ylabel('Energy')
        plt.title('Objective function')
        # smooth
        nrjz = self['log_energies'][1:]
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
