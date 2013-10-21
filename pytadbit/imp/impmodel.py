"""
06 Aug 2013


"""

from pytadbit.utils.extraviews import color_residues, chimera_view
from scipy.interpolate         import spline
from numpy                     import linspace
from warnings                  import warn

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
    def __str__(self):
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


    def write_cmm(self, directory, color=color_residues, rndname=True,
                  model_num=None):
        """
        Save a model in the cmm format, read by Chimera
        (http://www.cgl.ucsf.edu/chimera).

        **Note:** If none of model_num, models or cluster parameter are set,
        ALL the models will be written.

        :param directory: location where the file will be written (note: the
           name of the file will be model_1.cmm if model number is 1)
        :param None model_num: the number of the model to save
        :param True rndname: If True, file names will be formatted as:
           model.RND.cmm, where RND is the random number feed used by IMP to
           generate the corresponding model. If False, the format will be:
           model_NUM_RND.cmm where NUM is the rank of the model in terms of
           objective function value
        :param color_residues color: either a coloring function like
           :func:`pytadbit.imp.imp_model.color_residues` or a list of (r, g, b)
           tuples (as long as the number of particles)
        """
        if type(color) != list:
            color = color(len(self['x']))
        out = '<marker_set name=\"%s\">\n' % (self['rand_init'])
        form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
                ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
                'radius=\"' +
                str(self['radius']) +
                '\" note=\"%s\"/>\n')
        for n in xrange(len(self['x'])):
            out += form % (n + 1,
                           self['x'][n], self['y'][n], self['z'][n],
                           color[n][0], color[n][1], color[n][2], n + 1)
        form = ('<link id1=\"%s\" id2=\"%s\" r=\"1\" ' +
                'g=\"1\" b=\"1\" radius=\"' +
                str(self['radius']/2) +
                '\"/>\n')
        for n in xrange(1, len(self['x'])):
            out += form % (n, n + 1)
        out += '</marker_set>\n'

        if rndname:
            out_f = open('%s/model.%s.cmm' % (directory,
                                              self['rand_init']), 'w')
        else:
            out_f = open('%s/model_%s_rnd%s.cmm' % (
                directory, model_num, self['rand_init']), 'w')
        out_f.write(out)
        out_f.close()
        

    def write_xyz(self, directory, model_num=None, get_path=False, rndname=True):
        """
        Writes a xyz file containing the 3D coordinates of each particle in the
        model.

        **Note:** If none of model_num, models or cluster parameter are set,
        ALL the models will be written.

        :param directory: location where the file will be written (note: the
           file name will be model.1.xyz, if the model number is 1)
        :param None model_num: the number of the model to save
        :param True rndname: If True, file names will be formatted as:
           model.RND.xyz, where RND is the random number feed used by IMP to
           generate the corresponding model. If False, the format will be:
           model_NUM_RND.xyz where NUM is the rank of the model in terms of
           objective function value
        :param False get_path: whether to return, or not, the full path where
           the file has been written
        """
        if rndname:
            path_f = '%s/model.%s.xyz' % (directory, self['rand_init'])
        else:
            path_f = '%s/model_%s_rnd%s.xyz' % (directory, model_num,
                                                self['rand_init'])
        out = ''
        form = "%12s%12s%12.3f%12.3f%12.3f\n"
        for n in xrange(len(self['x'])):
            out += form % ('p' + str(n + 1), n + 1, round(self['x'][n], 3),
                           round(self['y'][n], 3), round(self['z'][n], 3))
        out_f = open(path_f, 'w')
        out_f.write(out)
        out_f.close()
        if get_path:
            return path_f
        else:
            return None


    def view_model(self, tool='chimera', savefig=None, cmd=None):
        """
        Visualize a selected model in the three dimensions.

        :param model_num: model to visualize
        :param 'chimera' tool: path to the external tool used to visualize the
           model
        :param None savefig: path to a file where to save the image OR movie
           generated (depending on the extension; accepted formats are png, mov
           and webm). if set to None, the image or movie will be shown using
           the default GUI.
        :param None cmd: list of commands to be passed to the viewer. The chimera list is:

           ::

             focus
             set bg_color white
             windowsize 800 600
             bonddisplay never #0
             represent wire
             shape tube #0 radius 5 bandLength 100 segmentSubdivisions 1 followBonds on
             clip yon -500
             ~label
             set subdivision 1
             set depth_cue
             set dc_color black
             set dc_start 0.5
             set dc_end 1
             scale 0.8

           Followed by the movie command to record movies:

           ::

             movie record supersample 1
             turn y 3 120
             wait 120
             movie stop
             movie encode output SAVEFIG

           Or the copy command for images:

           ::

             copy file SAVEFIG png

           Passing as the following list as 'cmd' parameter:
           ::

             cmd = ['focus', 'set bg_color white', 'windowsize 800 600', 'bonddisplay never #0', 'shape tube #0 radius 10 bandLength 200 segmentSubdivisions 100 followBonds on', 'clip yon -500', '~label', 'set subdivision 1', 'set depth_cue', 'set dc_color black', 'set dc_start 0.5', 'set dc_end 1', 'scale 0.8']

           will return the default image (other commands can be passed to
           modified the final image/movie).

        """
        self.write_cmm('/tmp/')
        chimera_view('/tmp/model.%s.cmm' % (self['rand_init']),
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd)



