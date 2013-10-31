"""
06 Aug 2013


"""

from pytadbit.utils.extraviews      import color_residues, chimera_view
from pytadbit.utils.three_dim_stats import generate_sphere_points
from pytadbit.utils.three_dim_stats import square_distance
from pytadbit.utils.three_dim_stats import generate_circle_points
from scipy.interpolate              import spline
from numpy                          import linspace
from warnings                       import warn
from re                             import findall, compile as compil
from math                           import sqrt
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
            rand_init = int(f_name.split('.')[-2])
        except:
            rand_init = None
    model = IMPmodel((('x', []), ('y', []), ('z', []), ('rand_init', rand_init),
                      ('objfun', None), ('radius', radius)))
    expr = compil(
        ' x="([0-9.-]+)" y="([0-9.-]+)" z="([0-9.-]+)".* radius="([0-9.]+)"')
    for x, y, z, radius in findall(expr, open(f_name).read()):
        model['x'].append(float(x))
        model['y'].append(float(y))
        model['z'].append(float(z))
    if not model['radius']:
        model['radius'] = float(radius)
    return model


def load_impmodel_from_xyz(f_name, rand_init=None, radius=None):
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
            rand_init = int(f_name.split('.')[-2])
        except:
            rand_init = None
    model = IMPmodel((('x', []), ('y', []), ('z', []), ('rand_init', rand_init),
                      ('objfun', None), ('radius', radius)))
    expr = compil('p[0-9]+\s+[0-9]+\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)')
    for x, y, z in findall(expr, open(f_name).read()):
        model['x'].append(float(x))
        model['y'].append(float(y))
        model['z'].append(float(z))
    return model


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


    def __len__(self):
        return len(self['x'])


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


    def distance(self, part1, part2):
        """
        """
        return sqrt((self['x'][part1] - self['x'][part2])**2 +
                    (self['y'][part1] - self['y'][part2])**2 +
                    (self['z'][part1] - self['z'][part2])**2)


    def square_distance_to(self, part1, part2):
        """
        """
        return ((self['x'][part1] - part2['x'])**2 +
                (self['y'][part1] - part2['y'])**2 +
                (self['z'][part1] - part2['z'])**2)


    def center_of_mass(self):
        """
        :returns: the center of mass of a given model
        """
        r_x = sum(self['x'])/len(self)
        r_y = sum(self['y'])/len(self)
        r_z = sum(self['z'])/len(self)
        return dict((('x', r_x), ('y', r_y), ('z', r_z)))


    def radius_of_gyration(self):
        """
        :returns: the radius of gyration for the components of the tensor
        """
        com = self.center_of_mass()
        rog = sqrt(sum([self.square_distance_to(i, com)
                        for i in xrange(len(self))]) / len(self))
        return rog


    def contour(self):
        """
        :returns: the totals length of the model
        """
        dist = 0
        for i in xrange(len(self)-1):
            dist += self.distance(i, i+1)
        return dist


    def longest_axe(self):
        """
        :returns: the maximum distance between two particles in the model
        """
        maxdist = 0
        for i in xrange(len(self) - 1):
            for j in xrange(i + 1, len(self)):
                dist = self.distance(i, j)
                if dist > maxdist:
                    maxdist = dist
        return maxdist


    def shortest_axe(self):
        """
        :returns: the minimum distance between two particles in the model
        """
        mindist = float('inf')
        for i in xrange(len(self) - 1):
            for j in xrange(i + 1, len(self)):
                dist = self.distance(i, j)
                if dist < mindist:
                    mindist = dist
        return mindist


    def min_max_by_axis(self):
        """
        :returns: the minimum and maximum coordinates for each x, y and z axis
        """
        return ((min(self['x']), max(self['x'])),
                (min(self['y']), max(self['y'])),
                (min(self['z']), max(self['z']))) 


    def cube_side(self):
        """
        :returns: the diagonal length of the cube containing the model
        """
        return sqrt((min(self['x']) - max(self['x']))**2 +
                    (min(self['y']) - max(self['y']))**2 +
                    (min(self['z']) - max(self['z']))**2)


    def cube_volume(self):
        """
        :returns: the volume of  the cube containing the model
        """
        return self.cube_side()**3


    def inaccessible_particles(self, radius):
        """
        Gives the number of loci/particles that are accessible to an object
        (i.e. a protein) of a given size.
        
        :param radius: radius of the object that we want to fit in the model

        :returns: a list of numbers, each being the ID of a particles that would
           never be reached by the given object
        
        """
        inaccessibles = []
        sphere = generate_sphere_points(100)
        for i in xrange(len(self)):
            impossibles = 0
            for x, y, z in sphere:
                thing = dict((('x', x * radius + self['x'][i]),
                              ('y', y * radius + self['y'][i]),
                              ('z', z * radius + self['z'][i])))
                # print form % (k+len(self), thing['x'], thing['y'], thing['z'], 0, 0, 0, k+len(self)),
                for j in xrange(len(self)):
                    if i == j: continue
                    # print self.square_distance_to(j, thing), radius
                    if self.square_distance_to(j, thing) < radius**2:
                        # print i, j
                        impossibles += 1
                        break
            if impossibles == 100:
                inaccessibles.append(i + 1)
        return inaccessibles


    def accessible_surface(self, radius):
        """
        """

        radius = 300
        between = 10 # number of cuts in edges
        points = []
        subpoints = []
        sphere = generate_sphere_points(100)
        for i in xrange(len(self)-1):
            points.append(dict((('x', self['x'][i]),
                                ('y', self['y'][i]),
                                ('z', self['z'][i]))))
            for x, y, z in sphere:
                thing = dict((('x', x * radius + self['x'][i]),
                              ('y', y * radius + self['y'][i]),
                              ('z', z * radius + self['z'][i])))
                subpoints.append(thing)
            difx = self['x'][i] - self['x'][i+1]
            dify = self['y'][i] - self['y'][i+1]
            difz = self['z'][i] - self['z'][i+1]
            normer = sqrt(difx**2 + dify**2 + difz**2)
            ortho = dict((('x', 1.), ('y', 1.), ('z', -(difx + dify) / difz)))
            normer = sqrt(ortho['x']**2 + ortho['y']**2 + ortho['z']**2)
            ortho['x'] = ortho['x'] / normer * radius
            ortho['y'] = ortho['y'] / normer * radius
            ortho['z'] = ortho['z'] / normer * radius
            for k in xrange(between-1, 0, -1):
                point = dict((('x', self['x'][i] - k * (difx / between)),
                              ('y', self['y'][i] - k * (dify / between)),
                              ('z', self['z'][i] - k * (difz / between))))
                points.append(point)
                for spoint in generate_circle_points(ortho['x'] + point['x'],
                                                     ortho['y'] + point['y'],
                                                     ortho['z'] + point['z'],
                                                     point['x']             ,
                                                     point['y']             ,
                                                     point['z']             ,
                                                     difx                  ,
                                                     dify                  ,
                                                     difz                  ,
                                                     20):
                    subpoints.append(dict((('x', spoint[0]),
                                           ('y', spoint[1]),
                                           ('z', spoint[2]))))
                    
            points.append(dict((('x', self['x'][i + 1]),
                                ('y', self['y'][i + 1]),
                                ('z', self['z'][i + 1]))))
        for x, y, z in sphere:
            thing = dict((('x', x * radius + self['x'][i+1]),
                          ('y', y * radius + self['y'][i+1]),
                          ('z', z * radius + self['z'][i+1])))
            subpoints.append(thing)

        impossibles = 0
        colors = []
        for spoint in subpoints:
            for point in points:
                if square_distance(point, spoint) < (radius-50)**2:
                    impossibles += 1
                    colors.append((100, 0, 0))
                    break
            else:
                colors.append((0, 100, 0))

        print '<marker_set name=\"1\">'
        for k, thing in enumerate(points):
            form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
                    ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
                    'radius=\"20\" note=\"%s\"/>\n')
            print form % (k+1+len(self), thing['x'], thing['y'], thing['z'], 0, 0, 0, k+1+len(self)),
        for k2, thing in enumerate(subpoints):
            form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
                    ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
                    'radius=\"10\" note=\"%s\"/>\n')
            print form % (2+k+k2+len(self), thing['x'], thing['y'], thing['z'], colors[k2][0], colors[k2][1], colors[k2][2], 2+k+k2+len(self)),
        print '</marker_set>'


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


    def view_model(self, tool='chimera', savefig=None, cmd=None, centroid=False,
                   gyradius=False):
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
        if gyradius:
            gyradius = self.radius_of_gyration()
            centroid=True
        self.write_cmm('/tmp/')
        chimera_view(['/tmp/model.%s.cmm' % (self['rand_init'])],
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd,
                     centroid=centroid, gyradius=gyradius)



