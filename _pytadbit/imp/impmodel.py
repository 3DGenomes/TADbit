"""
06 Aug 2013


"""

from pytadbit.utils.extraviews      import color_residues, chimera_view
from pytadbit.utils.extraviews      import tadbit_savefig, plot_3d_model
from pytadbit.utils.three_dim_stats import generate_sphere_points
from pytadbit.utils.three_dim_stats import fast_square_distance
from pytadbit.utils.three_dim_stats import build_mesh
from pytadbit.utils.extraviews      import tad_coloring
from pytadbit.utils.extraviews      import tad_border_coloring
from scipy.interpolate              import spline
from numpy                          import linspace
from warnings                       import warn
from re                             import findall, compile as compil
from math                           import sqrt, pi

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')


def model_header(model):
    """
    Defines the header to write in output files for a given model
    """
    if not 'description' in model:
        return ''
    outstr = ''
    for desc in sorted(model['description']):
        outstr += '# %-15s :  %s\n' % (desc.upper(), model['description'][desc])
    return outstr
    

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
          # RESOLUTION      : 20000
          # ASSEMBLY        : None
          # CHROMOSOME      : Test Chromosome
          # START           : 50
          # END             : 70
          1  19:1-10000        44.847     412.828    -162.673
          2  19:10001-20000   -55.574     396.869    -129.782

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


def load_impmodel_from_xyz_OLD(f_name, rand_init=None, radius=None):
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
    for xxx, yyy, zzz in findall(expr, open(f_name).read()):
        model['x'].append(float(xxx))
        model['y'].append(float(yyy))
        model['z'].append(float(zzz))
    return model


class IMPmodel(dict):
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


    def __len__(self):
        return len(self['x'])


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


    def distance(self, part1, part2):
        """
        Calculates the distance between one point of the model and an external
        coordinate
        
        :param part1: index of a particle in the model
        :param part2: index of a particle in the model

        :returns: distance between one point of the model and an external
           coordinate
        """
        return sqrt((self['x'][part1-1] - self['x'][part2-1])**2 +
                    (self['y'][part1-1] - self['y'][part2-1])**2 +
                    (self['z'][part1-1] - self['z'][part2-1])**2)


    def _square_distance_to(self, part1, part2):
        """
        :param part1: index of a particle in the model
        :param part2: external coordinate (dict format with x, y, z keys)

        :returns: square distance between one point of the model and an external
           coordinate
        """
        return ((self['x'][part1] - part2[0])**2 +
                (self['y'][part1] - part2[1])**2 +
                (self['z'][part1] - part2[2])**2)


    def center_of_mass(self):
        """
        Gives the center of mass of a model
        
        :returns: the center of mass of a given model
        """
        r_x = sum(self['x'])/len(self)
        r_y = sum(self['y'])/len(self)
        r_z = sum(self['z'])/len(self)
        return dict((('x', r_x), ('y', r_y), ('z', r_z)))


    def radius_of_gyration(self):
        """
        Calculates the radius of gyration or gyradius of the model
        
        Defined as:
        
        .. math::

          \sqrt{\\frac{\sum_{i=1}^{N} (x_i-x_{com})^2+(y_i-y_{com})^2+(z_i-z_{com})^2}{N}}

        with:
        
        * :math:`N` the number of particles in the model
        * :math:`com` the center of mass
        
        :returns: the radius of gyration for the components of the tensor
        """
        com = self.center_of_mass()
        rog = sqrt(sum([self._square_distance_to(i,
                                                 (com['x'], com['y'], com['z']))
                        for i in xrange(len(self))]) / len(self))
        return rog


    def contour(self):
        """
        Total length of the model
        
        :returns: the totals length of the model
        """
        dist = 0
        for i in xrange(1, len(self)):
            dist += self.distance(i, i+1)
        return dist


    def longest_axe(self):
        """
        Gives the distance between most distant particles of the model
        
        :returns: the maximum distance between two particles in the model
        """
        maxdist = 0
        for i in xrange(1, len(self)):
            for j in xrange(i + 1, len(self) + 1):
                dist = self.distance(i, j)
                if dist > maxdist:
                    maxdist = dist
        return maxdist


    def shortest_axe(self):
        """
        Minimum distance between two particles in the model
        
        :returns: the minimum distance between two particles in the model
        """
        mindist = float('inf')
        for i in xrange(1, len(self)):
            for j in xrange(i + 1, len(self) + 1):
                dist = self.distance(i, j)
                if dist < mindist:
                    mindist = dist
        return mindist


    def min_max_by_axis(self):
        """
        Calculates the minimum and maximum coordinates of the model
        
        :returns: the minimum and maximum coordinates for each x, y and z axis
        """
        return ((min(self['x']), max(self['x'])),
                (min(self['y']), max(self['y'])),
                (min(self['z']), max(self['z'])))


    def cube_side(self):
        """
        Calculates the side of a cube containing the model.
        
        :returns: the diagonal length of the cube containing the model
        """
        return sqrt((min(self['x']) - max(self['x']))**2 +
                    (min(self['y']) - max(self['y']))**2 +
                    (min(self['z']) - max(self['z']))**2)


    def cube_volume(self):
        """
        Calculates the volume of a cube containing the model.

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

        TODO: remove this function

        """
        inaccessibles = []
        sphere = generate_sphere_points(100)
        for i in xrange(len(self)):
            impossibles = 0
            for xxx, yyy, zzz in sphere:
                thing = (xxx * radius + self['x'][i],
                         yyy * radius + self['y'][i],
                         zzz * radius + self['z'][i])
                # print form % (k+len(self), thing['x'], thing['y'], thing['z'],
                # 0, 0, 0, k+len(self)),
                for j in xrange(len(self)):
                    if i == j:
                        continue
                    # print self._square_distance_to(j, thing), radius
                    if self._square_distance_to(j, thing) < radius**2:
                        # print i, j
                        impossibles += 1
                        break
            if impossibles == 100:
                inaccessibles.append(i + 1)
        return inaccessibles


    def accessible_surface(self, radius, nump=100, superradius=200,
                           include_edges=True, view_mesh=False, savefig=None,
                           write_cmm_file=None, verbose=False,
                           chimera_bin='chimera'):
        """
        Calculates a mesh surface around the model (distance equal to input
        **radius**) and checks if each point of this mesh could be replaced by
        an object (i.e. a protein) of a given **radius**

        Outer part of the model can be excluded from the estimation of
        accessible surface, as the occupancy outside the model is unkown (see
        superradius option).

        :param radius: radius of the object we want to fit in the model.
        :param None write_cmm_file: path to file in which to write cmm with the
           colored meshed (red inaccessible points, green accessible points)
        :param 100 nump: number of points to draw around a given particle. This
           number also sets the number of points drawn around edges, as each
           point occupies a given surface (see maths below). *Note that this
           number is considerably lowered by the occupancy of edges, depending
           of the angle formed by the edges surrounding a given particle, only
           10% to 50% of the ``nump`` will be drawn in fact.*
        :param True include_edges: if False, edges will not be included in the
           calculation of the accessible surface, only particles. Note that
           statistics on particles (like last item returned) will not change,
           and computation time will be significantly decreased.
        :param False view_mesh: launches chimera to display the mesh around the
           model
        :param None savefig: path where to save chimera image
        :param 'chimera' chimera_bin: path to chimera binary to use
        :param False verbose: prints stats about the surface
        :param 200 superradius: radius of an object used to exclude outer
           surface of the model. Superradius must be higher than radius.

        This function will first define a mesh around the chromatin,
        representing all possible position of the center of the object we want
        to fit. This mesh will be at a distance of *radius* from the chromatin
        strand. All dots in the mesh represents an equal area (*a*), the whole
        surface of the chromatin strand being: :math:`A=n \\times a` (*n* being
        the total number of dots in the mesh).

        The mesh consists of spheres around particles of the model, and
        cylinders around edges joining particles (no overlap is allowed between
        sphere and cylinders or cylinder and cylinder when they are
        consecutive).
        
        If we want that all dots of the mesh representing the surface of the
        chromatin, corresponds to an equal area (:math:`a`)

        .. math::

          a = \\frac{4\pi r^2}{s} = \\frac{2\pi r N_{(d)}}{c}

        with:

        * :math:`r` radius of the object to fit (as the input parameter **radius**)
        * :math:`s` number of points in sphere
        * :math:`c` number of points in circle (as the input parameter **nump**)
        * :math:`N_{(d)}` number of circles in an edge of length :math:`d`

        According to this, when the distance between two particles is equal
        to :math:`2r` (:math:`N=2r`), we would have :math:`s=c`.

        As :

        .. math::

          2\pi r = \sqrt{4\pi r^2} \\times \sqrt{\pi}
        
        It is fair to state the number of dots represented along a circle as:

        .. math::

          c = \sqrt{s} \\times \sqrt{\pi}

        Thus the number of circles in an edge of length :math:`d` must be:

        .. math::

          N_{(d)}=\\frac{s}{\sqrt{s}\sqrt{\pi}}\\times\\frac{d}{2r}

        :returns: a list of *1-* the number of dots in the mesh that could be
           occupied by an object of the given radius *2-* the total number of
           dots in the mesh *3-* the estimated area of the mesh (in square
           micrometers) *4-* the area of the mesh of a virtually straight strand
           of chromatin defined as 
           :math:`contour\\times 2\pi r + 4\pi r^2` (also in
           micrometers) *5-* a list of number of (accessibles, inaccessible) for
           each particle (percentage burried can be infered afterwards by
           accessible/(accessible+inaccessible) )

        """

        points, dots, superdots, points2dots = build_mesh(
            self['x'], self['y'], self['z'], len(self), nump, radius,
            superradius, include_edges)

        
        # calculates the number of inaccessible peaces of surface
        if superradius:
            radius2 = (superradius - 4)**2
            outdot  = []
            for x2, y2, z2 in superdots:
                for j, (x1, y1, z1) in enumerate(points):
                    if fast_square_distance(x1, y1, z1, x2, y2, z2) < radius2:
                        outdot.append(False)
                        break
                else:
                    outdot.append(True)
                    continue
                points.insert(0, points.pop(j))
        else:
            outdot = [False] * len(superdots)

        # calculates the number of inaccessible peaces of surface
        radius2 = (radius - 2)**2
        grey    = (0.6, 0.6, 0.6)
        red     = (1, 0, 0)
        green   = (0, 1, 0)
        colors  = []
        for i, (x2, y2, z2) in enumerate(dots):
            if outdot[i]:
                colors.append(grey)
                continue
            for j, (x1, y1, z1) in enumerate(points):
                if fast_square_distance(x1, y1, z1, x2, y2, z2) < radius2:
                    colors.append(red)
                    break
            else:
                colors.append(green)
                continue
            points.insert(0, points.pop(j))
        possibles = colors.count(green)

        acc_parts = []
        for p in sorted(points2dots.keys()):
            acc = 0
            ina = 0
            for dot in points2dots[p]:
                if colors[dot]==green:
                    acc += 1
                elif colors[dot]==red:
                    ina += 1
            acc_parts.append((p + 1, acc, ina))

        # some stats
        dot_area = 4 * pi * (float(radius) / 1000)**2 / nump
        area = (possibles * dot_area)
        total = (self.contour() / 1000 * 2 * pi * float(radius) / 1000 + 4 * pi
                 * (float(radius) / 1000)**2)
        if verbose:
            print (' Accessible surface: %s micrometers^2' +
                   '(%s accessible times %s micrometers)') % (
                round(area, 2), possibles, dot_area)
            print '    (%s accessible dots of %s total times %s micrometers)' % (
                possibles, outdot.count(False), round(dot_area, 5))
            print '  - %s%% of the contour mesh' % (
                round((float(possibles)/outdot.count(False))*100, 2))
            print '  - %s%% of a virtual straight chromatin (%s microm^2)' % (
                round((area/total)*100, 2), round(total, 2))

        # write cmm file
        if savefig:
            view_mesh = True
        if write_cmm_file or view_mesh:
            out = '<marker_set name=\"2\">\n'
            form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
                    ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
                    'radius=\"7\"/>\n')
            for k_2, thing in enumerate(dots):
                out += form % (1 + k_2, thing[0], thing[1], thing[2],
                                colors[k_2][0], colors[k_2][1], colors[k_2][2])
            if superradius:
                for k_3, thing in enumerate(superdots):
                    out += form % (1 + k_3 + k_2 + 1,
                                   thing[0], thing[1], thing[2],
                                    0.1, 0.1, 0.1)
            out += '</marker_set>\n'
            if view_mesh:
                out_f = open('/tmp/tmp_mesh.cmm', 'w')
                out_f.write(out)
                out_f.close()
            if write_cmm_file:
                out_f = open(write_cmm_file, 'w')
                out_f.write(out)
                out_f.close()
        if view_mesh:
            chimera_cmd = [
                'focus', 'set bg_color white', 'windowsize 800 600',
                'clip yon -500', 'set subdivision 1', 'set depth_cue',
                'set dc_color black', 'set dc_start 0.5', 'set dc_end 1',
                'scale 0.8']
            if savefig:
                if savefig.endswith('.png'):
                    chimera_cmd += ['copy file %s png' % (savefig)]
                elif savefig[-4:] in ('.mov', 'webm'):
                    chimera_cmd += [
                        'movie record supersample 1', 'turn y 3 120',
                        'wait 120', 'movie stop',
                        'movie encode output %s' % savefig]
            self.write_cmm('/tmp/')
            chimera_view(['/tmp/tmp_mesh.cmm',
                          '/tmp/model.%s.cmm' % (self['rand_init'])],
                         chimera_bin=chimera_bin,
                         savefig=savefig, chimera_cmd=chimera_cmd)

        return (possibles, outdot.count(False), area, total, acc_parts)


    def write_cmm(self, directory, color='index', rndname=True,
                  model_num=None, **kwargs):
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
        :param 'index' color: can be:

             * a string as:
                 * '**index**' to color particles according to their position in the
                   model (:func:`pytadbit.utils.extraviews.color_residues`)
                 * '**tad**' to color particles according to the TAD they belong to
                   (:func:`pytadbit.utils.extraviews.tad_coloring`)
                 * '**border**' to color particles marking borders. Color according to
                   their score (:func:`pytadbit.utils.extraviews.tad_border_coloring`)
                   coloring function like.
             * a function, that takes as argument a model and any other parameter
               passed through the kwargs.
             * a list of (r, g, b) tuples (as long as the number of particles).
               Each r, g, b between 0 and 1.
        :param kwargs: any extra argument will be passed to the coloring
           function
        """
        if type(color) is str:
            if color == 'index':
                color = color_residues(self, **kwargs)
            elif color == 'tad':
                if not 'tads' in kwargs:
                    raise Exception('ERROR: missing TADs\n   ' +
                                    'pass an Experiment.tads disctionary\n')
                color = tad_coloring(self, **kwargs)
            elif color == 'border':
                if not 'tads' in kwargs:
                    raise Exception('ERROR: missing TADs\n   ' +
                                    'pass an Experiment.tads disctionary\n')
                color = tad_border_coloring(self, **kwargs)
            else:
                raise NotImplementedError(('%s type of coloring is not yet ' +
                                           'implemeted\n') % color)
        elif hasattr(color, '__call__'): # its a function
            color = color(self, **kwargs)
        elif type(color) is not list:
            raise TypeError('one of function, list or string is required\n')
        out = '<marker_set name=\"%s\">\n' % (self['rand_init'])
        form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
                ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
                'radius=\"' + str(30) +
                # str(self['radius']) +
                '\" note=\"%s\"/>\n')
        for i in xrange(len(self['x'])):
            out += form % (i + 1,
                           self['x'][i], self['y'][i], self['z'][i],
                           color[i][0], color[i][1], color[i][2], i + 1)
        form = ('<link id1=\"%s\" id2=\"%s\" r=\"1\" ' +
                'g=\"1\" b=\"1\" radius=\"' + str(10) +
                # str(self['radius']/2) +
                '\"/>\n')
        for i in xrange(1, len(self['x'])):
            out += form % (i, i + 1)
        out += '</marker_set>\n'

        if rndname:
            out_f = open('%s/model.%s.cmm' % (directory,
                                              self['rand_init']), 'w')
        else:
            out_f = open('%s/model_%s_rnd%s.cmm' % (
                directory, model_num, self['rand_init']), 'w')
        out_f.write(out)
        out_f.close()


    def write_xyz_OLD(self, directory, model_num=None, get_path=False,
                  rndname=True):
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
        for i in xrange(len(self['x'])):
            out += form % ('p' + str(i + 1), i + 1, round(self['x'][i], 3),
                           round(self['y'][i], 3), round(self['z'][i], 3))
        out_f = open(path_f, 'w')
        out_f.write(out)
        out_f.close()
        if get_path:
            return path_f
        else:
            return None


    def write_xyz(self, directory, model_num=None, get_path=False,
                  rndname=True, header=True):
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
        :param True header: write a header describing the experiment from which
           the model was calculated.
        """
        if rndname:
            path_f = '%s/model.%s.xyz' % (directory, self['rand_init'])
        else:
            path_f = '%s/model_%s_rnd%s.xyz' % (directory, model_num,
                                                self['rand_init'])
        out = ''
        if header:
            out += model_header(self)
        form = "%4s %21s%12.3f%12.3f%12.3f\n"
        # TODO: do not use resolution directly -> specific to Hi-C
        for i in xrange(len(self['x'])):
            out += form % (
                i + 1,
                '%s:%s-%s' % (
                    self['description']['chromosome'],
                    int(self['description']['start']) + int(self['description']['resolution']) * i + 1,
                    int(self['description']['start']) + int(self['description']['resolution']) * (i + 1)),
                round(self['x'][i], 3),
                round(self['y'][i], 3), round(self['z'][i], 3))
        out_f = open(path_f, 'w')
        out_f.write(out)
        out_f.close()
        if get_path:
            return path_f
        else:
            return None


    def view_model(self, tool='chimera', savefig=None, cmd=None,
                   center_of_mass=False, gyradius=False, color='index',
                   **kwargs):
        """
        Visualize a selected model in the three dimensions. (either with Chimera
        or through matplotlib).

        :param model_num: model to visualize
        :param 'chimera' tool: path to the external tool used to visualize the
           model. Can also be 'plot', to use matplotlib.
        :param None savefig: path to a file where to save the image OR movie
           generated (depending on the extension; accepted formats are png, mov
           and webm). if set to None, the image or movie will be shown using
           the default GUI.
        :param 'index' color: can be:

             * a string as:
                 * '**index**' to color particles according to their position in the
                   model (:func:`pytadbit.utils.extraviews.color_residues`)
                 * '**tad**' to color particles according to the TAD they belong to
                   (:func:`pytadbit.utils.extraviews.tad_coloring`)
                 * '**border**' to color particles marking borders. Color according to
                   their score (:func:`pytadbit.utils.extraviews.tad_border_coloring`)
                   coloring function like.
             * a function, that takes as argument a model and any other parameter
               passed through the kwargs.
             * a list of (r, g, b) tuples (as long as the number of particles).
               Each r, g, b between 0 and 1.
        :param False center_of_mass: draws a dot representing the center of
           mass of the model
        :param False gyradius: draws the center of mass of the model as a sphere
           with radius equal to the radius of gyration of the model
        :param None cmd: list of commands to be passed to the viewer.
           The chimera list is:

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

             cmd = ['focus', 'set bg_color white', 'windowsize 800 600',
                    'bonddisplay never #0',
                    'shape tube #0 radius 10 bandLength 200 segmentSubdivisions 100 followBonds on',
                    'clip yon -500', '~label', 'set subdivision 1',
                    'set depth_cue', 'set dc_color black', 'set dc_start 0.5',
                    'set dc_end 1', 'scale 0.8']

           will return the default image (other commands can be passed to
           modified the final image/movie).

        :param kwargs: see :func:`pytadbit.utils.extraviews.plot_3d_model` or
           :func:`pytadbit.utils.extraviews.chimera_view` for other arguments
           to pass to this function

        """
        if gyradius:
            gyradius = self.radius_of_gyration()
            center_of_mass = True
        if tool=='plot':
            plot_3d_model(self, color=color, **kwargs)
            return
        self.write_cmm('/tmp/', color=color, **kwargs)
        chimera_view(['/tmp/model.%s.cmm' % (self['rand_init'])],
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd,
                     center_of_mass=center_of_mass, gyradius=gyradius)
