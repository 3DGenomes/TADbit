"""
25 Oct 2016


"""
from __future__ import print_function

from pytadbit.utils.extraviews      import color_residues, chimera_view
from pytadbit.utils.extraviews      import plot_3d_model
from pytadbit.utils.three_dim_stats import generate_sphere_points
from pytadbit.utils.three_dim_stats import fast_square_distance
from pytadbit.utils.three_dim_stats import build_mesh
from pytadbit.utils.extraviews      import tad_coloring
from pytadbit.utils.extraviews      import tad_border_coloring
from pytadbit.utils.tadmaths        import newton_raphson
from pytadbit                       import __version__ as version
from math                           import sqrt, pi
import hashlib

try:
    basestring
except NameError:
    basestring = str

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


class StructuralModel(dict):
    """
    A container for the structural modelling results. The container is a dictionary
    with the following keys:

    - rand_init: Random number generator feed (needed for model reproducibility)
    - x, y, z: 3D coordinates of each particles. Each represented as a list

    """
    def __len__(self):
        return len(self['x'])

    def distance(self, part1, part2):
        """
        Calculates the distance between one point of the model and an external
        coordinate

        :param part1: index of a particle in the model
        :param part2: index of a particle in the model

        :returns: distance between one point of the model and an external
           coordinate
        """
        if part1 == 0 or part2 == 0:
            raise Exception('Particle number must be strictly positive\n')
        return sqrt((self['x'][part1-1] - self['x'][part2-1])**2 +
                    (self['y'][part1-1] - self['y'][part2-1])**2 +
                    (self['z'][part1-1] - self['z'][part2-1])**2)


    def _square_distance(self, part1, part2):
        """
        Calculates the square instance between one point of the model and an
        external coordinate

        :param part1: index of a particle in the model
        :param part2: index of a particle in the model

        :returns: distance between one point of the model and an external
           coordinate
        """
        return ((self['x'][part1-1] - self['x'][part2-1])**2 +
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
        rog = sqrt(sum(self._square_distance_to(i,
                                                (com['x'], com['y'], com['z']))
                       for i in range(len(self))) / len(self))
        return rog


    def contour(self):
        """
        Total length of the model

        :returns: the totals length of the model
        """
        dist = 0
        for i in range(1, len(self)):
            dist += self.distance(i, i+1)
        return dist


    def longest_axe(self):
        """
        Gives the distance between most distant particles of the model

        :returns: the maximum distance between two particles in the model
        """
        maxdist = 0
        for i in range(1, len(self)):
            for j in range(i + 1, len(self) + 1):
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
        for i in range(1, len(self)):
            for j in range(i + 1, len(self) + 1):
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
        for i in range(len(self)):
            impossibles = 0
            for xxx, yyy, zzz in sphere:
                thing = (xxx * radius + self['x'][i],
                         yyy * radius + self['y'][i],
                         zzz * radius + self['z'][i])
                # print form % (k+len(self), thing['x'], thing['y'], thing['z'],
                # 0, 0, 0, k+len(self)),
                for j in range(len(self)):
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


    def persistence_length(self, start=1, end=None, return_guess=False):
        """
        Calculates the persistence length (Lp) of given section of the model.
        Persistence length is calculated according to [Bystricky2004]_ :

        .. math::

          <r^2> = 2 \\times Lp^2 \\times (\\frac{Lc}{Lp} - 1 + e^{\\frac{-Lc}{Lp}})

        with the contour length as :math:`Lc = \\frac{d}{c}` where :math:`d` is
        the genomic dstance in bp and :math:`c` the linear mass density of the
        chromatin (in bp/nm).

        :param 1 start: start particle from which to calculate the persistence
           length
        :param None end: end particle from which to calculate the persistence
           length. Uses the last particle by default
        :param False return_guess: Computes persistence length using the
           approximation :math:`Lp=\\frac{Lc}{Lp}`

        :returns: persistence length, or 2 times the Kuhn length
        """
        clength = float(self.contour())
        end = end or len(self)
        sq_length = float(self._square_distance(start, end))

        guess = sq_length / clength
        if return_guess:
            return guess # incredible!
        kuhn = newton_raphson(guess, clength, sq_length)
        return 2 * kuhn


    def accessible_surface(self, radius, nump=100, superradius=200,
                           include_edges=True, view_mesh=False, savefig=None,
                           write_cmm_file=None, verbose=False,
                           chimera_bin='chimera'):
        """
        Calculates a mesh surface around the model (distance equal to input
        **radius**) and checks if each point of this mesh could be replaced by
        an object (i.e. a protein) of a given **radius**

        Outer part of the model can be excluded from the estimation of
        accessible surface, as the occupancy outside the model is unknown (see
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
            print((' Accessible surface: %s micrometers^2' +
                   '(%s accessible times %s micrometers)') % (
                round(area, 2), possibles, dot_area))
            print('    (%s accessible dots of %s total times %s micrometers)' % (
                possibles, outdot.count(False), round(dot_area, 5)))
            print('  - %s%% of the contour mesh' % (
                round((float(possibles)/outdot.count(False))*100, 2)))
            print('  - %s%% of a virtual straight chromatin (%s microm^2)' % (
                round((area/total)*100, 2), round(total, 2)))

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
                'focus',
                'bonddisplay never #1',
                'shape tube #1 radius 15 bandLength 300 segmentSubdivisions 1 followBonds on',
                '~show #1',
                'set bg_color white', 'windowsize 800 600',
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
                         chimera_bin=chimera_bin, align=False,
                         savefig=savefig, chimera_cmd=chimera_cmd)

        return (possibles, outdot.count(False), area, total, acc_parts)


    def write_cmm(self, directory, color='index', rndname=True,
                  model_num=None, filename=None, **kwargs):
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
        :param None filename: overide the default file name writing
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
        if isinstance(color, basestring):
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
        elif hasattr(color, '__call__'): # it's a function
            color = color(self, **kwargs)
        elif not isinstance(color, list):
            raise TypeError('one of function, list or string is required\n')
        out = '<marker_set name=\"%s\">\n' % (self['rand_init'])
        form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
                ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
                'radius=\"' + #str(30) +
                str(self['radius']) +
                '\" note=\"%s\"/>\n')
        for i in range(len(self['x'])):
            out += form % (i + 1,
                           self['x'][i], self['y'][i], self['z'][i],
                           color[i][0], color[i][1], color[i][2], i + 1)
        form = ('<link id1=\"%s\" id2=\"%s\" r=\"1\" ' +
                'g=\"1\" b=\"1\" radius=\"' + str(self['radius']/10.) +
                # str(self['radius']/2) +
                '\"/>\n')
        break_chroms = [1]
        try:
            for beg, end in zip(self['description']['start'],self['description']['end']):
                break_chroms.append((end - beg)/self['description']['resolution']+break_chroms[-1])
        except:
            pass
        for i in range(1, len(self['x'])):
            if i in break_chroms[1:]:
                continue
            out += form % (i, i + 1)
        out += '</marker_set>\n'

        if filename:
                out_f = open('%s/%s' % (directory, filename), 'w')
        else:
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
        for i in range(len(self['x'])):
            out += form % ('p' + str(i + 1), i + 1, round(self['x'][i], 3),
                           round(self['y'][i], 3), round(self['z'][i], 3))
        out_f = open(path_f, 'w')
        out_f.write(out)
        out_f.close()
        if get_path:
            return path_f
        else:
            return None



    def write_json(self, directory, color='index', rndname=True,
                   model_num=None, title=None, filename=None, **kwargs):
        """
        Save a model in the json format, read by TADkit.

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
        :param None filename: overide the default file name writing
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
        if isinstance(color, basestring):
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
        elif hasattr(color, '__call__'): # it's a function
            color = color(self, **kwargs)
        elif not isinstance(color, list):
            raise TypeError('one of function, list or string is required\n')
        form = '''
{
    "chromatin" : {
        "id" : %(sha)s,
        "title" : "%(title)s",
        "source" : "TADbit %(version)s",
        "metadata": {%(descr)s},
        "type" : "tadbit",
            "data": {
            "models":     [
                { %(xyz)s },
            ],
                "clusters":[%(cluster)s],
                      "centroid":[%(centroid)s],
                "restraints": [[][]],
            "chromatinColor" : [ ]
        }
    }
}
'''
        fil = {}
        fil['sha']     = hashlib.new(fil['xyz']).hexdigest()
        fil['title']   = title or "Sample TADbit data"
        fil['version'] = version
        fil['descr']   = ''.join('\n', ',\n'.join([
            '"%s" : %s' % (k, ('"%s"' % (v)) if isinstance(v, basestring) else v)
            for k, v in list(self.get('description', {}).items())]), '\n')
        fil['xyz']     = ','.join(['[%.4f,%.4f,%.4f]' % (self['x'][i], self['y'][i],
                                                         self['z'][i])
                                   for i in range(len(self['x']))])


        if filename:
                out_f = open('%s/%s' % (directory, filename), 'w')
        else:
            if rndname:
                out_f = open('%s/model.%s.cmm' % (directory,
                                                  self['rand_init']), 'w')
            else:
                out_f = open('%s/model_%s_rnd%s.cmm' % (
                    directory, model_num, self['rand_init']), 'w')
        out_f.write(out)
        out_f.close()



    def write_xyz(self, directory, model_num=None, get_path=False,
                  rndname=True, filename=None, header=True):
        """
        Writes a xyz file containing the 3D coordinates of each particle in the
        model.
        Outfile is tab separated column with the bead number being the
        first column, then the genomic coordinate and finally the 3
        coordinates X, Y and Z

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
        :param None filename: overide the default file name writing
        :param True header: write a header describing the experiment from which
        """
        if filename:
            path_f = '%s/%s' % (directory, filename)
        else:
            if rndname:
                path_f = '%s/model.%s.xyz' % (directory, self['rand_init'])
            else:
                path_f = '%s/model_%s_rnd%s.xyz' % (directory, model_num,
                                                    self['rand_init'])
        out = ''
        if header:
            out += model_header(self)
        form = "%s\t%s\t%.3f\t%.3f\t%.3f\n"
        # TODO: do not use resolution directly -> specific to Hi-C
        chrom_list = self['description']['chromosome']
        chrom_start = self['description']['start']
        chrom_end = self['description']['end']
        if not isinstance(chrom_list, list):
            chrom_list = [chrom_list]
            chrom_start = [chrom_start]
            chrom_end = [chrom_end]

        chrom_start = [(int(c) // int(self['description']['resolution'])
                        if int(c) else 0)
                       for c in chrom_start]
        chrom_end = [(int(c) // int(self['description']['resolution'])
                      if int(c) else len(self['x']))
                     for c in chrom_end]

        offset = -chrom_start[0]
        for crm in range(len(chrom_list)):
            for i in range(chrom_start[crm] + offset, chrom_end[crm] + offset):
                out += form % (
                    i + 1,
                    '%s:%s-%s' % (
                        chrom_list[crm],
                        int(chrom_start[crm] or 0) +
                            int(self['description']['resolution']) * (i - offset) + 1,
                        int(chrom_start[crm] or 0) +
                            int(self['description']['resolution']) * (i + 1 - offset)),
                    round(self['x'][i], 3),
                    round(self['y'][i], 3), round(self['z'][i], 3))
            offset += (chrom_end[crm] - chrom_start[crm])
        out_f = open(path_f, 'w')
        out_f.write(out)
        out_f.close()
        if get_path:
            return path_f
        else:
            return None

    def write_xyz_babel(self, directory, model_num=None, get_path=False,
                        rndname=True, filename=None):
        """
        Writes a xyz file containing the 3D coordinates of each particle in the
        model using a file format compatible with babel
        (http://openbabel.org/wiki/XYZ_%28format%29).
        Outfile is tab separated column with the bead number being the
        first column, then the genomic coordinate and finally the 3
        coordinates X, Y and Z
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
        :param None filename: overide the default file name writing
        """
        if filename:
            path_f = '%s/%s' % (directory, filename)
        else:
            if rndname:
                path_f = '%s/model.%s.xyz' % (directory, self['rand_init'])
            else:
                path_f = '%s/model_%s_rnd%s.xyz' % (directory, model_num,
                                                    self['rand_init'])
        out = ''
        # Write header as number of atoms
        out += str(len(self['x']))
        # Write comment as type of molecule
        out += "\nDNA\n"

        form = "%s\t%.3f\t%.3f\t%.3f\n"
        # TODO: do not use resolution directly -> specific to Hi-C
        chrom_list = self['description']['chromosome']
        chrom_start = self['description']['start']
        chrom_end = self['description']['end']
        if not isinstance(chrom_list, list):
            chrom_list = [chrom_list]
            chrom_start = [chrom_start]
            chrom_end = [chrom_end]
        chrom_start = [int(c)/int(self['description']['resolution']) for c in chrom_start]
        chrom_end = [int(c)/int(self['description']['resolution']) for c in chrom_end]
        offset = 0
        for crm in range(len(chrom_list)):
            for i in range(chrom_start[crm]+offset,chrom_end[crm]+offset):
                out += form % (
                    i + 1,
                    '%s:%s-%s' % (
                        chrom_list[crm],
                        int(chrom_start[crm] or 0) +
                            int(self['description']['resolution']) * (i - offset) + 1,
                        int(chrom_start[crm] or 0) +
                            int(self['description']['resolution']) * (i + 1 - offset)),
                    round(self['x'][i], 3),
                    round(self['y'][i], 3), round(self['z'][i], 3))
            offset += (chrom_end[crm]-chrom_start[crm])
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
        if tool == 'plot':
            x, y, z = self['x'], self['y'], self['z']
            plot_3d_model(x, y, z, color=color, **kwargs)
            return
        self.write_cmm('/tmp/', color=color, **kwargs)
        chimera_view(['/tmp/model.%s.cmm' % (self['rand_init'])],
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd,
                     center_of_mass=center_of_mass, gyradius=gyradius)
