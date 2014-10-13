#! /usr/bin/python

"""
25 sept. 2014

script based on this publication:

Lesne, A., Riposo, J., Roger, P., Cournac, A., & Mozziconacci, J. (2014).
3D genome reconstruction from chromosomal contacts.
Nature Methods, 4, 10-13. doi:10.1038/nmeth.3104

"""

from argparse import ArgumentParser
from mpl_toolkits.mplot3d import Axes3D
from numpy import linalg, array, log2, copy, minimum, newaxis, percentile
from matplotlib import pyplot as plt


def floyd_warshall_numpy(mat):
    '''
    floyd_warshall_numpy(adjacency_matrix) -> shortest_path_distance_matrix
    A vectorized NumPy implementation of the Floyd-Warshall algorithm.
    Input
    An NxN NumPy array describing the directed distances between N nodes.
    adjacency_matrix[i,j] = distance to travel directly from node i to node j
    (without passing through other nodes)

    Notes:
    * If there is no edge connecting i->j then adjacency_matrix[i,j] should be
      equal to numpy.inf.
    * The diagonal of adjacency_matrix should be zero.

    Output
    An NxN NumPy array such that result[i,j] is the shortest distance to travel
    between node i and node j.
    If no such path exists then result[i,j] == numpy.inf
    '''
    for k in xrange(len(mat)):
        mat = minimum(mat, mat[newaxis, k, :] + mat[:, k, newaxis])
    return mat

def adj(g):
    """
    Convert a directed graph to an adjaceny matrix.
    >>> g = {1: {2: 3, 3: 8, 5: -4}, 2: {4: 1, 5: 7}, 3: {2: 4},
             4: {1: 2, 3: -5}, 5: {4: 6}}
    >>> adj(g)
    {1: {1: 0, 2: 3, 3: 8, 4: inf, 5: -4},
     2: {1: inf, 2: 0, 3: inf, 4: 1, 5: 7},
     3: {1: inf, 2: 4, 3: 0, 4: inf, 5: inf},
     4: {1: 2, 2: inf, 3: -5, 4: 0, 5: inf},
     5: {1: inf, 2: inf, 3: inf, 4: 6, 5: 0}}
    """
    vertices = range(len(g))

    dist = copy(g)  # copy g
    for i in vertices:
        dist[i][i] = 0.
    return dist


def main():
    """
    main function
    """
    opts          = get_options()

    fnam          = opts.fnam
    scale         = float(opts.scale)
    show          = opts.show
    cmm           = opts.cmm
    strict_julien = opts.strict

    if opts.verbose:
        print '- Getting data'
    try:
        data = [[float(v) for v in l.split()] for l in open(fnam)]
    except ValueError:
        data = [[float(v) for v in l.split(',')] for l in open(fnam)]

    to_remove = dict([(i, None) for i, d in enumerate(data) if sum(d)==10000])
    if opts.verbose:
        print to_remove
    N = len(data)
    clean_data = []
    for i in xrange(N):
        if i in to_remove:
            continue
        clean_data.append([])
        for j in xrange(N):
            if j in to_remove:
                continue
            clean_data[-1].append(data[i][j])
    data  = clean_data
    
    N = len(data)
    if opts.verbose:
        print '- Make the data symmetric'
    total = 0
    for i in xrange(N):
        for j in xrange(i + 1, N):
            data[i][j] = min(data[i][j], data[j][i])
            data[j][i] = data[i][j]
            total += data[j][i]

    if opts.verbose:
        print '- Replace zeros'
    if strict_julien:
        minv = total / N**2
    else:
        minv = min([data[i][j] for j in xrange(N)
                    for i in xrange(N) if data[i][j]]) /2
    if not opts.dist:
        data = array([array([1./(data[i][j] or minv) for j in xrange(N)])
                      for i in xrange(N)])
    else:
        data = array([array([data[i][j] for j in xrange(N)])
                      for i in xrange(N)])

    ############################################################################
    ### Here it starts
    ############################################################################
    if not opts.dist:
        if opts.verbose:
            print '- Finding shortest path'
        if opts.fw == 'numpy':
            D = floyd_warshall_numpy(adj(data))
        else:
            import pyximport
            pyximport.install(reload_support=True)
            import floyd_warshall
            if opts.fw == 'cythonP':
                D = floyd_warshall.floyd_warshall_parallelized(adj(data))
            elif opts.fw == 'cython':
                D = floyd_warshall.floyd_warshall_single_core(adj(data))
            else:
                raise NotImplementedError('choose between cython cythonP or numpy')
    else:
        D = data
    # get the square of D for performence
    D2 = [[(D[i][j])**2 for j in xrange(N)] for i in xrange(N)]

    if opts.verbose:
        print '- Distance to barycenter'
    uptri = sum([sum(D2[j][j + 1:N]) for j in xrange(N)]) / (N**2)
    dO2 = [sum(D2[i]) / N - uptri for i in xrange(N)]

    # metric matrix
    if opts.verbose:
        print '- Calculating metric matrix'
    M = [[(dO2[i] + dO2[j] - D2[i][j]) / 2 for j in xrange(N)]
         for i in xrange(N)]
    # Eigen vectors and values
    if opts.verbose:
        print '- Getting eigenvectors and eigenvalues'
    Eval, Evect = linalg.eigh(array(M))
    # sort eigen vectors
    idx = (-Eval).argsort()
    Eval = Eval[idx]
    Evect = Evect[:, idx]

    # print '-'*80
    # print ' '.join(['%f'%e for e in Eval])
    # print '-'*80
    # print (sum(Eval[:3]), sum([abs(e) for e in Eval]),
    #        sum(Eval[:3]) / sum([abs(e) for e in Eval]))
    if opts.verbose:
        print ('    contribution of the 3 first eigenvalues: %.1f%%'  %
               (sum(Eval[:3]) / sum([abs(e) for e in Eval]) * 100))

    # Coordinates
    if opts.verbose:
        print '- Estimating 3D coordinates'
    x = [Evect[i][0] * Eval[0]**.5 for i in xrange(N)]
    y = [Evect[i][1] * Eval[1]**.5 for i in xrange(N)]
    z = [Evect[i][2] * Eval[2]**.5 for i in xrange(N)]
    ############################################################################
    ### And that's it...
    ############################################################################
    if opts.verbose:
        print '- Generating outfiles'
    if scale:
        if opts.verbose:
            print '  * scaling'
        # scale distances
        distances = []
        for i in xrange(N - 1):
            distances.append(distance((x[i], y[i], z[i]),
                                      (x[i + 1], y[i + 1], z[i + 1])))
        meand = percentile(distances, 25)
        factor = scale / meand
        x = [i * factor for i in x]
        y = [i * factor for i in y]
        z = [i * factor for i in z]
        for p in sorted(to_remove.keys()):
            if p:
                x.insert(p, (x[p] + x[p-1])/2)
                y.insert(p, (y[p] + y[p-1])/2)
                z.insert(p, (z[p] + z[p-1])/2)
            else:
                x.insert(p, (x[p]))
                y.insert(p, (y[p]))
                z.insert(p, (z[p]))
    if show:
        plt.figure(figsize=(18, 7))
        plt.subplot(121)
        plt.imshow(log2([([data[k][l] for l in xrange(N)])
                         for k in xrange(N)]),
                   interpolation='none', cmap='jet')
        plt.colorbar()
        axe = plt.subplot(122, projection='3d')
        plot_3d_model(x, y, z, axe=axe)
        plt.show()
    if cmm:
        if opts.verbose:
            print '  * writing CMM and XYZ files'
        write_cmm(x, y, z, fnam + '.cmm', radius=30)
        write_xyz(x, y, z, fnam + '.xyz')
    if opts.mtx:
        if opts.verbose:
            print '  * writing distance matrix'
        out = open(fnam + '.tsv', 'w')
        for i in xrange(N):
            out.write('\t'.join(['%.4f' % distance((x[i], y[i], z[i]),
                                                   (x[j], y[j], z[j]))
                                 for j in xrange(N)]))
            out.write('\n')
        out.close()
    if opts.verbose:
        print '\nDone\n\n'


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH [options]")

    parser.add_argument('-i', '--infile', dest='fnam',
                        default=False, help='matrix file.')
    parser.add_argument('--cmm', dest='cmm', action='store_true',
                        default=False, help='generate cmm and xyz files')
    parser.add_argument('--mtx', dest='mtx', action='store_true',
                        default=False, help='generate distance matrix file')
    parser.add_argument('--show', dest='show', action='store_true',
                        default=False, help='display matplotlib 3D plot')
    parser.add_argument('--dist', dest='dist', action='store_true',
                        default=False, help='''input is already a distance
                        matrix''')
    parser.add_argument('--verbose', dest='verbose', action='store_true',
                        default=False, help='''print process info''')
    parser.add_argument('--strict', dest='strict', action='store_true',
                        default=False, help='''Replace zeros by sum of the
                        matrix divided by the square of the number of rows
                        (recommended for binary matrices). By default it uses
                        half the minimum value''')
    parser.add_argument('--scale', dest='scale',
                        default=0, help='''[%(default)s] Average distance (nm)
                        beetween two particles; by default no scalling is
                        applied''')
    parser.add_argument('--fw', dest='fw', default='numpy',
                        help='''[%(default)s] implementation to
                        search shortest path using Floyd-Warshall can be one of
                        "numpy", "cython" or "cythonP" (parallel version for
                        very large matrices)''')

    opts = parser.parse_args()
    if not opts.fnam:
        exit(parser.print_help())
    return opts


def distance(part1, part2):
    """
    Calculates the distance between two particles.

    :param part1: coordinate in list format (x, y, z)
    :param part2: coordinate in list format (x, y, z)

    :returns: distance between two points in space
    """
    return ((part1[0] - part2[0])**2 +
            (part1[1] - part2[1])**2 +
            (part1[2] - part2[2])**2)**.5


def color_residues(x):
    """
    Function to color residues from blue to red.

    :param x: list of x coordinates

    :returns: a list of rgb tuples (red, green, blue), each between 0 and 1.
    """
    result = []
    for n in xrange(len(x)):
        red = float(n + 1) / len(x)
        result.append((red, 0, 1 - red))
    return result


def plot_3d_model(x, y, z, label=False, axe=None, thin=False, savefig=None,
                  show_axe=False, azimuth=-90, elevation=0., color='index',
                  **kwargs):
    """
    Given a 3 lists of coordinates (x, y, z) plots a three-dimentional model
    using matplotlib

    :param x: list of x coordinates
    :param y: list of y coordinates
    :param z: list of z coordinates
    :param False label: show labels
    :param None axe: a matplotlib.axes.Axes object to define the plot
       appearance
    :param False thin: draw a thin black line instead of representing particles
       and edges
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param -90 azimuth: angle to rotate camera along the y axis
    :param 0 elevation: angle to rotate camera along the x axis
    """
    show = False
    if isinstance(color, str):
        if color == 'index':
            color = color_residues(x, **kwargs)
        else:
            raise NotImplementedError(('%s type of coloring is not yet ' +
                                       'implemeted\n') % color)
    elif hasattr(color, '__call__'):  # its a function
        color = color(x, **kwargs)
    elif not isinstance(color, list):
        raise TypeError('one of function, list or string is required\n')
    if not axe:
        fig = plt.figure(figsize=kwargs.get('figsize', (8, 8)))
        axe = fig.add_subplot(1, 1, 1, projection='3d')
        show = True
    if not show_axe:
        axe._axis3don = False
    axe.view_init(elev=elevation, azim=azimuth)
    if thin:
        axe.plot(x, y, z, color='black', lw=1, alpha=0.2)
    else:
        for i in xrange(len(x)-1):
            axe.plot(x[i:i+2], y[i:i+2], z[i:i+2],
                     color=color[i], lw=3)
            if label:
                axe.text(x[i], y[i], z[i], str(i), size=7)
        if label:
            axe.text(x[i + 1], y[i + 1], z[i + 1],str(i + 1), size=7)
        axe.scatter(x, y, z, color=color, s=50)
    axe.pbaspect = [1, 1, 1]
    if show:
        if savefig:
            nice_savefig(savefig)
        else:
            plt.show()


def nice_savefig(savefig):
    form = savefig[-4:].split('.')[1]
    if not form in ['png', 'pdf', 'ps', 'eps', 'svg']:
        raise NotImplementedError('File extension must be one of %s' %(
            ['png', 'pdf', 'ps', 'eps', 'svg']))
    plt.savefig(savefig, format=form)



def write_xyz(x, y, z, path):
    """
    :param x: list of x coordinates
    :param y: list of y coordinates
    :param z: list of z coordinates
    :param path: to outfile
    """
    out = ''
    form = "%s\t%s\t%.3f\t%.3f\t%.3f\n"
    for i in xrange(len(x)):
        out += form % (
            i + 1,
            '??:??-??',
            round(x[i], 3),
            round(y[i], 3), round(z[i], 3))
    out_f = open(path, 'w')
    out_f.write(out)
    out_f.close()

def write_cmm(x, y, z, path, radius=30, **kwargs):
    """
    Write cmm file for Chimera

    :param x: list of x coordinates
    :param y: list of y coordinates
    :param z: list of z coordinates
    :param path: to outfile
    :param radius: radius of each particle (nm)
    """

    color = color_residues(x, **kwargs)
    out = '<marker_set name=\"quick_model\">\n'
    form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
            ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
            'radius=\"' + #str(30) +
            str(radius) +
            '\" note=\"%s\"/>\n')
    for i in xrange(len(x)):
        out += form % (i + 1,
                       x[i], y[i], z[i],
                       color[i][0], color[i][1], color[i][2], i + 1)
    form = ('<link id1=\"%s\" id2=\"%s\" r=\"1\" ' +
            'g=\"1\" b=\"1\" radius=\"' + str(radius / 3) +
            # str(self['radius']/2) +
            '\"/>\n')
    for i in xrange(1, len(x)):
        out += form % (i, i + 1)
    out += '</marker_set>\n'

    out_f = open(path, 'w')
    out_f.write(out)
    out_f.close()


if __name__ == "__main__":
    exit(main())
