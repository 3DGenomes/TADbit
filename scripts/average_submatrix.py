import os
from argparse                  import ArgumentParser
from tarfile                   import open as taropen
from itertools                 import imap
from warnings                  import warn
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib                import pyplot as plt
from mpl_toolkits.mplot3d      import Axes3D
from scipy                     import ndimage
from pytadbit.utils.extraviews import tadbit_savefig


def write_mat(matrix, outfile, count):
    out = open(outfile, 'w')
    out.write('# average of %d sub-matrices\n' % (count))
    for i in matrix:
        out.write('\t'.join('%f' % v for v in i) + '\n')
    out.close()


def check_tar_index(tarfile):
    tarsize = os.path.getsize(tarfile)
    if os.path.exists(tarfile + 'i'):
        # check size of the associated tar file
        if int(open(tarfile + 'i').next().split()[-1]) != tarsize:
            raise Exception('ERROR: associated tar file changed sized.'
                            '       remove the *.tari file and index again.')
    else:
        return tarsize
    return


def index_tar(tarfile):
    tarh = taropen(tarfile, ignore_zeros=True)
    tarsize = check_tar_index(tarfile)
    if tarsize:  # index not created yet
        print '   * TAR index file not found, creating it...'
        out = open(tarfile + 'i', 'w')
        out.write('# TAR size: %d\n' % tarsize)
        for member in tarh:
            out.write('%s\t%d\t%d\n' % (
                member.name.rsplit('_', 1)[0], member.offset_data, member.size))
        out.close()


def load_index(tarfile):
    print ' - Loading TAR index file'
    def _trans(line):
        k, b, e = line.split()
        return k, (int(b), int(e))
    fh = open(tarfile + 'i')
    _ = fh.next()
    return dict(imap(_trans, fh))


def parse_tar(tarfile, listfile, index, outlist, reso, kind, size):
    if listfile:
        print ' - Parsing input list of coordinates'
        def divide(x):
            a, b = x.split()
            c, d, e = a.replace(':', '-').split('-')
            f, g, h = b.replace(':', '-').split('-')
            return ('%s:%d-%d' % (c, int(d) / reso, int(e) / reso),
                    '%s:%d-%d' % (f, int(g) / reso, int(h) / reso))
        list_coords = (divide(l) for l in open(listfile))
    else:
        list_coords = (k.split('_')[1:3] for k in index.keys() if k.startswith(kind))

    vmatrix = [[0 for _ in xrange(size)] for _ in xrange(size)]
    wmatrix = [[0 for _ in xrange(size)] for _ in xrange(size)]
    all_rows = set(range(size))
    all_cols = set(range(size))

    print ' - Extracting %s submatrices' % ('selected' if listfile else 'all')
    tarfh = open(tarfile, 'rb')
    outfh = open(outlist, 'w')
    count = 0
    for k1, k2 in list_coords:
        try:
            beg, end = index['%s_%s_%s' % (kind, k1, k2)]
        except KeyError:
            warn('WANRING: some elements not found. Writting them in:\n'
                 + ('%s\n' % outlist) +
                 '')
            outfh.write('%s\t%s\n' % (k1, k2))
            continue
        count += 1
        tarfh.seek(beg)
        lines = tarfh.read(end).strip().split('\n')
        badrows = lines[1].strip()[9:]
        try:
            rows = all_rows.difference(map(int, badrows.split(',')))
        except ValueError:
            rows = all_rows
        badcols = lines[2].strip()[9:]
        try:
            cols = all_cols.difference(map(int, badcols.split(',')))
        except ValueError:
            cols = all_cols
        for i in rows:
            for j in cols:
                wmatrix[i][j] += 1
        for line in lines[3:]:
            a, b, c = line.split()
            vmatrix[int(a)][int(b)] += float(c)

    outfh.close()
    if not os.path.getsize(outlist):
        os.system('rm -f ' + outlist)

    for i in xrange(size):
        for j in xrange(size):
            try:
                vmatrix[i][j] /= wmatrix[i][j]
            except ZeroDivisionError:
                vmatrix[i][j] = float('nan')
    return vmatrix, count


def main():
    opts     = get_options()
    tarfile  = opts.tarfile
    listfile = opts.listfile
    outfile  = opts.outfile
    outplot  = opts.outplot
    outlist  = opts.outfile + '.lst'
    size     = opts.size
    reso     = opts.reso
    kind     = opts.matrix
    minmax   = opts.minmax

    if not opts.plot_only:
        print ' - Checking TAR index file'
        index_tar(tarfile)

        index = load_index(tarfile)

        vmatrix, count = parse_tar(tarfile, listfile, index, outlist, reso, kind, size)

        write_mat(vmatrix, outfile, count)

    if outplot:
        print ' - plotting'
        do_3d_plot(outfile, outplot, size, count, minmax, sigma=opts.sigma, log=opts.log)

    print 'Done.'


def do_3d_plot(nam, outfile, size, count, minmax, sigma=0, log=False):
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    beg = -size / 2
    end =  size / 2
    X = np.arange(beg, end, 1)
    Y = np.arange(beg, end, 1)
    X, Y = np.meshgrid(X, Y)
    Z = np.array([np.array([float(i) for i in l.split()]) for l in open(nam) if not l.startswith('#')])
    plt.title(nam + '\nMean: %.3f, median: %.3f, standard-deviation: %.3f (N=%d)' % (np.mean(Z), np.median(Z), np.std(Z), count))
    if sigma:
        Z = ndimage.gaussian_filter(Z, sigma=sigma, order=0)
    if log:
        Z = np.log(Z)
        zspan = minmax if minmax else np.max(np.abs(Z))
        zmax =  zspan
        zmin = -zspan
    else:
        zspan = minmax if minmax else np.max(np.abs(Z - 1))
        zmin = -zspan + 1
        zmax =  zspan + 1
    cmap = 'coolwarm'  # 'coolwarm'
    _ = ax.contourf(X, Y, Z, zdir='z', offset=zmin,
                    cmap=cmap, vmin=zmin, vmax=zmax)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmap,
                           linewidth=0, antialiased=True, alpha=1,
                           vmin=zmin, vmax=zmax, shade=True)
    ax.set_zlim3d(zmin, zmax)
    ax.view_init(elev=15, azim=25)
    cb = fig.colorbar(surf, shrink=0.5, aspect=20)
    cb.set_label('%sverage normalized interactions%s' %
                 ('Log a' if log else 'A',
                  '\nSmoothed with $\sigma=%s$' % sigma))
    tadbit_savefig(outfile)


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -l PATH -o PATH -n INT [options]")

    parser.add_argument('-i', '--input_tar', dest='tarfile', metavar='',
                        required=False, default=False,
                        help='input TAR file with all matrices.')
    parser.add_argument('-l', '--input_list', dest='listfile', metavar='',
                        required=False, default=False,
                        help=('input file with the list of coordinates'
                              'corresponding to matrices to extract. '
                              'By default all TAR file is parsed.'
                              'If coordinates are not found, a new file with '
                              'the list of not-found coordinates will be '
                              'created, in order to compute the matrices.'))
    parser.add_argument('-n', '--size', dest='size', type=int, metavar='',
                        required=False, help='size of the matrices')
    parser.add_argument('-r', '--reso', dest='reso', type=int, metavar='',
                        required=False, help='resolution at which matrices were created')
    parser.add_argument('--matrix', dest='matrix', metavar='', type=str,
                        choices=['nrm', 'raw', 'dec'], default='dec',
                        help='''[%(default)s] which matrix to generate
                        (choices: %(choices)s)''')
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='',
                        required=True, default=False,
                        help=('output file to write matrix (or from where to '
                              'read it, with the plot_only option).'))
    parser.add_argument('--plot_only', dest='plot_only', action='store_true',
                        default=False, help=('in case average matrix was '
                                             'already generated'))
    parser.add_argument('-p', '--outplot', dest='outplot', metavar='',
                        default=False,
                        help='output file to plot matrix.')
    parser.add_argument('--minmax', dest='minmax', metavar='',
                        type=float,
                        help="to scale the Y axis.")
    parser.add_argument('--sigma', dest='sigma', type=float,
                        default=0.0,
                        help='[%(default)s] smoothing parameter for the plotting')
    parser.add_argument('--log', dest='log', action='store_true',
                        default=False,
                        help='apply log on interaction counts for plotting')

    opts = parser.parse_args()

    if not opts.plot_only:
        if not opts.size:
            raise Exception('ERROR: should input matrix size')
        if not opts.tarfile:
            raise Exception('ERROR: should input path to TAR file')
    else:
        if not opts.outplot:
            raise Exception('ERROR: should input path to output plot file')
        if not os.path.exists(opts.outfile):
            raise Exception('ERROR: "output" path to matrix file should exist '
                            'for plotting')
    return opts

if __name__ == "__main__":
    exit(main())
