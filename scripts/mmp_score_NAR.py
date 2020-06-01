#! /usr/bin/python
"""
12 janv. 2015

DESCRIPTION
-----------

Computes the modeling potential of a matrix (MMP).


DEPENDENCIES
------------

scipy
numpy
matplotlib *

(*) only to get figures


USAGE
-----

python mmp_score.py -i PATH -o PATH --plot


CITING
------
Assessing the limits of restraint-based 3D modeling of genomes and genomic domains.
Marie Trussart, Francois Serra, Davide Bau, Ivan Junier, Luis Serrano and Marc A. Marti-Renom
2015

"""

from argparse     import ArgumentParser
from numpy        import linalg, array, copy, log2, std, mean
from numpy.random import shuffle as np_shuffle
from scipy.stats  import skew, kurtosis, norm as sc_norm
import sys, os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    import matplotlib.gridspec as gridspec
except ImportError:
    sys.stderr.write('WARNING: Matplotlib not installed trying to plot ' +
             'something will raise uggly errors\n')

def randomize_matrix(data, savefig=None):
    size = len(data)
    rand_data = copy(data)
    for d in range(size):
        diag = list(zip(*[range(d, size), range(size - d)]))
        rdiag = diag[:]
        np_shuffle(rdiag)
        for v in range(len(diag)):
            val = data[diag[v][0]][diag[v][1]]
            a, b = rdiag[v][0], rdiag[v][1]
            rand_data[b][a] = rand_data[a][b] = val
    if savefig:
        plt.subplot(211)
        plt.imshow(log2(data), interpolation='none')
        plt.subplot(212)
        plt.imshow(log2(rand_data), interpolation='none')
        plt.savefig(savefig, format='pdf')
        plt.close('all')
    return rand_data

def read_pw_file(fnam):
    vals = {}
    mv = 0
    for line in open(fnam):
        if line.startswith('#'):
            continue
        a, b, v = line.split()
        a, b = int(a), int(b)
        mv = max(mv, max(a, b))
        vals[(a, b)] = v

    return array([array([float(vals.get((i, j), vals.get((j, i), 0)))
                         for j in range(mv)])
                  for i in range(mv)])

def main():
    """
    main function
    """
    opts          = get_options()
    fnam          = opts.fnam

    if opts.outdir == 'tmp/':
        opts.outdir = 'tmp_' + opts.fnam.split('/')[-1].replace('_xyz.txt', '').replace('list_files_', '').replace('.txt', '') + '/'

    sys.stdout.write('\n')
    sys.stdout.write('\n                     Processing\n')
    sys.stdout.write('                     ----------\n\n')

    sys.stdout.write('  - reading data\n')
    if opts.abc:
        data = read_pw_file(fnam)
    else:
        try:
            data = array([array([float(v) for v in l.split()]) for l in open(fnam)])
        except ValueError:
            data = array([array([float(v) for v in l.split(',')]) for l in open(fnam)])

    if opts.start != 1 or opts.end != -1:
        data = array([array([i for i in d[opts.start-1:opts.end]])
                      for d in data[opts.start-1:opts.end]])

    sys.stdout.write('  - getting EigenVectors\n')
    egval, _ = linalg.eigh(data)
    # sort eigenvalues/vectors
    idx = (-egval).argsort()
    egval = egval[idx]

    regvals = []

    sys.stdout.write('  - randomization\n')
    for i in range(int(opts.nrand)):
        sys.stdout.write('\r    ' + str(i + 1) + ' / ' + str(opts.nrand))
        sys.stdout.flush()
        regval, _ = linalg.eigh(randomize_matrix(data))
        regval = [abs(j) for j in regval]
        regval.sort(reverse=True)
        regvals.append( regval)
    sys.stdout.write('\n')
    regvals = list(zip(*regvals))
    rvmean = []
    for rv in regvals:
        rvmean.append(mean(rv))
    total = sum(rvmean)/100
    rvmean = [i/total for i in rvmean]

    err = []
    for rv in regvals:
        rvstd = std(rv/total)
        err.append(2 * rvstd)

    zdata = sorted(log2([data[i][j] for i in range(len(data))
                         for j in range(i, len(data)) if data[i][j]]))
    skewness = skew(zdata)
    kurtness = kurtosis(zdata)

    if opts.plot:
        os.system('mkdir -p %s' % opts.outdir)
        # matrix plot
        _ = plt.figure(figsize=(14, 8))
        gs = gridspec.GridSpec(7, 5, wspace=0.5, hspace=1.5)
        ax1 = plt.subplot(gs[:   , 0:3])
        ax2 = plt.subplot(gs[1:5 , 3: ])
        ax3 = plt.subplot(gs[5:7 , 3: ])
        img = ax2.imshow(log2(data), interpolation='none')
        plt.colorbar(img, ax=ax2)
        #plt.subplots_adjust(right=0.8, left=0.6, hspace=0.3)

    if opts.plot:
        ax2.set_title('Original matrix', size=12)
        ax2.tick_params(axis='both', which='major', labelsize=10)
        ax2.set_xlabel('Bin')
        ax2.set_ylabel('Bin')

        normfit = sc_norm.pdf(zdata, mean(zdata), std(zdata))
        normplot = ax3.plot(zdata, normfit, ':o', color='grey', alpha=.4,
                            markersize=.5)
        ax3.tick_params(axis='both', which='major', labelsize=10)

        ax3.hist(zdata, bins=20, density=True, alpha=0.7, color='r')
        ax3.set_xlabel('Z-score')
        ax3.set_ylabel('Frequency')
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        rcParams['axes.axisbelow']  = True
        #ax3.minorticks_on()
        #ax3.grid(ls=':', color='grey', alpha=.7, lw=.5, which='major')
        #plt.savefig(opts.outdir + '/matrix_small.png',
        #            format='png')
        #plt.close('all')


        #plt.imshow(log2(data), interpolation='none')
        #plt.title('Original matrix')
        #plt.colorbar()
        #plt.savefig(opts.outdir + '/matrix.png',
        #            format='png')
        #plt.close('all')

        # distribution plot
        #axe = plt.axes(axisbelow=True)
        rcParams['xtick.direction'] = 'out'
        rcParams['ytick.direction'] = 'out'
        rcParams['axes.axisbelow']  = True
        rcParams['axes.grid']       = True
        rcParams['grid.color']      = 'w'
        rcParams['grid.linestyle']  = '-'
        rcParams['grid.linewidth']  = 2
        # rcParams['grid.alpha']      = .3
        ax1.minorticks_on()
        ax1.grid(ls='-', color='w', alpha=.3, lw=2, which='major')
        ax1.grid(ls='-', b=True, color='w', alpha=.3, lw=1, which='minor')
        ax1.spines['top'].set_color('none')
        ax1.spines['right'].set_color('none')
        ax1.spines['bottom'].set_color('none')
        ax1.spines['left'].set_color('none')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.set_xscale('log')
        ax1.set_facecolor((.9,.9,.9))

        ax1.errorbar(range(1, 1 + len(rvmean)), rvmean, yerr=err, ecolor='red',
                     color='orange', lw=2,
        label='%s randomizations' % (opts.nrand))

    total = sum(abs(egval)) / 100
    egval = array(sorted([e/total for e in abs(egval)], reverse=True))

    for i in range(len(rvmean)):
        if rvmean[i] + err[i] > egval[i]:
            break
    signifidx = i
    signifsum = sum(abs(egval[:signifidx]))
    signifcontr = 100 * sum(abs(egval[:signifidx])) / sum(abs(egval))

    size = len(data)

    sev = sum(egval[:signifidx]-rvmean[:signifidx])

    if opts.plot:
        ax1.plot(range(1, 1 + len(rvmean)), egval,
                 color='green', lw=2, label='Observed data')

        ax1.fill_between(range(1, 1 + len(rvmean)), rvmean, egval,
             where=(array(rvmean) + array(err))<egval,
             facecolor='green', interpolate=True, alpha=0.2)
        ax1.fill_between(range(1, 1 + len(rvmean)), rvmean, egval,
             where=(array(rvmean) + array(err))>egval,
             facecolor='red'  , interpolate=True, alpha=0.2)
        ax1.set_xlim((1,len(rvmean)))
        ax1.set_ylim((0, max(max(rvmean), max(egval))))
        diff = float(ax1.get_ylim()[1] - ax1.get_ylim()[0]) / 15
        ax1.legend(frameon=False, loc='upper right', prop={'size': 10})
        ax1.set_xlabel('Log indexes of Eigenvalues')
        ax1.set_ylabel('Eigenvalues (percentage of total)')
        ax1.set_title(opts.fnam.split('/')[-1].replace('_xyz.txt', '').replace(
        'list_files_', '').replace('.txt', ''))
        #plt.subplots_adjust(right=0.6)

        #img = Image.open(opts.outdir + '/matrix_small.png')
        #fig.figimage(img, 640, -160)

    minv = float(min([i for d in data for i in d if i])) / 2
    if minv == 0.5:
        minv = 1./(len(data)**2)

    mmp = -0.0002 * size + 0.0335 * skewness - 0.0229 * kurtness + 0.0069 * sev + 0.8126

    sys.stdout.write('\n')
    sys.stdout.write('\n                       Results\n')
    sys.stdout.write('                       -------\n\n')


    sys.stdout.write('                  MMP score: %.4f\n\n' % mmp)

    ex_a1, ex_b1 = [0.6975926,  0.2548171]
    supa1, supb1 = [0.69300732000423904, 0.29858572176099613]
    lowa1, lowb1 = [0.70217788900976075, 0.211048473299004]

    scc     = (mmp - ex_b1 ) / ex_a1
    scc_up1 = (mmp - supb1 ) / supa1
    scc_lw1 = (mmp - lowb1 ) / lowa1

    sys.stdout.write('  predicted dSCC is %.3f (%.3f-%.3f 68%% confidence)\n' % (scc , scc_up1 , scc_lw1 ))

    ex_a75, ex_b75 = [0.6975926,  0.2548171]
    supa75, supb75 = [0.69230778430383244, 0.30526310790548261]
    lowa75, lowb75 = [0.70287742471016734, 0.20437108715451746]

    scc_up75 = (mmp - supb75 ) / supa75
    scc_lw75 = (mmp - lowb75 ) / lowa75

    sys.stdout.write('                        (%.3f-%.3f 75%% confidence)\n' % (scc_up75 , scc_lw75 ))

    ex_a2, ex_b2 = [0.6975926,  0.2548171]
    supa2, supb2 = [0.68855373600821357, 0.34109720480765293]
    lowa2, lowb2 = [0.70663147300578644, 0.16853699025234709]

    scc_up2 = (mmp - supb2 ) / supa2
    scc_lw2 = (mmp - lowb2 ) / lowa2

    sys.stdout.write('                        (%.3f-%.3f 95%% confidence)\n' % (scc_up2 , scc_lw2 ))

    if opts.plot:
        # write the log
        log = ''
        log +=  '    1- Matrix size (number of eigenvalues): %s\n' % (len(egval))
        log +=  "    2- Skewness of the distribution: %0.3f\n" % (skewness)
        log +=  "    3- Kurtosis of the distribution: %0.3f\n" % (kurtness)
        log +=  "    4- Sum of differences signif EV real-rand: %0.3f\n\n" % (sev)
        plt.figtext(0.62, 0.77, log, size='small')
        log =  "MMP score: %.3f\n" % (mmp)
        log +=  "Predicted dSCC: %.3f (%.3f-%.3f at 95%% conf)\n" % (scc, scc_up2, scc_lw2)
        plt.figtext(0.61, 0.87, log, size=12)
        plt.savefig(opts.outdir + '/distribution_eigenvalues_obs-rand_NAR.pdf',
                    format='pdf')
        plt.close('all')

def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH [options]")

    parser.add_argument('-i', dest='fnam', metavar='PATH', required=True,
                        default=False, help='input matrix file')
    parser.add_argument('-o', dest='outdir', metavar='PATH',
                        default='tmp/',
                        help='[%(default)s] directory to store results')
    parser.add_argument('--abc', dest='abc', action='store_true',
                        default=False,
                        help='''[%(default)s] in case the inputnfile is an
                        interaction pair list''')
    parser.add_argument('--plot', dest='plot', action='store_true',
                        default=False,
                        help='''[%(default)s] do the plot''')
    parser.add_argument('-s', dest='start',  default=1, type=int,
                        help='''[%(default)s] (inclusive) start position''')
    parser.add_argument('-e', dest='end',  default=-1, type=int,
                        help='''[%(default)s] (inclusive) end position --
                        -1 means last''')
    parser.add_argument('--nrand', dest='nrand',
                        default=10,metavar='INT',
                        help='''[%(default)s] number of randomizations''')

    opts = parser.parse_args()
    return opts


if __name__ == "__main__":
    exit(main())