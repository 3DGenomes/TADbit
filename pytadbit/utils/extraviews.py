"""
06 Aug 2013


"""
from warnings import warn
import numpy as np


try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')


def nicer(res):
    """
    writes resolution number for human beings.
    """
    if not res % 1000000000:
        return str(res)[:-9] + 'Gb'
    if not res % 1000000:
        return str(res)[:-6] + 'Mb'
    if not res % 1000:
        return str(res)[:-3] + 'Kb'
    return str(res) + 'b'


COLOR = {None: '\033[31m', # red
         0   : '\033[34m', # blue
         1   : '\033[34m', # blue
         2   : '\033[34m', # blue
         3   : '\033[36m', # cyan
         4   : '\033[0m' , # white
         5   : '\033[1m' , # bold white
         6   : '\033[33m', # yellow
         7   : '\033[33m', # yellow
         8   : '\033[35m', # purple
         9   : '\033[35m', # purple
         10  : '\033[31m'  # red
         }

COLORHTML = {None: '<span style="color:red;">'       , # red
             0   : '<span>'                          , # blue
             1   : '<span style="color:blue;">'      , # blue
             2   : '<span style="color:blue;">'      , # blue
             3   : '<span style="color:purple;">'    , # purple
             4   : '<span style="color:purple;">'    , # purple
             5   : '<span style="color:teal;">'      , # cyan
             6   : '<span style="color:teal;">'      , # cyan
             7   : '<span style="color:olive;">'     , # yellow
             8   : '<span style="color:olive;">'     , # yellow
             9   : '<span style="color:red;">'       , # red
             10  : '<span style="color:red;">'         # red
             }


def colorize(string, num, ftype='ansi'):
    """
    Colorize with ANSII colors a string for printing in shell. this acording to
    a given number between 0 and 10

    :param string: the string to colorize
    :param num: a number between 0 and 10 (if None, number will be equal to 10)

    :returns: the string 'decorated' with ANSII color code
    """
    color = COLOR if ftype=='ansi' else COLORHTML
    return '{}{}{}'.format(color[num], string,
                           '\033[m' if ftype=='ansi' else '</span>')


def color_residues(n_part):
    """
    :param n_part: number of particles
    
    :returns: a list of rgb tuples (red, green, blue)
    """
    result = []
    for n in xrange(n_part):
        red = float(n + 1) / n_part
        result.append((red, 0, 1 - red))
    return result


def augmented_dendrogram(clust_count=None, dads=None, objfun=None, color=False,
                         axe=None, savefig=None, *args, **kwargs):

    from scipy.cluster.hierarchy import dendrogram
    fig = plt.figure(figsize=(8, 8))
    if axe:
        ax = axe
        fig = axe.get_figure()
        ddata = dendrogram(*args, **kwargs)
        plt.clf()
    else:
        ddata = dendrogram(*args, **kwargs)
        plt.clf()
        ax = fig.add_subplot(111)
        ax.patch.set_facecolor('lightgrey')
        ax.patch.set_alpha(0.4)
        ax.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
        ax.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
        ax.set_axisbelow(True)
        # remove tick marks
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False)
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False, which='minor')
    
    # set dict to store data of each cluster (count and energy), depending on
    # x position in graph.
    leaves = {}
    dist = ddata['icoord'][0][2] - ddata['icoord'][0][1]
    for i, x in enumerate(ddata['leaves']):
        leaves[dist*i + dist/2] = x
    minnrj = min(objfun.values())-1
    maxnrj = max(objfun.values())-1
    difnrj = maxnrj - minnrj
    total = sum(clust_count.values())
    if not kwargs.get('no_plot', False):
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'],
                           ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            # plt.plot(x, y, 'ro')
            plt.hlines(y, i[1], i[2], lw=2, color='grey')
            # for eaxch branch
            for i1, d1, d2 in zip(i[1:3], [d[0], d[3]], [d[1], d[2]]):
                try:
                    lw = float(clust_count[leaves[i1]])/total*10*len(leaves)
                except KeyError:
                    lw = 1.0
                nrj = objfun[leaves[i1]] if leaves[i1] in objfun else maxnrj
                ax.vlines(i1, d1-(difnrj-(nrj-minnrj)), d2, lw=lw,
                          color=(c if color else 'grey'))
                if leaves[i1] in objfun:
                    ax.annotate("%.3g" % (leaves[i1]),
                                (i1, d1-(difnrj-(nrj-minnrj))),
                                xytext=(0, -8),
                                textcoords='offset points',
                                va='top', ha='center')
            leaves[(i[1] + i[2])/2] = dads[leaves[i[1]]]
    bot = -int(difnrj)/10000 * 10000
    plt.yticks([bot+i for i in xrange(0, -bot-bot/10, -bot/10)],
               ["{:,}".format(int(minnrj)/10000 * 10000  + i)
                for i in xrange(0, -bot-bot/10, -bot/10)], size='small')
    ax.set_ylabel('Minimum IMP objective function')
    ax.set_xticks([])
    ax.set_xlim((plt.xlim()[0] - 2, plt.xlim()[1]))
    ax.figure.suptitle("Dendogram of clusters of 3D models")
    ax.set_title("Branch length proportional to model's objective function " +
                 "final value\n" +
                 "Branch width to the number of models in the cluster",
                 size='small')
    if savefig:
        fig.savefig(savefig)
    elif not axe:
        plt.show()
    return ddata


def plot_hist_box(data, part1, part2, axe=None, savefig=None):
    # setup the figure and axes
    if axe:
        fig = axe.get_figure()
    else:
        fig = plt.figure(figsize=(6, 6))
    bpAx = fig.add_axes([0.2, 0.7, 0.7, 0.2])   # left, bottom, width, height:
                                                # (adjust as necessary)
    bpAx.patch.set_facecolor('lightgrey')
    bpAx.patch.set_alpha(0.4)
    bpAx.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
    bpAx.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
    bpAx.set_axisbelow(True)
    bpAx.minorticks_on() # always on, not only for log
    # remove tick marks
    bpAx.tick_params(axis='both', direction='out', top=False, right=False,
                   left=False, bottom=False)
    bpAx.tick_params(axis='both', direction='out', top=False, right=False,
                   left=False, bottom=False, which='minor')
    # plot stuff
    bp = bpAx.boxplot(data, vert=False)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['medians'], color='darkred')
    plt.setp(bp['fliers'], color='darkred', marker='+')
    bpAx.plot(sum(data)/len(data), 1, 
              color='w', marker='*', markeredgecolor='k')
    bpAx.annotate('{:.4}'.format(bp['boxes'][0].get_xdata()[0]),
                  (bp['boxes'][0].get_xdata()[0], bp['boxes'][0].get_ydata()[1]),
                  va='bottom', ha='center', xytext=(0, 2),
                  textcoords='offset points',
                  size='small')
    bpAx.annotate('{:.4}'.format(bp['boxes'][0].get_xdata()[2]),
                  (bp['boxes'][0].get_xdata()[2], bp['boxes'][0].get_ydata()[1]),
                  va='bottom', ha='center', xytext=(0, 2),
                  textcoords='offset points',
                  size='small')
    bpAx.annotate('{:.4}'.format(bp['medians'][0].get_xdata()[0]),
                  (bp['medians'][0].get_xdata()[0], bp['boxes'][0].get_ydata()[0]),
                  va='top', ha='center', xytext=(0, -2),
                  textcoords='offset points', color='darkred',
                  size='small')
    histAx = fig.add_axes([0.2, 0.2, 0.7, 0.5]) # left specs should match and
                                                # bottom + height on this line should
                                                # equal bottom on bpAx line
    histAx.patch.set_facecolor('lightgrey')
    histAx.patch.set_alpha(0.4)
    histAx.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
    histAx.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
    histAx.set_axisbelow(True)
    histAx.minorticks_on() # always on, not only for log
    # remove tick marks
    histAx.tick_params(axis='both', direction='out', top=False, right=False,
                   left=False, bottom=False)
    histAx.tick_params(axis='both', direction='out', top=False, right=False,
                   left=False, bottom=False, which='minor')
    h = histAx.hist(data, bins=20, alpha=0.5, color='darkgreen')
    # confirm that the axes line up 
    xlims = np.array([bpAx.get_xlim(), histAx.get_xlim()])
    for ax in [bpAx, histAx]:
        ax.set_xlim([xlims.min(), xlims.max()])
    bpAx.set_xticklabels([])  # clear out overlapping xlabels
    bpAx.set_yticks([])  # don't need that 1 tick mark
    plt.xlabel('Distance between particles (nm)')
    plt.ylabel('Number of observations')
    bpAx.set_title('Histogram and boxplot of distances between particles {} and {}'.format(part1, part2))
    if savefig:
        fig.savefig(savefig)
    elif not axe:
        plt.show()

