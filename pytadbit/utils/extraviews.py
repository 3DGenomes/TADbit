"""
06 Aug 2013


"""
from warnings import warn
import numpy as np
from subprocess import Popen, PIPE


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
    ax.set_xlim((plt.xlim()[0] - 2, plt.xlim()[1] + 2))
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



def chimera_view(cmm_file, chimera_bin='chimera',
                 shape='tube', chimera_cmd=None,
                 savefig=None):
    """
    """
    pref_f = '/tmp/tmp.cmd'
    out = open(pref_f, 'w')
    out.write('open {}\n'.format(cmm_file))
    if not chimera_cmd:
        out.write('''
focus
set bg_color white
windowsize 800 600
bonddisplay never #0
shape tube #0 radius 10 bandLength 200 segmentSubdivisions 100 followBonds on
clip yon -500
~label
set subdivision 1
set depth_cue
set dc_color black
set dc_start 0.5
set dc_end 1
scale 0.8
''')
        if savefig:
            if savefig.endswith('.png'):
                out.write('copy file {} png'.format(savefig))
            elif savefig[-4:] in ('.mov', 'webm'):
                out.write('''
movie record supersample 1
turn y 3 120
wait 120
movie stop
movie encode output {0}
'''.format(savefig))
            elif savefig:
                raise Exception('Not supported format, must be png, mov or webm\n')
    else:
        out.write('\n'.join(chimera_cmd) + '\n')
    out.close()
    
    Popen('{} {}'.format(chimera_bin, pref_f), shell=True).communicate()


def plot_3d_optimization_result(result, scale, scale_arange, max_dist_arange,
                                upfreq_arange, lowfreq_arange):
    """
    Displays a three dimensional scatter plot representing the result of the
    optimization.

    :param result: 3D numpy array contating correlation values
    :param scale: represent optimization result for the scale
    :param scale_arange: range of scale values used in the optimization
    :param max_dist_arange: range of max_dist values used in the optimization
    :param upfreq_arange: range of upfreq values used in the optimization
    :param lowfreq_arange: range of lowfreq values used in the optimization
    """
    # TOBEREMOVED
    try:
        scale_idx = scale_arange.index(scale)
    except IndexError:
        raise Exception('Scale not found, in scale_arange')
    result = result[scale_idx,:,:,:]

    x = [my_round(i, 3) for i in lowfreq_arange]
    y = [my_round(i, 3) for i in upfreq_arange]
    z = [my_round(i, 3) for i in max_dist_arange]
    sort_result =  sorted([(result[k, j, i], x[i], y[j], z[k])
                           for i in range(len(x))
                           for j in range(len(y))
                           for k in range(len(z))], key=lambda x: x[0],
                          reverse=True)[0]
    x = [i for i in max_dist_arange for j in upfreq_arange for k in lowfreq_arange]
    y = [j for i in max_dist_arange for j in upfreq_arange for k in lowfreq_arange]
    z = [k for i in max_dist_arange for j in upfreq_arange for k in lowfreq_arange]
    col = [result[i,j,k] for i in range(len(max_dist_arange)) for j in range(len(upfreq_arange)) for k in range(len(lowfreq_arange))]

    fig = plt.figure()
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('maxdist')
    ax.set_ylabel('upfreq')
    ax.set_zlabel('lowfreq')
    lol = ax.scatter(x, y, z, c=col, s=100, alpha=0.9)
    cbar = fig.colorbar(lol)
    cbar.ax.set_ylabel('Correlation value')
    plt.title(('Optimal IMP parameters (scale={})\n' +
               'Best for: lowfreq={}, upfreq={}, maxdist={}'
               ).format(round(scale,3), my_round(sort_result[1], 2),
                        my_round(sort_result[2], 2),
                        my_round(sort_result[3], 2)),
                 size='large')
    plt.show()


def my_round(num, val):
    num = round(num, val)
    return int(num) if num == int(num) else num


def plot_2d_optimization_result(result, axes=('scale', 'maxdist', 'upfreq', 'lowfreq'),
                                show_best=0, skip=None):
    """
    A grid of heatmaps representing the result of the optimization.

    :param result: 3D numpy array contating correlation values
    :param 'scale','maxdist','upfreq','lowfreq' axes: tuple of axes to
       represent. The order will define which parameter will be placed on the
       w, z, y or x axe.
    :param 0 show_best: number of best correlation value to identifie.
    :param None skip: a dict can be passed here in order to fix a given axe,
       e.g.: {'scale': 0.001, 'maxdist': 500}

    """

    from mpl_toolkits.axes_grid1 import AxesGrid
    import matplotlib.patches as patches

    ori_axes, axes_range, result = result
    trans = [ori_axes.index(a) for a in axes]
    axes_range = [axes_range[i] for i in trans]

    # transpose results
    result = result.transpose(trans)

    vmin = result.min()
    vmax = result.max()

    wax = [my_round(i, 3) for i in axes_range[0]]
    zax = [my_round(i, 3) for i in axes_range[1]]
    xax = [my_round(i, 3) for i in axes_range[3]]
    yax = [my_round(i, 3) for i in axes_range[2]]
    sort_result =  sorted([(result[i, j, k, l], wax[i], zax[j], xax[l], yax[k])
                           for i in range(len(wax))
                           for j in range(len(zax))
                           for k in range(len(yax))
                           for l in range(len(xax))
                           ], key=lambda x: x[0],
                          reverse=True)[:show_best+1]

    # skip
    wax_range = range(len(wax))[::-1]
    zax_range = range(len(zax))
    skip = {} if not skip else skip
    for i, k in enumerate(axes):
        if not k in skip:
            continue
        if i == 0:
            wax_range = [wax.index(skip[k])]
        elif i==1:
            zax_range = [zax.index(skip[k])]
        else:
            raise Exception('ERROR: skip keys must be one of the two first keywords passed as axes parameter')

    # best number of rows/columns
    ncols  = int(np.sqrt(len(zax_range)) + 0.999)
    nrows  = int(np.sqrt(len(zax_range)) + 0.5)
    nncols = int(np.sqrt(len(wax_range)) + 0.999)
    nnrows = int(np.sqrt(len(wax_range)) + 0.5)
    
    fig = plt.figure(figsize=((nncols+ncols)*3,(nrows+nrows)*3))
    
    for ii in wax_range:
        print ii
        grid = AxesGrid(fig, int(str(nnrows) + str(nncols) + str(ii)),
                        nrows_ncols = (nrows, ncols),
                        axes_pad = 0.0,
                        label_mode = "1",
                        share_all = True,
                        cbar_location="right",
                        cbar_mode="single" if not ii else None,
                        cbar_size="7%",
                        cbar_pad="2%",
                        )
        for i in zax_range:
            if not i:
                grid[i].text(len(xax) * ncols/2, len(yax) + 1.5,
                             axes[0] + ': ' + str(my_round(wax[ii], 3)),
                             {'ha':'center', 'va':'top'}, size='large')
            im = grid[i].imshow(result[ii, i, :, :], interpolation="nearest",
                                vmin=vmin, vmax=vmax)
            grid[i].tick_params(axis='both', direction='out', top=False,
                                right=False, left=False, bottom=False)
            rect = patches.Rectangle((-0.5, len(yax)-.5),len(xax), 1.5,
                                     facecolor='grey', alpha=0.5)
            grid[i].add_patch(rect)
            grid[i].text(np.mean(range(0, len(xax))),
                         max(range(0, len(yax))) + 1.25,
                         axes[1] + ': ' + str(my_round(zax[i], 3)),
                         {'ha':'center', 'va':'center'})
            for j, best  in enumerate(sort_result[:-1]):
                if best[2] == zax[i] and best[1] == wax[ii]:
                    grid[i].text(xax.index(best[3]), yax.index(best[4]), str(j),
                                 {'ha':'center', 'va':'center'})
        for i in range(len(zax), nrows * ncols):
            grid[i].set_visible(False)
        # This affects all axes because we set share_all = True.
        grid.axes_llc.set_ylim(-0.5, len(yax)+1)
        grid.axes_llc.set_xticks(range(0, len(xax), 2))
        grid.axes_llc.set_yticks(range(0, len(yax), 2))
        grid.axes_llc.set_xticklabels([my_round(i, 3) for i in xax][::2])
        grid.axes_llc.set_yticklabels([my_round(i, 3) for i in yax][::2])
        grid.axes_llc.set_ylabel(axes[2])
        grid.axes_llc.set_xlabel(axes[3])
    grid.cbar_axes[0].colorbar(im)
    grid.cbar_axes[0].set_ylabel('Correlation value')
    fig.suptitle(('Optimal IMP parameters\n' +
                  'Best for: {0}={4}, {1}={5}, {2}={6}, {3}={7}'
                  ).format(*(list(axes) + [my_round(i, 3)
                                           for i in sort_result[0][1:]])),
                 size='large')
    plt.show()
