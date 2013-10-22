"""
06 Aug 2013


"""
from warnings import warn
import numpy as np
from subprocess import Popen, PIPE
from itertools import combinations

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
    return '%s%s%s' % (color[num], string,
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
    minnrj = min(objfun.values())
    maxnrj = max(objfun.values())
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
                    lw = float(clust_count[leaves[i1] + 1])/total*10*len(leaves)
                except KeyError:
                    lw = 1.0
                nrj = objfun[leaves[i1] + 1] if (leaves[i1] + 1) in objfun else maxnrj
                ax.vlines(i1, d1-(difnrj-(nrj-minnrj)), d2, lw=lw,
                          color=(c if color else 'grey'))
                if leaves[i1] + 1 in objfun:
                    ax.annotate("%.3g" % (leaves[i1] + 1),
                                (i1, d1-(difnrj-(nrj-minnrj))),
                                xytext=(0, -8),
                                textcoords='offset points',
                                va='top', ha='center')
            leaves[(i[1] + i[2])/2] = dads[leaves[i[1]] + 1]
    try:
        cutter = 10**int(np.log10(difnrj))
    except OverflowError: # case that the two are exactly the same
        cutter = 1
    cut = 10 if cutter >= 10 else 1
    bot = (-int(difnrj)/cutter * cutter) or -1 # do not want this to be null
    # just to display nice numbers
    form = lambda x: ''.join([(s + ',') if not i%3 and i else s
                              for i, s in enumerate(str(x)[::-1])][::-1])
    plt.yticks([bot+i for i in xrange(0, -bot-bot/cut, -bot/cut)],
               # ["{:,}".format (int(minnrj)/cutter * cutter  + i)
               ["%s" % (form(int(minnrj)/cutter * cutter  + i))
                for i in xrange(0, -bot-bot/cut, -bot/cut)], size='small')
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
    bpAx.annotate('%.4f' % (bp['boxes'][0].get_xdata()[0]),
                  (bp['boxes'][0].get_xdata()[0], bp['boxes'][0].get_ydata()[1]),
                  va='bottom', ha='center', xytext=(0, 2),
                  textcoords='offset points',
                  size='small')
    bpAx.annotate('%.4f' % (bp['boxes'][0].get_xdata()[2]),
                  (bp['boxes'][0].get_xdata()[2], bp['boxes'][0].get_ydata()[1]),
                  va='bottom', ha='center', xytext=(0, 2),
                  textcoords='offset points',
                  size='small')
    bpAx.annotate('%.4f' % (bp['medians'][0].get_xdata()[0]),
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
    bpAx.set_title('Histogram and boxplot of distances between particles ' +
                   '%s and %s' % (part1, part2))
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
    out.write('open %s\n' % (cmm_file))
    if not chimera_cmd:
        out.write('''
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
''')
        if savefig:
            if savefig.endswith('.png'):
                out.write('copy file %s png' % (savefig))
            elif savefig[-4:] in ('.mov', 'webm'):
                out.write('''
movie record supersample 1
turn y 3 120
wait 120
movie stop
movie encode output %s
''' % (savefig))
            elif savefig:
                raise Exception('Not supported format, must be png, mov or webm\n')
    else:
        out.write('\n'.join(chimera_cmd) + '\n')
    out.close()
    
    Popen('%s %s' % (chimera_bin, pref_f), shell=True).communicate()


def plot_3d_optimization_result(result,
                                axes=('scale', 'maxdist', 'upfreq', 'lowfreq')):
    """
    Displays a three dimensional scatter plot representing the result of the
    optimization.

    :param result: 3D numpy array contating correlation values
    :param 'scale','maxdist','upfreq','lowfreq' axes: tuple of axes to
       represent. The order will define which parameter will be placed on the
       w, z, y or x axe.
    """

    ori_axes, axes_range, result = result
    trans = [ori_axes.index(a) for a in axes]
    axes_range = [axes_range[i] for i in trans]

    # transpose results
    result = result.transpose(trans)

    wax = [my_round(i, 3) for i in axes_range[0]]
    zax = [my_round(i, 3) for i in axes_range[1]]
    xax = [my_round(i, 3) for i in axes_range[3]]
    yax = [my_round(i, 3) for i in axes_range[2]]
    sort_result = sorted([(result[i, j, k, l], wax[i], zax[j], xax[l], yax[k])
                          for i in range(len(wax))
                          for j in range(len(zax))
                          for k in range(len(yax))
                          for l in range(len(xax))
                          if not np.isnan(result[i, j, k, l])
                          ], key=lambda x: x[0],
                         reverse=True)[0]
    x = [i for i in axes_range[1] for j in axes_range[2] for k in axes_range[3]]
    y = [j for i in axes_range[1] for j in axes_range[2] for k in axes_range[3]]
    z = [k for i in axes_range[1] for j in axes_range[2] for k in axes_range[3]]
    
    from mpl_toolkits.mplot3d import Axes3D

    ncols  = int(np.sqrt(len(wax)) + 0.999)
    nrows  = int(np.sqrt(len(wax)) + 0.5)
    fig = plt.figure(figsize=((ncols)*6,(nrows)*4.5))

    for i in xrange(len(wax)):
        col = [result[i, j, k, l] for j in range(len(axes_range[1]))
               for k in range(len(axes_range[2])) for l in range(len(axes_range[3]))]

        ax = fig.add_subplot(int(str(nrows) + str(ncols) + str(i)),
                             projection='3d')
        ax.set_xlabel(axes[1])
        ax.set_ylabel(axes[2])
        ax.set_zlabel(axes[3])
        lol = ax.scatter(x, y, z, c=col, s=100, alpha=0.9)
        cbar = fig.colorbar(lol)
        cbar.ax.set_ylabel('Correlation value')
        tit = 'Optimal IMP parameters (subplot %s=%s)\n' % (axes[0], wax[i])
        tit += 'Best: %s=%%s, %s=%%s, %s=%%s, %s=%%s' % (axes[0], axes[1],
                                                         axes[3], axes[4])
        plt.title(tit % tuple([my_round(r, 3) for r in sort_result[1:]]))
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
    from matplotlib.cm import jet

    ori_axes, axes_range, result = result
    trans = [ori_axes.index(a) for a in axes]
    axes_range = [axes_range[i] for i in trans]

    # transpose results
    result = result.transpose(trans)

    # set NaNs
    result = np.ma.array(result, mask=np.isnan(result))
    cmap = jet
    cmap.set_bad('w', 1.)

    # defines axes
    vmin = result.min()
    vmax = result.max()
    wax = [my_round(i, 3) for i in axes_range[0]]
    zax = [my_round(i, 3) for i in axes_range[1]]
    xax = [my_round(i, 3) for i in axes_range[3]]
    yax = [my_round(i, 3) for i in axes_range[2]]

    # get best correlations
    sort_result =  sorted([(result[i, j, k, l], wax[i], zax[j], xax[l], yax[k])
                           for i in range(len(wax))
                           for j in range(len(zax))
                           for k in range(len(yax))
                           for l in range(len(xax))
                           if str(result[i, j, k, l]) != '--'],
                          key=lambda x: x[0],
                          reverse=True)[:show_best+1]

    # skip axes?
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
            raise Exception(('ERROR: skip keys must be one of the two first' +
                             ' keywords passed as axes parameter'))

    # best number of rows/columns
    ncols  = len(zax_range)
    nrows  = len(wax_range)
    fig = plt.figure(figsize=(max(6, float(ncols) * len(xax) / 3),
                              max(6, float(nrows) * len(yax) / 3)))
    grid = AxesGrid(fig, [.1,.1,.9,.75],
                    nrows_ncols = (nrows+1, ncols+1),
                    axes_pad = 0.0,
                    label_mode = "1",
                    share_all = False,
                    cbar_location="right",
                    cbar_mode="single",
                    cbar_size="%s%%" % (7./(float(ncols) * len(xax) / 3)),
                    cbar_pad="10%",
                    )
    cell = ncols
    used = []
    for ii in wax_range:
        cell+=1
        for i in zax_range:
            used.append(cell)
            im = grid[cell].imshow(result[ii, i, :, :], interpolation="nearest",
                                   origin='lower', vmin=vmin, vmax=vmax,
                                   cmap=cmap)
            grid[cell].tick_params(axis='both', direction='out', top=False,
                                   right=False, left=False, bottom=False)
            for j, best  in enumerate(sort_result[:-1]):
                if best[2] == zax[i] and best[1] == wax[ii]:
                    grid[cell].text(xax.index(best[3]), yax.index(best[4]), str(j),
                                    {'ha':'center', 'va':'center'})
            if ii == wax_range[0]:
                rect = patches.Rectangle((-0.5, len(yax)-0.5),len(xax), 1.5,
                                         facecolor='grey', alpha=0.5)
                rect.set_clip_on(False)
                grid[cell].add_patch(rect)
                grid[cell].text(len(xax) / 2 - 0.5,
                                len(yax),
                                axes[1] + ' ' + str(my_round(zax[i], 3)),
                                {'ha':'center', 'va':'center'})
            cell += 1
        rect = patches.Rectangle((len(xax)-.5, -0.5), 1.5, len(yax),
                                 facecolor='grey', alpha=0.5)
        rect.set_clip_on(False)
        grid[cell-1].add_patch(rect)
        grid[cell-1].text(len(xax)+.5, len(yax)/2-.5, 
                          axes[0] + ' ' + str(my_round(wax[ii], 3)),
                          {'ha':'right', 'va':'center'}, 
                          rotation=90)
    for i in range(cell+1):
        if not i in used:
            grid[i].set_visible(False)
    # This affects all axes because we set share_all = True.
    # grid.axes_llc.set_ylim(-0.5, len(yax)+1)
    grid.axes_llc.set_xticks(range(0, len(xax), 2))
    grid.axes_llc.set_yticks(range(0, len(yax), 2))
    grid.axes_llc.set_xticklabels([my_round(i, 3) for i in xax][::2])
    grid.axes_llc.set_yticklabels([my_round(i, 3) for i in yax][::2])
    grid.axes_llc.set_ylabel(axes[2])
    grid.axes_llc.set_xlabel(axes[3])
    grid.cbar_axes[0].colorbar(im)
    grid.cbar_axes[0].set_ylabel('Correlation value')
    tit = 'Optimal IMP parameters\n'
    tit += 'Best: %s=%%s, %s=%%s, %s=%%s, %s=%%s' % (axes[0], axes[1],
                                                     axes[3], axes[2])
    fig.suptitle(tit % tuple([my_round(i, 3) for i in sort_result[0][1:]]),
                 size='large')
    plt.show()


def compare_models(sm1, sm2, cutoff=150,
                   models1=None, cluster1=None,
                   models2=None, cluster2=None):
    """
    Plots the difference of contact maps of two group of structural models.
    
    :param sm1: a StructuralModel
    :param sm2: a StructuralModel
    :param 150 dcutoff: distance threshold (nm) to determine if two
       particles are in contact
    :param None models: if None (default) the contact map will be computed
       using all the models. A list of numbers corresponding to a given set
       of models can be passed
    :param None cluster: compute the contact map only for the models in the
       cluster number 'cluster'       
    """
    mtx1 = sm1.get_contact_matrix(models=models1, cluster=cluster1, cutoff=cutoff)
    mtx2 = sm2.get_contact_matrix(models=models2, cluster=cluster2, cutoff=cutoff)
    mtx3 = [[mtx2[i][j] - mtx1[i][j]
             for j in xrange(len(mtx1))]
            for i in xrange(len(mtx1))]
    fig = plt.figure(figsize=(8, 6))
    axe = fig.add_subplot(111)
    im = axe.imshow(mtx3, origin='lower', interpolation="nearest")
    axe.set_ylabel('Particle')
    axe.set_xlabel('Particle')
    cbar = axe.figure.colorbar(im)
    cbar.ax.set_ylabel('Signed log difference between models')
    plt.show()
