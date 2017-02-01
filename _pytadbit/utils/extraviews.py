"""
06 Aug 2013


"""
from warnings import warn
import numpy as np
from subprocess import Popen, PIPE
from itertools import combinations

try:
    from matplotlib.ticker import MultipleLocator
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    warn('matplotlib not found\n')

NTH = {
    1 : "First",
    2 : "Second",
    3 : "Third",
    4 : "Fourth",
    5 : "Fifth",
    6 : "Sixth",
    7 : "Seventh",
    8 : "Eighth",
    9 : "Ninth",
    10: "Tenth",
    11: "Eleventh",
    12: "Twelfth"
}

def setup_plot(axe, figsize=None):
    if axe:
        ax = axe
        fig = ax.get_figure()
    else:
        fig = plt.figure(figsize=(11, 5) if not figsize else figsize)
        ax = fig.add_subplot(111)
        ax.patch.set_facecolor('lightgrey')
        ax.patch.set_alpha(0.4)
        ax.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
        ax.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
        ax.set_axisbelow(True)
        ax.minorticks_on() # always on, not only for log
        # remove tick marks
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False)
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False, which='minor')
    return ax

def tadbit_savefig(savefig):
    try:
        form = savefig[-4:].split('.')[1]
    except IndexError: # no dot in file name
        warn('WARNING: file extension not found saving in png')
        form = 'png'
    if not form in ['png', 'pdf', 'ps', 'eps', 'svg']:
        raise NotImplementedError('File extension must be one of %s' %(
            ['png', 'pdf', 'ps', 'eps', 'svg']))
    plt.savefig(savefig, format=form)

def nicer(res):
    """
    writes resolution number for human beings.
    """
    if not res % 1000000000:
        return str(res)[:-9] + ' Gb'
    if not res % 1000000:
        return str(res)[:-6] + ' Mb'
    if not res % 1000:
        return str(res)[:-3] + ' kb'
    if res == 1:
        return 'bin'
    return str(res) + ' b'


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

def color_residues(x, **kwargs):
    """
    Function to color residues from blue to red.
    
    :param model: a given :class:`pytadbit.imp.impmodel.IMPmodel`
    
    :returns: a list of rgb tuples (red, green, blue), each between 0 and 1.
    """
    result = []
    for n in xrange(len(x)):
        red = float(n + 1) / len(x)
        result.append((red, 0, 1 - red))
    return result


def tad_coloring(x, mstart=None, mend=None, tads=None, **kwargs):
    """
    Colors TADs from blue to red (first to last TAD). TAD borders are displayed
    in scale of grey, from light to dark grey (again first to last border)
    
    :param model: a given :class:`pytadbit.imp.impmodel.IMPmodel`
    :param tads: a dictionary of TADs, Experiments.tads can be directly passed
    
    :returns: a list of rgb tuples (red, green, blue), each between 0 and 1.
    """
    ltads = [t for t in tads if tads[t]['end'] > mstart
             and tads[t]['start'] < mend]
    ntads = len(ltads)
    grey = 0.95
    try:
        grey_step = (0.95 - 0.4) / ntads
    except ZeroDivisionError:
        raise Exception('ERROR: TAD borders are not predicted yet.')
    result = []
    for t in tads:
        if tads[t]['end'] < mstart or tads[t]['start'] > mend:
            continue
        red = float(t + 1 - min(ltads)) / ntads
        for _ in range(int(max(tads[t]['start'], mstart)),
                       int(min(tads[t]['end'], mend))):
            result.append((red, 0, 1 - red))
        if tads[t]['end'] <= mend:
               result.append((grey, grey, grey))
        grey -= grey_step
    return result


def tad_border_coloring(x, mstart=None, mend=None, tads=None, **kwargs):
    """
    Colors TAD borders from blue to red (bad to good score). TAD are displayed
    in scale of grey, from light to dark grey (first to last particle in the
    TAD)
    
    :param model: a given :class:`pytadbit.imp.impmodel.IMPmodel`
    :param tads: a dictionary of TADs, Experiments.tads can be directly passed
    
    :returns: a list of rgb tuples (red, green, blue), each between 0 and 1.
    """
    result = []
    if not tads:
        raise Exception('ERROR: TAD borders are not predicted yet.')
    for t in tads:
        if tads[t]['end'] < mstart or tads[t]['start'] > mend:
            continue
        grey = 0.95
        grey_step = (0.95 - 0.4) / (tads[t]['end'] - tads[t]['start'] + 1)
        red = float(tads[t]['score']) / 10
        for _ in range(int(max(tads[t]['start'], mstart)),
                       int(min(tads[t]['end'], mend))):
            result.append((grey, grey, grey))
            grey -= grey_step
        if tads[t]['end'] <= mend:
            result.append((red, 0, 1 - red))
    return result


def augmented_dendrogram(clust_count=None, dads=None, objfun=None, color=False,
                         axe=None, savefig=None, *args, **kwargs):
    """
    """
    from scipy.cluster.hierarchy import dendrogram
    fig = plt.figure(figsize=kwargs.get('figsize', (8,8)))
    if axe:
        ax = axe
        fig = axe.get_figure()
        ddata = dendrogram(*args)
        plt.clf()
    else:
        ddata = dendrogram(*args)
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
    debug = False
    leaves = {}
    dist = ddata['icoord'][0][2] - ddata['icoord'][0][1]
    for i, x in enumerate(ddata['leaves']):
        leaves[dist*i + dist/2] = x
    minnrj = min(objfun.values())
    maxnrj = max(objfun.values())
    difnrj = maxnrj - minnrj
    total = max(clust_count.values())
    if not kwargs.get('no_plot', False):
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'],
                           ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            # plt.plot(x, y, 'ro')
            plt.hlines(y, i[1], i[2], lw=2, color=(c if color else 'grey'))
            # for each branch
            for i1, d1, d2 in zip(i[1:3], [d[0], d[3]], [d[1], d[2]]):
                try:
                    lw = (fig.get_figwidth() *
                          kwargs.get('width_factor', 5.0) *
                          float(clust_count[leaves[i1] + 1]) / total)
                except KeyError:
                    lw = 1.0
                nrj = objfun[leaves[i1] + 1] if (leaves[i1] + 1) in objfun else maxnrj
                d1 = d1 - (difnrj - (nrj - minnrj))
                ax.vlines(i1, d1, d2, lw=lw, color=(c if color else 'grey'))
                if leaves[i1] + 1 in objfun or debug:
                    ax.annotate("%.3g" % (leaves[i1] + 1),
                                (i1, d1), size=kwargs.get('fontsize', 8),
                                xytext=(0, -8),
                                textcoords='offset points',
                                va='top', ha='center')
            leaves[(i[1] + i[2])/2] = dads[leaves[i[1]] + 1] - 1
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
               # "{:,}".format (int(minnrj)/cutter * cutter  + i)
               ["%s" % (form(int(minnrj)/cutter * cutter  + i))
                for i in xrange(0, -bot-bot/cut, -bot/cut)],
               size=kwargs.get('fontsize', 8))
    ax.set_ylabel('Minimum IMP objective function',
                  size=kwargs.get('fontsize', 8) + 1)
    ax.set_xticks([])
    ax.set_xlim((plt.xlim()[0] - 2, plt.xlim()[1] + 2))
    ax.figure.suptitle("Dendogram of clusters of 3D models",
                       size=kwargs.get('fontsize', 8) + 2)
    ax.set_title("Branch length proportional to model's objective function " +
                 "final value\n" +
                 "Branch width to the number of models in the cluster " +
                 "(relative to %s models)" % (total),
                 size=kwargs.get('fontsize', 8))
    if savefig:
        tadbit_savefig(savefig)
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
    bpAx.annotate('%.1f' % (bp['boxes'][0].get_xdata()[0]),
                  (bp['boxes'][0].get_xdata()[0], bp['boxes'][0].get_ydata()[1]),
                  va='bottom', ha='center', xytext=(0, 2),
                  textcoords='offset points',
                  size='small')
    bpAx.annotate('%.1f' % (bp['boxes'][0].get_xdata()[2]),
                  (bp['boxes'][0].get_xdata()[2], bp['boxes'][0].get_ydata()[1]),
                  va='bottom', ha='center', xytext=(0, 2),
                  textcoords='offset points',
                  size='small')
    bpAx.annotate('%.1f' % (bp['medians'][0].get_xdata()[0]),
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
        tadbit_savefig(savefig)
    elif not axe:
        plt.show()


def plot_3d_model(x, y, z, label=False, axe=None, thin=False, savefig=None,
                  show_axe=False, azimuth=-90, elevation=0., color='index',
                  **kwargs):
    """
    Given a 3 lists of coordinates (x, y, z) plots a three-dimentional model
    using matplotlib

    :param model: a :class:`pytadbit.imp.impmodel.IMPmodel`
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
    """
    show = False
    if isinstance(color, str):
        if color == 'index':
            color = color_residues(x, **kwargs)
        elif color == 'tad':
            if not 'tads' in kwargs:
                raise Exception('ERROR: missing TADs\n   ' +
                                'pass an Experiment.tads disctionary\n')
            color = tad_coloring(x, **kwargs)
        elif color == 'border':
            if not 'tads' in kwargs:
                raise Exception('ERROR: missing TADs\n   ' +
                                'pass an Experiment.tads disctionary\n')
            color = tad_border_coloring(x, **kwargs)
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
            axe.text(x[i+1], y[i+1], z[i+1],str(i+1), size=7)
        axe.scatter(x, y, z, color=color, s=50)
    axe.pbaspect = [1,1,1]
    if show:
        if savefig:
            tadbit_savefig(savefig)
        else:
            plt.show()
    
    
def chimera_view(cmm_files, chimera_bin='chimera', #shape='tube',
                 chimera_cmd=None, savefig=None, center_of_mass=False,
                 gyradius=False, align=True, grid=False, highlight='all',
                 **kwargs):
    """
    Open a list of .cmm files with Chimera (http://www.cgl.ucsf.edu/chimera)
    to view models.

    :param cmm_files: list of .cmm files
    :param 'chimera' chimera_bin: path to chimera binary
    :param None chimera_cmd: list of command lines for chimera
    :param None savefig: path to a file where to save generated image
    :param False center_of_mass: if True, draws the center of mass
    :param False gyradius: increment the radius of the center of mass by the
       value of the radius of girations given.
    :param False align: align models
    :param False grid: tile models
    :param 'all' higlight: higlights a given model, or group of models.
       Can be either 'all', or a model number
    """
    pref_f = '/tmp/tmp.cmd'
    out = open(pref_f, 'w')
    for cmm_file in cmm_files:
        out.write('open %s\n' % (cmm_file))
    nmodels = len(cmm_files)
    if nmodels > 1:
        for i in xrange(nmodels):
            if i == 0:
                continue
            if align:
                out.write('match #%s #%s\n' % (i, 0))
            if highlight != 'all':
                if highlight != i:
                    out.write('color black #%s\n' % (i))
    if not chimera_cmd:
        the_shape = '''bonddisplay never #%s
shape tube #%s radius 15 bandLength 300 segmentSubdivisions 1 followBonds on
~show #%s'''
        out.write(('''
focus
set bg_color white
windowsize 800 600
represent wire
%s
clip yon -500
~label
set subdivision 1
set depth_cue
set dc_color black
set dc_start 0.5
set dc_end 1
scale 0.8%s\n
''' % ('\n'.join([the_shape % (mdl, mdl, mdl) for mdl in (
                       [highlight] if highlight!='all' else range(nmodels))]),
       '\ntile' if grid else '')) +
                  ('define centroid radius %s color 1,0,0,0.2\n' % (
                      gyradius if gyradius else 10) if center_of_mass else '')
        + (kwargs.get('extra', '')))
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
    
    Popen('%s %s' % (chimera_bin, pref_f), shell=True)


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
                                                         axes[2], axes[3])
        plt.title(tit % tuple([my_round(r, 3) for r in sort_result[1:]]))
    plt.show()


def my_round(num, val):
    num = round(num, val)
    return int(num) if num == int(num) else num


def plot_2d_optimization_result(result, axes=('scale', 'maxdist', 'upfreq',
                                              'lowfreq'), dcutoff=None,
                                show_best=0, skip=None, savefig=None,clim=None):
    """
    A grid of heatmaps representing the result of the optimization.

    :param result: 3D numpy array contating correlation values
    :param 'scale','maxdist','upfreq','lowfreq' axes: tuple of axes to
       represent. The order will define which parameter will be placed on the
       w, z, y or x axe.
    :param 0 show_best: number of best correlation value to identifie.
    :param None skip: a dict can be passed here in order to fix a given axe,
       e.g.: {'scale': 0.001, 'maxdist': 500}
    :param None dcutoff: optimized cutoff
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param None clim: color scale. If None max and min values are used.
    
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
    
    if clim:
      vmin=clim[0]
      vmax=clim[1]
    
    wax = [my_round(i, 3) for i in axes_range[0]]
    zax = [my_round(i, 3) for i in axes_range[1]]
    xax = [my_round(i, 3) for i in axes_range[3]]
    yax = [my_round(i, 3) for i in axes_range[2]]
    # get best correlations
    sort_result =  sorted([(result[i, j, k, l],
                            wax[i], zax[j], xax[l], yax[k])
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
    width  = max(4, (float(ncols) * len(xax)) / 3)
    height = max(3, (float(nrows) * len(yax)) / 3)
    fig = plt.figure(figsize=(width, height))
    grid = AxesGrid(fig, [.2, .2, .6, .5],
                    nrows_ncols = (nrows + 1, ncols + 1),
                    axes_pad = 0.0,
                    label_mode = "1",
                    share_all = False,
                    cbar_location="right",
                    cbar_mode="single",
                    # cbar_size="%s%%" % (20./ width),
                    cbar_pad="15%",
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
                                    {'ha':'center', 'va':'center'}, size=8)
            if ii == wax_range[0]:
                rect = patches.Rectangle((-0.5, len(yax)-0.5), len(xax), 1.5,
                                         facecolor='grey', alpha=0.5)
                rect.set_clip_on(False)
                grid[cell].add_patch(rect)
                grid[cell].text(len(xax) / 2. - 0.5,
                                len(yax),
                                axes[1] + ' ' + str(my_round(zax[i], 3)),
                                {'ha':'center', 'va':'center'}, size=8)
            cell += 1
        rect = patches.Rectangle((len(xax)-.5, -0.5), 1.5, len(yax),
                                 facecolor='grey', alpha=0.5)
        rect.set_clip_on(False)
        grid[cell-1].add_patch(rect)
        grid[cell-1].text(len(xax)+.5, len(yax)/2.-.5, 
                          axes[0] + ' ' + str(my_round(wax[ii], 3)),
                          {'ha':'center', 'va':'center'}, 
                          rotation=90, size=8)
    for i in range(cell+1):
        if not i in used:
            grid[i].set_visible(False)
    # This affects all axes because we set share_all = True.
    # grid.axes_llc.set_ylim(-0.5, len(yax)+1)
    grid.axes_llc.set_xticks(range(0, len(xax), 2))
    grid.axes_llc.set_yticks(range(0, len(yax), 2))
    grid.axes_llc.set_xticklabels([my_round(i, 3) for i in xax][::2], size=9)
    grid.axes_llc.set_yticklabels([my_round(i, 3) for i in yax][::2], size=9)
    grid.axes_llc.set_ylabel(axes[2], size=9)
    grid.axes_llc.set_xlabel(axes[3], size=9)
    grid.cbar_axes[0].colorbar(im)
    grid.cbar_axes[0].set_ylabel('Correlation value', size=9)
    grid.cbar_axes[0].tick_params(labelsize=9)
    tit = 'Optimal IMP parameters\n'
    tit += 'Best: %s=%%s, %s=%%s,\n%s=%%s, %s=%%s %s=%%s' % (
        axes[0], axes[1], axes[3], axes[2], 'dcutoff')
    fig.suptitle(tit % tuple([my_round(i, 3) for i in sort_result[0][1:]] +
                             [str(dcutoff)]),
                 size=12)
    if savefig:
        tadbit_savefig(savefig)
    else:
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


def _tad_density_plot(xpr, maxys=None, fact_res=1., axe=None,
                     focus=None, extras=None, normalized=True,
                     savefig=None, shape='ellipse'):
    """
    """
    from matplotlib.cm import jet
    show=False
    if focus:
        siz = focus[1] - focus[0]
        figsiz = 4 + (focus[1] - focus[0])/30
        beg, end = focus
        tads = dict([(t, xpr.tads[t]) for t in xpr.tads
                     if (xpr.tads[t]['start'] + 1 >= beg
                         and xpr.tads[t]['end'] <= end)])
    else:
        siz = xpr.size
        figsiz = 4 + (siz)/30
        tads = xpr.tads

    if not axe:
        fig = plt.figure(figsize=(figsiz, 1 + 1 * 1.8))
        axe = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0)
        show=True

    zsin = np.sin(np.linspace(0, np.pi))
    
    shapes = {'ellipse'   : lambda h : [0] + list(h * zsin) + [0],
              'rectangle' : lambda h: [0] + [h] * 50 + [0],
              'triangle'  : lambda h: ([h/25 * i for i in xrange(26)] +
                                       [h/25 * i for i in xrange(25, -1, -1)])}
    
    try:
        shape = shapes[shape]
    except KeyError:
        import this
        table = ''.join([this.d.get(chr(i), chr(i)) for i in range(256)])
        if locals()['funcr'.translate(table)].translate(table) == ''.join(
            [this.s[i].upper() if this.s[i-1] is 'v' else this.s[i]
             for i in [24, 36, 163, 8, 6, 16, 36]]):
            shape = lambda h: (
                [h/25 * i for i in xrange(25)] + [h+0.2]*2 +
                [h/25 * i for i in xrange(24, -1, -1)])
        else:
            raise NotImplementedError(
                '%s not valid, use one of ellipse, rectangle or triangle')
    maxys = maxys if isinstance(maxys, list) else []
    zeros = xpr._zeros or {}
    if normalized and xpr.norm:
        norms = xpr.norm[0]
    elif xpr.hic_data:
        if normalized:
            warn("WARNING: weights not available, using raw data")
        norms = xpr.hic_data[0]
    else:
        warn("WARNING: raw Hi-C data not available, " +
             "TAD's height fixed to 1")
        norms = None
    if not 'height' in tads[tads.keys()[0]]:
        diags = []
        siz = xpr.size
        sp1 = siz + 1
        if norms:
            for k in xrange(1, siz):
                s_k = siz * k
                diags.append(sum([norms[i * sp1 + s_k]
                                 if not (i in zeros
                                         or (i + k) in zeros) else 0.
                                  for i in xrange(siz - k)]) / (siz - k))
        for tad in tads:
            start, end = (int(tads[tad]['start']) + 1,
                          int(tads[tad]['end']) + 1)
            if norms:
                matrix = sum([norms[i + siz * j]
                             if not (i in zeros
                                     or j in zeros) else 0.
                              for i in xrange(start - 1, end - 1)
                              for j in xrange(i + 1, end - 1)])
            try:
                if norms:
                    height = float(matrix) / sum(
                        [diags[i-1] * (end - start - i)
                         for i in xrange(1, end - start)])
                else:
                    height = 1.
            except ZeroDivisionError:
                height = 0.
            maxys.append(height)
            start = float(start) / fact_res  # facts[iex]
            end   = float(end) / fact_res  # facts[iex]
            axe.fill([start] + list(np.linspace(start, end)) + [end], shape(height),
                     alpha=.8 if height > 1 else 0.4,
                     facecolor='grey', edgecolor='grey')
    else:
        for tad in tads:
            start, end = (int(tads[tad]['start']) + 1,
                          int(tads[tad]['end']) + 1)
            height = float(tads[tad]['height'])
            maxys.append(height)
            axe.fill([start] + list(np.linspace(start, end)) + [end],
                     shape(height),
                     alpha=.8 if height > 1 else 0.4,
                     facecolor='grey', edgecolor='grey')
    if extras:
        axe.plot(extras, [.5 for _ in xrange(len(extras))], 'rx')
    axe.grid()
    axe.patch.set_visible(False)
    axe.set_ylabel('Relative\nHi-C count')
    #
    for tad in tads:
        if not tads[tad]['end']:
            continue
        tad = tads[tad]
        axe.plot(((tad['end'] + 1.) / fact_res, ), (0., ),
                 color=jet(tad['score'] / 10) if tad['score'] else 'w',
                 mec=jet(tad['score'] / 10) if tad['score'] else 'k',
                 marker=6, ms=9, alpha=1, clip_on=False)
    axe.set_xticks([1] + range(100, int(tad['end'] + 1), 50))
    axe.minorticks_on()
    axe.xaxis.set_minor_locator(MultipleLocator(10))
    axe.hlines(1, tads[tads.keys()[0]]['start'], end, 'k', lw=1.5)
    if show:
        tit1 = fig.suptitle("TAD borders", size='x-large')
        plt.subplots_adjust(top=0.76)
        fig.set_facecolor('white')
        plots = []
        for scr in xrange(1, 11):
            plots += plt.plot((100,),(100,), marker=6, ms=9,
                              color=jet(float(scr) / 10), mec='none')
        try:
            axe.legend(plots,
                       [str(scr) for scr in xrange(1, 11)],
                       numpoints=1, title='Border scores',
                       fontsize='small', loc='lower left',
                       bbox_to_anchor=(1, 0.1))
        except TypeError:
            axe.legend(plots,
                       [str(scr) for scr in xrange(1, 11)],
                       numpoints=1, title='Border scores',
                       loc='lower left',
                       bbox_to_anchor=(1, 0.1))
        axe.set_ylim((0, max(maxys) + 0.4))
        if savefig:
            tadbit_savefig(savefig)
        else:
            plt.show()
        
def plot_compartments(crm, first, cmprts, matrix, show, savefig,
                      vmin=-1, vmax=1, whichpc=1):
    heights = []
    val = 0
    for i in xrange(len(matrix)):
        try:
            val = [c['dens'] for c in cmprts[crm] if c['start']==i][0]
        except IndexError:
            pass
        heights.append(val-1)
    
    maxheights = max([abs(f) for f in heights])
    try:
        heights = [f / maxheights for f in heights]
    except ZeroDivisionError:
        warn('WARNING: no able to plot chromosome %s' % crm)
        return

    # definitions for the axes
    left, width = 0.1, 0.75
    bottom, height = 0.1, 0.75
    bottom_h = left + height + 0.02
    
    rect_scatter = [left               , bottom  , width, height]
    rect_histx   = [left               , bottom_h, width, 0.08 ]
    rect_histy   = [left + width + 0.02, bottom  , 0.02 , 0.75  ]
    
    # start with a rectangular Figure
    plt.figure(1, figsize=(14, 14))
    
    axim = plt.axes(rect_scatter)
    axex = plt.axes(rect_histx, sharex=axim)
    axey = plt.axes(rect_histy)
    axey.set_ylabel('density (orange)')

    im = axim.imshow(matrix, interpolation='nearest', cmap='coolwarm',
                     vmin=vmin, vmax=vmax)
    axim.minorticks_on()
    cb = plt.colorbar(im, cax=axey)
    cb.set_label('Pearson product-moment correlation coefficients')
    axim.grid()

    # scale first PC
    mfirst = max(max(first), abs(min(first)))
    first = [f/mfirst for f in first]
    
    axex.plot(first, color='green', alpha=0.5)
    axex.plot(heights, color='black', alpha=1, linewidth=2)
    axex.plot(heights, color='orange', alpha=1, linewidth=1)
    axex.fill_between(range(len(matrix)), [0] * len(matrix), first,
                      where=np.array(first) > 0, color='olive', alpha=0.5)
    axex.fill_between(range(len(matrix)), [0] * len(matrix), first,
                      where=np.array(first) < 0, color='darkgreen', alpha=0.5)
    axex.set_yticks([0])
    axex.set_ylabel('%s PC (green)\ndensity (orange)' % (NTH[whichpc]))
    breaks = [0] + [i + 0.5 for i, (a, b) in
                    enumerate(zip(first[1:], first[:-1]))
                    if a * b < 0] + [len(first)]
    # COMPARTMENTS A/B
    a_comp = []
    b_comp = []
    breaks = []
    for cmprt in cmprts[crm]:
        breaks.append(cmprt['start'])
        try:
            if cmprt['type'] == 'A':
                a_comp.append((cmprt['start'], cmprt['end']))
            elif cmprt['type'] == 'B':
                b_comp.append((cmprt['start'], cmprt['end']))
        except KeyError:
            if cmprt['dens'] > 1:
                a_comp.append((cmprt['start'], cmprt['end']))
            else:
                b_comp.append((cmprt['start'], cmprt['end']))            
    a_comp.sort()
    b_comp.sort()
    axex.hlines([0.05]*len(a_comp), [a[0] for a in a_comp],
                [a[1] for a in a_comp], color='red' , linewidth=6)
    axex.hlines([-0.05]*len(b_comp), [b[0] for b in b_comp],
                [b[1] for b in b_comp], color='blue' , linewidth=6)

    # axex.hlines([0]*(len(breaks)/2), breaks[ :-1:2], breaks[1::2],
    #             color='red' , linewidth=4, alpha=0.7)
    # axex.hlines([0]*((len(breaks) - 1)/2), breaks[1:-1:2], breaks[2::2],
    #             color='blue', linewidth=4, alpha=0.7)

    axex.grid()
    axex.minorticks_on()
    axex.grid(b=True, which='minor')
    plt.setp(axex.get_xticklabels(), visible=False)
    
    axex.set_xlim(0,len(matrix))
    axim.set_ylim(0,len(matrix))
    if show:
        plt.show()
    if savefig:
        tadbit_savefig(savefig)
        plt.close('all')

def plot_compartments_summary(crm, cmprts, show, savefig, title=None):

    plt.close('all')
    # start with a rectangular Figure
    a_comp = []
    b_comp = []
    breaks = []
    for cmprt in cmprts[crm]:
        breaks.append(cmprt['start'])
        try:
            if cmprt['type'] == 'A':
                a_comp.append((cmprt['start'], cmprt['end']))
            elif cmprt['type'] == 'B':
                b_comp.append((cmprt['start'], cmprt['end']))
        except KeyError:
            if cmprt['dens'] > 1:
                a_comp.append((cmprt['start'], cmprt['end']))
            else:
                b_comp.append((cmprt['start'], cmprt['end']))            
    a_comp.sort()
    b_comp.sort()
    fig, ax = plt.subplots(figsize=(3 + cmprt['end'] / 100., 2))
    plt.subplots_adjust(top=0.7, bottom=0.25)
    plt.hlines([0.05]*len(a_comp), [a[0] for a in a_comp],
               [a[1] for a in a_comp], color='red' , linewidth=6)
    plt.hlines([-0.05]*len(b_comp), [b[0] for b in b_comp],
               [b[1] for b in b_comp], color='blue' , linewidth=6)
    plt.vlines(breaks, [-0.3]*len(breaks), [0.3]*len(breaks), color='black',
               linestyle=':')
    plt.title(title if title else ('Chromosome %s' % crm))
    plt.text(1, 0.55, 'A compartments')
    plt.text(1, -0.75, 'B compartments')
    plt.xlabel('Genomic bin')
    plt.ylim((-1, 1))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.set_yticks([])
    plt.xlim((0, cmprts[crm][-1]['end']))
    for item in [fig, ax]:
        item.patch.set_visible(False)
    if show:
        plt.show()
    if savefig:
        tadbit_savefig(savefig)
        plt.close('all')
