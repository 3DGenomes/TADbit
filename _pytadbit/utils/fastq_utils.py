"""
12 nov. 2014
"""

from warnings import warn
from gzip import open as gopen
import numpy as np
from pytadbit.utils.extraviews import tadbit_savefig

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')

def quality_plot(fnam, nreads=None, axe=None, savefig=None):
    """
    Plot the qualities

    :param fnam: path to FASTQ file
    :param None nreads: max number of reads to read, not necesary to read all
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    """
    phred = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
    quals = []
    try:
        fhandler = open(fnam)
    except IOError:
        fhandler = gopen(fnam)
    if nreads:
        while True:
            try:
                next(fhandler)
            except EOFError:
                break
            next(fhandler)
            next(fhandler)
            line = next(fhandler)
            quals.append([phred.index(i) for i in line.strip()])
            if len(quals) > nreads:
                break
    else: # do this because it's faster
        while True:
            try:
                next(fhandler)
            except EOFError:
                break
            next(fhandler)
            next(fhandler)
            next(fhandler)
            line = next(fhandler)
            quals.append([phred.index(i) for i in line.strip()])
    fhandler.close()

    quals = zip(*quals)
    meanquals = [np.mean(q) for q in quals]
    errorquals = [np.std(q) for q in quals]

    if axe:
        ax = axe
        fig = axe.get_figure()
        plt.clf()
    else:
        fig = plt.figure()
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
    plt.figure(figsize=(15, 7))
    plt.errorbar(range(len(line.strip())), meanquals,
                 yerr=errorquals, ecolor='orange')

    plt.xlim((0, len(line)))
    plt.xlabel('Sequence')
    plt.ylabel('PHRED score')
    if savefig:
        tadbit_savefig(savefig)
    elif not axe:
        plt.show()

