"""
12 nov. 2014
"""

from warnings import warn
from gzip import open as gopen
import numpy as np
from pytadbit.utils.extraviews import tadbit_savefig
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, religated, repaired
from os import SEEK_END
from numpy import std, mean
from random import random
import re

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')


def quality_plot(fnam, r_enz=None, nreads=None, axe=None, savefig=None, paired=False):
    """
    Plots the sequencing quality of a given FASTQ file. If a restrinction enzyme
    (RE) name is provided, can also represent the distribution of digested and
    undigested RE sites and estimate an expected proportion of dangling-ends.

    Proportion of dangling-ends is inferred by counting the number of times a
    dangling-end site, is found at the beginning of any of the reads (divided by
    the number of reads).

    :param fnam: path to FASTQ file
    :param None nreads: max number of reads to read, not necesary to read all
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param False paired: is input FASTQ contains both ends

    :returns: the percentage of dangling-ends (sensu stricto) and the percentage of
       reads with at least a ligation site.
    """
    phred = dict([(c, i) for i, c in enumerate(
        '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~')])
    quals = []
    henes = []
    sites = []
    fixes = []
    liges = []
    ligep = 0
    tkw = dict(size=4, width=1.5)
    if fnam.endswith('.gz'):
        fhandler = gopen(fnam)
    else:
        fhandler = open(fnam)
    if not r_enz:
        if nreads:
            while True:
                try:
                    next(fhandler)
                except EOFError:
                    break
                seq = next(fhandler)
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
                if len(quals) > nreads:
                    break
        else: # do this because it's faster
            while True:
                try:
                    next(fhandler)
                except EOFError:
                    break
                seq = next(fhandler)
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
    else:
        r_site = RESTRICTION_ENZYMES[r_enz].replace('|', '')
        l_site = religated(r_enz)
        d_site = repaired(r_enz)
        if r_site*2 == l_site:
            # in case the religated site equals 2 restriction sites (like DnpII)
            site = re.compile('(?<!%s)' % r_site + r_site + '(?!%s)' % r_site)
            fixe = re.compile('(?<!%s)' % d_site + d_site + '(?!%s)' % d_site)
        else:
            site = re.compile(r_site)
            fixe = re.compile(d_site)
        lige = re.compile(l_site)
        if nreads:
            while True:
                try:
                    next(fhandler)
                except StopIteration:
                    break
                seq = next(fhandler)
                sites.extend([m.start() for m in site.finditer(seq)])
                fixes.extend([m.start() for m in fixe.finditer(seq)])
                liges.extend([m.start() for m in lige.finditer(seq)])
                ligep += l_site in seq
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
                if len(quals) > nreads:
                    break
        else: # do this because it's faster
            while True:
                try:
                    next(fhandler)
                except StopIteration:
                    break
                seq = next(fhandler)
                sites.extend([m.start() for m in site.finditer(seq)])
                fixes.extend([m.start() for m in fixe.finditer(seq)])
                liges.extend([m.start() for m in lige.finditer(seq)])
                ligep += l_site in seq
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
    fhandler.close()
    if not nreads:
        nreads = len(quals)
    quals = zip(*quals)
    meanquals = [np.mean(q) for q in quals]
    errorquals = [np.std(q) for q in quals]

    if axe:
        ax = axe
        fig = axe.get_figure()
        ax2 = fig.add_subplot(212)
    else:
        if r_enz:
            _, (ax, ax2) = plt.subplots(2,1, figsize=(15, 12))
        else:
            _, ax = plt.subplots(1,1, figsize=(15, 6))
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
    ax.errorbar(range(len(line.strip())), meanquals,
                linewidth=1, elinewidth=1, color='darkblue',
                yerr=errorquals, ecolor='orange')

    ax.set_xlim((0, len(line)))
    ax.set_xlabel('Nucleotidic position')
    ax.set_ylabel('PHRED score')
    ax.set_title('Sequencing Quality (%d reads)' % (nreads))
    ax.yaxis.label.set_color('darkblue')
    ax.tick_params(axis='y', colors='darkblue', **tkw)
    axb = ax.twinx()
    axb.plot([henes.count(i) for i in xrange(len(line))], linewidth=1,
             color='black', linestyle='--')
    axb.yaxis.label.set_color('black')
    axb.tick_params(axis='y', colors='black', **tkw)
    axb.set_ylabel('Number of "N" per position')
    try: # no Ns found (yes... it happens)
        axb.set_yscale('log')
        axb.set_ylim((0, axb.get_ylim()[1] * 1000))
    except ValueError:
        axb.set_yscale('linear')
    ax.set_ylim((0, ax.get_ylim()[1]))
    ax.set_xlim((0, len(line)))

    if r_enz:
        ax.set_title('Sequencing Quality and deconvolution (%s %d reads)' % (
            r_enz, nreads))
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)
        ax2.patch.set_facecolor('lightgrey')
        ax2.patch.set_alpha(0.4)
        ax2.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
        ax2.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
        ax2.set_axisbelow(True)
        ax2.set_xlabel('Nucleotidic position')
        seq_len = len(line) - max((len(r_site), len(l_site), len(d_site)))
        sites = [sites.count(k) for k in xrange(seq_len)] # Undigested
        liges = [liges.count(k) for k in xrange(seq_len)] # OK
        fixes = [fixes.count(k) for k in xrange(seq_len)] # DE
        if d_site in r_site:
            pos = r_site.find(d_site)
            fixes = (fixes[:pos] +
                     [fixes[k] - sites[k-pos] for k in xrange(pos, seq_len)])
        if d_site in l_site:
            pos = l_site.find(d_site)
            fixes = (fixes[:pos] +
                     [fixes[k] - liges[k-pos] for k in xrange(pos, seq_len)])
        site_len = max((len(r_site), len(l_site), len(d_site)))
        if paired:
            sites[len(line) / 2 - site_len:
                  len(line) / 2] = [float('nan')] * site_len
            liges[len(line) / 2 - site_len:
                  len(line) / 2] = [float('nan')] * site_len
            fixes[len(line) / 2 - site_len:
                  len(line) / 2] = [float('nan')] * site_len
        ax2.plot(sites, linewidth=2, color='darkred')
        ax2.set_ylabel('Undigested RE site (%s)' % r_site)
        ax2.yaxis.label.set_color('darkred')
        ax2.tick_params(axis='y', colors='darkred', **tkw)
        ax3 = ax2.twinx()
        ax3.plot(liges, linewidth=2, color='darkblue')
        ax3.yaxis.label.set_color('darkblue')
        ax3.tick_params(axis='y', colors='darkblue', **tkw)
        ax3.set_ylabel('Religated (%s)' % l_site)
        if any([f > 0 for f in fixes]):
            ax4 = ax2.twinx()
            ax4.spines["right"].set_position(("axes", 1.07))
            make_patch_spines_invisible(ax4)
            ax4.spines["right"].set_visible(True)        
            ax4.plot(fixes, linewidth=2, color='darkorange')
            ax4.yaxis.label.set_color('darkorange')
            ax4.tick_params(axis='y', colors='darkorange', **tkw)
            ax4.set_ylabel('Dangling-ends (%s)' % d_site)
        else:
            ax2.set_ylabel('RE site & Dangling-ends  (%s)' % r_site)
        ax2.set_xlim((0, len(line)))
        lig_cnt = (np.nansum(liges) - liges[0] - liges[len(line) / 2])
        sit_cnt = (np.nansum(sites) - sites[0] - sites[len(line) / 2])
        des = ((100. * (fixes[0] + (fixes[(len(line) / 2)]
                                            if paired else 0)))
                       / nreads) if any([f > 0 for f in fixes]) else (
            100. * (sites[0] + (sites[(len(line) / 2)] if paired else 0))) / nreads
        plt.title(('Percentage of digested sites: %.0f%%, of dangling-ends: %.0f%%\n' +
                   'Percentage of reads with ligation site: %.0f%%') %(
                      (100. * lig_cnt) / (lig_cnt + sit_cnt),
                      des,
                      (ligep * 100.) / nreads))
        plt.subplots_adjust(right=0.85)
    if savefig:
        tadbit_savefig(savefig)
        plt.close('all')
    elif not axe:
        plt.show()
    return des, (ligep * 100.) / nreads


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)


def count_reads(fnam):
    """
    Count the number of reads in a FASTQ file (can be slow on big files, try
    count_reads_approx)

    :param fnam: path to file

    :returns: the number of reads (number of lines divided by four)
    """
    nlines = sum(1 for _ in open(fnam))
    if nlines % 4:
        raise IOError('ERROR: Number of lines not multiple of four\n')
    return nlines / 4


def count_reads_approx(fnam, samples=1000, verbose=True):
    """
    Get the approximate number of reads in a FASTQ file. By averaging the sizes
    of a given sample od randomly selected reads, and relating this mean to the
    size of the file.

    :param fnam: path to FASTQ file
    :param 1000 samples: number of reads to sample. 1000 generally gives an
       accuracy bellow 0.1%
    :param True verbose: prints Number of reads and accuracy (based on standard
       error of the mean)

    :returns: number of reads estimated
    """
    fhandler = open(fnam)
    fhandler.seek(0, SEEK_END)
    flen = fhandler.tell()
    values = []
    def _read_size(rnd):
        fhandler.seek(rnd)
        while True:
            line = fhandler.next()
            if line.startswith('@'):
                line2 = fhandler.next()
                if not line2.startswith('@'):
                    break
        return len(line) + 2 * len(line2) +  len(fhandler.next())
    for _ in xrange(samples):
        rnd = int(random() * flen)
        try:
            values.append(_read_size(rnd))
        except StopIteration:
            samples-=1
    mean_len = float(mean(values))
    nreads = flen / mean_len
    
    if verbose:
        dev = std(values) / samples**.5 * 2
        nreads_sup = flen / (mean_len - dev)
        nreads_bel = flen / (mean_len + dev)
        # print nreads_sup > 186168911 > nreads_bel, ':',
        print ' %d reads -- 95%% between %d and %d (%f %% accuracy)' % (
            int(nreads), int(nreads_sup), int(nreads_bel),
            (nreads_sup - nreads_bel) / nreads * 100)
        # print nreads, '- 186168911 = ',
        # print int(nreads) - 186168911, '(',
        # print abs(nreads - 186168911.00000) / nreads * 100, '% )'
    return int(nreads)
