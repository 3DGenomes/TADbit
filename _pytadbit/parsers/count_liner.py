"""
18 may 2015
"""

from random import random
from os import SEEK_END
from numpy import std, mean

def wc(fnam):
    return sum(1 for _ in open(fnam))

def wc_approx(fnam, samples=1000, verbose=True):
    """
    Get the approximate number of reads in a FASTQ file. By averaging the sizes
    of a given sample od randomly selected reads, and relating this mean to the
    size of the file.

    :param fnam: path to FASTQ file
    :param 1000 samples: number of reads to sample. 1000 generally gives an
       accuracy bellow 0.1%
    :param True verbose: prints Number of reads and accuracy

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
