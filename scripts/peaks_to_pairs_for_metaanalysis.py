#!/usr/bin/env python

from argparse    import ArgumentParser
from itertools import combinations

import os

def main():
    opts   = get_options()
    reso   = opts.reso
    nbins  = opts.nbins
    left   = nbins
    right  = nbins + 1
    infile = os.path.realpath(opts.inbed)
    out = open(opts.outpairs, 'w')
    if opts.extra:
        split_elts = lambda w, x, y, z: (w, x, y, z)
        ncols  = 4
        ex1, ex2 = opts.extra
        test_extras = lambda x, y: x==ex1 and y== ex2
    else:
        split_elts = lambda w, x, y: (w, x, y, None)
        ncols  = 3
        test_extras = lambda x, y: True
    for elts1, elts2 in combinations((l.split()[:ncols]
                                      for l in open(infile)), 2):
        c1, b1, e1, s1 = split_elts(*elts1)
        c2, b2, e2, s2 = split_elts(*elts2)
        if c1 != c2:
            continue
        pos1 = (int(e1) + int(b1)) / 2
        pos2 = (int(e2) + int(b2)) / 2
        pos1, pos2 = sorted((pos1, pos2))
        beg1 = (pos1 / reso - left ) * reso
        if beg1 < 0:
            continue
        end1 = (pos1 / reso + right) * reso
        beg2 = (pos2 / reso - left ) * reso
        if beg2 < 0:
            continue
        end2 = (pos2 / reso + right) * reso
        if end1 > beg2:
            continue
        if test_extras(s1, s2):
            out.write('%s:%d-%d\t%s:%d-%d\n' % (c1, beg1, end1, c2, beg2, end2))
    out.close()


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbed', metavar='',
                        required=True, default=False,
                        help="input bed-like file with list of coordinates (e.g. peaks).")
    parser.add_argument('-o', '--outfile', dest='outpairs', metavar='',
                        required=True, default=False,
                        help="output list of pairs")
    parser.add_argument('-n','--nbins', dest='nbins', metavar='',
                        default=10, type=int,
                        help="Number of bins in left and right to consider")
    parser.add_argument('-r','--reso', dest='reso', metavar='',
                        required=True, type=int,
                        help="Number of bins in left and right to consider")
    parser.add_argument('--extra', dest='extra', metavar='',
                        type=str, nargs='+',
                        help="Use fourth column of the bed file as filter, each pair should match input")
    opts = parser.parse_args()

    return opts

if __name__ == '__main__':
    exit(main())
