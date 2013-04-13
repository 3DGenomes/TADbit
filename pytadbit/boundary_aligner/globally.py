"""
02 Dec 2012

global aligner for Topologically Associated Domains
"""
from math import log


def needleman_wunsch(tads1, tads2, penalty=-6., ext_pen=-5.6,
                     max_dist=500000, verbose=False):
    """
    Align two lists of TAD's boundaries.
    
    :param tads1: list of boundaries for one chromosome under one condition
    :param tads2: list of boundaries for the same chromosome under other
        conditions
    :param -0.1 penalty: penalty to open a gap in the alignment of boundaries
    :param 500000 max_dist: distance from which match are denied. A bin_size
        of 20Kb the number of bins corresponding to 0.5Mb is 25.
    :param False verbose: print the Needleman-Wunsch score matrix, and the
        alignment of boundaries

    :returns: the max score in the Needleman-Wunsch score matrix.
    """
    ##############
    tads1 = [0.0] + tads1
    tads2 = [0.0] + tads2
    l_tads1  = len(tads1)
    l_tads2  = len(tads2)
    dister = lambda x, y: log(1. / (abs(x - y) + 1))
    scores = virgin_score(penalty, l_tads1, l_tads2)
    pen = penalty
    for i in xrange(1, l_tads1):
        for j in xrange(1, l_tads2):
            d_dist = dister(tads2[j], tads1[i])
            match  = d_dist + scores[i-1][j-1]
            insert = scores[i-1][j] + pen
            delete = scores[i][j-1] + pen
            if d_dist > max_dist:
                scores[i][j] = max((insert, delete))
                pen = ext_pen
            else:
                if match >= max(insert, delete):
                    pen = ext_pen
                else:
                    pen = ext_pen
                scores[i][j] = max((match, insert, delete))
    align1 = []
    align2 = []
    i = l_tads1 -1
    j = l_tads2 -1
    max_score = None
    while i and j:
        score      = scores[i][j]
        if score > max_score:
            max_score = score
        d_dist     = dister(tads2[j], tads1[i])
        value      = scores[i-1][j-1] + d_dist
        if equal(score, value):
            align1.insert(0, tads1[i])
            align2.insert(0, tads2[j])
            i -= 1
            j -= 1
        elif equal(score, scores[i-1][j] + penalty):
            align1.insert(0, tads1[i])
            align2.insert(0, '-')
            i -= 1
        elif equal(score, scores[i-1][j] + ext_pen):
            align1.insert(0, tads1[i])
            align2.insert(0, '-')
            i -= 1
        elif equal(score, scores[i][j-1] + penalty):
            align1.insert(0, '-')
            align2.insert(0, tads2[j])
            j -= 1
        elif equal(score, scores[i][j-1] + ext_pen):
            align1.insert(0, '-')
            align2.insert(0, tads2[j])
            j -= 1
        else:
            for scr in scores: 
                print ' '.join(['%6s' % (round(y, 2)) for y in scr])
            raise Exception('Something  is failing and it is my fault...',
                            i, j, tads1[i], tads2[j])
    while i:
        align1.insert(0, tads1[i])
        align2.insert(0, '-')
        i -= 1
    while j:
        align1.insert(0, '-')
        align2.insert(0, tads2[j])
        j -= 1
        
    if verbose:
        print '\n Alignment:'
        print 'TADS 1: '+'|'.join(['%9s' % (str(int(x)) if x!='-' else '-'*3) \
                                   for x in align1])
        print 'TADS 2: '+'|'.join(['%9s' % (str(int(x)) if x!='-' else '-'*3) \
                                   for x in align2])
    return [align1, align2], max_score


def virgin_score(penalty, l_tads1, l_tads2):
    """
    create empty matrix
    """
    zeros    = [0.0 for _ in xrange(l_tads2)]
    return [[penalty * j for j in xrange(l_tads2)]] + \
           [[penalty * i] + zeros for i in xrange(1, l_tads1)]


def equal(a, b, cut_off=1e-9):
    """
    is equal?
    """
    return abs(a-b) < cut_off
