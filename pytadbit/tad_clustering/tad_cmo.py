"""
23 Jan 2013


"""

import numpy as np
from itertools import product

def core_nw(p_scores, penalty, l_p1, l_p2):
    scores = virgin_score(penalty, l_p1 + 1, l_p2 + 1)
    for i in xrange(1, l_p1 + 1):
        for j in xrange(1, l_p2 + 1):
            match  = p_scores[i - 1][j - 1] + scores[i-1][j-1]
            insert = scores[i-1][j] + penalty
            delete = scores[i][j-1] + penalty
            scores[i][j] = max((match, insert, delete))
    align1 = []
    align2 = []
    i = l_p1 
    j = l_p2 
    while i and j:
        score      = scores[i][j]
        value = scores[i-1][j-1] + p_scores[i - 1][j - 1]
        if equal(score, value):
            align1.insert(0, i - 1)
            align2.insert(0, j - 1)
            i -= 1
            j -= 1
        elif equal(score, scores[i-1][j] + penalty):
            align1.insert(0, i - 1)
            align2.insert(0, '-')
            i -= 1
        elif equal(score, scores[i][j-1] + penalty):
            align1.insert(0, '-')
            align2.insert(0, j - 1)
            j -= 1
        else:
            for x in scores:
                print ' '.join(['%7s' % (round(y, 4)) for y in x])
            raise Exception('Something  is failling and it is my fault...')
    return align1, align2

def optimal_cmo(tad1, tad2, num_v=None):
    """
    :argument p1: first matrix to align
    :argument p2: second matrix to align
    :argument None num_v: number of eigen vectors to consider, max is:\
                          max(min(len(p1), len(p2)))

    :return: 2 lists, one per aligned matrix
        
    NOTE: np.linalg.eig returns eigenvalues/eigenvectors sorted by\
          their absolute values
    """

    l_p1 = len(tad1)
    l_p2 = len(tad2)
    val1, vec1 = np.linalg.eig(tad1)
    val2, vec2 = np.linalg.eig(tad2)
    #
    val1 = np.array([np.sqrt(abs(v)) for v in val1])
    val2 = np.array([np.sqrt(abs(v)) for v in val2])
    idx = val1.argsort()[::-1]
    val1 = val1[idx]
    vec1 = vec1[idx]
    idx = val2.argsort()[::-1]
    val2 = val2[idx]
    vec2 = vec2[idx]
    #
    vec1 = np.array([val1[i] * vec1[:,i] for i in xrange(num_v)]).transpose()
    vec2 = np.array([val2[i] * vec2[:,i] for i in xrange(num_v)]).transpose()
    best_sc = 0
    best_alis = []
    for num in xrange(1, num_v + 1):
        for factors in product([1,-1], repeat=num):
            vec1p = factors * vec1[:,:num]
            vec2p = vec2[:,:num]
            p_scores = prescoring(vec1p, vec2p, num, l_p1, l_p2)
            penalty = min([min(s) for s in p_scores] + [0.0])
            align1, align2 = core_nw(p_scores, penalty, l_p1, l_p2)
            score = get_score(align1, align2, tad1, tad2)
            if score > best_sc:
                best_sc = score
                best_alis = [align1, align2]
    try:
        align1, align2 = best_alis
    except ValueError:
        pass
    print_score(align1, align2, tad1, tad2)
    print '\n Alignment (score = {}):'.format(best_sc)
    print 'TADS 1: '+'|'.join(['%4s' % (str(int(x)) \
                                        if x!='-' else '-'*3) for x in align1])
    print 'TADS 2: '+'|'.join(['%4s' % (str(int(x)) \
                                        if x!='-' else '-'*3) for x in align2])

    return align1, align2
    

def virgin_score(penalty, l_p1, l_p2):
    zeros    = [0.0 for _ in xrange(l_p2)]
    return [[penalty * j for j in xrange(l_p2)]] + \
           [[penalty * i] + zeros for i in xrange(1, l_p1)]


def equal(a, b, cut_off=1e-9):
    return abs(a-b) < cut_off


def prescoring(vec1p, vec2p, num, l_p1, l_p2):
    """
    yes... this is the bottle neck, half of the time spent here
    """
    return [[sum(vec1p[i] * vec2p[j]) for j in xrange(l_p2)] for i in xrange(l_p1)]


def get_score(align1, align2, tad1, tad2):
    """
    using Spearman Rho value
    """
    map1 = []
    map2 = []
    for i, j in zip(align1, align2):
        if j != '-' and i != '-':
            map1.append(i)
            map2.append(j)
    pp1 = [tad1[i][j] for i in map1 for j in map1]
    pp2 = [tad2[i][j] for i in map2 for j in map2]
    return np.corrcoef(pp1, pp2, rowvar=0)[1, 0]


def get_OLD_score(align1, align2, p1, p2):
    map1 = []
    map2 = []
    for i, j in zip(align1, align2):
        if j != '-' and i != '-':
            map1.append(i)
            map2.append(j)
    com = len(map1)
    pp1 = [[p1[i][j] for j in map1] for i in map1]
    pp2 = [[p2[i][j] for j in map2] for i in map2]
    cm1 = sum([pp1[i][j] for i in xrange(com) for j in xrange(i + 1, com)])
    cm2 = sum([pp2[i][j] for i in xrange(com) for j in xrange(i + 1, com)])
    cmo = sum([1 - abs(pp2[i][j] - pp1[i][j]) \
               for i in xrange(com) for j in xrange(i + 1, com)])
    return float(2 * cmo)/(cm1+cm2)


def prescoring(vec1p, vec2p, num, l_p1, l_p2):
    """
    yes... this is the bottle neck, half of the time spent here
    """
    return [[sum(vec1p[i] * vec2p[j]) for j in xrange(l_p2)] for i in xrange(l_p1)]


def main():
    """
    main function
    """
    pass


if __name__ == "__main__":
    exit(main())
