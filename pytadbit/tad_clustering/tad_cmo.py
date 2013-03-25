"""
23 Jan 2013

functions to align contact maps.
Algorithm based on:
Di Lena, P., Fariselli, P., Margara, L., Vassura, M., & Casadio, R. (2010). 
Fast overlapping of protein contact maps by alignment of eigenvectors. 
Bioinformatics (Oxford, England), 26(18), 2250-8. doi:10.1093/bioinformatics/btq402

"""

from numpy import array, sqrt#, corrcoef
from numpy import min as npmin
from numpy import sum as npsum
from numpy.linalg import eig
from scipy.stats import spearmanr
from itertools import product, combinations

# for aleigen:
from numpy import median

import re
from subprocess import Popen, PIPE


def core_nw(p_scores, penalty, l_p1, l_p2):
    """
    Core of the Needleman-Wunsch algorithm 
    """
    scores = virgin_score(penalty, l_p1 + 1, l_p2 + 1)
    for i in xrange(1, l_p1 + 1):
        for j in xrange(1, l_p2 + 1):
            match  = p_scores[i - 1][j - 1] + scores[i - 1][j - 1]
            insert = scores[i - 1][j] + penalty
            delete = scores[i][j - 1] + penalty
            scores[i][j] = max((match, insert, delete))
    align1 = []
    align2 = []
    i = l_p1 
    j = l_p2 
    while i and j:
        score = scores[i][j]
        value = scores[i - 1][j - 1] + p_scores[i - 1][j - 1]
        if equal(score, value):
            i -= 1
            j -= 1
            align1.insert(0, i)
            align2.insert(0, j)
        elif equal(score, scores[i - 1][j] + penalty):
            i -= 1
            align1.insert(0, i)
            align2.insert(0, '-')
        elif equal(score, scores[i][j - 1] + penalty):
            j -= 1
            align1.insert(0, '-')
            align2.insert(0, j)
        else:
            for scr in scores:
                print ' '.join(['%7s' % (round(y, 4)) for y in scr])
            raise Exception('Something  is failing and it is my fault...')
    return align1, align2, score


def equal(a, b, cut_off=1e-9):
    """
    Equality for floats
    """
    return abs(a - b) < cut_off


def optimal_cmo(hic1, hic2, num_v=None, max_num_v=None, verbose=False,
                method='score'):
    """

    Note: penalty is defined as the minimum value of the pre-scoring matrix
    
    :param hic1: first matrix to align
    :param hic2: second matrix to align
    :param None num_v: number of eigen vectors to consider, max is:
        max(min(len(hic1), len(hic2)))
    :param None max_num_v: maximum number of eigen vectors to consider.
    :param score method: distance function to use as alignment score. if 'score'
       distance will be the result of the last value of the Needleman-Wunsch
       algorithm. If 'frobenius' a modification of the Frobenius distance will
       be used.

    :returns: 2 lists, one per aligned matrix, plus a dict summarizing the
        goodness of the alignment with the distance between matrices, their 
        Spearman correlation Rho value and pvalue.
    """

    l_p1 = len(hic1)
    l_p2 = len(hic2)
    num_v = num_v or min(l_p1, l_p2)
    if max_num_v:
        num_v = min(max_num_v, num_v)
    if num_v > l_p1 or num_v > l_p2:
        raise Exception('\nnum_v should be at most {}\n'.format(min(l_p1,
                                                                    l_p2)))
    val1, vec1 = eig(hic1)
    if npsum(vec1).imag:
        raise Exception("ERROR: Hi-C data is not symmetric.\n" +
                        '{}\n\n{}'.format(hic1, vec1))
    val2, vec2 = eig(hic2)
    if npsum(vec2).imag:
        raise Exception("ERROR: Hi-C data is not symmetric.\n" +
                        '{}\n\n{}'.format(hic2, vec2))
    #
    val1 = array([sqrt(abs(v)) for v in val1])
    val2 = array([sqrt(abs(v)) for v in val2])
    idx = val1.argsort()[::-1]
    val1 = val1[idx]
    vec1 = vec1[idx]
    idx = val2.argsort()[::-1]
    val2 = val2[idx]
    vec2 = vec2[idx]
    #
    vec1 = array([val1[i] * vec1[:, i] for i in xrange(num_v)]).transpose()
    vec2 = array([val2[i] * vec2[:, i] for i in xrange(num_v)]).transpose()
    nearest = 100000000000
    best_alis = []
    for num in xrange(1, num_v + 1):
        for factors in product([1, -1], repeat=num):
            vec1p = factors * vec1[:, :num]
            vec2p = vec2[:, :num]
            p_scores = prescoring(vec1p, vec2p, l_p1, l_p2)
            penalty = min([npmin(p_scores)] + [0.0])
            align1, align2, dist = core_nw(p_scores, penalty, l_p1, l_p2)
            try:
                if method == 'frobenius':
                    dist = get_dist(align1, align2, hic1, hic2)
                else:
                    dist = -dist
                if dist < nearest:
                    nearest = dist
                    best_alis = [align1, align2]
                    best_pen = penalty
            except IndexError:
                pass
    try:
        align1, align2 = best_alis
    except ValueError:
        pass
    if verbose:
        print '\n Alignment (score = {}):'.format(nearest)
        print 'TADS 1: '+'|'.join(['%4s' % (str(int(x)) \
                                            if x!='-' else '-'*3) for x in align1])
        print 'TADS 2: '+'|'.join(['%4s' % (str(int(x)) \
                                            if x!='-' else '-'*3) for x in align2])
    rho, pval = get_score(align1, align2, hic1, hic2)
    print best_pen
    if not best_pen:
        print 'WARNING: penalty NULL!!!\n\n'
    return align1, align2, {'dist': nearest, 'rho': rho, 'pval': pval}
    

def virgin_score(penalty, l_p1, l_p2):
    """
    Fill a matrix with zeros, except first row and first column filled with \
    multiple values of penalty.
    """
    zeros    = [0.0 for _ in xrange(l_p2)]
    return [[penalty * j for j in xrange(l_p2)]] + \
           [[penalty * i] + zeros for i in xrange(1, l_p1)]


def prescoring(vc1, vc2, l_p1, l_p2):
    """
    yes... this is the bottle neck, almost 2/3 of the time spent here
    """
    p_score = []
    for i in xrange(l_p1):
        vci = vc1[i].__mul__
        p_score.append([vci(vc2[j]).sum() for j in xrange(l_p2)])
    return p_score
    #return [[sum(vc1[i] * vc2[j]) for j in xrange(l_p2)] for i in xrange(l_p1)]


def get_dist(align1, align2, tad1, tad2):
    """
    Frobenius norm
    """
    map1 = []
    map2 = []
    ext = 0
    for i, j in zip(align1, align2):
        if i != '-' and j != '-':
            map1.append(i)
            map2.append(j)
            ext = 0
        elif i == '-':
            if ext > 1:
                map1.append(0)
                map2.append(j)
            ext += 1
        else:
            if ext > 1:
                map1.append(i)
                map2.append(0)
            ext += 1
        #     ext += 1
        #     continue
        # if ext < 3:
        #     while ext:
        #         map1.pop(-1)
        #         map2.pop(-1)
        #         ext -= 1
        # else:
        #     ext = 0
    pp1 = [tad1[i][j] for i, j in combinations(map1, 2)]
    pp2 = [tad2[i][j] for i, j in combinations(map2, 2)]
    return sum([(pp - pp2[p])**2 for p, pp in enumerate(pp1)])/(len(pp1)+1)
    #return corrcoef(pp1, pp2, rowvar=0)[1, 0]


    
def get_score(align1, align2, tad1, tad2):
    """
    using Spearman Rho value
    TODO: perhaps use the gapped matrices instead.
    """
    map1 = []
    map2 = []
    for i, j in zip(align1, align2):
        if j != '-' and i != '-':
            map1.append(i)
            map2.append(j)
    pp1 = [[tad1[i][j] for j in map1] for i in map1]
    pp2 = [[tad2[i][j] for j in map2] for i in map2]
    return spearmanr(pp1, pp2, axis=None)


def get_OLD_score(align1, align2, p1, p2):
    """
    Original scoring function, based on contact map overlap.
    """
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


###
# Following is for aleigen


def matrix2binnary_contacts(tad1, tad2):
    cutoff = median(tad1)
    contacts1 = []
    contacts2 = []
    for i in xrange(len(tad1)):
        for j in xrange(i, len(tad1)):
            if tad1[i][j] > cutoff:
                contacts1.append((i, j))
    cutoff = median(tad2)
    for i in xrange(len(tad2)):
        for j in xrange(i, len(tad2)):
            if tad2[i][j] > cutoff:
                contacts2.append((i, j))
    return contacts1, contacts2


def run_aleigen(contacts1, contacts2, num_v):
    """

    * c1, c2 = number of contacts of the first and seconf contact map
               (after removing non-matching columns/rows)
    * cmo = total number of matching contacts (above the first diagonal)
      of the computed overlap
    * score = 2*CMO/(C1+C2)

    """
    f_string = '/tmp/lala{}.txt'
    f_name1 = f_string.format(1)
    f_name2 = f_string.format(2)
    write_contacts(contacts1, contacts2, f_string)
    sc_str = re.compile('Score\s+C1\s+C2\s+CMO\n([0-9.]+)\s+[0-9]+\s+.*')
    out = Popen('aleigen {} {} {}'.format(f_name1, f_name2, num_v),
                shell=True, stdout=PIPE).communicate()[0]
    score = [float(c) for c in re.findall(sc_str, out)]
    print out
    align1 = []
    align2 = []
    for line in out.split('\n')[2:]:
        if not re.match('[0-9]+\s+[0-9]+', line):
            continue
        el1, el2 = [int(c) for c in line.split()]
        align1.append(el1)
        align2.append(el2)
    return align1, align2, score


def write_contacts(contacts1, contacts2, f_string):
    for i, contacts in enumerate([contacts1, contacts2]):
        out = open(f_string.format(i+1), 'w')
        out.write(str(max([max(c) for c in contacts])+1) + '\n')
        out.write('\n'.join([str(c1) + ' ' + str(c2) for c1, c2 in contacts]))
        out.write('\n')
        out.close()

