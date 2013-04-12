"""
02 Dec 2012

Aligner based on reciprocal closest hits for Topologically Associated Domains
"""


def find_closest(num, tads1):
    closest = 0
    diff = inf = float('inf')
    for n in tads1:
        if abs(n - num) < diff:
            diff = abs(n - num)
            closest = n
        elif diff != inf:
            break
    return closest


def find_closest_reciprocal(num, tads1, tads2, penalty, start=0):
    """
    function to check the needleman_wunsch algorithm.
    """
    closest = None
    diff = inf = float('inf')
    gap = 0
    for n in tads2:
        if n < start: continue
        if abs(n - num) < diff:
            other_num = find_closest(n, tads1)
            if num == other_num:
                diff = abs(n - num)
                if closest > -1:
                    gap += 1
                closest = n
        elif diff != inf:
            break
    if diff == inf:
        return '-', None, penalty
    #return None, None
    return closest, gap, diff


def reciprocal(tads1, tads2, penalty=None, verbose=False):
    """

    :argument None penalty: if 
    
    tads1 = [1, 5, 6, 9, 18, 22, 33, 34, 36]
    tads2 = [0, 1, 7, 9, 17, 21, 26, 29, 35]
    """
    if not penalty:
        # set penalty to the average length of a TAD
        penalty = float(reduce(lambda x, y: abs(x - y), tads1))/len(tads1)
        penalty += float(reduce(lambda x, y: abs(x - y), tads2))/len(tads2)
    start = 0
    diffs = []
    align1 = []
    align2 = []
    i = 0
    for t in tads1:
        closest, gap, diff = find_closest_reciprocal(t, tads1, tads2,
                                                     penalty, start=start)
        if closest != '-':
            start=closest
        diffs.append(diff)
        while gap:
            align1.append('-')
            align2.append(tads2[i])
            i += 1
            gap -= 1
        align1.append(t)
        align2.append(closest)
        i += 1
    if verbose:
        print 'TAD1: ' + '|'.join(['{:>3}'.format(a) for a in align1])
        print 'TAD2: ' + '|'.join(['{:>3}'.format(a) for a in align2])
    return align1, align2, float(sum(diffs))/len(align1)
        
