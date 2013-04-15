"""
02 Dec 2012

Aligner based on reciprocal closest hits for Topologically Associated Domains
"""


def find_closest(num, tads1, start=0):
    closest = 0
    diff = inf = float('inf')
    for n in tads1:
        if n < start: continue
        if abs(n - num) < diff:
            diff = abs(n - num)
            closest = n
        elif diff != inf:
            break
    return closest


def find_closest_reciprocal(t1, tads1, tads2, start=0):
    """
    function to check the needleman_wunsch algorithm.
    """
    closest = None
    diff = inf = float('inf')
    gap = 0
    for t2 in tads2:
        if t2 <= start: continue
        if abs(t2 - t1) < diff:
            t1prim = find_closest(t2, tads1, start=t1)
            if t1 == t1prim:
                if diff != inf:
                    gap += 1
                diff = abs(t2 - t1)
                closest = t2
        elif diff != inf:
            break
    else:
        return '-', 0
    return closest, gap


def reciprocal(tads1, tads2, penalty=None, verbose=False, max_dist=None):
    """

    :argument None penalty: if 
    
    tads1 = [1, 5, 6, 9, 18, 22, 33, 34, 36]
    tads2 = [0, 1, 7, 9, 17, 21, 26, 29, 35]
    tads1 = [15560, 17240, 18680, 18960, 19100, 19220]
    tads2 = [15600, 17200, 18680, 18980, 19340]
    penalty = None
    """
    if not penalty:
        # set penalty to the average length of a TAD
        penalty  = float(reduce(lambda x, y: abs(x - y),
                                tads1)) / len(tads1)
        penalty += float(reduce(lambda x, y: abs(x - y),
                                tads2)) / len(tads2)
    start = 0
    diffs = []
    align1 = []
    align2 = []
    i = 0
    for t in tads1:
        closest, gap = find_closest_reciprocal(t, tads1, tads2,
                                               start=start)
        diff = penalty
        if closest != '-':
            diff = abs(t - closest)
            if diff > max_dist:
                print 'MAAAAX: ', t, closest, diff, max_dist
                closest = '-'
                diff = penalty
            else:
                start = closest
        diffs.append(diff)
        while gap > 0:
            align1.append('-')
            align2.append(tads2[i])
            i += 1
            gap -= 1
        align1.append(t)
        align2.append(closest)
        if closest != '-':
            i += 1
    if verbose:
        print '\n Alignment:'
        print 'TADS 1: '+'|'.join(['%9s' % (str(int(x)) if x!='-' else '-'*3) \
                                   for x in align1])
        print 'TADS 2: '+'|'.join(['%9s' % (str(int(x)) if x!='-' else '-'*3) \
                                   for x in align2])
    return [align1, align2], float(sum(diffs))/len(align1)
        
