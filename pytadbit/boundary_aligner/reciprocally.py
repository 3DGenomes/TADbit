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


def find_closest_reciprocal(t1, tads1, tads2, penalty, start=0):
    """
    function to check the needleman_wunsch algorithm.
    """
    closest = None
    diff = inf = float('inf')
    gap = 0
    for t2 in tads2:
        if t2 < start: continue
        if abs(t2 - t1) < diff:
            t1prim = find_closest(t2, tads1)
            if t1 == t1prim:
                diff = abs(t2 - t1)
                if closest > -1:
                    gap += 1
                closest = t2
        elif diff != inf:
            break
    if diff == inf:
        return '-', None, penalty
    #return None, None
    return closest, gap, diff


def reciprocal(tads1, tads2, penalty=None, verbose=False, max_dist=None):
    """

    :argument None penalty: if 
    
    tads1 = [1, 5, 6, 9, 18, 22, 33, 34, 36]
    tads2 = [0, 1, 7, 9, 17, 21, 26, 29, 35]
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
        closest, gap, diff = find_closest_reciprocal(t, tads1, tads2,
                                                     penalty, start=start)
        if closest != '-':
            start = closest + 1 # FIXME!!!!!
            if diff > max_dist:
                # FIXME
                pass
        diffs.append(diff)
        while gap:
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
        
