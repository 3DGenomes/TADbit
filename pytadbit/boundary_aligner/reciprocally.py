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


def find_closest_reciprocal(num, tads1, tads2, start=0):
    """
    function to check the needleman_wunsch algorithm.
    """
    closest = 0
    diff = inf = float('inf')
    for n in tads2:
        if n < start: continue
        if abs(n - num) < diff:
            other_num = find_closest(n, tads1)
            print n, num, other_num
            if num == other_num:
                diff = abs(n - num)
                closest = n
            else:
                print 'gap', num, other_num
        elif diff != inf:
            break
    if diff == inf:
        return '---', None
    #return None, None
    return closest, diff


def reciprocal(tads1, tads2):
    start = 0
    diffs = []
    align1 = []
    align2 = []
    for t in tads1:
        closest, diff = find_closest_reciprocal(t, tads1, tads2, start=start)
        if closest != '---':
            start=closest
        print t, closest
        diffs.append(diff)
        align1.append(t)
        align2.append(closest)
        
