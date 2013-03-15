"""
02 Dec 2012

this is not working... Local aligner for Topologically Associated Domains
"""
from numpy import zeros, log

def smith_waterma(tads1, tads2):
    """
    WARNING: not working!
    """
    tads1 = [1.0, 3.0, 5.0, 6.0, 12.0, 15.0, 17.0, 19.0]
    tads2 = [1.0, 5.0, 6.0, 7.0, 10.0, 11.0, 13.0, 17.0, 19.0, 20.0]
    bin_size = 1.
    chr_len  = 20.
    for penalty in xrange(1,10):
        penalty  = -float(penalty)/10
        for penalty2 in xrange(1,10):
            penalty2 = -float(penalty2)/10
            penalty  = -0.3
            penalty2 = -0.1
            print '-'*50
            print ' ->', penalty, penalty2
            l_tads1= len(tads1)
            l_tads2= len(tads2)
            pointers = zeros((l_tads1+1, l_tads2+1))
            scores   = zeros((l_tads1+1, l_tads2+1))
            fieled   = zeros((l_tads1+1, l_tads2+1))
            #penalty  = -.5
            for i in xrange(1, l_tads1+1):
                for j in xrange(1, l_tads2+1):
                    value = 1-log((1+bin_size*abs(tads2[j-1]-tads1[i-1]))**2)/log((1+chr_len)**2)
                    penalty_i = penalty_j = penalty
                    penalty_j = penalty2 if pointers[i, j-1] == 1 else penalty
                    penalty_i = penalty2 if pointers[i-1, j] == 2 else penalty
                    vals = {'diag': value + (scores[i-1, j-1]),
                            'left': value + (scores[i, j-1]  ) + penalty_j,
                            'up'  : value + (scores[i-1, j]  ) + penalty_i,
                            'end' : 0}
                    if bin_size*abs(tads2[j-1]-tads1[i-1]) > 3000000:
                        vals['diag'] = scores[i-1, j-1] + penalty
                    best = max(vals, key=lambda x: vals[x])
                    scores[i, j] = vals[best]
                    fieled[i, j] = value
                    if   best == 'end' :
                        pointers[i, j] = 0
                    elif best == 'left':
                        pointers[i, j] = 1
                    elif best == 'up'  :
                        pointers[i, j] = 2
                    else:
                        pointers[i, j] = 3
            align1 = []
            align2 = []
            val = scores.argmax()
            i = val / (len(tads2)+1)
            j = val % (len(tads2)+1)
            while pointers[i, j]:
                if pointers[i, j] == 3:
                    align1.insert(0,tads1[i-1])
                    align2.insert(0,tads2[j-1])
                    i -= 1
                    j -= 1
                elif pointers[i, j] == 2:
                    align1.insert(0,tads1[i-1])
                    align2.insert(0,'-')
                    i -= 1
                elif pointers[i, j] == 1:
                    align1.insert(0,'-')
                    align2.insert(0,tads2[j-1])
                    j -= 1
            print ' '.join(['%6s' % (int(y)) for y in [0.0]+tads2])
            for x in scores: print ' '.join(['%6s' % (round(y, 2)) for y in x])
            for x in pointers: print ' '.join(['%s' % (int(y)) for y in x])
            print '|'.join(['%4s' % (str(x)[:-2] if x!='-' else '-'*3) for x in align1])
            print '|'.join(['%4s' % (str(x)[:-2] if x!='-' else '-'*3) for x in align2])
            #for x in fieled: print ' '.join(['%6s' % (round(y, 2)) for y in x])
