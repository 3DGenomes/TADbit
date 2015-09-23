"""
20 juin 2014

Implementation of iterative correction Imakaev 2012

Schematic flow chart for five iterations (it = 0->5) starting from the symmetric
matrix W of size N:

   +---------->  Wij
   |    
   |              |
   |              v
   |             __
   |             \
   |        Si = /_    Wij
   |            j=0->N
   |    
   |              |
   |              v
   |    
   |                   _
   |        DBi = Si / S
   |    
   |              |
   |              v
   |    
   |        Bi = Bi x DBi             ---> keep track, used as expected value
   |    
   |              |
   |              v
   |    
   |                 Wij
   |       Wij = -----------
   |              DBi x DBj  
   |    
   |              |
   |        it<5 / \ it=5
   |____________/   \_________   TADbit          _           
    it++                     \`----------> Wij / S    meaning that: Si = O(1)
                             |                        ('Si' tends towards one
                             |                         when 'it' -> infinite)
                             |Strict Imakaev
                             |
                             v
                             
                            Wij
                          -------  meaning that: Si = 1
                           ___
                           \
                           /__ Wi

"""

def _update_S(W):
    S = {}
    meanS = 0.0
    for bin1 in W:
        S[bin1] = sum(W[bin1].values())
        meanS += S[bin1]
    meanS /= len(W)
    return S, meanS

def _updateDB(S, meanS, B):
    DB = {}
    for bin1 in S:
        DB[bin1] = float(S[bin1]) / meanS
        B[bin1] *= DB[bin1]
    return DB

def _update_W(W, DB):
    for bin1 in W:
        DBbin1 = DB[bin1]
        W1 = W[bin1]
        for bin2 in W1:
            try:
                W1[bin2] /= DBbin1 * DB[bin2]
            except ZeroDivisionError: # whole row is empty
                continue

def iterative(hic_data, bads=None, iterations=0, max_dev=0.00001,
              verbose=False):
    """
    Implementation of iterative correction Imakaev 2012
    
    :param hic_data: dictionary containing the interaction data
    :param None remove: columns not to consider
    :param 0 iterations: number of iterations to do (99 if a fully smoothed
       matrix with no visibility differences between columns is desired)
    :param 0.00001 max_dev: maximum difference allowed between a row and the
       mean value of all raws
    :returns: a vector of biases (length equal to the size of the matrix)
    """
    if verbose:
        print 'iterative correction'
    size = len(hic_data)
    if not bads:
        bads = {}
    remove = [i in bads for i in xrange(size)]
    remove = remove or tuple([int(hic_data[i+i*size]==0) for i in xrange(size)])
    W = {}
    for i in xrange(size):
        if remove[i]:
            continue
        W[i] = {}
        for j in xrange(size):
            if remove[j]:
                continue
            if hic_data[i, j]:
                W[i][j] = hic_data[i, j]
    B = dict([(b, 1.) for b in W])
    if len(W) == 0:
        raise ZeroDivisionError('ERROR: normalization failed, all bad columns')
    for it in xrange(iterations + 1):
        S, meanS = _update_S(W)
        DB = _updateDB(S, meanS, B)
        _update_W(W, DB)
        S = sorted(S.values())
        dev = max(abs(S[0]  / meanS - 1), abs(S[-1] / meanS - 1))
        if verbose:
            print '   %15.3f %15.3f %15.3f %4s %9.5f' % (S[0], meanS, S[-1], it, dev)
        if dev < max_dev:
            break
    for i in xrange(size):
        try:
            if B[i]:
                B[i] *= meanS**.5
            else:
                B[i] = 1.
        except KeyError:
            B[i] = 1.
    return B

