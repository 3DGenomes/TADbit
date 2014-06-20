"""
20 juin 2014

Implementation of iterative correction Imakaev 2012
"""

def _update_S(W):
    S = {}
    meanS = 0.0
    for crm1 in W:
        S[crm1] = {}
        for bin1 in W[crm1]:
            W1 = W[crm1][bin1]
            S[crm1][bin1] = sum([sum(W1[crm2].values()) for crm2 in W1])
            meanS += S[crm1][bin1]
    meanS /= len(W[crm1])
    return S, meanS

def _updateDB(S, meanS, B):
    DB = {}
    for crm in S:
        DB[crm] = {}
        for bin1 in S[crm]:
            DB[crm][bin1] = float(S[crm][bin1]) / meanS
            B[crm][bin1] *= DB[crm][bin1]
    return DB

def _update_W(W, DB):
    for crm1 in W:
        Wc1 = W[crm1]
        for bin1 in Wc1:
            DBbin1 = DB[crm1][bin1]
            W1 = Wc1[bin1]
            for crm2 in W1:
                DBc2 = DB[crm2]
                W1c2 = W1[crm2]
                for bin2 in W1c2:
                    try:
                        W1c2[bin2] /= DBbin1 * DBc2[bin2]
                    except ZeroDivisionError:
                        # whole row is empty
                        continue

def iterative(W, iterations=0):
    """
    """
    B = dict([(c, dict([(b, 1.) for b in W[c]])) for c in W])
    for it in xrange(iterations + 1):
        S, meanS = _update_S(W)
        DB = _updateDB(S, meanS, B)
        _update_W(W, DB)
    for c in B:
        for i in B[c]:
            B[c][i] *= meanS**.5
    return B
