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
from __future__ import print_function
#===============================================================================
# try:
#     import rpy2.robjects as robjects
#     from rpy2.robjects.packages   import importr
#     from rpy2.rinterface          import RRuntimeError
#     try:
#         dryhic = importr('dryhic')
#         from numpy import float64
#     except RRuntimeError:
#         pass
# except ImportError:
#     pass
#===============================================================================

from subprocess import Popen, PIPE
from os import path

from numpy import genfromtxt

from pytadbit.utils.file_handling import which


def oneD(tmp_dir='.', form='tot ~ s(map) + s(cg) + s(res)', p_fit=None,
         seed=1, **kwargs):
    """
    Normalizes according to oneD normalization that takes into account the GC
    content, mappability and the number of restriction sites per bin.

    Vidal, E., le Dily, F., Quilez, J., Stadhouders, R., Cuartero, Y., Graf, T., Marti-Renom, Marc A., Beato, M., Filion, G. (2017).
    OneD: increasing reproducibility of Hi-C Samples with abnormal karyotypes.
    bioRxiv. http://doi.org/10.1101/148254

    :param form: string representing an R Formulae
    :param None p_fit: proportion of data to be used in fitting (for very
       large datasets). Number between 0 and 1
    :param kwargs: dictionary with keys present in the formula and values being
       lists of equal length.
       for example:
           oneD(tot=[1,2,3...],
                map=[1,2,3...],
                res=[1,2,3...],
                cg =[1,2,3...])

    :returns: list of biases to use to normalize the raw matrix of interactions
    """

    script_path = which('normalize_oneD.R')
    proc_par = ["Rscript", "--vanilla", script_path]

    in_csv = path.join(tmp_dir, 'tot.csv')
    proc_par.append(in_csv)

    csvfile = open(in_csv, 'w')
    headers = sorted(kwargs.keys())
    csvfile.write(','.join(headers) + '\n')
    csvfile.write('\n'.join(','.join(str(kwargs[k][i]) for k in headers)
                            for i in range(len(kwargs['tot']))) + '\n')
    csvfile.close()

    out_csv = path.join(tmp_dir, 'biases.csv')
    proc_par.append(out_csv)

    proc_par.append('"%s"' % (form))

    if p_fit:
        proc_par.append(str(p_fit))

    if seed > 1:
        proc_par.append(str(seed))
    elif seed < 1:
        raise Exception(('ERROR: seed number (currently: %d) should be an '
                         'interger greater than 1 (because of R)') % (seed))

    proc = Popen(proc_par, stderr=PIPE)
    err = proc.stderr.readlines()
    print('\n'.join(err))

    biases_oneD = genfromtxt(out_csv, delimiter=',', dtype=float)

    return biases_oneD


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


def copy_matrix(hic_data, bads):
    W = {}
    N = len(hic_data)
    for k, v in hic_data.items():
        i, j = divmod(k, N)
        if i in bads or j in bads:
            continue
        try:
            W[i][j] = v
        except KeyError:
            W[i] = {}
            W[i][j] = v
    return W


def iterative(hic_data, bads=None, iterations=0, max_dev=0.00001,
              verbose=False, **kwargs):
    """
    Implementation of iterative correction Imakaev 2012

    :param hic_data: dictionary containing the interaction data
    :param None bads: dictionary with column not to be considered
    :param None remove: columns not to consider
    :param 0 iterations: number of iterations to do (99 if a fully smoothed
       matrix with no visibility differences between columns is desired)
    :param 0.00001 max_dev: maximum difference allowed between a row and the
       mean value of all raws
    :returns: a vector of biases (length equal to the size of the matrix)
    """
    if verbose:
        print('iterative correction')
    size = len(hic_data)
    if not bads:
        bads = {}
    remove = [i in bads for i in range(size)]
    remove = remove or tuple([int(hic_data[i+i*size]==0) for i in range(size)])

    if verbose:
        print("  - copying matrix")

    W = copy_matrix(hic_data, bads)
    B = dict([(b, 1.) for b in W])
    if len(B) == 0:
        raise ZeroDivisionError('ERROR: normalization failed, all bad columns')
    if verbose:
        print("  - computing biases")
    for it in range(iterations + 1):
        S, meanS = _update_S(W)
        DB = _updateDB(S, meanS, B)
        if iterations == 0: # exit before, we do not need to update W
            break
        _update_W(W, DB)
        S = sorted(S.values())
        dev = max(abs(S[0]  / meanS - 1), abs(S[-1] / meanS - 1))
        if verbose:
            print('   %15.3f %15.3f %15.3f %4s %9.5f' % (S[0], meanS, S[-1], it, dev))
        if dev < max_dev:
            break
    for i in range(size):
        try:
            if B[i]:
                B[i] *= meanS**.5
            else:
                B[i] = 1.
        except KeyError:
            B[i] = 1.
    return B


def expected(hic_data, bads=None, signal_to_noise=0.05, inter_chrom=False, **kwargs):
    """
    Computes the expected values by averaging observed interactions at a given
    distance in a given HiC matrix.

    :param hic_data: dictionary containing the interaction data
    :param None bads: dictionary with column not to be considered
    :param 0.05 signal_to_noise: to calculate expected interaction counts,
       if not enough reads are observed at a given distance the observations
       of the distance+1 are summed. a signal to noise ratio of < 0.05
       corresponds to > 400 reads.

    :returns: a vector of biases (length equal to the size of the matrix)
    """
    min_n = signal_to_noise ** -2. # equals 400 when default

    size = len(hic_data)
    try:
        if not inter_chrom:
            size = max(hic_data.chromosomes.values())
    except AttributeError:
        pass

    if size > 1200:
        import sys
        sys.setrecursionlimit(size + 100)

    expc = {}
    dist = 0
    while dist < size:
        diag = []
        new_dist, val = _meandiag(hic_data, dist, diag, min_n, size, bads)
        for dist in range(dist, new_dist + 1):
            expc[dist] = val
    return expc


def _meandiag(hic_data, dist, diag, min_n, size, bads):
    if hic_data.section_pos:
        for crm in hic_data.section_pos:
            for i in range(hic_data.section_pos[crm][0],
                            hic_data.section_pos[crm][1] - dist):
                if i in bads:
                    continue
                diag.append(hic_data[i, i + dist])
    else:
        for i in range(size - dist):
            if i in bads:
                continue
            diag.append(hic_data[i, i + dist])
    sum_diag = sum(diag)
    if not diag:
        return dist + 1, 0.
    if sum_diag > min_n:
        return dist + 1, float(sum_diag) / len(diag)
    if dist >= size:
        return dist + 1, float(sum_diag) / len(diag)
    return _meandiag(hic_data, dist + 1, diag, min_n, size, bads)
