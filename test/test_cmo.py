"""
25 Mar 2013


"""
from random import random
from copy import deepcopy
from numpy import array, log2
from scipy.interpolate import interp1d
from cPickle import load
from pytadbit.tad_clustering.tad_cmo import optimal_cmo
from matplotlib import pyplot as plt
from scipy.stats import zscore

inf = open('tad_distr.pik')
WINA, CNTA, WIND, CNTD = load(inf)
inf.close()

DISTRA = interp1d(WINA, CNTA)
DISTRD = interp1d(WIND, CNTD)


def generate_1_interaction(diag=False):
    if diag:
        return float(DISTRD(random()))
    if random() > .8: # prob of having a 0.0
        return 0.0
    return float(DISTRA(random()))


def generate_interaction(size):
    rnd_inter = {}
    for i in xrange(size):
        for j in xrange(i, size):
            if i == j:
                rnd_inter[(i, j)] = float(DISTRD(random()))
                continue
            if random() > .8:
                rnd_inter[(i, j)] = 0.1
                rnd_inter[(j, i)] = 0.1
                continue
            rnd_inter[(i, j)] = float(DISTRA(random()))
            rnd_inter[(j, i)] = rnd_inter[(i, j)]
    return rnd_inter


def generate_random_contacts(size=None, tad1=None, prob=0.5,
                             indel=2, ext=5, p_insert=0.5):
    """
    returns 2 random matrices correposnding to TADs
    and 2 contacts lists.
    """
    if not size:
        size1 = size2  = 10# + int(random()*5)
    else:
        size1 = size2  = size# + int(random()*5)
    if tad1:
        contacts1 = {}
        for i in xrange(len(tad1m)):
            for j in xrange(len(tad1m)):
                contacts1[(i, j)] = tad1m[i][j]
    else:
        contacts1 = generate_interaction(size1)
    contacts2 = deepcopy(contacts1)
    # randomly change some contacts
    for i in xrange(size1):
        for j in xrange(i, size1):
            if random() < prob:
                contacts2[(i, j)] = generate_1_interaction(i==j)
                contacts2[(j, i)] = contacts2[(i, j)]
    # a bit removing and inserting columns
    changes = 0
    p_delete = 1 - p_insert
    indels = []
    while changes < indel:
        rnd = random()
        my_j = int(random()*size2)
        if rnd < p_insert:
            for _ in xrange(ext):
                indels.append(my_j + 1)
                size2 += 1
                for i, d in enumerate(indels):
                    if abs(d) > (my_j + 1):
                        indels[i] = (d + 1) if d > 0 else -(abs(d) + 1)
                for i in xrange(size2 - 2, -1, -1):
                    for j in xrange(size2 - 2, i - 1, -1):
                        if i >= my_j and j >= my_j:
                            contacts2[(i + 1, j + 1)] = contacts2[(j, i)]
                        if j >= my_j:
                            contacts2[(i, j + 1)]     = contacts2[(j, i)]
                    if i == my_j:
                        contacts2[(i, j + 1)]     = contacts2[(j, i)]
                for i in xrange(size2):
                    for j in xrange(i, size2):
                        contacts2[(j, i)] = contacts2[(i, j)]
                for i in xrange(size2):
                    contacts2[(i, my_j)] = generate_1_interaction(i==my_j)
                    contacts2[(my_j, i)] = contacts2[(i, my_j)]
        elif rnd > 1 - p_delete:
            my_j = min(my_j, size2 - ext)
            for _ in xrange(ext):
                size2 -= 1
                indels.append(-(my_j + 1))
                for i, d in enumerate(indels):
                    if abs(d) > (my_j + 1):
                        indels[i] = d - 1  if d > 0 else -(abs(d) - 1)
                for i in xrange(size2):
                    for j in xrange(i, size2):
                        if i >= my_j and i >= my_j:
                            contacts2[(i, j)] = contacts2[(i + 1, j + 1)]
                        elif j >= my_j:
                            contacts2[(i, j)] = contacts2[(i , j + 1)]
                for i in xrange(size2):
                    for j in xrange(i + 1, size2):
                        contacts2[(j, i)] = contacts2[(i, j)]
                for i in xrange(size2):
                    del(contacts2[(i, size2)])
                    del(contacts2[(size2, i)])
                del(contacts2[(size2, size2)])
        else:
            continue
        changes += 1
    print 'inserted/deleted', [str(i-1) if i>0 else '-' + str((abs(i)-1)) \
                               for i in indels]
    for i in range(size1):
        if not (i, i) in contacts1:
            contacts1[(i, i)] = generate_1_interaction(True)
    for i in range(size2):
        if not (i, i) in contacts2:
            contacts2[(i, i)] = generate_1_interaction(True)
    tad1 = contact2matrix(contacts1, size1)
    tad2 = contact2matrix(contacts2, size2)
    for t in tad1: print '\t'.join(['{:.1f}'.format(i) for i in t])# + '\n'
    print '-'*20
    for t in tad2: print '\t'.join(['{:.1f}'.format(i) for i in t])# + '\n'
    return tad1, tad2, indels#, contacts1, contacts2


def contact2matrix(contacts, size):
    matrix = [[0  for _ in xrange(size)] for _ in xrange(size)]
    for i, j in contacts:
        try:
            matrix[i][j] = contacts[(i, j)]
            matrix[j][i] = contacts[(j, i)]
        except IndexError:
            pass
    return matrix


def merge_tads(tad1, tad2, ali):
    ali1, ali2, _ = deepcopy(ali)
    for i in xrange(max(ali1[0], ali2[0]) - 1, -1, -1):
        if ali[0][0]:
            ali1.insert(0, i)
            ali2.insert(0, '-')
        elif ali[1][0]:
            ali1.insert(0, '-')
            ali2.insert(0, i)
    size = len(ali1)
    matrix1 = [[0.1 for _ in xrange(size)] for _ in xrange(size)]
    matrix2 = [[0.1 for _ in xrange(size)] for _ in xrange(size)]
    matrixm = [[0.1 for _ in xrange(size)] for _ in xrange(size)]
    for i in xrange(size):
        if ali1[i] == '-' or ali2[i] == '-':
            matrixm[i] = [float('nan') for _ in xrange(size)]
            if ali1[i] == '-':
                matrix1[i]  = [float('nan') for _ in xrange(size)]
                matrix2[i]  = [tad2[ali2[i]][ali2[j]] \
                               if ali2[j] != '-' else float('nan')\
                               for j in xrange(size)]
            elif ali2[i] == '-':
                matrix2[i]  = [float('nan') for _ in xrange(size)]
                matrix1[i]  = [tad1[ali1[i]][ali1[j]] \
                               if ali1[j] != '-' else float('nan')\
                               for j in xrange(size)]
            continue
        for j in xrange(size):
            if ali1[j] == '-' or ali2[j] == '-':
                matrixm[i][j] = float('nan')
                if ali1[j] == '-':
                    matrix1[i][j] = float('nan')
                    matrix2[i][j]  = tad2[ali2[i]][ali2[j]]
                elif ali2[j] == '-':
                    matrix2[i][j] = float('nan')
                    matrix1[i][j] = tad1[ali1[i]][ali1[j]]
                continue
            matrix1[i][j]  = tad1[ali1[i]][ali1[j]]
            matrix2[i][j]  = tad2[ali2[i]][ali2[j]]
            matrixm[i][j]  = tad1[ali1[i]][ali1[j]]
            matrixm[i][j] += tad2[ali2[i]][ali2[j]]
            matrixm[i][j] /= 2
    return matrix1, matrix2, matrixm


def to_binary(tad1, tad2):
    cells = []
    for i in range(len(tad1)):
        for j in range(i, len(tad1)):
            cells.append(tad1[i][j])
    for i in range(len(tad2)):
        for j in range(i, len(tad2)):
            cells.append(tad1[i][j])
    cells = zscore(cells)
            

def main():
    """
    main function
    """

    results = {'nw-ext frob-ext': [],
               'nw-ext frob-sam': [],
               'nw-sam frob-ext': [],
               'nw-sam frob-sam': [],
               'binary frob'}
    for II in xrange(100):

        print II,'='*80
        size = 11
        ext = 1
        tad1, tad2, indels = generate_random_contacts(prob=0.0, indel=3,
                                                      size=size, ext=ext,
                                                      p_insert=0.5)


        align1 = range(size)
        align2 = range(size)
        ali1 = 0
        ali2 = 0
        for i in sorted(indels, key=abs):
            if i > 0:
                align1.insert(i-1+ali2, '-')
                ali2 += 1
            else:
                align2.insert(abs(i)-1 + ali1, '-')
                ali1 += 1
                # ali2 -= 1
        plus = len(align2) - align2.count('-')
        if ali2 > 0:
            align2 += range(plus, plus+ali2)
        else:
            align2 = align2[:ali2]
        # align1 = []
        # align2 = []
        # j1 = 0
        # j2 = 0
        # for i in range(size + 1):
        #     if i in [a - 1 for a in indels if a > 0]:
        #         align1 += ['-']
        #         j2 -= 1
        #     if size > i:
        #         align1.append(i)
        #     if i + j2 in [abs(a) - 1 for a in indels if a < 0]:
        #         align2 += ['-']
        #         j1 -= 1
        #     if size + (j1-j2) > i:
        #         align2.append(i)
        print '|'.join(['{:^3}'.format(a) for a in align1])
        print '|'.join(['{:^3}'.format(a) for a in align2])
        
       # tad1 = tad1[::-1]
        # tad2 = tad2[::-1]
        aliF1 = optimal_cmo(tad1, tad2, max_num_v=12, long_nw=True,
                            long_dist=True, method='frobenius', verbose=True)
        aliF2 = optimal_cmo(tad2, tad1, max_num_v=12, long_nw=True,
                            long_dist=True, method='frobenius', verbose=True)
        aliF = list(min((aliF1, aliF2), key=lambda x: x[2]['dist']))
        if aliF == list(aliF2):
            aliF[0], aliF[1] = aliF[1], aliF[0]
        print align1
        print aliF[0]
        print align2
        print aliF[1]
        results['nw-ext frob-ext'].append(align1==aliF[0] and align2==aliF[1])
        aliF = optimal_cmo(tad1, tad2, max_num_v=12, long_nw=True,
                           long_dist=False, method='frobenius', verbose=True)
        results['nw-ext frob-sam'].append(align1==aliF[0] and align2==aliF[1])
        aliF = optimal_cmo(tad1, tad2, max_num_v=12, long_nw=False,
                           long_dist=True, method='frobenius', verbose=True)
        results['nw-sam frob-ext'].append(align1==aliF[0] and align2==aliF[1])
        aliF = optimal_cmo(tad1, tad2, max_num_v=12, long_nw=False,
                           long_dist=False, method='frobenius', verbose=True)
        results['nw-sam frob-sam'].append(align1==aliF[0] and align2==aliF[1])
    for r in results:
        print r, results[r].count(True), '/', len(results[r])


    matrix1, matrix2, matrixF = merge_tads(tad1, tad2, aliF)
    cmap = plt.get_cmap()
    cmap.set_bad(color = 'k', alpha = .7)
    plt.subplot(2, 2, 2)
    plt.imshow(log2(matrix1), origin='lower', interpolation="nearest")
    plt.title('Original matrix')
    plt.subplot(2, 2, 3)
    plt.imshow(log2(matrix2), origin='lower', interpolation="nearest")
    plt.title('Changed Matrix (indels: {})'.format(indels))
    plt.subplot(2, 2, 1)
    plt.imshow(log2(matrixF), origin='lower', interpolation="nearest")
    plt.title('Merged Frobenius')
    plt.subplot(2, 2, 4)
    diff = len(tad1) - len(tad2)
    for i in xrange(len(tad2)):
        tad2[i] += [0.]*diff
    for i in xrange(diff):
        tad2.append([0.]*len(tad1))
    plt.imshow(log2(tad2), origin='lower', interpolation="nearest")
    plt.title('Changed Matrix')
    plt.show()


# align1 = aliF[0]
align1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
# align1 = [0, 1, 2, 3, 4, 5, 6, 7, '-', 8, 9]
# # align1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# align2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
align2 = [0, '-', 1, 2, 3, 4, 5, 6, 7, 8]
# align2 = [0, 1, '-', 2, 3, 4, 5, 6, 7, 8]
# align2 = [0, 1, 2, '-', 3, 4, 5, 6, 7, 8]
# align2 = [0, 1, 2, 3, '-', 4, 5, 6, 7, 8]
# align2 = [0, 1, 2, 3, 4, '-', 5, 6, 7, 8]
# align2 = [0, 1, 2, 3, 4, 5, 6, '-', 7, 8]
# align2 = [0, 1, 2, 3, 4, 5, 6, 7, '-', 8]
# print get_dist(align1, align2, tad1, tad2)

# out = open('/home/fransua/Projects/tad-ledily/scripts/lala.pik')
# tad1m = load(out)
# out.close()




if __name__ == "__main__":
    exit(main())

for p, pp in enumerate(pp1):
    print 
