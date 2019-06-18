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
# from scipy.stats import zscore
import multiprocessing as mu
from pytadbit import Experiment
from pytadbit.utils.tadmaths import zscore


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
                             indel=2, ext=5, p_insert=0.5, verbose=False):
    """
    returns 2 random matrices correposnding to TADs
    and 2 contacts lists.
    """
    if not size:
        size1 = size2  = 10# + int(random()*5)
    else:
        size1 = size2  = size# + int(random()*5)
    if tad1:
        size1 = size2 = size = len(tad1)
        contacts1 = {}
        for i in xrange(len(tad1)):
            for j in xrange(len(tad1)):
                contacts1[(i, j)] = tad1[i][j]
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
        elif rnd > 1 - p_delete and size2 > ext:
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
    if verbose:
        print '\t'.join(['\033[0;31m{}\033[0m'.format(i) for i in range(len(tad1))])# + '\n'
        for t in tad1: print '\t'.join(['{:.1f}'.format(i) for i in t])# + '\n'
        print '-'*20
        print '\t'.join(['\033[0;31m{}\033[0m'.format(i) for i in range(len(tad2))])# + '\n'
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


def to_binary(exp, method='bin'):
    exp.normalize_hic(method='bytot')
    hic = exp.hic_data[0][:]
    nhic = [(float(hic[i]) / exp.wght[0][i]) if exp.wght[0][i] else 0.0 \
            for i in xrange(exp.size**2)]
    zscore(nhic, exp.size)
    if method == 'trin':
        return [1 if z.__str__() == 'nan' \
                else 0 if -1 < z < 1 else 1 if z > 0 else -1 for z in nhic]
    elif method == 'fivn':
        return [1 if z.__str__() == 'nan' \
                else 0 if z < -2 \
                else 1 if z < -1 \
                else 2 if z < 1 \
                else 3 if z < 2 \
                else 4 for z in nhic]
    return [0 if z < 0 else 1 for z in nhic]
            

def get_original(size, indels, verbose=False):
    align1 = range(size)
    align2 = range(size)
    for i in sorted(indels, key=abs):
        if i > 0:
            i -= 1
            try:
                j = align2.index(i)
            except ValueError:
                j = align2.count('-') + i
            align1.insert(j, '-')
            last = [a for a in align2[:j] if a != '-']
            last = last[-1] if last else 0
            align2 = (align2[:j] +
                      [last + 1] +
                      [a + 1 if a != '-' else a for a in align2[j:]])
        else:
            i = abs(i) - 1
            try:
                j = align2.index(i)
            except ValueError:
                j = align2.count('-') + i
            align2.insert(j, '-')
        for i in xrange(min(len(align1), len(align2))):
            if align1[i] == align2[i] == '-':
                align1.pop(i)
                align2.pop(i)
                break
    align1 = align1[:min(len(align1), len(align2))]
    align2 = align2[:min(len(align1), len(align2))]
    if verbose:
        print 'TADS 1: ' + '|'.join(['{:>4}'.format(a if a != '-' else '---') for a in align1])
        print 'TADS 2: ' + '|'.join(['{:>4}'.format(a if a != '-' else '---') for a in align2])
    return align1, align2


def print_ali(ali1, ali2, align1, align2):
    colsc = zip(ali1, ali2)
    colsr = zip(align1, align2)
    real = ['\033[0;32m' if col in colsr else '' for col in colsc]
    end = ['\033[0m' if col in colsr else '' for col in colsc]
    print 'TADS 1: ' + '|'.join([real[i] + '{:>4}'.format(a if a != '-' else '---') + end[i] for i, a in enumerate(ali1)])
    print 'TADS 2: ' + '|'.join([real[i] + '{:>4}'.format(a if a != '-' else '---') + end[i] for i, a in enumerate(ali2)])
    return float(real.count('\033[0;32m'))/len(colsc)
    

def main():
    """
    main function
    """

    out = open('lala.pik')
    tad1m = load(out)
    out.close()
    
    use_real = True
    kind = 'raw'
    results = {'nw-ext frob-ext': [],
               'nw-ext frob-sam': [],
               'nw-sam frob-ext': [],
               'nw-sam frob-sam': [],
               'binary frob'    : []}
    pool = mu.Pool(4)
    jobs = {}
    for II in xrange(10):
        
        size   = len(tad1m) if use_real else 21
        ext    = 1 + int(random()*4)
        max_ev = 8
        tad1, tad2, indels = generate_random_contacts(tad1=tad1m if use_real else None,
                                                      prob=0.05, indel=3,
                                                      size=size, ext=ext,
                                                      p_insert=0.5)
        if kind in ['bin', 'trin']:
            exp1 = Experiment('test1', resolution=20000, xp_handler=[tad1])
            exp2 = Experiment('test2', resolution=20000, xp_handler=[tad2])
            exp1.hic_data[0] = to_binary(exp1)
            exp2.hic_data[0] = to_binary(exp2)
            tad1 = exp1.get_hic_matrix()
            tad2 = exp2.get_hic_matrix()
        # both sides
        jobs[II] = {}
        jobs[II][0] = indels
        jobs[II][1] = pool.apply_async(optimal_cmo, args=(tad1, tad2),
                                       kwds={'max_num_v' : max_ev,
                                             'long_nw'   : True,
                                             'long_dist' : True,
                                             'method'    : 'frobenius',
                                             'verbose'   : False})
        jobs[II][2] = pool.apply_async(optimal_cmo, args=(tad1, tad2),
                                       kwds={'max_num_v' : max_ev,
                                             'long_nw'   : True,
                                             'long_dist' : True,
                                             'method'    : 'frobenius',
                                             'verbose'   : False})
        # jobs[II][3] = pool.apply_async(optimal_cmo, args=(tad1, tad2),
        #                                kwds={'max_num_v' : max_ev,
        #                                      'long_nw'   : False,
        #                                      'long_dist' : True,
        #                                      'method'    : 'frobenius',
        #                                      'verbose'   : False})
        # jobs[II][4] = pool.apply_async(optimal_cmo, args=(tad1, tad2),
        #                                kwds={'max_num_v' : max_ev,
        #                                      'long_nw'   : False,
        #                                      'long_dist' : True,
        #                                      'method'    : 'frobenius',
        #                                      'verbose'   : False})
        # aliF1 = optimal_cmo(tad1, tad2, max_num_v=max_ev, long_nw=True,
        #                     long_dist=True, method='frobenius', verbose=False)
        # aliF2 = optimal_cmo(tad2, tad1, max_num_v=max_ev, long_nw=True,
        #                     long_dist=True, method='frobenius', verbose=False)
    pool.close()
    pool.join()
    
    for II in xrange(10):
        print II, '=' * 80
        indels = jobs[II][0]
        print 'inserted/deleted', II, [str(i-1) if i>0 else '-' + str((abs(i)-1)) \
                                       for i in indels]
        align1, align2 = get_original(size, indels, verbose=True)
        del(jobs[II][0])
        #aliF1, aliF2, aliF3, aliF4 = [jobs[II][j].get() for j in jobs[II]]
        aliF1, aliF2 = [jobs[II][j].get() for j in jobs[II]]
        aliF = list(min((aliF1, aliF2), key=lambda x: x[2]['dist']))
        if aliF == list(aliF2):
            aliF[0], aliF[1] = aliF[1], aliF[0]
        # print_ali(aliF1[0], aliF1[1], align1, align2)
        # print_ali(aliF2[1], aliF2[0], align1, align2)
        results['nw-ext frob-ext'].append(print_ali(aliF[1], aliF[0],
                                                    align1, align2))
        # 1 side no gap extension penalty in frobenius calculus
        # aliF1 = optimal_cmo(tad1, tad2, max_num_v=max_ev, long_nw=True,
        #                     long_dist=False, method='frobenius', verbose=False)
        # aliF2 = optimal_cmo(tad2, tad1, max_num_v=max_ev, long_nw=True,
        #                     long_dist=False, method='frobenius', verbose=False)
        # aliF = list(min((aliF3, aliF4), key=lambda x: x[2]['dist']))
        # if aliF == list(aliF4):
        #     aliF[0], aliF[1] = aliF[1], aliF[0]
        # results['nw-ext frob-sam'].append(print_ali(aliF[0], aliF[1],
        #                                             align1, align2))
        # aliF = optimal_cmo(tad1, tad2, max_num_v=max_ev, long_nw=False,
        #                    long_dist=True, method='frobenius', verbose=True)
        # results['nw-sam frob-ext'].append(align1==aliF[0] and align2==aliF[1])
        # aliF = optimal_cmo(tad1, tad2, max_num_v=12, long_nw=False,
        #                    long_dist=False, method='frobenius', verbose=True)
        # results['nw-sam frob-sam'].append(align1==aliF[0] and align2==aliF[1])
    for r in results:
        print r, sum(results[r]), '/', len(results[r])


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
# align1 = [0, 1, 2, 3, 4, 5, 6, 7, '-', 8, 9]
# # align1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# align2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# align2 = [0, '-', 1, 2, 3, 4, 5, 6, 7, 8]
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
