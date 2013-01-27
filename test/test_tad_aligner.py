"""
21 Jan 2013


"""

import numpy as np
from random import random
from copy import deepcopy
from pytadbit.tad_clustering.tad_cmo import optimal_cmo

def main():
    """
    main function
    """

    tad1, tad2, contacts1, contacts2 = generate_random_contacts(size=13)
    size1 = len(tad1)
    size2 = len(tad2)#max([max(c) for c in contacts2])
    num_v = min(size1, size2)
    #align1a, align2a, scorea = run_aleigen(contacts1,
    #                                       contacts2, num_v)
    align1, align2 = optimal_cmo(tad1, tad2, num_v)
    align1b = []
    align2b = []
    for c in xrange(len(align1)):
        if align1[c] != '-' and align2[c] != '-':
            align1b.append(align1[c])
            align2b.append(align2[c])
    print '='*80
    print ' '+''.join(['{:>2}'.format(i) for i in range(len(tad1))])
    print tad1
    print '-'*80
    print ' '+''.join(['{:>2}'.format(i) for i in range(len(tad2))])
    print tad2
    print '='*80
    #print ''.join(['{:>4}'.format(a) for a in align1a])
    #print ''.join(['{:>4}'.format(a) for a in align2a])
    print '-'*80
    print ''.join(['{:>4}'.format(a) for a in align1b])
    print ''.join(['{:>4}'.format(a) for a in align2b])
    print '='*80


def generate_random_contacts(size=10, prob_contact=0.5, 
                             prob_mut=0.1, num_changes=2):
    """
    returns 2 random matrices corresponding to TADs
    and 2 contacts lists.
    """
    size1 = size2  = size# + int(random()*5)
    contacts1 = []
    for i in xrange(size1):
        for j in xrange(i, size1):
            if random() < prob_contact:
                contacts1.append((i, j))
    contacts2 = deepcopy(contacts1)
    # randomly change some contacts
    for c in xrange(len(contacts2) -1, -1, -1):
        if random() < prob_mut:
            del(contacts2[c])
    for i in xrange(size1):
        for j in xrange(i, size1):
            if random() < prob_mut:
                contacts2.append((i, j))
    # a bit removing and inserting columns
    changes = 0
    p_insert = 0.5
    p_delete = 0.5
    inserted = []
    deleted = []
    while changes < num_changes:
        rnd = random()
        my_j = int(random()*size2)
        if rnd < p_insert:
            inserted.append(my_j)
            size2 += 1
            for i, d in enumerate(inserted):
                if d >= inserted:
                    inserted[i] = d + 1
            for i, (c1, c2) in enumerate(contacts2):
                if c1 >= my_j:
                    c1 += 1
                    contacts2[i] = (c1, c2)
                if c2 >= my_j:
                    contacts2[i] = (c1, c2 + 1)
            for i in xrange(j):
                if random() < prob_contact:
                    contacts2.append((i, j))
        elif rnd > 1 - p_delete:
            size2 -= 1
            deleted.append(my_j)
            for i, d in enumerate(deleted):
                if d >= deleted:
                    deleted[i] = d - 1
            for i, (c1, c2) in enumerate(contacts2):
                if c1 >= my_j:
                    c1 -= 1
                    contacts2[i] = (c1, c2)
                if c2 >= my_j:
                    contacts2[i] = (c1, c2 - 1)
            to_del = []
            for c, (i, j) in enumerate(contacts2):
                if j == my_j:
                    to_del.append(c)
            for c in reversed(to_del):
                del(contacts2[c])
        else:
            continue
        changes += 1
    print 'inserted', inserted
    print 'deleted', deleted
    for i in range(size1+1):
        if not (i, i) in contacts1:
            contacts1.append((i, i))
    for i in range(size2+1):
        if not (i, i) in contacts2:
            contacts2.append((i, i))
    tad1 = contact2matrix(contacts1, size1)
    tad2 = contact2matrix(contacts2, size2)
    return tad1, tad2, contacts1, contacts2


def contact2matrix(contacts, size):
    matrix = [[0  for _ in xrange(size+1)] for _ in xrange(size+1)]
    for i, j in contacts:
        matrix[i][j] = 1
        matrix[j][i] = 1
    return np.array(matrix)


if __name__ == "__main__":
    exit(main())
