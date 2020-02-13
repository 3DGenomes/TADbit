"""
21 Jan 2013


"""
from __future__ import print_function

import numpy as np
from random import random
from copy import deepcopy
from pytadbit.tad_clustering.tad_cmo import optimal_cmo

def main():
    """
    main function
    """

    tad1, tad2, contacts1, contacts2 = generate_random_contacts(size=12,
                                                                prob_mut=0.5,
                                                                num_changes=2)
    size1 = len(tad1)
    size2 = len(tad2)#max([max(c) for c in contacts2])
    num_v = min(size1, size2)
    #align1a, align2a, scorea = run_aleigen(contacts1,
    #                                       contacts2, num_v)
    align1, align2, score = optimal_cmo(tad1, tad2, num_v)
    align1b = []
    align2b = []
    for c in range(len(align1)):
        if align1[c] != '-' and align2[c] != '-':
            align1b.append(align1[c])
            align2b.append(align2[c])
    print('='*80)
    print(''.join(['{:>5}'.format(i) for i in range(len(tad1))]))
    print('\n'.join([''.join(['{:>5}'.format(round(j, 2)) for j in i]) for i in tad1]))
    print('-'*80)
    print(''.join(['{:>5}'.format(i) for i in range(len(tad2))]))
    print('\n'.join([''.join(['{:>5}'.format(round(j, 2)) for j in i]) for i in tad2]))
    print('='*80)
    #print ''.join(['{:>4}'.format(a) for a in align1a])
    #print ''.join(['{:>4}'.format(a) for a in align2a])
    print('-'*80)
    print(''.join(['{:>4}'.format(a) for a in align1b]))
    print(''.join(['{:>4}'.format(a) for a in align2b]))
    print('='*80)


def generate_random_contacts(size=10, prob_contact=0.7, 
                             prob_mut=0.1, num_changes=2):
    """
    returns 2 random matrices corresponding to TADs
    and 2 contacts lists.
    """
    size1 = size2  = size# + int(random()*5)
    contacts1 = []
    for i in range(size1):
        for j in range(i+1, size1):
            if random() < prob_contact:
                contacts1.append((i, j, random()))
            else:
                contacts1.append((i, j, 0))
    contacts2 = deepcopy(contacts1)
    # randomly change some contacts
    for c in range(len(contacts2) -1, -1, -1):
        if random() < prob_mut:
            contacts2[c] = (contacts2[c][0], contacts2[c][1],
                            contacts2[c][2] + random()/2 - contacts2[c][2]/2)
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
            size2 += 1
            for i, d in enumerate(inserted):
                if d >= my_j:
                    inserted[i] = d + 1
            inserted.append(my_j)
            deleted.append(-1)
            for i, (c1, c2, c3) in enumerate(contacts2):
                if c1 >= my_j:
                    c1 += 1
                    contacts2[i] = (c1, c2, c3)
                if c2 >= my_j:
                    contacts2[i] = (c1, c2 + 1, c3)
            for i in range(size2):
                if random() < prob_contact:
                    contacts2.append((i, my_j, random()))
                else:
                    contacts2.append((i, my_j, 0))
        elif rnd > 1 - p_delete:
            size2 -= 1
            for i, d in enumerate(deleted):
                if d > my_j:
                    deleted[i] = d - 1
            deleted.append(my_j)
            inserted.append(-1)
            to_del = []
            for c, (_, j, _) in enumerate(contacts2):
                if j == my_j:
                    to_del.append(c)
            for c in reversed(to_del):
                del(contacts2[c])
            for i, (c1, c2, c3) in enumerate(contacts2):
                if c1 >= my_j:
                    c1 -= 1
                    contacts2[i] = (c1, c2, c3)
                if c2 >= my_j:
                    contacts2[i] = (c1, c2 - 1, c3)
        else:
            continue
        changes += 1
    print('inserted', inserted)
    print('deleted ', deleted)
    for i in range(size1):
        if not (i, i) in contacts1:
            contacts1.append((i, i, 1))
    for i in range(size2):
        if not (i, i) in contacts2:
            contacts2.append((i, i, 1))
    tad1 = contact2matrix(contacts1, size1)
    tad2 = contact2matrix(contacts2, size2)
    return tad1, tad2, contacts1, contacts2


def contact2matrix(contacts, size):
    matrix = [[0  for _ in range(size)] for _ in range(size)]
    for i, j, k in contacts:
        matrix[i][j] = k
        matrix[j][i] = k
    return np.array(matrix)


if __name__ == "__main__":
    exit(main())
