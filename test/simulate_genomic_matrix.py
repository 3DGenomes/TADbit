"""
23 abr 2015

This script generates a false genome for testing purposes
"""
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES
from random import random
import re

# VARIABLES
num_crms      = 9
mean_crm_size = 1000000

r_enz         = 'DpnII'

selfcircle    = 1000
dangling      = 1000
error         = 1000
extradangling = 1000
duplicates    = 1000

re_seq = RESTRICTION_ENZYMES[r_enz].replace('|', '')
enz_cut = RESTRICTION_ENZYMES[r_enz].index('|')

nts = ('ATGC' * 100) +'N'

# RANDOM GENOME GENERATION
genome = {}
for crm in xrange(1, num_crms + 1):
    crm_len = int(mean_crm_size * random())
    genome['chr' + str(crm)] = ''.join([nts[int(401 * random())]
                                        for _ in xrange(crm_len)])

out = open('test.fa~', 'w')
for crm in xrange(1, num_crms + 1):
    out.write('>chr%d\n' % crm)
    crm = 'chr' + str(crm)
    for p in xrange(0, len(genome[crm]), 60):
        out.write(genome[crm][p:p+60] + '\n')
out.close()

from pytadbit.parsers.genome_parser import parse_fasta

genome_bis = parse_fasta('test.fa~')

if genome_bis == genome:
    genome = genome_bis
else:
    raise Exception('problem with genome parser')

# RE FRAGMENTS
frags = {}
for crm in genome:
    frags[crm] = {}
    beg = 0
    for pos in re.finditer(re_seq, genome[crm]):
        end = pos.start() + 1 + enz_cut
        if beg == end:
            continue
        frags[crm][beg] = [beg, end]
        beg = end
    if beg != end:
        frags[crm][beg] = [beg, len(genome[crm])]

# Prepare files
sam_crm = '@SQ\tSN:%s\tLN:%d\n'
sam_head = """@HD\tVN:1.3
%s@RG\tID:0\tPG:GEM\tPL:ILLUMINA\tSM:0
@PG\tID:GEM\tPN:gem-2-sam\tVN:1.847
"""
crm_heads = ''
for crm in genome:
    crm_heads += sam_crm % (crm, len(genome[crm]))

sam_head = sam_head % crm_heads

sam_read = '{id}\t{flag}\t{crm}\t{pos}\t254\t3M\t*\t0\t0\tAAA\tHHH\tRG:Z:0\tNH:i:1\tNM:i:0\tXT:A:U\tmd:Z:3\n'
map_read = "{id}\tAAA\tHHH\t1\t{crm}:{flag}:{pos}:3\n"

out1 = open('test_read1.sam~', 'w')
out1.write(sam_head)
out2 = open('test_read2.sam~', 'w')
out2.write(sam_head)

out11 = open('test_read1.map~', 'w')
out22 = open('test_read2.map~', 'w')

map_flags = ['+', '-']
sam_flags = [66 , 82 ]

# SELF-CIRCLES
for i in xrange(selfcircle):
    crm  = genome.keys()    [int(random() * len(genome))]
    while True:
        frag = frags[crm].keys()[int(random() * len(frags[crm]))]
        if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
            break
    while True:
        pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                   + frags[crm][frag][0] + 1)
        pos2 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                   + frags[crm][frag][0] + 1)
        if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
            continue
        if pos2 > pos1:
            sd1 = 1
            sd2 = 0
            pos1 -= 2
            break
        elif pos2 < pos1:
            sd1 = 0
            sd2 = 1
            pos2 -= 2
            break
    read1 = {'crm': crm, 'pos': pos1, 'flag': sam_flags[sd1], 'id': 'lala01.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': sam_flags[sd2], 'id': 'lala01.%012d' % (i)}
    out1.write(sam_read.format(**read1))
    out2.write(sam_read.format(**read2))
    read1 = {'crm': crm, 'pos': pos1, 'flag': map_flags[sd1], 'id': 'lala01.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': map_flags[sd2], 'id': 'lala01.%012d' % (i)}
    out11.write(map_read.format(**read1))
    out22.write(map_read.format(**read2))

# DANGLING-ENDS
for i in xrange(dangling):
    crm  = genome.keys()    [int(random() * len(genome))]
    while True:
        frag = frags[crm].keys()[int(random() * len(frags[crm]))]
        if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
            break
    while True:
        pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                   + frags[crm][frag][0] + 1)
        pos2 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                   + frags[crm][frag][0] + 1)
        if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
            continue
        if pos2 > pos1:
            sd1 = 0
            sd2 = 1
            pos2 -= 2
            break
        elif pos2 < pos1:
            sd1 = 1
            sd2 = 0
            pos1 -= 2
            break
    read1 = {'crm': crm, 'pos': pos1, 'flag': sam_flags[sd1], 'id': 'lala02.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': sam_flags[sd2], 'id': 'lala02.%012d' % (i)}
    out1.write(sam_read.format(**read1))
    out2.write(sam_read.format(**read2))
    read1 = {'crm': crm, 'pos': pos1, 'flag': map_flags[sd1], 'id': 'lala02.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': map_flags[sd2], 'id': 'lala02.%012d' % (i)}
    out11.write(map_read.format(**read1))
    out22.write(map_read.format(**read2))

# ERRORS
for i in xrange(error):
    crm  = genome.keys()    [int(random() * len(genome))]
    while True:
        frag = frags[crm].keys()[int(random() * len(frags[crm]))]
        if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
            break
    while True:
        pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                   + frags[crm][frag][0] + 1)
        pos2 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                   + frags[crm][frag][0] + 1)
        sd1 = sd2 = int(random()*2)
        if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
            continue
        if pos1 != pos2:
            if sd1 == 1:
                pos1 -= 2
                pos2 -= 2
            break
    read1 = {'crm': crm, 'pos': pos1, 'flag': sam_flags[sd1], 'id': 'lala03.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': sam_flags[sd2], 'id': 'lala03.%012d' % (i)}
    out1.write(sam_read.format(**read1))
    out2.write(sam_read.format(**read2))
    read1 = {'crm': crm, 'pos': pos1, 'flag': map_flags[sd1], 'id': 'lala03.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': map_flags[sd2], 'id': 'lala03.%012d' % (i)}
    out11.write(map_read.format(**read1))
    out22.write(map_read.format(**read2))

# EXTRA-DANGLING-ENDS
for i in xrange(extradangling):
    crm  = genome.keys()    [int(random() * len(genome))]
    while True:
        frag = frags[crm].keys()[int(random() * len(frags[crm]))]
        if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
            break
    while True:
        while True:
            pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                       + frags[crm][frag][0] + 1)
            if frags[crm][frag][1] - pos1 < 490 or pos1 - frags[crm][frag][0] < 490:
                break
        while True:
            pos2 = int(random() * 999) - 499 + pos1
            if not frags[crm][frag][0] <= pos2 <= frags[crm][frag][1]:
                break
        if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
            continue
        if pos2 > pos1:
            sd1 = 0
            sd2 = 1
            pos2 -= 2
            break
        elif pos2 < pos1:
            sd1 = 1
            sd2 = 0
            pos1 -= 2
            break
    read1 = {'crm': crm, 'pos': pos1, 'flag': sam_flags[sd1], 'id': 'lala04.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': sam_flags[sd2], 'id': 'lala04.%012d' % (i)}
    out1.write(sam_read.format(**read1))
    out2.write(sam_read.format(**read2))
    read1 = {'crm': crm, 'pos': pos1, 'flag': map_flags[sd1], 'id': 'lala04.%012d' % (i)}
    read2 = {'crm': crm, 'pos': pos2, 'flag': map_flags[sd2], 'id': 'lala04.%012d' % (i)}
    out11.write(map_read.format(**read1))
    out22.write(map_read.format(**read2))

# TOO CLOSE FROM RE

# VAID PAIRS

# DUPPLICATES
i = 0
while i < duplicates * 2:
    crm1  = genome.keys()    [int(random() * len(genome))]
    crm2  = genome.keys()    [int(random() * len(genome))]
    while True:
        frag1 = frags[crm1].keys()[int(random() * len(frags[crm1]))]
        frag2 = frags[crm2].keys()[int(random() * len(frags[crm2]))]
        if frags[crm2][frag2][1] - frags[crm2][frag2][0] < 6:
            continue
        if frags[crm1][frag1][1] - frags[crm1][frag1][0] < 6:
            continue
        if frag1 != frag2:
            break
    while True:
        pos1 = int(random() * (frags[crm1][frag1][1] - frags[crm1][frag1][0])
                   + frags[crm1][frag1][0])
        pos2 = int(random() * (frags[crm2][frag2][1] - frags[crm2][frag2][0])
                   + frags[crm2][frag2][0])
        if not (3 < pos1 < (len(genome[crm1]) - 3) and 3 < pos2 < (len(genome[crm2]) - 3)):
            continue
        break
    if pos2 > pos1:
        sd1 = 0
        sd2 = 1
        pos2 -= 2
    elif pos2 < pos1:
        sd1 = 1
        sd2 = 0
        pos1 -= 2
    read1 = {'crm': crm1, 'pos': pos1, 'flag': sam_flags[sd1], 'id': 'lala05.%012d' % (i)}
    read2 = {'crm': crm2, 'pos': pos2, 'flag': sam_flags[sd2], 'id': 'lala05.%012d' % (i)}
    out1.write(sam_read.format(**read1))
    out2.write(sam_read.format(**read2))
    read1 = {'crm': crm1, 'pos': pos1, 'flag': map_flags[sd1], 'id': 'lala05.%012d' % (i)}
    read2 = {'crm': crm2, 'pos': pos2, 'flag': map_flags[sd2], 'id': 'lala05.%012d' % (i)}
    out11.write(map_read.format(**read1))
    out22.write(map_read.format(**read2))
    i += 1
    if random() > 0.5:
        read2 = {'crm': crm1, 'pos': pos1, 'flag': sam_flags[sd1], 'id': 'lala05.1%011d' % (i)}
        read1 = {'crm': crm2, 'pos': pos2, 'flag': sam_flags[sd2], 'id': 'lala05.1%011d' % (i)}
        out1.write(sam_read.format(**read1))
        out2.write(sam_read.format(**read2))
        read2 = {'crm': crm1, 'pos': pos1, 'flag': map_flags[sd1], 'id': 'lala05.1%011d' % (i)}
        read1 = {'crm': crm2, 'pos': pos2, 'flag': map_flags[sd2], 'id': 'lala05.1%011d' % (i)}
        out11.write(map_read.format(**read1))
        out22.write(map_read.format(**read2))
    else:
        read2 = {'crm': crm1, 'pos': pos1, 'flag': sam_flags[sd1], 'id': 'lala05.1%011d' % (i)}
        read1 = {'crm': crm2, 'pos': pos2, 'flag': sam_flags[sd2], 'id': 'lala05.1%011d' % (i)}
        out1.write(sam_read.format(**read1))
        out2.write(sam_read.format(**read2))
        read2 = {'crm': crm1, 'pos': pos1, 'flag': map_flags[sd1], 'id': 'lala05.1%011d' % (i)}
        read1 = {'crm': crm2, 'pos': pos2, 'flag': map_flags[sd2], 'id': 'lala05.1%011d' % (i)}
        out11.write(map_read.format(**read1))
        out22.write(map_read.format(**read2))
    i += 1

out1.close()
out2.close()
out11.close()
out22.close()

# PARSE SAM
from pytadbit.parsers.map_parser import parse_map
from pytadbit.parsers.sam_parser import parse_sam

parse_map(['test_read1.map~'], ['test_read2.map~'],
          './lala1-map~', './lala2-map~', genome, re_name='DPNII', mapper='GEM')

parse_sam(['test_read1.sam~'], ['test_read2.sam~'],
          './lala1-sam~', './lala2-sam~', genome, re_name='DPNII', mapper='GEM')


# GET INTERSECTION
from pytadbit.mapping.mapper import get_intersection

get_intersection('lala1-map~', 'lala2-map~', 'lala-map~')
get_intersection('lala1-sam~', 'lala2-sam~', 'lala-sam~')

# FILTER
from pytadbit.mapping.filter import filter_reads

masked1 = filter_reads('lala-map~')

masked2 = filter_reads('lala-sam~')


