"""
10 jun 2015
"""

import sys

mapper = int(sys.argv[1])
win = int(sys.argv[2])

if win == 1:
    r_beg1 = [1  ] * len(range(25 , 105, 5))
    r_end1  = range(25 , 105, 5)
    r_beg2 = [101] * len(range(125, 205, 5))
    r_end2  = range(125, 205, 5)
if win == 2:
    r_beg1 = [1  ] * len(range(30 , 105, 5))
    r_end1  = range(30 , 105, 5)
    r_beg2 = [101] * len(range(130, 205, 5))
    r_end2  = range(130, 205, 5)
elif win == 3:
    r_beg1 = [1  ] * len(range(25 , 105, 10))
    r_end1  = range(25 , 105, 10)
    r_beg2 = [101] * len(range(125, 205, 10))
    r_end2  = range(125, 205, 10)
elif win == 4:
    r_beg1 = [1, 1]
    r_end1  = [50, 100]
    r_beg2 = [101, 101]
    r_end2  = [150, 200]
elif win == 5:
    r_beg1 = [1]
    r_end1  = [100]
    r_beg2 = [101]
    r_end2  = [200]

fastq          = '/scratch/db/FASTQs/hsap/dixon_2012/dixon-2012_200bp.fastq'
fastq          = 'mid_dixon-2012_200bp.fastq'
# fastq          = 'short_dixon-2012_200bp.fastq'
# fastq        = '/scratch/test/sample_dataset/FASTQs/sample_hsap_HindIII.fastq'
gem_index_path = '/scratch/db/index_files/Homo_sapiens-79/Homo_sapiens.gem'
out_map_dir1   = 'read1_%s-%s/' % (mapper, win)
out_map_dir2   = 'read2_%s-%s/' % (mapper, win)
temp_dir1      = 'tmp1_%s-%s/' % (mapper, win)
temp_dir2      = 'tmp2_%s-%s/' % (mapper, win)

# print 'read 1'
# mapping(gem_index_path, fastq, out_map_dir1, 'HindIII',
#         temp_dir=temp_dir1, windows=((0,100),))
# print 'read 2'
# mapping(gem_index_path, fastq, out_map_dir2, 'HindIII',
#         temp_dir=temp_dir2, windows=((100, 200),))

from pytadbit.mapping.mapper import iterative_mapping  
from pytadbit.mapping.full_mapper import full_mapping 
from pytadbit.parsers.map_parser import parse_map
from pytadbit.parsers.sam_parser import parse_sam
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.mapping.mapper import get_intersection
from pytadbit.mapping.filter import filter_reads, apply_filter

if mapper == 1:
    print 'read 1'
    outfiles1 = iterative_mapping(gem_index_path, fastq, out_map_dir1,
                                  r_beg1, [e + 2 for e in r_end1],
                                  temp_dir=temp_dir1)
    print 'read 2'
    outfiles2 = iterative_mapping(gem_index_path, fastq, out_map_dir2,
                                  r_beg2, [e + 2 for e in r_end2],
                                  temp_dir=temp_dir2)
    parse_thing = parse_sam
elif mapper == 2:
    print 'read 1'
    outfiles1 = full_mapping(gem_index_path, fastq, out_map_dir1, 'HindIII',
                             temp_dir=temp_dir1, frag_map=False,
                             windows=(zip(*(r_beg1, r_end1))))
    print 'read 2'
    outfiles2 = full_mapping(gem_index_path, fastq, out_map_dir2, 'HindIII',
                             temp_dir=temp_dir2, frag_map=False,
                             windows=(zip(*(r_beg2, r_end2))))
    parse_thing = parse_map
elif mapper == 3:
    print 'read 1'
    outfiles1 = full_mapping(gem_index_path, fastq, out_map_dir1, 'HindIII',
                             temp_dir=temp_dir1,
                             windows=(zip(*(r_beg1, r_end1))))
    print 'read 2'
    outfiles2 = full_mapping(gem_index_path, fastq, out_map_dir2, 'HindIII',
                             temp_dir=temp_dir2,
                             windows=(zip(*(r_beg2, r_end2))))
    parse_thing = parse_map

read1, read2 = 'read1.tsv_%s-%s' % (mapper, win), 'read2.tsv_%s-%s' % (mapper, win)
parse_thing(outfiles1, outfiles2, out_file1=read1, out_file2=read2,
            genome_seq=parse_fasta('/scratch/db/index_files/Homo_sapiens-79/Homo_sapiens.fa'),
            re_name='HindIII', verbose=True)

reads = 'both_reads.tsv_%s-%s' % (mapper, win)
get_intersection(read1, read2, reads)

masked = filter_reads(reads)
freads = 'filtered_reads.tsv_%s-%s' % (mapper, win)
apply_filter(reads, freads, masked)
