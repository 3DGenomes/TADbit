from pytadbit.mapping.mapper import iterative_mapping
from pytadbit.parsers.sam_parser import parse_sam
import sys, os
from pytadbit.mapping import map_re_sites
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.mapping.mapper import get_intersection

chunk = 10000000
chunk = int(sys.argv[1])


PATH   = '/home/fransua/Documents/Courses/given/2014_CSDM/notebooks/'
INFILE = '/home/fransua/Documents/Courses/given/2014_CSDM/notebooks/fastq/%s.fastq'
rep = 'SRR_test'
INFILE = INFILE % rep

OUTPATH = PATH + rep + '_' + str(chunk) + '/'

chr_names = ['2L', '2R', '3L', '3R', '4', 'X']
genome_seq = parse_fasta([PATH + 'dmel_reference/chr%s.fa' % crm for crm in chr_names], chr_names)

frags = map_re_sites('HindIII', genome_seq, verbose=True)

sams1 = iterative_mapping(
            gem_index_path       = PATH + 'dmel_reference/dm3.genome.gem',
            fastq_path           = INFILE,
            out_sam_path         = OUTPATH + '%s_r1.txt' % rep,
            temp_dir             = PATH + 'tmp_dir/',
            range_start          = [10] * 5, # starts with a flag sequence
            range_stop           = range(30, 55, 5),
            nthreads             = 8,  # on intel corei7 CPUs 4 threads are as fast as
                                       # 8, but leave some room for you other applications
            max_reads_per_chunk  = chunk,
            single_end           = True)
print 'created thes SAM files:', sams1

sams2 = iterative_mapping(
            gem_index_path       = PATH + 'dmel_reference/dm3.genome.gem',
            fastq_path           = INFILE,
            out_sam_path         = OUTPATH + '%s_r2.txt' % rep,
            temp_dir             = PATH + 'tmp_dir/',
            range_start          = range(80, 55, -5), # starts with a flag sequence
            range_stop           = [100] * 5,
            nthreads             = 8,  # on intel corei7 CPUs 4 threads are as fast as
                                       # 8, but leave some room for you other applications
            max_reads_per_chunk  = chunk,
            single_end           = True)

sams1 = [OUTPATH + fnam for fnam in os.listdir(OUTPATH) if fnam.rsplit('.', 2)[0].endswith('_r1.txt')]
sams2 = [OUTPATH + fnam for fnam in os.listdir(OUTPATH) if fnam.rsplit('.', 2)[0].endswith('_r2.txt')]

print 'created thes SAM files:', sams2

parse_sam(sams1, sams2, frags, 
          OUTPATH + 'reads1_%s.tsv' % rep, OUTPATH + 'reads2_%s.tsv' % rep,
          genome_seq, 'HindIII', verbose=True)

reads1 = OUTPATH + 'reads1_%s.tsv' % rep
reads2 = OUTPATH + 'reads2_%s.tsv' % rep
reads  = OUTPATH + 'reads12_%s.tsv' % rep

get_intersection(reads1, reads2, reads, verbose=True)

from pytadbit.mapping.analyze import hic_map
hic_map(reads, genome_seq, resolution=100000, savedata='lala')

from pytadbit.mapping.analyze import plot_genomic_distribution

plot_genomic_distribution(reads, resolution=50000, genome_seq=genome_seq) # because I know it

#  690691 SRR_test_10000000/reads1_SRR_test.tsv
#  690691 SRR_test_100000/reads1_SRR_test.tsv
# 1242927 SRR_test_200000/reads1_SRR_test.tsv
# 1035866 SRR_test_500000/reads1_SRR_test.tsv
