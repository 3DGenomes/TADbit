#! /usr/bin/env python
"""
Takes TADbit bam and converts it to HICUP-like BAM format (standard BAM)
"""

import os

from argparse     import ArgumentParser

from pysam                        import AlignmentFile

from pytadbit.mapping.filter          import MASKED
from pytadbit.parsers.hic_bam_parser  import filters_to_bin
from pytadbit.utils                   import printime
from pytadbit.utils.file_handling     import magic_open


# def get_mapped_chunk(fhandler, nreads, ncpus):
#     seqs = dict((i, {}) for i in range(ncpus))
#     pos_file = 0
#     for line in fhandler:
#         pos_file += 1
#         rid, seq, qal, _, pos = line.split()
#         pos = int(pos.split(':')[2])
#         seqs[pos_file % ncpus][rid, pos] = (seq, qal)
#         if pos_file >= nreads:
#             yield seqs
#             seqs = dict((i, {}) for i in range(ncpus))
#     return seqs

def get_mapped_chunk(map_folder, nreads):
    seqs = {}
    printime(' - loading chunk')
    pos_file = 0
    for fname in os.listdir(map_folder):
        printime('    - ' + fname)
        fhandler = magic_open(os.path.join(map_folder, fname))
        for line in fhandler:
            pos_file += 1
            rid, seq, qal, _, pos = line.split()
            pos = int(pos.split(':')[2])
            rid = rid.split('~')[0]
            seqs[rid, pos] = (seq, qal)
            if pos_file >= nreads:
                yield seqs
                printime(' - loading chunk')
                seqs = {}
                pos_file = 0
    yield seqs


def main():
    """
    main function
    """
    opts           = get_options()
    filter_exclude = filters_to_bin(opts.filter)
    tadbit_bam     = opts.tadbit_bam
    hicup_bam      = opts.hicup_bam
    map_folder     = opts.map_folder
    nreads         = opts.nreads * 1_000_000

    tag_dict = {
        (1, 1): ( 67, 131),
        (0, 0): (115, 179),
        (1, 0): ( 99, 147),
        (0, 1): ( 83, 163),
    }

    out = open(hicup_bam, 'w')
    for seqs in get_mapped_chunk(map_folder, nreads):
        bamfile = AlignmentFile(tadbit_bam, 'rb')
        refs = bamfile.references
        printime(f' - processing BAM (for {len(seqs) / 1_000_000}M reads)')
        for r in bamfile.fetch(multiple_iterators=False):
            if r.flag & filter_exclude:
                continue
            rid  = r.qname
            ridname = rid.split('#')[0]
            pos1 = r.reference_start + 1
            which, len1 = r.cigar[0]
            tags = dict(r.tags)
            if which == 6:  # first read-end
                s1, s2 = tags['S1'], tags['S2']
            else:
                s2, s1 = tags['S1'], tags['S2']
            if s1 == 0:
                pos1 = pos1 - len1 + 1
            try:
                seq, qal = seqs[ridname, pos1]
            except KeyError:
                continue
            crm1 = r.reference_name
            crm2 = refs[r.mrnm]
            pos2 = r.mpos + 1
            len2 = r.tlen

            dist = 0 if crm1 != crm2 else abs(pos2 - pos1)
            tags = dict(r.tags)

            if s2 == 0:
                pos2 = pos2 - len2 + 1

            flag = tag_dict[s1, s2][0]

            out.write(
                (
                    f'{r.qname}\t{flag}\t{crm1}\t{pos1}\t{len1}\t'
                    f'{len(seq)}M\t{crm2}\t{pos2}\t{dist}\t{seq}\t'
                    f'{qal}\tMD:Z:{len1}\tPG:Z:MarkDuplicates\tNM:i:0\t'
                    f'AS:i:{len1}\tXS:i:1\n')
            )
        bamfile.close()
        seqs.clear()
    out.close()



def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH [options]")

    parser.add_argument('-i', dest='tadbit_bam', metavar='PATH', required=True,
                        default=False, help='input TADbit BAM file')
    parser.add_argument('-M', dest='map_folder', metavar='PATH', required=True,
                        default=False, help='path to FOLDER with mapped reads.')
    parser.add_argument('-o', dest='hicup_bam', metavar='PATH', required=True,
                        help='output HICUP BAM file')
    parser.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                        choices = list(range(1, 11)),
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))
    parser.add_argument('-c', dest='cpus',  default=1, type=int,
                        help='[%(default)s] number of CPUs')
    parser.add_argument('-n', dest='nreads',  default=1, type=float,
                        help='''[%(default)s] number of million reads to process
                         in each batch (the more reads the more RAM used)''')
    opts = parser.parse_args()
    return opts


if __name__ == "__main__":
    exit(main())