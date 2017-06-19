
from argparse                        import ArgumentParser
from cPickle                         import load
from pytadbit.parsers.hic_bam_parser import filters_to_bin, write_matrix
from pytadbit.parsers.hic_bam_parser import _write_small_matrix
from pytadbit.mapping.filter         import MASKED



def main():
    opts          = get_options()
    biases = load(open(opts.biases))
    resolution = biases['resolution']

    first_line = opts.first_line
    last_line = opts.last_line

    for count, line in enumerate(open(opts.list_coords)):
        if not first_line <= count < last_line:
            continue
        crm1, beg1, end1, crm2, beg2, end2 = line.split()
        beg1, end1, beg2, end2 = map(int, (beg1, end1, beg2, end2))

        # _write_small_matrix(opts.inbam, resolution, biases,
        #                     opts.outdir,
        #                     filter_exclude=opts.filter,
        #                     normalizations=opts.matrices,
        #                     region1=crm1, start1=beg1, end1=end1,
        #                     region2=crm2, start2=beg2, end2=end2,
        #                     append_to_tar=opts.tarfile,
        #                     ncpus=opts.cpus, tmpdir=opts.tmpdir,
        #                     verbose=not opts.quiet)

        write_matrix(opts.inbam, resolution, biases,
                     opts.outdir,
                     filter_exclude=opts.filter,
                     normalizations=opts.matrices,
                     region1=crm1, start1=beg1, end1=end1,
                     region2=crm2, start2=beg2, end2=end2,
                     append_to_tar=opts.tarfile,
                     ncpus=opts.cpus, tmpdir=opts.tmpdir,
                     verbose=not opts.quiet)


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbam', metavar='',
                        required=True, default=False, help='input HiC-BAM file.')
    parser.add_argument('-l', '--list_coords', dest='list_coords', metavar='',
                        required=True, default=False,
                        help='''input file, with list of coords (in the form:
                        "chr1 34000 44000 chr1 77000 87000").''')
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='',
                        default=True, help='output directory.')
    parser.add_argument('-t', '--tarfile', dest='tarfile', metavar='',
                        default=False, help='''skip the generation of files, directly
                        append them to a tar file
                        (does not need to be created).''')
    parser.add_argument('--first_line', dest='first_line', type=int,
                        default=0,
                        help='[%(default)s] truncate input list to start at given line')
    parser.add_argument('--last_line', dest='last_line', type=int,
                        default=float('inf'),
                        help='[%(default)s] truncate input list to stop at given line')
    parser.add_argument('--tmp', dest='tmpdir', metavar='',
                        default=False, help='''path where to store temporary
                        files (by default outdir is used).''')
    parser.add_argument('-b', '--biases', dest='biases', metavar='',
                        help='''path to pickle file with array of biases''')
    parser.add_argument('-C', '--cpus', dest='cpus', metavar='', type=int,
                        default=8, help='''[%(default)s] number of cpus to be
                        used for parsing the HiC-BAM file''')
    parser.add_argument('--matrices', dest='matrices', metavar='', type=str,
                        nargs='+', default=['norm', 'raw', 'decay'],
                        help='''[%(default)s] which matrix to generate''')
    parser.add_argument('-f', '--format', dest='format', default='abc',
                        choices=['abc', 'mat'], required=False, help='''[%(default)s]
                        format in which to write the output matrix (choose from %(choices)s)''')
    parser.add_argument('-q', '--quiet', dest='quiet', default=False, action='store_true',
                        help='display no running information')
    parser.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 8, 9, 10],
                        choices = range(1, 11),
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))

    opts = parser.parse_args()
    # convert filters to binary for samtools
    opts.filter = filters_to_bin(opts.filter)
    if not opts.biases and ('norm' in opts.matrices or
                            'decay' in opts.matrices):
        raise Exception('ERROR: should provide path to bias file.')
    if not opts.tmpdir:
        opts.tmpdir = opts.outdir

    return opts


if __name__=='__main__':
    exit(main())
