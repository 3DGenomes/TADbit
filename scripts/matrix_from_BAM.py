#! /usr/bin/python

"""
extract subset-matrix from a BAM file, and eventually normalizes it using
 precomputed biases
 From list of sub-matrices to extract, use index list to write in a file of nGb
 Check that not pass character limit
 Not output raw matrix and normalized with decay
"""

from argparse                        import ArgumentParser
from cPickle                         import load
from sys                             import stdout

from pytadbit.utils.file_handling    import mkdir
from pytadbit.mapping.filter         import MASKED
from pytadbit.parsers.hic_bam_parser import filters_to_bin, write_matrix


def main():
    opts          = get_options()
    inbam          = opts.inbam
    resolution     = opts.reso
    filter_exclude = opts.filter
    ncpus          = opts.cpus
    if opts.biases:
        biases     = load(open(opts.biases))
    else:
        biases     = {}
    outdir         = opts.outdir
    tmpdir         = opts.tmpdir
    coord1         = opts.coord1
    coord2         = opts.coord2

    if biases and biases['resolution'] != resolution:
        raise Exception('ERROR: different resolution in bias file (you want %d,'
                        ' there is %d).\n' % (resolution, biases['resolution']))
    if coord2 and not coord1:
        coord1, coord2 = coord2, coord1

    if not coord1:
        region1 = None
        start1  = None
        end1    = None
        region2 = None
        start2  = None
        end2    = None
    else:
        try:
            crm1, pos1   = coord1.split(':')
            start1, end1 = pos1.split('-')
            region1 = crm1
            start1  = int(start1)
            end1    = int(end1)
        except ValueError:
            region1 = coord1
            start1  = None
            end1    = None
        if coord2:
            try:
                crm2, pos2   = coord2.split(':')
                start2, end2 = pos2.split('-')
                region2 = crm2
                start2  = int(start2)
                end2    = int(end2)
            except ValueError:
                region2 = coord2
                start2  = None
                end2    = None
        else:
            region2 = None
            start2  = None
            end2    = None

    mkdir(outdir)
    mkdir(tmpdir)
    if region1:
        if region1:
            if not opts.quiet:
                stdout.write('\nExtraction of %s' % (region1))
            if start1:
                if not opts.quiet:
                    stdout.write(':%s-%s' % (start1, end1))
            else:
                if not opts.quiet:
                    stdout.write(' (full chromosome)')
            if region2:
                if not opts.quiet:
                    stdout.write(' intersection with %s' % (region2))
                if start2:
                    if not opts.quiet:
                        stdout.write(':%s-%s\n' % (start2, end2))
                else:
                    if not opts.quiet:
                        stdout.write(' (full chromosome)\n')
            else:
                if not opts.quiet:
                    stdout.write('\n')
    else:
        if not opts.quiet:
            stdout.write('\nExtraction of full genome\n')

    write_matrix(inbam, resolution, biases, outdir,
                 filter_exclude=filter_exclude,
                 normalizations=opts.matrices,
                 region1=region1, start1=start1, end1=end1,
                 region2=region2, start2=start2, end2=end2,
                 append_to_tar=opts.tarfile,
                 ncpus=ncpus, tmpdir=tmpdir, verbose=not opts.quiet)


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbam', metavar='',
                        required=True, default=False, help='input HiC-BAM file.')
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='',
                        default=True, help='output directory.')
    parser.add_argument('-t', '--tarfile', dest='tarfile', metavar='',
                        default=False, help='''skip the generation of files, directly
                        append them to a tar file
                        (does not need to be created).''')
    parser.add_argument('--tmp', dest='tmpdir', metavar='',
                        default=False, help='''path where to store temporary
                        files (by default outdir is used).''')
    parser.add_argument('-r', '--resolution', dest='reso', type=int, metavar='',
                        required=True, help='''wanted resolution form the
                        generated matrix''')
    parser.add_argument('-b', '--biases', dest='biases', metavar='',
                        help='''path to pickle file with array of biases''')
    parser.add_argument('-c', '--coord', dest='coord1',  metavar='',
                        default=None, help='''Coordinate of the region to
                        retrieve. By default all genome, arguments can be
                        either one chromosome name, or the coordinate in
                        the form: "-c chr3:110000000-120000000"''')
    parser.add_argument('-c2', '--coord2', dest='coord2',  metavar='',
                        default=None, help='''Coordinate of a second region to
                        retrieve the matrix in the intersection with the first
                        region.''')
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
