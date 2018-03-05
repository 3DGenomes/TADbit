"""
Computes the insulation index using different windows.

Insulation index assigned to bins at a given resolution, not between bins.
"""

from argparse                    import ArgumentParser
from cPickle                     import load

from pytadbit.parsers.hic_parser import load_hic_data_from_bam
from pytadbit.mapping.filter     import MASKED


def main():
    opts          = get_options()
    biases = load(open(opts.biases))
    resolution = biases['resolution']

    hic_data = load_hic_data_from_bam(opts.inbam, resolution, biases=biases,
                                      filter_exclude=opts.filter,
                                      verbose=not opts.quiet,
                                      tmpdir=opts.tmpdir, ncpus=opts.cpus)

    insidx = {}
    bias = hic_data.bias
    bads = hic_data.bads
    decay = hic_data.expected
    for dist, end in opts.dists:
        if not opts.quiet:
            print ' - computing insulation in band %d-%d' % (dist, end)
        insidx[(dist, end)] = {}
        for crm in hic_data.chromosomes:
            for pos in range(hic_data.section_pos[crm][0] + end,
                             hic_data.section_pos[crm][1] - end):
                insidx[(dist, end)][pos] = sum(
                    hic_data[i, j] / bias[i] / bias[j] / decay[crm][abs(j-i)]
                    for i in range(pos - end, pos - dist + 1)
                    if not i in bads
                    for j in range(pos + dist, pos + end + 1)
                    if not j in bads)
    out = open(opts.outfile, 'w')
    out.write('# CRM\tCOORD\t' + '\t'.join(['%d-%d' % (d1, d2)
                                            for d1, d2 in opts.dists]) +
              '\n')

    for crm in hic_data.section_pos:
        for pos in range(*hic_data.section_pos[crm]):
            beg = (pos - hic_data.section_pos[crm][0]) * resolution
            out.write('{}\t{}-{}\t{}\n'.format(
                crm, beg + 1, beg + resolution,
                '\t'.join([str(insidx[dist].get(pos, 'NaN'))
                           for dist in opts.dists])))
    out.close()


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbam', metavar='',
                        required=True, default=False, help='input HiC-BAM file')
    parser.add_argument('-l', '--list_dists', dest='dists', type=str,
                        default=['2,2', '4,5', '8,11', '14,19'], nargs='+',
                        help='''[%(default)s] list of pairs of distances between
                        which to compute the insulation index. E.g. 4,5 means
                        that for a given bin B(i), all interactions between
                        B(i-4) to B(i-5) and B(i+4) to B(i+5) will be summed and
                        used to compute the insulation index''')
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='',
                        required=True, default=True, help='path to output file')
    parser.add_argument('--tmp', dest='tmpdir', metavar='',
                        default='.', help='''path where to store temporary
                        files.''')
    parser.add_argument('-b', '--biases', dest='biases', metavar='',
                        required=True,
                        help='''path to pickle file with array of biases''')
    parser.add_argument('-C', '--cpus', dest='cpus', metavar='', type=int,
                        default=8, help='''[%(default)s] number of cpus to be
                        used for parsing the HiC-BAM file''')
    parser.add_argument('-q', '--quiet', dest='quiet', action='store_true',
                        default=False, help='display no running information')
    parser.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT',choices = range(1, 11),
                        default=[1, 2, 3, 4, 6, 7, 9, 10],
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))

    opts = parser.parse_args()
    opts.dists = [tuple(map(int, d.split(','))) for d in opts.dists]
    if not all([len(d) == 2 for d in opts.dists]):
        raise Exception('ERROR: distance should be input by pairs.')
    return opts

if __name__ == '__main__':
    exit(main())
