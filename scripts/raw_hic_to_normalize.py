"""
17 Apr 2014


"""
from pytadbit import Chromosome
from argparse import ArgumentParser
from sys import stdout
from warnings import warn


def main():
    """
    main function
    """

    opts = get_options()
    crm = Chromosome(':P')

    for i, data in enumerate(opts.data):
        crm.add_experiment('exp' + str(i), resolution=int(opts.resolution[i]),
                           hic_data=data)
        crm.experiments['exp' + str(i)].normalize_hic()

    if len(opts.data) > 1:
        exp = crm.experiments[0] + crm.experiments[1]
        for i in range(2, len(opts.data)):
            exp += crm.experiments[i]
    else:
        exp = crm.experiments[0]

    if opts.abc:
        exp.write_interaction_pairs(opts.output, normalized=opts.norm,
                                    zscored=False)
    else:
        if type(opts.output) == file:
            out = opts.output
        else:
            out = open(opts.output, 'w')
        out.write(exp.print_hic_matrix(print_it=False,
                                       normalized=opts.norm))

def get_options():
    '''
    parse option from call
    '''
    parser = ArgumentParser(
        usage=("%(prog)s [options] file [options] file [options] " +
               "file [options [file ...]]"))
    parser.add_argument('-d', '--data', dest='data', metavar="PATH",
                      action='append', default=None, nargs='+',
                      help='''path to a file containing Hi-C data in matrix
                      format. If used several times, experiments will be
                      summed up. I.e.: --data hic_replicate_1.txt
                      hic_replicate_2.txt''')
    parser.add_argument('-r', '--resolution', dest='resolution', metavar="int",
                      action='store', default=None, nargs='+',
                      help='''resolution of Hi-C experiments passed using the
                      hic_files option. Use same order as "data" argument.''')
    parser.add_argument('--output', dest='output', action="store", metavar="PATH", 
                      default=stdout, help=
                      '''[stdout] path to out-file where to store result''')
    parser.add_argument('--abc', dest='abc', action='store_true',
                      help='Result writen in column format')
    parser.add_argument('--raw', dest='raw', action='store_true', default=False,
                      help='[%(default)s] do not normalize the data.')
    opts = parser.parse_args()
    if not opts.data:
        exit(parser.print_help())
    if not opts.resolution:
        warn('ERROR: should provide resolution')
        exit(parser.print_help())
    opts.norm = not opts.raw
    if len(opts.resolution) == 1 and len(opts.data) > 1:
        opts.resolution *= len(opts.data)
    return opts


if __name__ == "__main__":
    exit(main())
