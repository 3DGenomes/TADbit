#! /usr/bin/python

"""
05 May 2015


"""
# MatPlotLib not asking for X11
import matplotlib as mpl
mpl.use('Agg')

from argparse import ArgumentParser, HelpFormatter

from pytadbit.mapping.mapper import iterative_mapping
from pytadbit.utils.fastq_utils import quality_plot
from pytadbit.parsers.genome_parser import parse_fasta


def main():
    opts = get_options()
    iterative_mapping()

    ## PARSE FASTA
    genome = parse_fasta(opts.fasta if len(opts.fasta) <= 1 else opts.fasta[0],
                         chr_names=opts.chr_name, verbose=True)

    ## CREATE INDEX


def get_options():
    """
    parse option from call

    """
    parser = ArgumentParser(
        usage="%(prog)s [options] [--cfg CONFIG_PATH]",
        formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                   max_help_position=27))
    glopts = parser.add_argument_group('General arguments')
    taddet = parser.add_argument_group('TAD detection arguments')
    optimo = parser.add_argument_group('Optimization of IMP arguments')
    modelo = parser.add_argument_group('Modeling with optimal IMP arguments')
    descro = parser.add_argument_group('Descriptive, optional arguments')
    analyz = parser.add_argument_group('Output arguments')

    parser.add_argument('--cfg', dest='cfg', metavar="PATH", action='store',
                      default=None, type=str,
                      help='path to a configuration file with predefined ' +
                      'parameters')
    parser.add_argument('--analyze_only', dest='analyze_only',
                        action='store_true', default=False,
                        help=('load precomputed models in outdir, ' +
                              'skip optimization, modeling'))
    parser.add_argument('--optimize_only', dest='optimize_only', default=False,
                        action='store_true',
                        help='do the optimization of the region and exit')

    #########################################
    # GENERAL
    glopts.add_argument(
        '--root_path', dest='root_path', metavar="PATH",
        default='', type=str,
        help=('path to search for data files (just pass file name' +
              'in "data")'))
    glopts.add_argument('--data', dest='data', metavar="PATH", nargs='+',
                        type=str,
                        help='''paths to file(s) with FASTQ files. If many,
                        files will be summed up. I.e.: --data
                        replicate_1.txt replicate_2.txt''')
    glopts.add_argument('--fasta', dest='fasta', metavar="PATH", nargs='+',
                        type=str,
                        help='''paths to file(s) with FASTA files of the
                        reference genome. If many, files will be concatenated.
                        I.e.: --fasta chr_1.fa chr_2.fa
                        In this last case, order is important or the rest of the
                        analysis.''')
    glopts.add_argument('--chr_name', dest='chr_name', metavar="STR", nargs='+',
                        default=[], type=str,
                        help='''[fasta header] chromsome name(s). Used in the
                        same order as data.''')

    # update fasta paths adding root directory
    if opts.root_path and opts.fasta[0]:
        for i in xrange(len(opts.fasta)):
            logging.info(os.path.join(opts.root_path, opts.fasta[i]))
            opts.data[i] = os.path.join(opts.root_path, opts.fasta[i])

    return opts


if __name__ == "__main__":
    exit(main())
