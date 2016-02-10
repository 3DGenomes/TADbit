"""

information needed

 - path working directory with parsed reads

"""

from argparse import HelpFormatter
from pytadbit import get_dependencies_version
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.parsers.map_parser import parse_map
from os import system, path
import logging
import fcntl
from cPickle import load, UnpicklingError

DESC = "Parse mapped Hi-C reads and get the intersection"

def get_out_files(workdir):
    fnames1 = []
    fnames2 = []
    for line in open(path.join(workdir, 'trace.log')):
        if   line.startswith('# MAPPED READ1 '):
            fnames1.append(line.split()[-1])
        elif line.startswith('# MAPPED READ2 '):
            fnames2.append(line.split()[-1])
    return fnames1, fnames2

def run(opts):
    check_options(opts)

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    # glopts.add_argument('--qc_plot', dest='quality_plot', action='store_true',
    #                   default=False,
    #                   help='generate a quality plot of FASTQ and exits')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')


    parser.add_argument_group(glopts)

def check_options(opts):

    if not opts.workdir: raise Exception('ERROR: output option required.')
