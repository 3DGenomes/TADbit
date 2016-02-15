"""

information needed

 - path working directory with parsed reads

"""

from argparse import HelpFormatter
from pytadbit.mapping.mapper     import get_intersection
from os import path
import logging
import fcntl
import sqlite3 as lite
import time

DESC = "Parse mapped Hi-C reads and get the intersection"

def run(opts):
    check_options(opts)
    launch_time = time.localtime()

    fname1, fname2 = load_parameters_fromdb(opts.workdir)

    mreads = opts.workdir, '02_parsed_reads', ''

    get_intersection(fname1, fname2, mreads)

def load_parameters_fromdb(workdir):
    con = lite.connect(path.join(workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        # get the JOBid of the parsing job
        cur.execute("""
        select distinct Id from JOBs
        where Type = 'Parse'
        """)
        jobids = cur.fetchall()
        if len(jobids) > 1:
            raise NotImplementedError('ERROR: only one parsing per working '
                                      'directory supported')
        parse_jobid = jobids[0][0]
        # fetch path to parsed BED files
        cur.execute("""
        select distinct Path from PATHs
        where JOBid = %d and Type = 'BED'
        """ % parse_jobid)
        fname1, fname2 = [path.join(workdir, e[0]) for e in cur.fetchall()]

    return fname1, fname2


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
