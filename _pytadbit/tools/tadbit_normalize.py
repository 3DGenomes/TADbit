"""

information needed

 - path working directory with parsed reads

"""
from argparse                    import HelpFormatter
from os import path
import sqlite3 as lite
import time

DESC = 'normalize Hi-C data and write results to file as matrices'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()

    finish_time = time.localtime()

    save_to_db(opts, launch_time, finish_time)

def save_to_db(opts,
               launch_time, finish_time):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()

def load_parameters_fromdb(opts):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        # get the JOBid of the parsing job
        cur.execute("""
        select distinct Id from JOBs
        where Type = 'Parse'
        """)
        jobids = cur.fetchall()
        if len(jobids) > 1:
            cur.execute("delete from JOBs where Id = jobids")
        parse_jobid = jobids[0][0]
        # fetch path to parsed BED files
        cur.execute("""
        select distinct Path from PATHs
        where JOBid = %d and Type = 'BED'
        """ % parse_jobid)
        fname1, fname2 = [path.join(opts.workdir, e[0]) for e in cur.fetchall()]


def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    parser.add_argument_group(glopts)

def check_options(opts):

    if not opts.workdir: raise Exception('ERROR: output option required.')

    # check resume
    if not path.exists(opts.workdir) and opts.resume:
        print ('WARNING: can use output files, found, not resuming...')
        opts.resume = False

    # sort filters
    if opts.apply:
        opts.apply.sort()
