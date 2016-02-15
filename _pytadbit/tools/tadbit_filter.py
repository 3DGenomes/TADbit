"""

information needed

 - path working directory with parsed reads

"""

from argparse                    import HelpFormatter
from pytadbit.mapping.mapper     import get_intersection
from os                          import path, system
from pytadbit.utils.sqlite_utils import get_jobid, add_path, get_path_id, print_db
from hashlib                     import md5
import logging
import fcntl
import sqlite3 as lite
import time

DESC = "Parse mapped Hi-C reads and get the intersection"

def run(opts):
    check_options(opts)
    launch_time = time.localtime()

    fname1, fname2 = load_parameters_fromdb(opts.workdir)

    jobid = get_jobid(workdir=opts.workdir)

    system('mkdir -p ' + path.join(opts.workdir, '%03d_filtered_reads' % jobid))

    mreads = path.join(opts.workdir, '%03d_filtered_reads' % jobid,
                       'all_r1-r2_intersection.tsv')

    count, multiples = get_intersection(fname1, fname2, mreads)

    finish_time = time.localtime()

    # save all job information to sqlite DB
    save_to_db(opts, count, multiples, mreads,
               launch_time, finish_time)

def save_to_db(opts, count, multiples, mreads,
               launch_time, finish_time):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='INTERSECTIONs'""")
        if not cur.fetchall():
            cur.execute("""
        create table INTERSECTIONs
           (Id integer primary key,
            PATHid int,
            Total_interactions int,
            Multiple_interactions text,
            unique (PATHid))""")
        try:
            parameters = ' '.join(
                ['%s:%s' % (k, v) for k, v in opts.__dict__.iteritems()
                 if not k in ['fastq', 'index', 'renz', 'iterative', 'workdir',
                              'func', 'tmp'] and not v is None])
            parameters = parameters.replace("'", '"')
            param_hash = md5(' '.join(
                ['%s:%s' % (k, v) for k, v in sorted(opts.__dict__.iteritems())
                 if not k in ['workdir', 'func', 'tmp']])).hexdigest()
            cur.execute("""
    insert into JOBs
     (Id  , Parameters, Launch_time, Finish_time,    Type, Parameters_md5)
    values
     (NULL,       '%s',        '%s',        '%s', 'Parse',           '%s')
     """ % (parameters,
            time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
            time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass

        jobid = get_jobid(cur)
        
        add_path(cur, mreads, '2D_BED', jobid, opts.workdir)
        try:
            cur.execute("""
            insert into INTERSECTIONs
            (Id  , PATHid, Total_interactions, Multiple_interactions)
            values
            (NULL,    %d,     %d,      '%s')
            """ % (get_path_id(cur, mreads, opts.workdir),
                   count, ' '.join(['%s:%d' % (k, multiples[k])
                                    for k in sorted(multiples)])))
        except lite.IntegrityError:
            print 'WARNING: already parsed'
        print_db(cur, 'FASTQs')
        print_db(cur, 'PATHs')
        print_db(cur, 'SAMs')
        print_db(cur, 'BEDs')
        print_db(cur, 'JOBs')
        print_db(cur, 'INTERSECTIONs')
        

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
