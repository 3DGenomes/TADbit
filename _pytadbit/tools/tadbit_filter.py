"""

information needed

 - path working directory with parsed reads

"""

from argparse                    import HelpFormatter
from pytadbit.mapping.mapper     import get_intersection
from os                          import path, system
from pytadbit.utils.sqlite_utils import get_jobid, add_path, get_path_id, print_db
from pytadbit.mapping.analyze    import plot_iterative_mapping, insert_sizes
from pytadbit.mapping.filter     import filter_reads, apply_filter
from hashlib                     import md5
import logging
import fcntl
import sqlite3 as lite
import time

DESC = "Parse mapped Hi-C reads and get the intersection"

def run(opts):
    check_options(opts)
    launch_time = time.localtime()

    fname1, fname2 = load_parameters_fromdb(opts)

    jobid = get_jobid(workdir=opts.workdir) + 1

    reads = path.join(opts.workdir, '%03d_filtered_reads' % jobid,
                      'all_r1-r2_intersection.tsv')
    mreads = path.join(opts.workdir, '%03d_filtered_reads' % jobid,
                       'valid_r1-r2_intersection.tsv')

    if not opts.resume:
        system('mkdir -p ' + path.join(opts.workdir, '%03d_filtered_reads' % jobid))

        # compute the intersection of the two read ends
        count, multiples = get_intersection(fname1, fname2, reads)

        # compute insert size
        print 'Get insert size...'
        median, ori_max_mol = insert_sizes(reads, nreads=1000000)
        
        print '  - median insert size =', median
        print '  - max (99.9%) insert size =', ori_max_mol
    
        max_mol = 4 * median
        print ('   Using 4 times median insert size (%d bp) to check '
               'for random breaks') % max_mol
    
        print "identify pairs to filter..."
        masked = filter_reads(reads, max_molecule_length=max_mol,
                              over_represented=0.001, max_frag_size=100000,
                              min_frag_size=50, re_proximity=5,
                              min_dist_to_re=max_mol, fast=True)

    n_valid_pairs = apply_filter(reads, mreads, masked,
                                 filters=opts.apply)

    finish_time = time.localtime()

    # save all job information to sqlite DB
    save_to_db(opts, count, multiples, mreads, n_valid_pairs, masked,
               launch_time, finish_time)

def save_to_db(opts, count, multiples, mreads, n_valid_pairs, masked,
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
            cur.execute("""
        create table FILTERs
           (Id integer primary key,
            PATHid int,
            Name text,
            Count int,
            JOBid int,
            unique (PATHid))""")
        try:
            parameters = ' '.join(
                ['%s:%s' % (k, v) for k, v in opts.__dict__.iteritems()
                 if not k in ['fastq', 'index', 'renz', 'iterative', 'workdir',
                              'func', 'tmp'] and not v is None])
            parameters = parameters.replace("'", '"')
            param_hash = md5(' '.join(
                ['%s:%s' % (k, v) for k, v in sorted(opts.__dict__.iteritems())
                 if not k in ['force', 'workdir', 'func', 'tmp']])).hexdigest()
            cur.execute("""
    insert into JOBs
     (Id  , Parameters, Launch_time, Finish_time,    Type, Parameters_md5)
    values
     (NULL,       '%s',        '%s',        '%s', 'Filter',           '%s')
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
        for f in masked:
            add_path(cur, masked[f]['fnam'], 'FILTER', jobid, opts.workdir)
            try:
                cur.execute("""
            insert into FILTERs
            (Id  , PATHid, Name, Count, JOBid)
            values
            (NULL,    %d,     '%s',      '%s', %d)
                """ % (get_path_id(cur, masked[f]['fnam'], opts.workdir),
                       masked[f]['name'], masked[f]['reads'], jobid))
            except lite.IntegrityError:
                print 'WARNING: already filtered'
        try:
            cur.execute("""
        insert into FILTERs
        (Id  , PATHid, Name, Count, JOBid)
        values
        (NULL,    %d,     '%s',      '%s', %d)
            """ % (get_path_id(cur, mreads, opts.workdir),
                   'Valid-pairs', n_valid_pairs, jobid))
        except lite.IntegrityError:
            print 'WARNING: already filtered'
        print_db(cur, 'FASTQs')
        print_db(cur, 'PATHs')
        print_db(cur, 'SAMs')
        print_db(cur, 'BEDs')
        print_db(cur, 'JOBs')
        print_db(cur, 'INTERSECTIONs')        
        print_db(cur, 'FILTERs')

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

    return fname1, fname2

def populate_args(parser):
    """
    parse option from call
    """
    masked = {1 : {'name': 'self-circle'       }, 
              2 : {'name': 'dangling-end'      },
              3 : {'name': 'error'             },
              4 : {'name': 'extra dangling-end'},
              5 : {'name': 'too close from RES'},
              6 : {'name': 'too short'         },
              7 : {'name': 'too large'         },
              8 : {'name': 'over-represented'  },
              9 : {'name': 'duplicated'        },
              10: {'name': 'random breaks'     }}
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    glopts.add_argument('--force', dest='force', action='store_true',
                      default=False,
                      help='overwrite previously run job')

    glopts.add_argument('--resume', dest='resume', action='store_true',
                      default=False,
                      help='use filters of previously run job')

    glopts.add_argument('--apply', dest='apply', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 8, 9, 10],
                        help=("""[%(default)s] Use filters to define a set os valid pair of reads
                        e.g.: '--apply 1 2 3 4 6 7 8 9'. Where these numbers""" + 
                        "correspond to: %s" % (', '.join(
                            ['%2d: %15s' % (k, masked[k]['name']) for k in masked]))))

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
