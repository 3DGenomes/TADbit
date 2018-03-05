"""

information needed

 - path working directory with parsed reads

"""

from argparse                        import HelpFormatter
from string                          import ascii_letters
from random                          import random
from os                              import path, remove
from shutil                          import copyfile
from multiprocessing                 import cpu_count
from warnings                        import warn

import sqlite3 as lite
import time

from pytadbit.mapping                import get_intersection
from pytadbit.utils.file_handling    import mkdir
from pytadbit.utils.sqlite_utils     import get_jobid, add_path, get_path_id, print_db
from pytadbit.utils.sqlite_utils     import already_run, digest_parameters
from pytadbit.mapping.analyze        import insert_sizes
from pytadbit.mapping.filter         import filter_reads, apply_filter
from pytadbit.parsers.hic_bam_parser import bed2D_to_BAMhic


DESC = "Filter parsed Hi-C reads and get valid pair of reads to work with"


def run(opts):
    check_options(opts)
    launch_time = time.localtime()

    fname1, fname2 = load_parameters_fromdb(opts)

    param_hash = digest_parameters(opts)

    reads = path.join(opts.workdir, '03_filtered_reads',
                      'all_r1-r2_intersection_%s.tsv' % param_hash)
    mreads = path.join(opts.workdir, '03_filtered_reads',
                       'valid_r1-r2_intersection_%s.tsv' % param_hash)

    if not opts.resume:
        mkdir(path.join(opts.workdir, '03_filtered_reads'))

        # compute the intersection of the two read ends
        print 'Getting intersection between read 1 and read 2'
        count, multiples = get_intersection(fname1, fname2, reads)

        # compute insert size
        print 'Get insert size...'
        hist_path = path.join(opts.workdir,
                              'histogram_fragment_sizes_%s.pdf' % param_hash)
        median, max_f, mad = insert_sizes(
            reads, nreads=1000000, stats=('median', 'first_decay', 'MAD'),
            savefig=hist_path)

        print '  - median insert size =', median
        print '  - double median absolution of insert size =', mad
        print '  - max insert size (when a gap in continuity of > 10 bp is found in fragment lengths) =', max_f

        max_mole = max_f # pseudo DEs
        min_dist = max_f + mad # random breaks
        print ('   Using the maximum continuous fragment size'
               '(%d bp) to check '
               'for pseudo-dangling ends') % max_mole
        print ('   Using maximum continuous fragment size plus the MAD '
               '(%d bp) to check for random breaks') % min_dist

        print "identify pairs to filter..."
        masked = filter_reads(reads, max_molecule_length=max_mole,
                              over_represented=opts.over_represented,
                              max_frag_size=opts.max_frag_size,
                              min_frag_size=opts.min_frag_size,
                              re_proximity=opts.re_proximity,
                              min_dist_to_re=min_dist, fast=True)

    n_valid_pairs = apply_filter(reads, mreads, masked, filters=opts.apply)

    outbam = path.join(opts.workdir, '03_filtered_reads',
                       'intersection_%s' % param_hash)

    if opts.valid:
        infile = mreads
    else:
        infile = reads
    bed2D_to_BAMhic(infile, opts.valid, opts.cpus, outbam, opts.format, masked,
                    samtools=opts.samtools)

    finish_time = time.localtime()
    print median, max_f, mad
    # save all job information to sqlite DB
    save_to_db(opts, count, multiples, reads, mreads, n_valid_pairs, masked,
               outbam + '.bam', hist_path, median, max_f, mad, launch_time, finish_time)


def save_to_db(opts, count, multiples, reads, mreads, n_valid_pairs, masked,
               outbam, hist_path, median, max_f, mad, launch_time, finish_time):
    if 'tmpdb' in opts and opts.tmpdb:
        # check lock
        while path.exists(path.join(opts.workdir, '__lock_db')):
            time.sleep(0.5)
        # close lock
        open(path.join(opts.workdir, '__lock_db'), 'a').close()
        # tmp file
        dbfile = opts.tmpdb
        try: # to copy in case read1 was already mapped for example
            copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
        except IOError:
            pass
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='INTERSECTION_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
        create table INTERSECTION_OUTPUTs
           (Id integer primary key,
            PATHid int,
            Total_interactions int,
            Multiple_interactions text,
            Median_fragment_length,
            MAD_fragment_length,
            Max_fragment_length,
            unique (PATHid))""")
            cur.execute("""
        create table FILTER_OUTPUTs
           (Id integer primary key,
            PATHid int,
            Name text,
            Count int,
            Applied text,
            JOBid int,
            unique (PATHid))""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            param_hash = digest_parameters(opts, get_md5=True )
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
        add_path(cur, outbam, 'HIC_BAM', jobid, opts.workdir)
        add_path(cur, outbam + '.bai', 'HIC_BAI', jobid, opts.workdir)
        add_path(cur,  reads, '2D_BED', jobid, opts.workdir)
        add_path(cur, hist_path, 'FIGURE', jobid, opts.workdir)
        try:
            cur.execute("""
            insert into INTERSECTION_OUTPUTs
            (Id  , PATHid, Total_interactions, Multiple_interactions, Median_fragment_length, MAD_fragment_length, Max_fragment_length)
            values
            (NULL,    %d,                  %d,                  '%s',                     %d,                  %d,                  %d)
            """ % (get_path_id(cur, mreads, opts.workdir),
                   count, ' '.join(['%s:%d' % (k, multiples[k])
                                    for k in sorted(multiples)]),
                   median, mad, max_f))
        except lite.IntegrityError:
            print 'WARNING: already filtered'
            if opts.force:
                cur.execute(
                    'delete from INTERSECTION_OUTPUTs where PATHid = %d' % (
                        get_path_id(cur, mreads, opts.workdir)))
                cur.execute("""
                insert into INTERSECTION_OUTPUTs
                (Id  , PATHid, Total_interactions, Multiple_interactions, Median_fragment_length, MAD_fragment_length, Max_fragment_length)
                values
                (NULL,    %d,                  %d,                  '%s',                     %d,                  %d,                  %d)
                """ % (get_path_id(cur, mreads, opts.workdir),
                       count, ' '.join(['%s:%d' % (k, multiples[k])
                                        for k in sorted(multiples)]),
                       median, mad, max_f))
        for nf, f in enumerate(masked, 1):
            try:
                add_path(cur, masked[f]['fnam'], 'FILTER', jobid, opts.workdir)
            except KeyError:
                continue
            try:
                cur.execute("""
            insert into FILTER_OUTPUTs
                (Id  , PATHid, Name, Count, Applied, JOBid)
            values
                (NULL,     %d, '%s',  '%s',    '%s',    %d)
                """ % (get_path_id(cur, masked[f]['fnam'], opts.workdir),
                       masked[f]['name'], masked[f]['reads'],
                       'True' if nf in opts.apply else 'False', jobid))
            except lite.IntegrityError:
                print 'WARNING: already filtered'
                if opts.force:
                    cur.execute(
                        'delete from FILTER_OUTPUTs where PATHid = %d' % (
                            get_path_id(cur, masked[f]['fnam'], opts.workdir)))
                    cur.execute("""
                insert into FILTER_OUTPUTs
                    (Id  , PATHid, Name, Count, Applied, JOBid)
                values
                    (NULL,     %d, '%s',  '%s',    '%s',    %d)
                    """ % (get_path_id(cur, masked[f]['fnam'], opts.workdir),
                           masked[f]['name'], masked[f]['reads'],
                           'True' if nf in opts.apply else 'False', jobid))
        try:
            cur.execute("""
        insert into FILTER_OUTPUTs
            (Id  , PATHid, Name, Count, Applied, JOBid)
        values
            (NULL,     %d, '%s',  '%s',    '%s',    %d)
            """ % (get_path_id(cur, mreads, opts.workdir),
                   'valid-pairs', n_valid_pairs, '', jobid))
        except lite.IntegrityError:
            print 'WARNING: already filtered'
            if opts.force:
                cur.execute(
                    'delete from FILTER_OUTPUTs where PATHid = %d' % (
                        get_path_id(cur, mreads, opts.workdir)))
                cur.execute("""
                insert into FILTER_OUTPUTs
                (Id  , PATHid, Name, Count, Applied, JOBid)
                values
                (NULL,     %d, '%s',  '%s',    '%s',    %d)
                """ % (get_path_id(cur, mreads, opts.workdir),
                       'valid-pairs', n_valid_pairs, '', jobid))
        print_db(cur, 'MAPPED_INPUTs')
        print_db(cur, 'PATHs')
        print_db(cur, 'MAPPED_OUTPUTs')
        print_db(cur, 'PARSED_OUTPUTs')
        print_db(cur, 'JOBs')
        print_db(cur, 'INTERSECTION_OUTPUTs')
        print_db(cur, 'FILTER_OUTPUTs')
    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
    # release lock
    try:
        remove(path.join(opts.workdir, '__lock_db'))
    except OSError:
        pass


def load_parameters_fromdb(opts):
    if 'tmpdb' in opts and opts.tmpdb:
        dbfile = opts.tmpdb
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        # get the JOBid of the parsing job
        if not opts.pathids:
            cur.execute("""
            select distinct Id from PATHs
            where Type = 'BED'
            """)
            pathids = cur.fetchall()
            if len(pathids) > 2:
                raise Exception('ERROR: more than one possible input found'
                                '(PATHids: %s), use "tadbit describe" and '
                                'select corresponding PATHid with --pathids' % (
                                    ', '.join([str(j[0]) for j in pathids])))
            pathids = [p[0] for p in pathids]
        else:
            pathids = opts.pathids
        # fetch path to parsed BED files
        cur.execute("""
        select distinct Path from PATHs
        where Id = %d or Id = %d
        """ % (pathids[0], pathids[1]))
        fname1, fname2 = [path.join(opts.workdir, e[0]) for e in cur.fetchall()]
        if not path.exists(fname1):
            if path.exists(fname1 + '.gz') and path.exists(fname2 + '.gz'):
                fname1 += '.gz'
                fname2 += '.gz'
            else:
                raise IOError('ERROR: unput file_handling does not exist')

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

    glopts  = parser.add_argument_group('General options')
    output  = parser.add_argument_group('Storage options')
    filter_ = parser.add_argument_group('Filtering options')

    glopts.add_argument('--force', dest='force', action='store_true',
                        default=False,
                        help='overwrite previously run job')

    glopts.add_argument('--resume', dest='resume', action='store_true',
                        default=False,
                        help='use filters of previously run job')

    filter_.add_argument('--apply', dest='apply', nargs='+',
                         type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                         choices = range(1, 11),
                         help=("""[%(default)s] Use filters to define a set os valid pair of reads
                         e.g.: '--apply 1 2 3 4 6 7 8 9'. Where these numbers""" +
                               "correspond to: %s" % (', '.join(
                                   ['%2d: %15s' % (k, masked[k]['name']) for k in masked]))))

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    glopts.add_argument("-C", "--cpus", dest="cpus", type=int,
                        default=cpu_count(), help='''[%(default)s] Maximum
                        number of CPU cores  available in the execution host.
                         If higher than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    filter_.add_argument('--over_represented', dest='over_represented', metavar="NUM",
                         action='store', default=0.001, type=float,
                         help='''[%(default)s%%] percentage of restriction-enzyme
                         (RE) genomic fragments with more coverage to exclude
                         (possible PCR artifact).''')

    filter_.add_argument('--max_frag_size', dest='max_frag_size', metavar="NUM",
                         action='store', default=100000, type=int,
                         help='''[%(default)s] to exclude large genomic RE
                         fragments (probably resulting from gaps in the reference
                         genome)''')

    filter_.add_argument('--min_frag_size', dest='min_frag_size', metavar="NUM",
                         action='store', default=50, type=int,
                         help='''[%(default)s] to exclude small genomic RE
                         fragments (smaller than sequenced reads)''')

    filter_.add_argument('--re_proximity', dest='re_proximity', metavar="NUM",
                         action='store', default=5, type=int,
                         help ='''[%(default)s] to exclude read-ends falling too
                         close from RE site (pseudo-dangling-ends)''')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument('--pathids', dest='pathids', metavar="INT",
                        action='store', default=None, nargs='+', type=int,
                        help='''Use as input data generated by a job under a given
                        pathids. Use tadbit describe to find out which.
                        Needs one PATHid per read, first for read 1,
                        second for read 2.''')
    output.add_argument('--format', dest='format', default='mid',
                        choices=['short', 'mid', 'long'],
                        help='''[%(default)s] for compression into pseudo-BAM
                        format. Short contains only positions of reads mapped,
                        mid everything but restriction sites.''')
    output.add_argument('--valid', dest='valid', default=False,
                        action='store_true',
                        help='''stores only valid-pairs discards filtered out
                        reads.''')
    output.add_argument('--clean', dest='clean', default=False,
                        action='store_true',
                        help='''remove intermediate files. WARNING: together
                        with format "short" or valid options, this may results
                        in losing data''')
    glopts.add_argument('--samtools', dest='samtools', metavar="PATH",
                        action='store', default='samtools', type=str,
                        help='''path samtools binary''')

    parser.add_argument_group(glopts)


def check_options(opts):

    if not opts.workdir:
        raise Exception('ERROR: output option required.')

    # check resume
    if not path.exists(opts.workdir) and opts.resume:
        print 'WARNING: can use output files, found, not resuming...'
        opts.resume = False

    # sort filters
    if opts.apply:
        opts.apply.sort()

    # for lustre file system....
    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)
        try:
            copyfile(path.join(opts.workdir, 'trace.db'), opts.tmpdb)
        except IOError:
            pass

    # check if job already run using md5 digestion of parameters
    if already_run(opts):
        if not opts.force:
            if 'tmpdb' in opts and opts.tmpdb:
                remove(path.join(dbdir, dbfile))
            exit('WARNING: exact same job already computed, see JOBs table above')
        else:
            warn('WARNING: exact same job already computed, overwritting...')
