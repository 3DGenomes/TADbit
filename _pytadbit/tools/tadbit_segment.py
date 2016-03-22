"""

information needed

 - path working directory with parsed reads

"""
from argparse                     import HelpFormatter
from pytadbit                     import load_hic_data_from_reads
from pytadbit                     import tadbit
from pytadbit.utils.sqlite_utils  import already_run, digest_parameters
from pytadbit.utils.sqlite_utils  import add_path, get_jobid, print_db
from pytadbit.utils.file_handling import mkdir
from pytadbit.parsers.tad_parser  import parse_tads
from os                           import path, remove
from time                         import sleep
from shutil                       import copyfile
from string                       import ascii_letters
from random                       import random
import sqlite3 as lite
import time

DESC = 'Finds TAD or compartment segmentation in Hi-C data.'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()
    param_hash = digest_parameters(opts)

    if not opts.nosql:
        (bad_co, bad_co_id, biases, biases_id,
         mreads, mreads_id, reso) = load_parameters_fromdb(opts)
        # store path ids to be saved in database
        inputs = bad_co_id, biases_id, mreads_id
    else:
        bad_co = opts.bad_co
        biases = opts.biases
        mreads = opts.mreads
        reso   = opts.reso
        

    mreads = path.join(opts.workdir, mreads)
    bad_co = path.join(opts.workdir, bad_co)
    biases = path.join(opts.workdir, biases)

    mkdir(path.join(opts.workdir, '05_segmentation'))

    print 'loading %s at resolution %s' % (mreads, nice(reso))
    hic_data = load_hic_data_from_reads(mreads, reso)
    hic_data.bads = dict((int(l.strip()), True) for l in open(bad_co))
    hic_data.bias = dict((int(l.split()[0]), float(l.split()[1]))
                         for l in open(biases))

    # compartments
    cmp_result = {}
    if not opts.only_tads:
        print 'Searching compartments'
        hic_data.find_compartments(crms=opts.crms)

        cmprt_dir = path.join(opts.workdir, '05_segmentation',
                              'compartments_%s_%s' % (
            nice(reso), param_hash))
        mkdir(cmprt_dir)
        for crm in hic_data.chromosomes:
            cmprt_file = path.join(cmprt_dir, '%s.tsv' % crm)
            hic_data.write_compartments(cmprt_file,
                                        chrom=crm)
            cmp_result[crm] = {'path': cmprt_file,
                               'num' : len(hic_data.compartments[crm])}

    # TADs
    tad_result = {}
    if not opts.only_compartments:
        print 'Searching TADs'
        tad_dir = path.join(opts.workdir, '05_segmentation',
                             'tads_%s_%s' % (
                                 nice(reso), param_hash))
        mkdir(tad_dir)
        for crm in hic_data.chromosomes:
            if opts.crms and not crm in opts.crms:
                continue
            print '  - %s' % crm
            matrix = hic_data.get_matrix(focus=crm)
            beg, end = hic_data.section_pos[crm]
            size = len(matrix)
            if size < 10:
                print "     Chromosome too short (%d bins), skipping..." % size
                continue
            # transform bad column in chromosome referential
            remove = tuple([1 if i in hic_data.bads else 0 for i in xrange(beg, end)])
            # maximum size of a TAD
            max_tad_size = size if opts.max_tad_size is None else opts.max_tad_size
            result = tadbit([matrix], remove=remove,
                            n_cpus=opts.cpus, verbose=True,
                            max_tad_size=max_tad_size,
                            no_heuristic=True)
            tads = load_tad_height(result, size, beg, end, hic_data)
            table = ''
            table += '%s\t%s\t%s\t%s%s\n' % ('#', 'start', 'end', 'score', 'density')
            for tad in tads:
                table += '%s\t%s\t%s\t%s%s\n' % (
                    tad, int(tads[tad]['start'] + 1), int(tads[tad]['end'] + 1),
                    abs(tads[tad]['score']), '\t%s' % (round(
                        float(tads[tad]['height']), 3)))
            out_tad = path.join(tad_dir, '%s.tsv' % (crm))
            out = open(out_tad, 'w')
            out.write(table)
            out.close()
            tad_result[crm] = {'path' : out_tad,
                               'num': len(tads)}

    finish_time = time.localtime()

    if not opts.nosql:
        save_to_db(opts, cmp_result, tad_result, reso, inputs, 
                   launch_time, finish_time)

def save_to_db(opts, cmp_result, tad_result, reso, inputs,
               launch_time, finish_time):
    if 'tmp' in opts and opts.tmp:
        # check lock
        while path.exists(path.join(opts.workdir, '__lock_db')):
            sleep(0.5)
        # close lock
        open(path.join(opts.workdir, '__lock_db'), 'wa').close()
        # tmp file
        dbfile = opts.tmp
        copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='SEGMENT_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
            create table SEGMENT_OUTPUTs
               (Id integer primary key,
                JOBid int,
                Inputs text,
                TADs int,
                Compartments int,
                Chromosome text,
                Resolution int)""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            param_hash = digest_parameters(opts, get_md5=True )
            cur.execute("""
            insert into JOBs
            (Id  , Parameters, Launch_time, Finish_time, Type , Parameters_md5)
            values
            (NULL,       '%s',        '%s',        '%s', 'Segment',           '%s')
            """ % (parameters,
                   time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
                   time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass
        jobid = get_jobid(cur)
        for crm in max(cmp_result.keys(), tad_result.keys(),
                       key=lambda x: len(x)):
            if crm in cmp_result:
                add_path(cur, cmp_result[crm]['path'], 'COMPARTMENT',
                         jobid, opts.workdir)
            if crm in tad_result:
                add_path(cur, tad_result[crm]['path'], 'TAD', jobid, opts.workdir)
            cur.execute("""
            insert into SEGMENT_OUTPUTs
            (Id  , JOBid, Inputs, TADs, Compartments, Chromosome, Resolution)
            values
            (NULL,    %d,   '%s',   %d,           %d,       '%s',         %d)
            """ % (jobid,
                   ','.join([str(i) for i in inputs]),
                   cmp_result[crm]['num'] if crm in tad_result else 0,
                   tad_result[crm]['num'] if crm in cmp_result else 0,
                   crm,
                   reso))
            print_db(cur, 'PATHs')
            print_db(cur, 'JOBs')
            print_db(cur, 'SEGMENT_OUTPUTs')
    if 'tmp' in opts and opts.tmp:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
        # release lock
        remove(path.join(opts.workdir, '__lock_db'))

def load_parameters_fromdb(opts):
    if 'tmp' in opts and opts.tmp:
        dbfile = opts.tmp
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        if not opts.jobid:
            # get the JOBid of the parsing job
            cur.execute("""
            select distinct Id from JOBs
            where Type = 'Normalize'
            """)
            jobids = cur.fetchall()
            if len(jobids) > 1:
                raise Exception('ERROR: more than one possible input found, use'
                                '"tadbit describe" and select corresponding '
                                'jobid with --jobid')
            parse_jobid = jobids[0][0]
        else:
            parse_jobid = opts.jobid
        # fetch path to parsed BED files
        cur.execute("""
        select distinct Path, PATHs.Id from PATHs
        where paths.jobid = %s and paths.Type = 'BAD_COLUMNS'
        """ % parse_jobid)
        bad_co, bad_co_id  = cur.fetchall()[0]
        cur.execute("""
        select distinct Path, PATHs.Id from PATHs
        where paths.jobid = %s and paths.Type = 'BIASES'
        """ % parse_jobid)
        biases, biases_id = cur.fetchall()[0]
        cur.execute("""
        select distinct Path, PATHs.Id from PATHs
        inner join NORMALIZE_OUTPUTs on PATHs.Id = NORMALIZE_OUTPUTs.Input
        where NORMALIZE_OUTPUTs.JOBid = %d;
        """ % parse_jobid)
        mreads, mreads_id = cur.fetchall()[0]
        cur.execute("""
        select distinct Resolution from NORMALIZE_OUTPUTs
        where NORMALIZE_OUTPUTs.JOBid = %d;
        """ % parse_jobid)
        reso = int(cur.fetchall()[0][0])
        return (bad_co, bad_co_id, biases, biases_id,
                mreads, mreads_id, reso)

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    glopts.add_argument('--tmp', dest='tmp', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument('--nosql', dest='nosql',
                        action='store_true', default=False, 
                        help='do not load/store data from/in sqlite database')

    glopts.add_argument('--mreads', dest='mreads', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path valid-pairs file''')

    glopts.add_argument('--bad_cols', dest='bad_co', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to file with bad columns''')

    glopts.add_argument('--biases',   dest='biases', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to file with precalculated biases by
                        columns''')

    glopts.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int,
                        help='''resolution at which to output matrices''')

    glopts.add_argument('--norm_matrix', dest='norm_matrix', metavar="PATH",
                        action='store', default=None, type=str, 
                        help='''path to a matrix file with normalized read
                        counts''')
    
    glopts.add_argument('--raw_matrix', dest='raw_matrix', metavar="PATH",
                        action='store', default=None, type=str, 
                        help='''path to a matrix file with raw read
                        counts''')

    glopts.add_argument('--only_compartments', dest='only_compartments',
                        action='store_true', default=False, 
                        help='''search A/B compartments using first eigen vector
                        of the correlation matrix''')

    glopts.add_argument('--only_tads', dest='only_tads',
                        action='store_true', default=False, 
                        help='''search TAD boundaries break-point detection
                        algorithm''')

    glopts.add_argument('--perc_zeros', dest='perc_zeros', metavar="FLOAT",
                        action='store', default=95, type=float, 
                        help='maximum percentage of zeroes allowed per column')

    glopts.add_argument('-j', '--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')    

    glopts.add_argument('-c', '--chromosomes', dest='crms', metavar="STR",
                        action='store', default=None, type=str, nargs='+',
                        help='''Name of the chromosomes on which to search
                        for TADs or compartments.''')    

    glopts.add_argument('--max_tad_size', dest='max_tad_size', metavar="INT",
                        action='store', default=None, type=int,
                        help='''an integer defining the maximum size of TAD. Default
                        defines it as the number of rows/columns''')

    glopts.add_argument("-C", "--cpu", dest="cpus", type=int,
                        default=0, help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    glopts.add_argument('--force', dest='force', action='store_true',
                      default=False,
                      help='overwrite previously run job')

    parser.add_argument_group(glopts)

def check_options(opts):

    # check resume
    if not path.exists(opts.workdir) and opts.resume:
        print ('WARNING: can use output files, found, not resuming...')
        opts.resume = False

    if 'tmp' in opts and opts.tmp:
        dbdir = opts.tmp
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmp = path.join(dbdir, dbfile)
        copyfile(path.join(opts.workdir, 'trace.db'), opts.tmp)

    if already_run(opts) and not opts.force:
        remove(path.join(dbdir, dbfile))
        exit('WARNING: exact same job already computed, see JOBs table above')

def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)

def load_tad_height(tad_def, size, beg, end, hic_data):
    bias, zeros = hic_data.bias, hic_data.bads
    tads, _ = parse_tads(tad_def)
    diags = []
    for k in xrange(1, size):
        try:
            diags.append(sum(
                hic_data[i, i + k] / bias[i + k] / bias[i]
                for i in xrange(beg, end - k) if not i in zeros and not i + k in zeros
                ) / float(sum(1 for i in range(beg, end - k)
                              if not i in zeros and not i + k in zeros)))
        except ZeroDivisionError:
            diags.append(0.)
    for tad in tads:
        start, final = (int(tads[tad]['start']) + 1,
                        int(tads[tad]['end']) + 1)
        matrix = sum(
            hic_data[i, j] / bias[j] / bias[i]
            for i in xrange(beg + start - 1, beg + final - 1) if not i in zeros
            for j in xrange(i + 1          , beg + final - 1) if not j in zeros)
        try:
            height = float(matrix) / sum(
                [diags[i - 1] * (final - start - i)
                 for i in xrange(1, final - start)])
        except ZeroDivisionError:
            height = 0.
        tads[tad]['height'] = height
    return tads
