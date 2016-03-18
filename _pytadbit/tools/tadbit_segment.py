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
from os                           import path
import sqlite3 as lite
import time

DESC = 'Finds TAD or compartment segmentation in Hi-C data.'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()
    param_hash = digest_parameters(opts)

    bad_co, biases, mreads, reso = load_parameters_fromdb(opts)

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
    if not opts.only_tads:
        print 'Searching compartments'
        hic_data.find_compartments()

        cmprt_file = path.join(opts.workdir, '05_segmentation',
                               'compartments_%s_%s.tsv' % (
            nice(reso), param_hash))
        hic_data.write_compartments(cmprt_file)

    if not opts.only_compartments:
        print 'Searching TADs'
        out_tads = []
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
            out_tad = path.join(opts.workdir, '05_segmentation',
                                'tads_%s_%s_%s.tsv' % (
                                    crm, nice(reso), param_hash))
            out = open(out_tad, 'w')
            out.write(table)
            out.close()
            out_tads.append(out_tad)

    finish_time = time.localtime()

    save_to_db(opts, cmprt_file, out_tads, reso, launch_time, finish_time)


def save_to_db(opts, cmprt_file, out_tads, reso, launch_time, finish_time):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='SEGMENT_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
            create table SEGMENT_OUTPUTs
               (Id integer primary key,
                JOBid int,
                Resolution int,
                unique (JOBid))""")
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
        add_path(cur, cmprt_file, 'COMPARTMENT', jobid, opts.workdir)
        for tad_file in out_tads:
            add_path(cur, tad_file, 'TAD', jobid, opts.workdir)

        cur.execute("""
        insert into SEGMENT_OUTPUTs
        (Id  , JOBid, Resolution)
        values
        (NULL,    %d,        %d)
        """ % (jobid, reso))
        print_db(cur, 'PATHs')
        print_db(cur, 'JOBs')
        print_db(cur, 'SEGMENT_OUTPUTs')

def load_parameters_fromdb(opts):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
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
        select distinct path from paths
        where paths.jobid = %s and paths.Type = 'BAD_COLUMNS'
        """ % parse_jobid)
        bad_columns = cur.fetchall()[0][0]
        cur.execute("""
        select distinct path from paths
        where paths.jobid = %s and paths.Type = 'BIASES'
        """ % parse_jobid)
        biases = cur.fetchall()[0][0]
        cur.execute("""
        select distinct Path from PATHs
        inner join NORMALIZE_OUTPUTs on PATHs.Id = NORMALIZE_OUTPUTs.Input
        where NORMALIZE_OUTPUTs.JOBid = %d;
        """ % parse_jobid)
        mreads = cur.fetchall()[0][0]
        cur.execute("""
        select distinct Resolution from NORMALIZE_OUTPUTs
        where NORMALIZE_OUTPUTs.JOBid = %d;
        """ % parse_jobid)
        reso = int(cur.fetchall()[0][0])
        return bad_columns, biases, mreads, reso

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

    if already_run(opts) and not opts.force:
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
