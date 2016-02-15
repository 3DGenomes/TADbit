"""

information needed

 - path to FASTQ
 - path to reference genome
 - path to indexed reference genome
 - read number (1/2)
 - restriction enzyme used
 - species name
 - chromosome names (optional)
 - descriptive fields (optional, e.g. --descr=flowcell:C68AEACXX,lane:4,index:24nf)

mapping strategy

 - iterative/fragment
 - mapper

"""

from argparse                             import HelpFormatter
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES
from pytadbit.utils.fastq_utils           import quality_plot
from pytadbit.mapping.full_mapper         import full_mapping
from pytadbit.utils.sqlite_utils          import get_path_id, add_path, print_db
from pytadbit.utils.sqlite_utils          import get_jobid
from pytadbit                             import get_dependencies_version
from os                                   import system, path
from hashlib                              import md5
from multiprocessing                      import cpu_count
import logging
import fcntl
import sqlite3 as lite
import time

DESC = "Map Hi-C reads and organize results in an output working directory"

def run(opts):
    check_options(opts)

    launch_time = time.localtime()

    if opts.quality_plot:
        logging.info('Generating Hi-C QC plot at:\n  ' +
               path.join(opts.workdir, path.split(opts.fastq)[-1] + '.pdf'))
        quality_plot(opts.fastq, r_enz=opts.renz,
                     nreads=100000, paired=False,
                     savefig=path.join(opts.workdir,
                                       path.split(opts.fastq)[-1] + '.pdf'))
        return

    jobid = get_jobid(workdir=opts.workdir) + 1

    logging.info('mapping %s read %s to %s', opts.fastq, opts.read, opts.workdir)
    outfiles = full_mapping(opts.index, opts.fastq,
                            path.join(opts.workdir,
                                      '%03d_mapped_r%d' % (jobid, opts.read)),
                            opts.renz, temp_dir=opts.tmp, nthreads=opts.cpus,
                            frag_map=not opts.iterative, clean=True,
                            windows=opts.windows, get_nread=True, skip=opts.skip)

    # adjust line count
    if opts.skip:
        for i, (out, _) in enumerate(outfiles[1:], 1):
            outfiles[i] = out, outfiles[i-1][1] - sum(1 for _ in open(outfiles[i-1][0]))
    
    finish_time = time.localtime()

    # save all job information to sqlite DB
    save_to_db(opts, outfiles, launch_time, finish_time)
    
    # write machine log
    with open(path.join(opts.workdir, 'trace.log'), "a") as mlog:
        fcntl.flock(mlog, fcntl.LOCK_EX)
        mlog.write('\n'.join([
            ('# MAPPED READ%s\t%d\t%s' % (opts.read, num, out))
            for out, num in outfiles]) + '\n')
        fcntl.flock(mlog, fcntl.LOCK_UN)

    logging.info('cleaning temporary files')
    # clean
    system('rm -rf ' + opts.tmp)

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')
    mapper = parser.add_argument_group('Mapping options')
    descro = parser.add_argument_group('Descriptive, optional arguments')

    glopts.add_argument('--cfg', dest='cfg', metavar="PATH", action='store',
                      default=None, type=str,
                      help='path to a configuration file with predefined ' +
                      'parameters')

    glopts.add_argument('--qc_plot', dest='quality_plot', action='store_true',
                      default=False,
                      help='generate a quality plot of FASTQ and exits')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='path to an output folder.')

    glopts.add_argument('--fastq', dest='fastq', metavar="PATH", action='store',
                      default=None, type=str,
                      help='path to a FASTQ files (can be compressed files)')

    glopts.add_argument('--index', dest='index', metavar="PATH",
                        type=str,
                        help='''paths to file(s) with indexed FASTA files of the
                        reference genome.''')

    glopts.add_argument('--read', dest='read', metavar="INT", 
                        type=int,
                        help='read number')

    glopts.add_argument('--renz', dest='renz', metavar="STR", 
                        type=str,
                        help='restriction enzyme name')

    glopts.add_argument('--chr_name', dest='chr_name', metavar="STR", nargs='+',
                        default=[], type=str,
                        help='''[fasta header] chromosome name(s). Used in the
                        same order as data.''')

    glopts.add_argument('--tmp', dest='tmp', metavar="PATH", action='store',
                      default=None, type=str,
                      help='''path to a temporary directory (default next to
                      "workdir" directory)''')

    mapper.add_argument('--iterative', dest='iterative', default=False,
                        action='store_true',
                        help='''default mapping strategy is fragment based
                        use this flag for iterative mapping''')

    mapper.add_argument('--windows', dest='windows', default=None,
                        nargs='+',
                        help='''defines windows to be used to trim the input
                        FASTQ reads, for example an iterative mapping can be defined
                        as: "--windows 1:20 1:25 1:30 1:35 1:40 1:45 1:50". But
                        this parameter can also be used for fragment based mapping
                        if for example pair-end reads are both in the same FASTQ,
                        for example: "--windows 1:50" (if the length of the reads
                        is 100). Note: that the numbers are both inclusive.''')

    # mapper.add_argument('--mapping_only', dest='mapping_only', action='store_true',
    #                     help='only do the mapping does not parse results')

    descro.add_argument('--species', dest='species', metavar="STR", 
                        type=str,
                        help='species name')

    descro.add_argument('--descr', dest='description', metavar="LIST", nargs='+',
                        type=str,
                        help='''extra descriptive fields each filed separated by
                        coma, and inside each, name and value separated by column: 
                        --descr=cell:lymphoblast,flowcell:C68AEACXX,index:24nf''')

    glopts.add_argument('--skip', dest='skip', action='store_true',
                      default=False,
                      help='[DEBUG] in case already mapped.')

    mapper.add_argument("-C", "--cpu", dest="cpus", type=int,
                        default=0, help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

def check_options(opts):
    if opts.cfg:
        get_options_from_cfg(opts.cfg, opts)

    # check RE name
    try:
        _ = RESTRICTION_ENZYMES[opts.renz]
    except KeyError:
        print ('\n\nERROR: restriction enzyme not found. Use one of:\n\n'
               + ' '.join(sorted(RESTRICTION_ENZYMES)) + '\n\n')
        raise KeyError()
    except AttributeError:
        pass

    # check skip
    if not path.exists(opts.workdir) and opts.skip:
        print ('WARNING: can use output files, found, not skipping...')
        opts.skip = False

    # number of cpus
    if opts.cpus == 0:
        opts.cpus = cpu_count()
    else:
        opts.cpus = min(opts.cpus, cpu_count())

    # check compulsory options
    if not opts.quality_plot:
        if not opts.index: raise Exception('ERROR: index  parameter required.')
    if not opts.workdir:   raise Exception('ERROR: workdir parameter required.')
    if not opts.fastq  :   raise Exception('ERROR: fastq  parameter required.')
    if not opts.renz   :   raise Exception('ERROR: renz   parameter required.')

    # create tmp directory
    if not opts.tmp:
        opts.tmp = opts.workdir + '_tmp_r%d' % opts.read

    try:
        opts.windows = [[int(i) for i in win.split(':')]
                        for win in opts.windows]
    except TypeError:
        pass
        
    system('mkdir -p ' + opts.workdir)
    # write log
    # if opts.mapping_only:
    log_format = '[MAPPING {} READ{}]   %(message)s'.format(opts.fastq, opts.read)
    # else:
    #     log_format = '[DEFAULT]   %(message)s'

    # reset logging
    logging.getLogger().handlers = []

    try:
        print 'Writting log to ' + path.join(opts.workdir, 'process.log')
        logging.basicConfig(level=logging.INFO,
                            format=log_format,
                            filename=path.join(opts.workdir, 'process.log'),
                            filemode='aw')
    except IOError:
        logging.basicConfig(level=logging.DEBUG,
                            format=log_format,
                            filename=path.join(opts.workdir, 'process.log2'),
                            filemode='aw')

    # to display log on stdout also
    logging.getLogger().addHandler(logging.StreamHandler())

    # write version log
    vlog_path = path.join(opts.workdir, 'TADbit_and_dependencies_versions.log')
    dependencies = get_dependencies_version()
    if not path.exists(vlog_path) or open(vlog_path).readlines() != dependencies:
        logging.info('Writting versions of TADbit and dependencies')
        vlog = open(vlog_path, 'w')
        vlog.write(dependencies)
        vlog.close()

def save_to_db(opts, outfiles, launch_time, finish_time):
    # write little DB to keep track of processes and options
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        # check if table exists
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='FASTQs'""")
        if not cur.fetchall():
            cur.execute("""
            create table PATHs
               (Id integer primary key,
                JOBid int, Path text, Type text,
                unique (Path))""")
            cur.execute("""
            create table JOBs
               (Id integer primary key,
                Parameters text,
                Launch_time text,
                Finish_time text,
                Type text,
                Parameters_md5 text,
                unique (Parameters_md5))""")
            cur.execute("""
            create table FASTQs
               (Id integer primary key,
                PATHid int,
                Entries int,
                Trim text,
                Frag text,
                Read int,
                Enzyme text,
                WRKDIRid int,
                SAMid int,
                INDEXid int,
                unique (PATHid,Entries,Read,Enzyme,WRKDIRid,SAMid,INDEXid))""")

        try:
            parameters = ' '.join(
                ['%s:%s' % (k, v) for k, v in opts.__dict__.iteritems()
                 if not k in ['fastq', 'index', 'renz', 'iterative', 'workdir',
                              'func', 'tmp'] and not v is None])
            param_hash = md5(' '.join(
                ['%s:%s' % (k, v) for k, v in sorted(opts.__dict__.iteritems())
                 if not k in ['workdir', 'func', 'tmp']])).hexdigest()
            cur.execute("""
    insert into JOBs
     (Id  , Parameters, Launch_time, Finish_time, Type , Parameters_md5)
    values
     (NULL,       '%s',        '%s',        '%s', 'Map',           '%s')
     """ % (parameters,
            time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
            time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass
        jobid = get_jobid(cur)
        add_path(cur, opts.workdir, 'WORKDIR', jobid)
        add_path(cur, opts.fastq  ,  'FASTQ' , jobid, opts.workdir)
        add_path(cur, opts.index  , 'INDEX'  , jobid, opts.workdir)
        for i, (out, num) in enumerate(outfiles):
            try:
                window = opts.windows[i]
            except IndexError:
                window = opts.windows[-1]
            except TypeError:
                window = 'None'
            add_path(cur, out, 'SAM/MAP', jobid, opts.workdir)
            frag = ('none' if opts.iterative else 'frag' if i==len(outfiles) - 1
                    else 'full')
            try:
                cur.execute("""
    insert into FASTQs
     (Id  , PATHid, Entries, Trim, Frag, Read, Enzyme, WRKDIRid, SAMid, INDEXid)
    values
     (NULL,      %d,     %d, '%s', '%s',   %d,   '%s',       %d,    %d,      %d)
     """ % (get_path_id(cur, opts.fastq, opts.workdir), num, window, frag,
            opts.read, opts.renz, get_path_id(cur, opts.workdir),
            get_path_id(cur, out, opts.workdir),
            get_path_id(cur, opts.index, opts.workdir)))
            except lite.IntegrityError:
                pass
        print_db(cur, 'FASTQs')
        print_db(cur, 'PATHs' )
        print_db(cur, 'JOBs'  )

def get_options_from_cfg(cfg_file, opts):
    raise NotImplementedError()

    
