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

from argparse import HelpFormatter
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES
from pytadbit.utils.fastq_utils   import quality_plot
from pytadbit.mapping.full_mapper import full_mapping
from pytadbit import get_dependencies_version
from os import system, path
from multiprocessing import cpu_count
import logging
import fcntl
import sqlite3 as lite

DESC = "Map Hi-C reads and organize results in an output working directory"

def run(opts):
    check_options(opts)

    if opts.quality_plot:
        logging.info('Generating Hi-C QC plot at:\n  ' +
               path.join(opts.workdir, path.split(opts.fastq)[-1] + '.pdf'))
        quality_plot(opts.fastq, r_enz=opts.renz,
                     nreads=100000, paired=False,
                     savefig=path.join(opts.workdir,
                                       path.split(opts.fastq)[-1] + '.pdf'))
        return

    logging.info('mapping %s read %s to %s', opts.fastq, opts.read, opts.workdir)
    outfiles = full_mapping(opts.index, opts.fastq,
                            path.join(opts.workdir, '01_mapped_r%d' % opts.read),
                            opts.renz, temp_dir=opts.tmp, nthreads=opts.cpus,
                            frag_map=not opts.iterative, clean=True,
                            windows=opts.windows, get_nread=True, skip=opts.skip)

    # adjust line count
    if opts.skip:
        for i, (out, _) in enumerate(outfiles[1:], 1):
            outfiles[i] = out, outfiles[i-1][1] - sum(1 for _ in open(outfiles[i-1][0]))
    
    # save all job information to sqlite DB
    save_to_db(opts, outfiles)
    
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

    glopts.add_argument('-o', '--workdir', dest='workdir', metavar="PATH",
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
        raise Exception('ERROR: can use output files, workdir missing')

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

def save_to_db(opts, outfiles):
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
                Path text, Type text,
                unique (Path))""")
            cur.execute("""
            create table FASTQs
               (Id integer primary key,
                FASTQid int,
                Entries int,
                Trim text,
                Frag int,
                Read int,
                Enzyme text,
                WRKDIRid int,
                SAMid int,
                INDEXid int,
                unique (FASTQid,Entries,Read,Enzyme,WRKDIRid,SAMid,INDEXid))""")
        add_path(cur, opts.workdir, 'WORKDIR')
        add_path(cur, opts.fastq,  'FASTQ')
        add_path(cur, opts.index, 'INDEX')
        for i, (out, num) in enumerate(outfiles):
            try:
                window = opts.windows[i]
            except IndexError:
                window = opts.windows[-1]
            except TypeError:
                window = 'None'
            add_path(cur, out, 'SAM/MAP')
            try:
                cur.execute("""
    insert into FASTQs
     (Id  , FASTQid, Entries, Trim, Frag, Read, Enzyme, WRKDIRid, SAMid, INDEXid)
    values
     (NULL,      %d,      %d, '%s',   %d,   %d,   '%s',       %d,    %d,      %d)
     """ % (get_id(cur, opts.fastq), num, window, not opts.iterative,
            opts.read, opts.renz, get_id(cur, opts.workdir), get_id(cur, out),
            get_id(cur, opts.index)))
            except lite.IntegrityError:
                pass
        print_db(cur, 'FASTQs')
        print_db(cur, 'PATHs')

def get_options_from_cfg(cfg_file, opts):
    raise NotImplementedError()

def add_path(cur, path, type):
    try:
        cur.execute("insert into PATHs (Id  , Path, Type) values (NULL, '%s', '%s')" % (
            path, type))
    except lite.IntegrityError:
        pass

def get_id(cur, name):
    cur.execute('SELECT Id from PATHS where Path="%s"' % name)
    return cur.fetchall()[0][0]

def print_db(cur, name):
    cur.execute('select * from %s' % name)
    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    cols = [max(vals) for vals in zip(*[[len(str(v)) for v in row]
                                        for row in rows + [names]])]
    print ',-' + '-.-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-.'
    print '| ' + ' | '.join([('%{}s'.format(cols[i])) % str(v) for i, v in enumerate(names)]) + ' |'
    print '|-' + '-+-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-|'
    print '| ' + '\n| '.join([' | '.join([('%{}s'.format(cols[i])) % str(v)
                                        for i, v in enumerate(row)]) + ' |'  for row in rows])
    print "'-" + '-^-'.join(['-' * cols[i] for i, v in enumerate(names)]) + "-'"
    
