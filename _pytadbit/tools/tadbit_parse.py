"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""

from argparse                       import HelpFormatter
from pytadbit                       import get_dependencies_version
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.parsers.new_map_parser    import parse_map
from os                             import system, path
from pytadbit.utils.sqlite_utils    import get_path_id, add_path, print_db
from hashlib                        import md5
import time
import logging
import fcntl
from cPickle import load, UnpicklingError
import sqlite3 as lite

DESC = "Parse mapped Hi-C reads and get the intersection"

def run(opts):
    check_options(opts)

    launch_time = time.localtime()
    
    f_names1, f_names2, renz = load_parameters_fromdb(opts.workdir)

    name = path.split(opts.workdir)[-1]
    
    system('mkdir -p ' + path.join(opts.workdir, '02_parsed_reads'))
    if not opts.read:
        out_file1 = path.join(opts.workdir, '02_parsed_reads', '%s_r1.tsv' % name)
        out_file2 = path.join(opts.workdir, '02_parsed_reads', '%s_r2.tsv' % name)
    elif opts.read == 1:
        out_file1 = path.join(opts.workdir, '02_parsed_reads', '%s_r1.tsv' % name)
        out_file2 = None
        f_names2  = None
    elif opts.read == 2:
        out_file2 = None
        f_names1  = f_names2
        f_names2  = None
        out_file1 = path.join(opts.workdir, '02_parsed_reads', '%s_r2.tsv' % name)
        
    logging.info('parsing genomic sequence')
    try:
        # allows the use of cPickle genome to make it faster
        genome = load(open(opts.genome[0]))
    except UnpicklingError:
        genome = parse_fasta(opts.genome)

    if not opts.skip:
        logging.info('parsing reads in %s project', name)
        counts, multis = parse_map(f_names1, f_names2, out_file1=out_file1,
                                   out_file2=out_file2, re_name=renz, verbose=True,
                                   genome_seq=genome)
    else:
        counts = {}
        counts[0] = {}
        fhandler = open(out_file1)
        for line in fhandler:
            if line.startswith('# MAPPED '):
                _, _, item, value = line.split()
                counts[0][item] = int(value)
            elif not line.startswith('#'):
                break
        multis = {}
        multis[0] = 0
        for line in fhandler:
            if '|||' in line:
                multis[0] += line.count('|||')
        if out_file2:
            counts[1] = {}
            fhandler = open(out_file2)
            for line in fhandler:
                if line.startswith('# MAPPED '):
                    _, _, item, value = line.split()
                    counts[1][item] = int(value)
                elif not line.startswith('#'):
                    break
            multis[1] = 0
            for line in fhandler:
                if '|||' in line:
                    multis[1] += line.count('|||')                

    # write machine log
    with open(path.join(opts.workdir, 'trace.log'), "a") as mlog:
        fcntl.flock(mlog, fcntl.LOCK_EX)
        for read in counts:
            for item in counts[read]:
                mlog.write('# PARSED READ%s PATH\t%d\t%s\n' % (
                    read, counts[read][item],
                    out_file1 if read == 1 else out_file2))
        fcntl.flock(mlog, fcntl.LOCK_UN)

    finish_time = time.localtime()

    # save all job information to sqlite DB
    save_to_db(opts, counts, multis, f_names1, f_names2, out_file1, out_file2,
               launch_time, finish_time)

def save_to_db(opts, counts, multis, f_names1, f_names2, out_file1, out_file2,
               launch_time, finish_time):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='SAMs'""")
        if not cur.fetchall():
            cur.execute("""
        create table SAMs
           (Id integer primary key,
            PATHid int,
            BEDid int,
            Entries int,
            unique (PATHid, BEDid, Entries))""")
            cur.execute("""
        create table BEDs
           (Id integer primary key,
            PATHid int,
            Entries int,
            Multiples int,
            unique (PATHid, Entries))""")
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
        cur.execute("select Id from JOBs where Id = (select max(id)  from JOBs)")
        jobid = cur.fetchall()[0][0]        
        add_path(cur, out_file1, 'BED', jobid, opts.workdir)
        for genome in opts.genome:
            add_path(cur, genome, 'FASTA', jobid, opts.workdir)
        if out_file2:
            add_path(cur, out_file2, 'BED', jobid, opts.workdir)
        if not opts.read == 2:
            try:
                sum_reads = 0
                for i, item in enumerate(counts[0]):
                    cur.execute("""
                    insert into SAMs
                    (Id  , PATHid, BEDid, Entries)
                    values
                    (NULL,    %d,     %d,      %d)
                    """ % (get_path_id(cur, f_names1[i], opts.workdir),
                           get_path_id(cur, out_file1, opts.workdir),
                           counts[0][item]))
                    sum_reads += counts[0][item]
                cur.execute("""
                insert into BEDs
                (Id  , PATHid, Entries, Multiples)
                values
                (NULL,     %d,      %d,        %d)
                """ % (get_path_id(cur, out_file1, opts.workdir),
                       sum_reads, multis[0]))
            except lite.IntegrityError:
                print 'WARNING: already parsed'
        if not opts.read == 1:
            try:
                sum_reads = 0
                for i, item in enumerate(counts[1]):
                    cur.execute("""
                    insert into SAMs
                    (Id  , PATHid, BEDid, Entries)
                    values
                    (NULL,     %d,    %d,      %d)
                    """ % (get_path_id(cur, f_names2[i], opts.workdir),
                           get_path_id(cur, out_file2, opts.workdir),
                           counts[1][item]))
                    sum_reads += counts[1][item]
                cur.execute("""
                insert into BEDs
                (Id  , PATHid, Entries, Multiples)
                values
                (NULL,     %d,      %d,        %d)
                """ % (get_path_id(cur, out_file2, opts.workdir),
                       sum_reads, multis[1]))
            except lite.IntegrityError:
                print 'WARNING: already parsed'
        print_db(cur, 'FASTQs')
        print_db(cur, 'PATHs')
        print_db(cur, 'SAMs')
        print_db(cur, 'BEDs')
        print_db(cur, 'JOBs')

def load_parameters_fromdb(workdir):
    con = lite.connect(path.join(workdir, 'trace.db'))
    fnames1 = []
    fnames2 = []
    ids = []
    with con:
        cur = con.cursor()
        # fetch file names to parse
        cur.execute("""
        select distinct PATHs.Id,PATHs.Path from PATHs
        inner join FASTQs on PATHs.Id = FASTQs.SAMid
        where FASTQs.Read = 1
        """)
        for fname in cur.fetchall():
            ids.append(fname[0])
            fnames1.append(path.join(workdir, fname[1]))
        cur.execute("""
        select distinct PATHs.Id,PATHs.Path from PATHs
        inner join FASTQs on PATHs.Id = FASTQs.SAMid
        where FASTQs.Read = 2
        """)
        for fname in cur.fetchall():
            ids.append(fname[0])
            fnames2.append(path.join(workdir, fname[1]))
        # GET enzyme name
        enzymes = []
        for fid in ids:
            cur.execute("""
            select distinct FASTQs.Enzyme from FASTQs
            where FASTQs.SAMid=%d
            """ % fid)
            enzymes.extend(cur.fetchall())
        if len(set(reduce(lambda x, y: x+ y, enzymes))) != 1:
            raise Exception(
                'ERROR: different enzymes used to generate these files')
        renz = enzymes[0][0]
    return fnames1, fnames2, renz
        

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

    glopts.add_argument('--type', dest='type', metavar="STR", 
                        type=str, default='map', choices=['map', 'sam', 'bam'], 
                        help='''[%(default)s]file type to be parser, map
                        (GEM-mapper), sam or bam''')

    glopts.add_argument('--read', dest='read', metavar="INT",
                        type=int, default=None, 
                        help='In case only one of the reads needs to be parsed')

    glopts.add_argument('--skip', dest='skip', action='store_true',
                      default=False,
                      help='[DEBUG] in case already mapped.')

    glopts.add_argument('--genome', dest='genome', metavar="PATH", nargs='+',
                        type=str,
                        help='''paths to file(s) with FASTA files of the
                        reference genome. If many, files will be concatenated.
                        I.e.: --fasta chr_1.fa chr_2.fa
                        In this last case, order is important or the rest of the
                        analysis. Note: it can also be the path to a previously
                        parsed genome in pickle format.''')

    parser.add_argument_group(glopts)

def check_options(opts):

    if not opts.workdir: raise Exception('ERROR: output option required.')
    if opts.type != 'map':
        raise NotImplementedError('ERROR: not yet there')

    if not opts.genome: raise Exception('ERROR: genome parameter required.')
    if not opts.workdir: raise Exception('ERROR: workdir parameter required.')

    # check skip
    if not path.exists(opts.workdir) and opts.skip:
        print ('WARNING: can use output files, found, not skipping...')
        opts.skip = False

    if opts.workdir.endswith('/'):
        opts.workdir = opts.workdir[:-1]

    # write log
    log_format = '[PARSING]   %(message)s'

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
