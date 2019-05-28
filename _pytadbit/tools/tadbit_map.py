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

from os                                   import path, remove, system
from string                               import ascii_letters
from random                               import random
from shutil                               import copyfile
from multiprocessing                      import cpu_count
from argparse                             import HelpFormatter
from traceback                            import print_exc
import logging
import sqlite3 as lite
import time

from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, identify_re
from pytadbit.utils.fastq_utils           import quality_plot
from pytadbit.utils.file_handling         import which, mkdir, is_fastq
from pytadbit.mapping.full_mapper         import full_mapping
from pytadbit.utils.sqlite_utils          import get_path_id, add_path, print_db
from pytadbit.utils.sqlite_utils          import get_jobid, already_run, digest_parameters
from pytadbit                             import get_dependencies_version

DESC = "Map Hi-C reads and organize results in an output working directory"

def run(opts):
    check_options(opts)

    launch_time = time.localtime()

    # hash that gonna be append to output file names
    param_hash = digest_parameters(opts, get_md5=True)

    # create tmp directory
    if not opts.tmp:
        temp_dir = opts.workdir + '_tmp_r%d_%s' % (opts.read, param_hash)
    else:
        temp_dir = path.join(opts.tmp,
                             'TADbit_tmp_r%d_%s' % (opts.read, param_hash))

    # QC plot
    fig_path = path.join(opts.workdir,
                            '%s_%s_%s.png' % (path.split(opts.fastq)[-1],
                                              '-'.join(map(str, opts.renz)), param_hash))
    logging.info('Generating Hi-C QC plot')

    dangling_ends, ligated = quality_plot(opts.fastq, r_enz=opts.renz,
                                            nreads=100000, paired=False,
                                            savefig=fig_path)
    for renz in dangling_ends:
        logging.info('  - Dangling-ends (sensu-stricto): %.3f%%', dangling_ends[renz])
    for renz in ligated:
        logging.info('  - Ligation sites: %.3f%%', ligated[renz])
    if opts.skip_mapping:
        save_to_db(opts, dangling_ends, ligated, fig_path, [], launch_time, time.localtime())
        return

    # Mapping
    logging.info('mapping %s read %s to %s', opts.fastq, opts.read, opts.workdir)

    outfiles = full_mapping(opts.index, opts.fastq,
                            path.join(opts.workdir,
                                      '01_mapped_r%d' % (opts.read)), mapper=opts.mapper,
                            r_enz=opts.renz, temp_dir=temp_dir, nthreads=opts.cpus,
                            frag_map=not opts.iterative, clean=not opts.keep_tmp,
                            windows=opts.windows, get_nread=True, skip=opts.skip,
                            suffix=param_hash, mapper_params=opts.mapper_param)

    # adjust line count
    if opts.skip:
        for i, (out, _) in enumerate(outfiles[1:], 1):
            outfiles[i] = out, outfiles[i-1][1] - sum(1 for _ in open(outfiles[i-1][0]))

    finish_time = time.localtime()

    # save all job information to sqlite DB
    save_to_db(opts, dangling_ends, ligated, fig_path, outfiles, launch_time, finish_time)
    try:
        save_to_db(opts, dangling_ends, ligated, fig_path, outfiles, launch_time, finish_time)
    except Exception as e:
        # release lock
        remove(path.join(opts.workdir, '__lock_db'))
        print_exc()
        exit(1)

    # write machine log
    try:
        while path.exists(path.join(opts.workdir, '__lock_log')):
            time.sleep(0.5)
            open(path.join(opts.workdir, '__lock_log'), 'a').close()
        with open(path.join(opts.workdir, 'trace.log'), "a") as mlog:
            mlog.write('\n'.join([
                ('# MAPPED READ%s\t%d\t%s' % (opts.read, num, out))
                for out, num in outfiles]) + '\n')
            # release lock
        try:
            remove(path.join(opts.workdir, '__lock_log'))
        except OSError:
            pass
    except Exception as e:
        # release lock
        remove(path.join(opts.workdir, '__lock_db'))
        print_exc()
        exit(1)

    # clean
    if not opts.keep_tmp:
        logging.info('cleaning temporary files')
        system('rm -rf ' + temp_dir)


def check_options(opts):

    opts.mapper_binary = which(opts.mapper_binary)
    if not opts.mapper_binary:
        raise Exception('\n\nERROR: Mapper binary not found, for GEM install it from:'
                        '\nhttps://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%202/'
                        '\n - Download the GEM-binaries-Linux-x86_64-core_i3 if'
                        'have a recent computer, the '
                        'GEM-binaries-Linux-x86_64-core_2 otherwise\n - '
                        'Uncompress with "tar xjvf GEM-binaries-xxx.tbz2"\n - '
                        'Copy the binary gem-mapper to /usr/local/bin/ for '
                        'example (somewhere in your PATH).\n\nNOTE: GEM does '
                        'not provide any binary for MAC-OS.')

    # check RE name
    if opts.renz == ['CHECK']:
        print '\nSearching for most probable restriction enzyme in file: %s' % (opts.fastq)
        try:
            pat, enz, pv = identify_re(opts.fastq, nreads=100000)
            print ' -> Most probable digested site: %s (pv: %f)' % (pat, pv)
            print ' -> Enzymes matching: %s' % (', '.join(enz))
        except ValueError:
            print ' -> Nothing found...'
        exit()
    for n, renz in enumerate(opts.renz):
        if renz == 'NONE':
            opts.renz[n] = None
            continue
        try:
            _ = RESTRICTION_ENZYMES[renz]
        except KeyError:
            print ('\n\nERROR: restriction enzyme %s not found.' % (renz)
                   +  'Use one of:\n\n'
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

    # check paths
    if opts.mapper == 'gem' and not path.exists(opts.index):
        raise IOError('ERROR: index file not found at ' + opts.index)

    if not path.exists(opts.fastq):
        raise IOError('ERROR: FASTQ file not found at ' + opts.fastq)

    if not is_fastq(opts.fastq):
        raise IOError(('ERROR: FASTQ file %s wrong format, check') % (opts.fastq))

    try:
        opts.windows = [[int(i) for i in win.split(':')]
                        for win in opts.windows]
    except TypeError:
        pass

    mkdir(opts.workdir)
    # write log
    # if opts.mapping_only:
    log_format = '[MAPPING {} READ{}]   %(message)s'.format(opts.fastq, opts.read)
    # else:
    #     log_format = '[DEFAULT]   %(message)s'

    # reset logging
    logging.getLogger().handlers = []

    try:
        print 'Writing log to ' + path.join(opts.workdir, 'process.log')
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
        logging.info('Writing versions of TADbit and dependencies')
        vlog = open(vlog_path, 'w')
        vlog.write(dependencies)
        vlog.close()

    # check GEM mapper extra options
    if opts.mapper_param:
        opts.mapper_param = dict([o.split(':') for o in opts.mapper_param])
    else:
        opts.mapper_param = {}
    if opts.mapper == 'gem':
        gem_valid_option = set(["granularity", "q", "quality-format",
                                "gem-quality-threshold", "mismatch-alphabet",
                                "m", "e", "min-matched-bases",
                                "max-big-indel-length", "s", "strata-after-best",
                                "fast-mapping", "unique-mapping", "d", "D",
                                "allow-incomplete-strata", "max-decoded-matches",
                                "min-decoded-strata", "p", "paired-end-alignment",
                                "b", "map-both-ends", "min-insert-size",
                                "max-insert-size", "E", "max-extendable-matches",
                                "max-extensions-per-match", "unique-pairing"])
        for k in opts.mapper_param:
            if not k in gem_valid_option:
                raise NotImplementedError(('ERROR: option "%s" not a valid GEM option'
                                           'or not suported by this tool.') % k)

    # create empty DB if don't exists
    dbpath = path.join(opts.workdir, 'trace.db')
    open(dbpath, 'a').close()

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
        if 'tmpdb' in opts and opts.tmpdb:
            remove(path.join(dbdir, dbfile))
        exit('WARNING: exact same job already computed, see JOBs table above')


def save_to_db(opts, dangling_ends, ligated, fig_path, outfiles, launch_time, finish_time):
    """
    write little DB to keep track of processes and options
    """
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
        # check if table exists
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='MAPPED_INPUTs'""")
        if not cur.fetchall():
            try:
                cur.execute("""
                create table PATHs
                (Id integer primary key,
                JOBid int, Path text, Type text,
                unique (Path))""")
            except lite.OperationalError:
                pass  # may append when mapped files cleaned
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
            create table MAPPED_INPUTs
               (Id integer primary key,
                PATHid int,
                Entries int,
                Trim text,
                Frag text,
                Read int,
                Enzyme text,
                Dangling_Ends text,
                Ligation_Sites text,
                WRKDIRid int,
                MAPPED_OUTPUTid int,
                INDEXid int,
                unique (PATHid,Entries,Read,Enzyme,WRKDIRid,MAPPED_OUTPUTid,INDEXid))""")

        try:
            parameters = digest_parameters(opts, get_md5=False)
            param_hash = digest_parameters(opts, get_md5=True)
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
        add_path(cur, opts.fastq  ,  'MAPPED_FASTQ' , jobid, opts.workdir)
        add_path(cur, opts.index  , 'INDEX'  , jobid, opts.workdir)
        add_path(cur, fig_path  , 'FIGURE'  , jobid, opts.workdir)
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
    insert into MAPPED_INPUTs
     (Id  , PATHid, Entries, Trim, Frag, Read, Enzyme, Dangling_Ends, Ligation_Sites, WRKDIRid, MAPPED_OUTPUTid, INDEXid)
    values
     (NULL,      %d,     %d, '%s', '%s',   %d,   '%s',         '%s',          '%s',       %d,              %d,      %d)
     """ % (get_path_id(cur, opts.fastq, opts.workdir), num, window, frag,
            opts.read, '-'.join(map(str, opts.renz)),
            ' '.join('%s:%.3f%%' % (r, dangling_ends.get(r, float('nan'))) for r in opts.renz),
            ' '.join('%s:%.3f%%' % ('-'.join(r), ligated.get(r, float('nan'))) for r in ligated),
            get_path_id(cur, opts.workdir),
            get_path_id(cur, out, opts.workdir),
            get_path_id(cur, opts.index, opts.workdir)))
            except lite.IntegrityError:
                pass
        print_db(cur, 'MAPPED_INPUTs')
        print_db(cur, 'PATHs' )
        print_db(cur, 'JOBs'  )
    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
    # release lock
    try:
        remove(path.join(opts.workdir, '__lock_db'))
    except OSError:
        pass

def get_options_from_cfg(cfg_file, opts):
    raise NotImplementedError()



def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')
    mapper = parser.add_argument_group('Mapping options')
    descro = parser.add_argument_group('Descriptive, optional arguments')

    # glopts.add_argument('--cfg', dest='cfg', metavar="PATH", action='store',
    #                     default=None, type=str,
    #                     help='path to a configuration file with predefined ' +
    #                     'parameters')

    glopts.add_argument('--skip_mapping', dest='skip_mapping', action='store_true',
                        default=False,
                        help='generate a Hi-C specific quality plot from FASTQ and exits')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='path to an output folder.')

    glopts.add_argument('--fastq', dest='fastq', metavar="PATH", action='store',
                        default=None, type=str, required=True,
                        help='path to a FASTQ files (can be compressed files)')

    glopts.add_argument('--index', dest='index', metavar="PATH",
                        type=str, required=True,
                        help='''paths to file(s) with indexed FASTA files of the
                        reference genome.''')

    glopts.add_argument('--read', dest='read', metavar="INT",
                        type=int, required=True,
                        help='read number')

    glopts.add_argument('--renz', dest='renz', metavar="STR",
                        type=str, required=True, nargs='+',
                        help='''restriction enzyme name(s). Use "--renz CHECK"
                        to search for most probable and exit; and use
                        "--renz NONE" to avoid using RE site information.''')

    glopts.add_argument('--chr_name', dest='chr_name', metavar="STR", nargs='+',
                        default=[], type=str,
                        help='''[fasta header] chromosome name(s). Used in the
                        same order as data.''')

    glopts.add_argument('--tmp', dest='tmp', metavar="PATH", action='store',
                        default=None, type=str,
                        help='''path to a temporary directory (default next to
                        "workdir" directory)''')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument('--noX', action='store_true', help='no display server (X screen)')

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

    glopts.add_argument('--keep_tmp', dest='keep_tmp', action='store_true',
                        default=False,
                        help='[DEBUG] keep temporary files.')

    mapper.add_argument("-C", "--cpus", dest="cpus", type=int,
                        default=cpu_count(), help='''[%(default)s] Maximum
                        number of CPU cores  available in the execution host.
                         If higher than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    mapper.add_argument('--mapper', dest='mapper', metavar="STR",
                        type=str, default='gem',
                        help='[%(default)s] mapper used, options are gem or bowtie2')

    mapper.add_argument('--mapper_binary', dest='mapper_binary', metavar="STR",
                        type=str, default='gem-mapper',
                        help='[%(default)s] path to mapper binary')

    mapper.add_argument('--mapper_param', dest="mapper_param", type=str, default=0,
                        nargs='+',
                        help='''any parameter that could be passed to the GEM or BOWTIE2
                        mapper. e.g. if we want to set the proportion of
                        mismatches to 0.05 and the maximum indel length to 10,
                        (in GEM it would be: -e 0.05 --max-big-indel-length 10),
                        here we could write: "--gem_param e:0.05
                        max-big-indel-length:10". IMPORTANT: some options are
                        incompatible with 3C-derived experiments.''')
