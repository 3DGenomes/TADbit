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


DESC = "Map Hi-C reads and organize results in an output working directory"

def run(opts):
    check_options(opts)

    if opts.quality_plot:
        logging.info('Generating Hi-C QC plot at:\n  ' +
               path.join(opts.output, path.split(opts.fastq)[-1] + '.pdf'))
        quality_plot(opts.fastq, r_enz=opts.renz,
                     nreads=100000, paired=False,
                     savefig=path.join(opts.output,
                                       path.split(opts.fastq)[-1] + '.pdf'))
        return

    logging.info('mapping %s read %s to %s', opts.fastq, opts.read, opts.output)
    outfiles = full_mapping(opts.index, opts.fastq,
                            path.join(opts.output, '01_mapped_r' + opts.read),
                            opts.renz, temp_dir=opts.tmp, nthreads=opts.cpus,
                            frag_map=opts.strategy=='frag', clean=True,
                            windows=opts.windows, get_nread=True)

    # write machine log
    with open(path.join(opts.output, 'trace.log'), "a") as mlog:
        fcntl.flock(mlog, fcntl.LOCK_EX)
        mlog.write('\n'.join([('# MAPPED READ%s PATH\t%d\t' % (opts.read, num)) + out
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

    glopts.add_argument('-o', '--output', dest='output', metavar="PATH",
                        action='store', default=None, type=str,
                        help='path to output folder')

    glopts.add_argument('--fastq', dest='fastq', metavar="PATH", action='store',
                      default=None, type=str,
                      help='path to a FASTQ files (can be compressed files)')

    glopts.add_argument('--index', dest='index', metavar="PATH",
                        type=str,
                        help='''paths to file(s) with indexed FASTA files of the
                        reference genome.''')

    glopts.add_argument('--read', dest='read', metavar="INT", 
                        type=str,
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
                      output directory)''')

    mapper.add_argument('--strategy', dest='strategy', default='frag',
                        choices=['frag', 'iter'],
                        help='''mapping strategy, can be "frag" for fragment
                        based mapping or "iter" for iterative mapping''')

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

    mapper.add_argument("-C", "--cpu", dest="cpus", type=int,
                        default=0, help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    parser.add_argument_group(glopts)
    parser.add_argument_group(descro)
    parser.add_argument_group(mapper)

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

    # number of cpus
    if opts.cpus == 0:
        opts.cpus = cpu_count()
    else:
        opts.cpus = min(opts.cpus, cpu_count())

    # check compulsory options
    if not opts.quality_plot:
        if not opts.index : raise Exception('ERROR: index  parameter required.')
    if not opts.output:     raise Exception('ERROR: output parameter required.')
    if not opts.fastq :     raise Exception('ERROR: fastq  parameter required.')
    if not opts.renz  :     raise Exception('ERROR: renz   parameter required.')
    if not opts.tmp:
        opts.tmp = opts.output + '_tmp_r' + opts.read

    if opts.strategy == 'frag':
        opts.windows = None
        
    if opts.strategy == 'iter':
        raise NotImplementedError()

    system('mkdir -p ' + opts.output)
    # write log
    # if opts.mapping_only:
    log_format = '[MAPPING {} READ{}]   %(message)s'.format(opts.fastq, opts.read)
    # else:
    #     log_format = '[DEFAULT]   %(message)s'

    # reset logging
    logging.getLogger().handlers = []

    try:
        print 'Writting log to ' + path.join(opts.output, 'process.log')
        logging.basicConfig(level=logging.INFO,
                            format=log_format,
                            filename=path.join(opts.output, 'process.log'),
                            filemode='aw')
    except IOError:
        logging.basicConfig(level=logging.DEBUG,
                            format=log_format,
                            filename=path.join(opts.output, 'process.log2'),
                            filemode='aw')

    # to display log on stdout also
    logging.getLogger().addHandler(logging.StreamHandler())

    # write version log
    vlog_path = path.join(opts.output, 'TADbit_and_dependencies_versions.log')
    dependencies = get_dependencies_version()
    if not path.exists(vlog_path) or open(vlog_path).readlines() != dependencies:
        logging.info('Writting versions of TADbit and dependencies')
        vlog = open(vlog_path, 'w')
        vlog.write(dependencies)
        vlog.close()

def get_options_from_cfg(cfg_file, opts):
    raise NotImplementedError()

