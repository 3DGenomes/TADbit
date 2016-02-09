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

from argparse import ArgumentParser, HelpFormatter
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES
from pytadbit.utils.fastq_utils   import quality_plot
from pytadbit.mapping.full_mapper import full_mapping
from pytadbit import get_dependencies_version
from os import system, path
import logging
import fcntl

def main():
    opts = get_options()

    if opts.quality_plot:
        logging.info('Generating Hi-C QC plot at:\n  ' +
               path.join(opts.output, path.split(opts.fastq)[-1] + '.pdf'))
        quality_plot(opts.fastq, r_enz=opts.renz,
                     nreads=100000, paired=False,
                     savefig=path.join(opts.output,
                                       path.split(opts.fastq)[-1] + '.pdf'))
        return

    windows = opts.windows

    logging.info('mapping %s read %s to %s', opts.fastq, opts.read, opts.output)
    outfiles = full_mapping(opts.index, opts.fastq,
                            path.join(opts.output, '01_mapped_r' + opts.read),
                            opts.renz, temp_dir=opts.tmp,
                            frag_map=opts.strategy=='frag', clean=True,
                            windows=windows)

    # write machine log
    with open(path.join(opts.output, 'trace.log'), "a") as mlog:
        fcntl.flock(mlog, fcntl.LOCK_EX)
        mlog.write('\n'.join([('# MAPPED READ%s PATH\t%d\t' % (opts.read, num)) + out
                              for out, num in outfiles]))
        fcntl.flock(mlog, fcntl.LOCK_UN)

    logging.info('cleaning temporary files')
    # clean
    system('rm -rf ' + opts.tmp)

def get_options():
    """
    parse option from call
    """
    parser = ArgumentParser(
        usage="%(prog)s [options] [--cfg CONFIG_PATH]",
        formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                   max_help_position=27))

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

    glopts.add_argument('--genome', dest='genome', metavar="PATH", nargs='+',
                        type=str,
                        help='''paths to file(s) with FASTA files of the
                        reference genome. If many, files will be concatenated.
                        I.e.: --fasta chr_1.fa chr_2.fa
                        In this last case, order is important or the rest of the
                        analysis.''')

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

    mapper.add_argument('--windows', dest='windows', default='auto',
                        nargs='+',
                        help='''for iterative mapping, defines windows. e.g.
                        --windows 20 25 30 35 40 45 50''')

    mapper.add_argument('--read_length', dest='read_length',
                        type=int,
                        help='''read length, compulsory in iterative mapping with
                        --windows auto''')

    mapper.add_argument('--mapping_only', dest='mapping_only', action='store_true',
                        help='only do the mapping does not parse results')

    descro.add_argument('--species', dest='species', metavar="STR", 
                        type=str,
                        help='species name')

    descro.add_argument('--descr', dest='description', metavar="LIST", nargs='+',
                        type=str,
                        help='''extra descriptive fields each filed separated by
                        coma, and inside each, name and value separated by column: 
                        --descr=cell:lymphoblast,flowcell:C68AEACXX,index:24nf''')

    parser.add_argument_group(glopts)
    parser.add_argument_group(descro)
    parser.add_argument_group(mapper)
    opts = parser.parse_args()

    if opts.cfg:
        get_options_from_cfg(opts.cfg, opts)

    if (opts.strategy == 'iter' and opts.window == 'auto'
        and not opts.read_length):
        raise Exception('ERROR: need to input read_length')
    # check RE name
    try:
        _ = RESTRICTION_ENZYMES[opts.renz]
    except KeyError:
        print ('\n\nERROR: restriction enzyme not found. Use one of:\n\n'
               + ' '.join(sorted(RESTRICTION_ENZYMES)) + '\n\n')
        raise KeyError()

    # check compulsory options
    if not opts.quality_plot:
        if not opts.genome: raise Exception('ERROR: genome option required.')
        if not opts.index : raise Exception('ERROR: index  option required.')
    if not opts.output: raise Exception('ERROR: output option required.')
    if not opts.fastq : raise Exception('ERROR: fastq  option required.')
    if not opts.renz  : raise Exception('ERROR: renz   option required.')
    if not opts.tmp:
        opts.tmp = opts.output + '_tmp_r' + opts.read

    if opts.strategy == 'frag':
        opts.windows = None
        
    if opts.strategy == 'iter':
        raise NotImplementedError()

    # write log
    if opts.mapping_only:
        log_format = '[MAPPING {} READ{}]   %(message)s'.format(opts.fastq, opts.read)
    else:
        log_format = '[DEFAULT]   %(message)s'
    try:
        logging.basicConfig(filename=path.join(opts.output, 'process.log'),
                            level=logging.INFO, format=log_format)
    except IOError:
        logging.basicConfig(filename=path.join(opts.outdir, 'process.log2'),
                            level=logging.INFO, format=log_format)
    logging.getLogger().addHandler(logging.StreamHandler())

    # write version log
    system('mkdir -p ' + opts.output)
    if not path.exists(path.join(opts.output,
                                 'TADbit_and_dependencies_versions.log')):
        vlog = path.join(opts.output, 'TADbit_and_dependencies_versions.log')
        vlog = open(vlog, 'w')
        vlog.write(get_dependencies_version())
        vlog.close()

    return opts

def get_options_from_cfg(cfg_file, opts):
    raise NotImplementedError()


if __name__ == "__main__":
    exit(main())
