"""
22 Jan 2014

This script uses ensembl REST service to convert genomic coordinates of one
species to another species.
As borders from one chromosome in one species might be distribute along several
chromosome in an other species, this script needs to be used over whole genomes.

WARNING: Ensembl REST server is beta and only provides for latest genome builts
WARNING: needs internet :)
"""

from pytadbit.utils.remap_tads_ensembl import load_genome, remap_genome
from pytadbit.utils.remap_tads_ensembl import save_new_genome, decode_resolution
from optparse                          import OptionParser
import os


def main():
    """
    main function
    """
    opts = get_options()

    # load all chromosomes of original genome
    genome = {}
    if opts.verbose:
        print '\nLoad %s chromosomes from \n%s:' % (
            opts.original_species.capitalize().replace('_', ' ') + "'s",
            ' ' * 10 + opts.genome)
    genome = load_genome(opts.genome, res=opts.res, verbose=opts.verbose)

    # remap TADs
    if opts.verbose:
        print '\nCreating new TADbit chromosomes with new coordinates\n'

    # create directory for output
    rootdir = os.path.abspath(opts.out_path)
    if not os.path.exists(rootdir):
        os.mkdir(rootdir)

    # remap TAD coordinates of all chromosomes
    new_genome, trace = remap_genome(genome, opts.original_assembly,
                                     opts.target_assembly,
                                     opts.original_species,
                                     opts.target_species, opts.log)

    # save new chromosomes
    save_new_genome(genome if opts.check else new_genome, trace, opts.check,
                    opts.target_species, rootdir)


def get_options():
    """
    parse option from call

    sys.argv = ['', '--genome', '/home/fransua/db/hi-c/dixon_hsap-mmus/mmus/100Kb/',
                '--species', 'mus_musculus', '--to_species', 'homo_sapiens',
                '--from_map', 'NCBIM37', '--to_map', 'GRCm38', '--res', '100000',
                '-o', '/home/fransua/Box/tadbit/scripts/hsap2mmus/']

    """

    parser = OptionParser(
        usage=("%prog [options] file [options] file [options] " +
               "file [options [file ...]]"))
    parser.add_option('--genome', dest='genome', metavar="PATH",
                      action='store', default=None,
                      help='''path to a directory with a list of
                      chromosomes saved through tadbit (required if not
                      passing chromosomes)''')
    # parser.add_option('--crm', dest='crm', metavar="PATH",
    #                   help='''path to input file, a chromosome saved through
    #                   tadbit (required if not passing genome)''')
    parser.add_option('--species', dest='original_species', metavar="STRING",
                      help='''Species name (no spaces in name, use underscore)
                      of the input chromosome(s) (i.e.: homo_sapiens''')
    parser.add_option('--to_species', dest='target_species', metavar="STRING",
                      default=None,
                      help='''Name of the species name (no spaces in name, use
                      underscore) to which we want to remap the chromosome(s)
                      (i.e.: mus_musculus''')
    parser.add_option('--from_map', dest='original_assembly', metavar="STRING",
                      default=None,
                      help='''NCBI ID of the original assembly (i.e.: NCBIM37
                      for mouse or NCBI36 for human)''')
    parser.add_option('--to_map', dest='target_assembly', metavar="STRING",
                      default=None,
                      help='''NCBI ID of the assembly we want to map to (i.e.:
                      GRCm38 for mouse or GRCh37 for human) -- Needed with
                      "from_map"''')
    parser.add_option('--check',
                      dest='check', default=False, action='store_true',
                      help='''[%default] do not convert coordinates, just check
                      if convertion is possible. If not border is removed.''')
    parser.add_option('-o', dest='out_path', metavar="PATH",
                      default='./',
                      help='''[%default] path to out file where converted TADbit
                      chromosome(s) will be stored''')
    parser.add_option('--log', dest='log', metavar="PATH",
                      default='',
                      help='''[%default] path to log file.''')
    parser.add_option('--res', dest='res',
                      default=None,
                      help='''Wanted resolution for the detection of TADs (i.e.:
                      100Kb)''')
    parser.add_option('--verbose',
                      dest='verbose', default=True, action='store_false',
                      help='''[%default] verbosity''')
    parser.add_option('--continue',
                      dest='skip', default=True, action='store_false',
                      help='''[%default] skips chromosomes already saved in the
                      output directory''')
    opts = parser.parse_args()[0]

    if not opts.genome:
        exit(parser.print_help())
    if not opts.original_species:
        print '   Needs an from_species argument\n\n'
        exit(parser.print_help())
    if opts.original_assembly and not opts.target_assembly:
        print '   Needs an to_map argument if original assembly is passed\n\n'
        exit(parser.print_help())
    if opts.res:
        try:
            opts.res = int(opts.res)
        except ValueError:
            opts.res = decode_resolution(opts.res)

    return opts



if __name__ == "__main__":
    exit(main())
