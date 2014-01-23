"""
17 May 2013


Liftover (1) wrapper applied to the comparison of topologically associated
domains.


This script allows to compare Hi-C experiments (mainly align TAD boundaries)
done with different assemblies (e.g.: NCBI36 and GRCh37 for human genome), or
in different species.


INSTALL:

 - liftover tool needs to be downloaded from
   (http://hgdownload.cse.ucsc.edu/admin/exe/), installed, and appended to the
   path.
  - depending on the data a 'chain' file may also be downloaded. For example
    from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/


(1) Fujita, P. A., Rhead, B., Zweig, A. S., Hinrichs, A. S., Karolchik, D.,
    Cline, M. S., Goldman, M., et al. (2011).
    The UCSC Genome Browser database: update 2011.
    Nucleic Acids Research, 39(Database issue), D876-82. doi:10.1093/nar/gkq963
    
"""

from os                        import system, listdir
from os.path                   import isdir
from pytadbit                  import load_chromosome
from pytadbit.utils.remap_tads import remap_chr, reorder
from optparse                  import OptionParser


def check_pik(path):
    with open(path, "r") as f:
        f.seek (0, 2)           # Seek @ EOF
        fsize = f.tell()        # Get Size
        f.seek (max (fsize-2, 0), 0) # Set pos @ last n chars
        key = f.read()       # Read to end
    return key == 's.'


def main():
    """
    main function
    """
    opts = get_options()
    res = opts.res
    if opts.genomes:
        # load all chromosomes of reference genomes
        ref_genome = {}
        for crm in listdir(opts.ref_genome):
            crm_path = opts.ref_genome + crm + '/'
            if not isdir(crm_path):
                continue
            for crm_fh in listdir(crm_path):
                crm_pik = crm_path + crm_fh
                if not check_pik(crm_pik):
                    continue
                ref_genome[crm] = load_chromosome(crm_pik)
        if not opts.res:
            resolutions = []
            for crm in ref_genome:
                for exp in ref_genome[crm].experiments:
                    resolutions.append(exp.resolution)
            if not all([r == resolutions[0] for r in resolutions]):
                raise AssertionError('Not all Experiments have the ' +
                                     'same resolution\n')
            res = resolutions[0]
        alt_genomes = {}
        for i, genome in enumerate(opts.genomes):
            alt_genomes[i] = {}
            for crm in listdir(genome):
                crm_path = genome + crm + '/'
                if not isdir(crm_path):
                    continue
                for crm_fh in listdir(crm_path):
                    crm_pik = crm_path + crm_fh
                    if not check_pik(crm_pik):
                        continue
                    try:
                        alt_genomes[i][crm] = load_chromosome(crm_pik)
                    except:
                        print ('SKIPPING: {} \n not a valid ' +
                               'chromosome').format(crm_pik)
            genome = {}
            for crm in alt_genomes[i]:
                genome = remap_chr(alt_genomes[i][crm], crm, '/tmp/',
                                   opts.lft_path, opts.chain_path,
                                   genome=genome)
            reorder(genome)
            for exp in genome:
                for crm in genome[exp]:
                    try:
                        ref_genome[crm].add_experiment(
                            exp, res, tad_handler=genome[exp][crm])
                    except KeyError:
                        print ('Chromosome {} skipped, not in reference ' +
                               'genome').format(crm)

    system('mkdir -p ' + opts.out_path)
    for crm in ref_genome:
        system('mkdir -p ' + opts.out_path + '/' + crm)
        out_f = opts.out_path + '/' + crm + '/chr' + crm + '.tdb'
        ref_genome[crm].save_chromosome(out_f, force=True)

    # TODO: the same for 1 chromosome
    

def get_options():
    '''
    parse option from call
    '''
    def vararg_callback(option, _, value, parser):
        assert value is None
        value = []
        rargs = parser.rargs
        while rargs:
            arg = rargs[0]
            if ((arg[:2] == "--" and len(arg) > 2) or
                (arg[:1] == "-" and len(arg) > 1 and arg[1] != "-")):
                break
            else:
                value.append(arg)
                del rargs[0]
        setattr(parser.values, option.dest, value)
    #
    parser = OptionParser(
        usage=("%prog [options] file [options] file [options] " +
               "file [options [file ...]]"))
    parser.add_option('--genomes', dest='genomes', metavar="PATH",
                      action='callback', default=None,
                      callback=vararg_callback, 
                      help='''path(s) to a directory/ies with a list of
                      chromosomes saved through tadbit (required if not
                      passing chromosomes)''')
    parser.add_option('--ref_genome', dest='ref_genome', metavar="PATH", 
                      help='''path to a directory with a list of chromosomes
                      saved through tadbit (required with genomes option)''')
    parser.add_option('--crm', dest='crm', metavar="PATH", 
                      help='''path to input file, a chromosome saved through
                      tadbit (required if not passing genomes)''')
    parser.add_option('--ref_crm', dest='ref_crm', metavar="PATH", 
                      help='''path to second input file, a reference chromosome
                      saved through tadbit (required)''')
    parser.add_option('--chain', dest='chain_path', action="store", \
                      help=
                      '''path to UCSC chain file (required)''')
    parser.add_option('-o', dest='out_path', metavar="PATH",
                      default='./',
                      help='''path to out file where merged tadbit chromosome
                      will be stored''')
    parser.add_option('--res', dest='res',
                      default=None,
                      help='''Wanted resolution for the detection of TADs (i.e.:
                      100Kb)''')
    parser.add_option('--crm_name', dest='crm_name',
                      default=None,
                      help='''Chromosome name for crm1 (e.g. 21).''')
    parser.add_option('--tmp', dest='tmp_path', metavar="PATH",
                      default='./',
                      help='''path to temporary directory to store liftover
                      outfiles''')
    parser.add_option('--liftover', 
                      dest='lft_path', default='/usr/local/bin/',\
                      help='''[%default] path to liftover binary''')
    opts = parser.parse_args()[0]
    if not opts.crm or not opts.ref_crm or not opts.chain_path:
        if not opts.genomes or not opts.ref_genome or not opts.chain_path:
            exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
