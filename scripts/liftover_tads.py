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

from subprocess import Popen, PIPE
from os         import system, listdir
from os.path    import isdir
from pytadbit   import load_chromosome
from optparse   import OptionParser


def liftover(coords, tmp_path, lft_path, chn_path):
    tmp = open(tmp_path + '/tmp', 'w')
    for coord in coords:
        tmp.write('chr{}\t{}\t{}\n'.format(coord[1], coord[2], coord[2]+1))
    tmp.close()
    _, err = Popen((lft_path + 'liftOver ' +
                    tmp_path + '/tmp ' +
                    chn_path + ' ' +
                    tmp_path + '/out ' +
                    tmp_path + '/out.log'
                    ), shell =True, stdout=PIPE, stderr=PIPE).communicate()
    if (not 'Reading' in err) or (not 'Mapping' in err):
        raise Exception(err)
    founds = [[l.split()[0],
               int(l.split()[1])] for l in open(tmp_path + '/out').readlines()]
    missed = [int(l.split()[1]) for l in open(tmp_path + '/out.log').readlines()\
              if not l.startswith('#')]
    system('rm -f {}/tmp'.format(tmp_path))
    mapped = {}
    j = 0
    for k, (i, _, end) in enumerate(coords):
        mapped[i] = None
        if end in missed:
            j += 1
            continue
        mapped[i] = (founds[k - j][0][3:], founds[k - j][1])
    return mapped
    
    
def remap_chr(crm_obj, crm, tmp, lft_path, chain_path,
              wnt_exp=None, genome=None):
    missed = 0
    found  = 0
    genome = {} or genome
    for exp in crm_obj.experiments:
        if wnt_exp:
            if not exp.name in wnt_exp:
                continue
        print exp
        coords = []
        res = exp.resolution
        for t in exp.tads:
            coords.append((t, crm,
                           int(exp.tads[t]['end']   * res)))
        mapped = liftover(coords, tmp, lft_path, chain_path)
        for t in mapped:
            try:
                crm2, end = mapped[t]
                genome.setdefault(exp.name, {})
                if crm2 not in genome[exp.name]:
                    genome[exp.name][crm2] = {'end': [],
                                             'score': []}
                print ' -> crm{:<2}:{:>10} | crm{:<2}:{:>10}'.format(
                    crm, int(exp.tads[t]['end']) * 100000,
                    crm2, int(end / res + 0.5) * 100000)
                genome[exp.name][crm2]['end'  ].append(float(int(end / res + 0.5)))
                genome[exp.name][crm2]['score'].append(exp.tads[t]['score'])
                found += 1
            except TypeError:
                missed += 1
    print 'missed: {}, found: {}'.format(missed, found)
    return genome


def reorder(genome):
    for exp in genome:
        for crm in genome[exp]:
            srt = sorted(genome[exp][crm]['end'])
            order = [genome[exp][crm]['end'].index(elt) for elt in srt]
            genome[exp][crm]['end']   = [genome[exp][crm]['end'  ][i] for i in order]
            genome[exp][crm]['score'] = [genome[exp][crm]['score'][i] for i in order]
            todel = []
            for brk in xrange(1, len(genome[exp][crm]['end'])):
                if genome[exp][crm]['end'][brk] == genome[exp][crm]['end'][brk-1]:
                    todel.append(brk-1)
                    genome[exp][crm]['score'][brk] = max(
                        genome[exp][crm]['score'][brk], genome[exp][crm]['score'][brk-1])
            for brk in todel[::-1]:
                genome[exp][crm]['end'].pop(brk)
                genome[exp][crm]['score'].pop(brk)
            genome[exp][crm]['start'] = []
            for brk in xrange(len(genome[exp][crm]['end']) + 1):
                genome[exp][crm]['start'].append(genome[exp][crm]['end'][brk - 1] + 1 \
                                                 if brk > 0 else 0)
                    

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
