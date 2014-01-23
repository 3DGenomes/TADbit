"""
22 Jan 2014

This script uses ensembl REST service to convert genomic coordinates of one
species to another species.
As borders from one chromosome in one species might be distribute along several
chromosome in an other species, this script needs to be used over whole genomes.

WARNING: Ensembl REST server is beta and only provides for latest genome builts
WARNING: needs internet :)
"""

import httplib2, sys
import json
from os                           import system, listdir
from os.path                      import isdir
from optparse                     import OptionParser
from pytadbit.utils.file_handling import check_pik
from pytadbit                     import load_chromosome

def get_best_alignment_coord(ali_list):
    """
    returns coordinates of the longest ungapped alignment
    """
    best_sumlen = 0
    best_ali = None
    for i, ali in enumerate(ali_list):
        sumlen = len([c for c in zip(*[ali['alignments'][0]['seq'],
                                       ali['alignments'][1]['seq']])
                      if not '-' in c])
        if sumlen > best_sumlen:
            best_ali = i
            best_sumlen = sumlen
    return {'chr'    : ali_list[best_ali]['alignments'][1]['seq_region'],
            'start'  : ali_list[best_ali]['alignments'][1]['start'],
            'end'    : ali_list[best_ali]['alignments'][1]['end'],
            'strand' : ali_list[best_ali]['alignments'][1]['strand']}


def rest_query(crm, beg, end, seed_spe='homo_sapiens', targ_spe='mus_musculus',
               alignment='LASTZ_NET', server="http://beta.rest.ensembl.org"):
    """
    Given a seed genomic segment find the coordinate of its syntenis segment on
    a given target species.
    
    :param crm: chromosome name
    :param beg: region start position
    :param end: region end position
    :param seed_spe: name of the seed species (from which com the original
       coordinates). Scientific name with underscore instead of spaces
    :param targ_spe: name of the target species
    :param 'LASTZ_NET' alignment: see http://beta.rest.ensembl.org/documentation/info/genomic_alignment_block_region
       for more explanations
    :param 'http://beta.rest.ensembl.org' server: see http://beta.rest.ensembl.org/documentation/info/genomic_alignment_block_region
       for more explanations

    :returns: a dict with the coordinates of the syntenic segment found (chr,
       start, end, strand)
       
    """
    http = httplib2.Http(".cache")
    output = 'json'
    strand = 1
    ext = ("/alignment/block/region/%s/%s:%s-%s:%s?" +
           "species_set=%s&species_set=%s;method=%s")
    ext = ext % (seed_spe, crm, beg, end, strand, seed_spe, targ_spe, alignment)
    resp, content = http.request(
        server + ext, method="GET",
        headers={"Content-Type":"application/%s" % output})

    if content.startswith('Content-Type application/json had a problem'):
        print 'Region %s:%s-%s not found' % (crm, beg, end)
        return
    
    if not resp.status == 200:
        print "Invalid response: ", resp.status
        sys.exit()

    return get_best_alignment_coord(json.loads(content))


def load_genome(genome_path, res=None, verbose=False):
    ref_genome = {}
    for crm in listdir(genome_path):
        crm_path = genome_path + crm + '/'
        if not isdir(crm_path):
            continue
        for crm_fh in listdir(crm_path):
            crm_pik = crm_path + crm_fh
            if not check_pik(crm_pik):
                continue
            if crm in ref_genome:
                raise Exception('More than 1 TADbit chromosome file found\n')
            if verbose:
                print '  Chromosome:', crm
            ref_genome[crm] = load_chromosome(crm_pik)
    # check resolution
    resolutions = [] if not res else [res]
    for crm in ref_genome:
        for exp in ref_genome[crm].experiments:
            resolutions.append(exp.resolution)
    if not all([r == resolutions[0] for r in resolutions]):
        raise AssertionError('Not all Experiments have the ' +
                             'same resolution\n')
    res = resolutions[0]
    return ref_genome, res
  

def main():
    """
    main function
    """
    opts = get_options()
    res = opts.res

    # load all chromosomes of reference genome
    if opts.verbose:
        print 'Load seed chromosomes from %s' % (opts.ref_genome)
    ref_genome, res_ref = load_genome(opts.ref_genome, res=opts.res,
                                      verbose=opts.verbose)
    
    # load all chromosomes of target genome
    alt_genomes = {}
    for i, genome in enumerate(opts.genomes):
        if opts.verbose:
            print 'Load target chromosomes from %s' % (genome)
        alt_genomes[i], res_alt = load_genome(genome, res=opts.res,
                                              verbose=opts.verbose)

    ## create a new genome with TAD borders from each genome
    # first map alt_genome TAD borders to reference genome
    new_genome = {}
    for crm in alt_genome:
        for exp in alt_genome[crm].experiments:
            new_genome.setdefault(exp.name, {})
            for tad in exp.tads.values():
                coords = rest_query(crm.name, tad['start'], tad['end'],
                                    seed_spe='homo_sapiens',
                                    targ_spe='mus_musculus')
                if coords:
                    new_genome[exp.name].setdefault(coords['chr'],
                                                    {'end': [], 'score': []})
                    new_genome[exp.name][coords['chr']]['end'].append(coords['end'])
                    new_genome[exp.name][coords['chr']]['score'].append(tad['score'])

    ### HERE!!
    
    # then check if borders in reference genome are found in alt_genome
    # otherwise remove them
    reorder(new_genome)
    for exp in new_genome:
        for crm in new_genome[exp]:
            try:
                ref_genome[crm].add_experiment(
                    exp, res, tad_handler = new_genome[exp][crm])
            except KeyError:
                print ('Chromosome {} skipped, not in reference ' +
                       'genome').format(crm)

    system('mkdir -p ' + opts.out_path)
    for crm in ref_genome:
        system('mkdir -p ' + opts.out_path + '/' + crm)
        out_f = opts.out_path + '/' + crm + '/chr' + crm + '.tdb'
        ref_genome[crm].save_chromosome(out_f, force=True)


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
