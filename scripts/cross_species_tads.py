"""
22 Jan 2014

This script uses ensembl REST service to convert genomic coordinates of one
species to another species.
As borders from one chromosome in one species might be distribute along several
chromosome in an other species, this script needs to be used over whole genomes.

WARNING: Ensembl REST server is beta and only provides for latest genome builts
WARNING: needs internet :)
"""

import httplib2, sys, re, os
import json
from os                           import listdir
from os.path                      import isdir
from optparse                     import OptionParser
from pytadbit.utils.file_handling import check_pik
from pytadbit                     import load_chromosome, Chromosome
from time                         import sleep, time
from pytadbit.parsers.tad_parser  import parse_tads


GOOD_CRM = re.compile("^[A-Za-z]?[0-9]{0,3}[XVI]{0,3}(?:ito)?[A-Z-a-z]?$")


def decode_resolution(res):
    if 'Kb' in res:
        mult = 1000
    elif 'Mb' in res:
        mult = 1000000
    else:
        raise NotImplementedError('%s not know' % (res[-2:]))
    return int(res[:-2]) * mult


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
    return ref_genome


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


def get_best_map_coord(map_list):
    """
    returns coordinates of the longest mapped region
    """
    best_length = 0
    best_mapped = None
    for i, mapped in enumerate(map_list):
        length = mapped['mapped']['end'] - mapped['mapped']['start']
        if length > best_length:
            best_mapped = i
            best_length = length
    return {'chr'    : map_list[best_mapped]['mapped']['seq_region_name'],
            'start'  : map_list[best_mapped]['mapped']['start'],
            'end'    : map_list[best_mapped]['mapped']['end'],
            'strand' : map_list[best_mapped]['mapped']['strand']}


def syntenic_segment(crm, beg, end, from_species='homo_sapiens',
                     to_species='mus_musculus', alignment='LASTZ_NET',
                     server="http://beta.rest.ensembl.org", **kwargs):
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
    output = 'json'
    strand = 1
    ext = ("/alignment/block/region/%s/%s:%s-%s:%s?" +
           "species_set=%s&species_set=%s;method=%s")
    ext = ext % (from_species, crm, beg, end, strand, from_species,
                 to_species, alignment)
    resp, content = HTTP.request(
        server + ext, method="GET",
        headers={"Content-Type":"application/%s" % output})

    if content.startswith('Content-Type application/json had a problem'):
        return 'Region %s:%s-%s not found' % (crm, beg, end)
    
    if not resp.status == 200:
        try:
            jsonel = json.loads(content)['error']
            if re.findall('greater than ([0-9]+) for ',
                          jsonel):
                end = int(re.findall('greater than ([0-9]+) for ',
                                     jsonel)[0])
                return end
        except Exception, e:
            print str(e)
            print content
            print "Invalid response: ", resp.status
        raise Exception

    return get_best_alignment_coord(json.loads(content))


def remap_segment(crm, beg, end, species, from_map=None, to_map=None,
                  server="http://beta.rest.ensembl.org", **kwargs):
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
    output = 'json'
    strand = 1
    ext = ('/map/%s/%s/%s:%s..%s:%s/%s')
    ext = ext % (species, from_map, crm, beg, end, strand, to_map)
    resp, content = HTTP.request(
        server + ext, method="GET",
        headers={"Content-Type":"application/%s" % output})
    if content.startswith('Content-Type application/json had a problem'):
        return 'Region %s:%s-%s not found' % (crm, beg, end)
    
    if not resp.status == 200:
        try:
            jsonel = json.loads(content)['error']
            if re.findall('greater than ([0-9]+) for ',
                          jsonel):
                end = int(re.findall('greater than ([0-9]+) for ',
                                     jsonel)[0])
                return end
        except Exception, e:
            print str(e)
            print content
            print "Invalid response: ", resp.status
        raise Exception
    try:
        map_list = json.loads(content)['mappings']
    except KeyError:
        return
    if not map_list:
        return 'Region %s:%s-%s not found' % (crm, beg, end)
    return get_best_map_coord(map_list)


def map_tad(i, tad, crm, resolution, from_species, synteny=True, mapping=True,
            trace=None, **kwargs):
    """
    TODO: do synteny search right after mapping, in order to keep the trace
    """
    beg = int(tad['end']       * resolution)
    end = int((tad['end'] + 1) * resolution)
    ori_crm = crm
    ## keep trace
    trace.setdefault(crm, {})
    if not tad['end'] in trace[crm]:
        trace[crm][tad['end']] = {'from': {
            'crm'  : crm,
            'start': int(tad['end']       * resolution),
            'end'  : int((tad['end'] + 1) * resolution)}}

    coords = {}
    if mapping:
        while True:
            try:
                coords = remap_segment(crm, beg, end, from_species, **kwargs)
                if type(coords) is int:
                    if coords > beg:
                        beg = int(tad['end']       * resolution)
                        end = coords
                    else:
                        beg = int((tad['end'] - 1) * resolution)
                        end = coords
                else:
                    if not 'mapped to' in trace[crm][tad['end']]:
                        if type(coords) is dict:
                            trace[crm][tad['end']]['mapped to'] = coords
                        else:
                            trace[crm][tad['end']]['syntenic at'] = {
                                'chr': None, 'start':None, 'end': None}
                    break
            except Exception, e:
                print str(e)
                print '\n... reconnecting (mapping)...'
                print ' ' * ((i%50) + 9 + (i%50)/10),
                sleep(1)
                global HTTP
                HTTP = httplib2.Http(".cache")

    if synteny and type(coords) is dict:
        crm, beg, end = coords['chr'], coords['start'], coords['end']
        while True:
            try:
                coords = syntenic_segment(crm, beg, end, from_species, **kwargs)
                if type(coords) is int:
                    if coords > beg:
                        beg = int(tad['end']       * resolution)
                        end = coords
                    else:
                        beg = int((tad['end'] - 1) * resolution)
                        end = coords
                else:
                    if not 'syntenic at' in trace[ori_crm][tad['end']]:
                        if type(coords) is dict:
                            trace[ori_crm][tad['end']]['syntenic at'] = coords
                        else:
                            trace[ori_crm][tad['end']]['syntenic at'] = {
                                'chr': None, 'start':None, 'end': None}
                    break
            except Exception, e:
                print str(e)
                print '\n... reconnecting (synteny)...'
                print ' ' * ((i%50) + 9 + (i%50)/10),
                sleep(1)
                global HTTP
                HTTP = httplib2.Http(".cache")
    return coords


def convert_chromosome(crm, new_genome, from_species, synteny=True,
                       mapping=True, trace=None, **kwargs):
    new_crm = Chromosome(crm.name, species=crm.species, assembly=crm.assembly)
    log = []
    crm_name = crm.name.replace('chr', '').replace('crm', '')
    for exp in crm.experiments:
        print '      Experiment %10s (%s TADs)' % (exp.name, len(exp.tads))
        sys.stdout.write('         ')
        connections = 0 # variable to avoid reaching the limit of 6 connection
                        # per second allowed by Ensembl.
        t0 = time()
        for i, tad in enumerate(exp.tads.values()):
            if not tad['end'] in trace[crm]:
                # MAP
                coords = map_tad(i, tad, crm_name, exp.resolution,
                                 from_species, synteny=synteny, mapping=mapping,
                                 trace=trace, **kwargs)
                connections += synteny + mapping
                if connections >= 4 :
                    to_sleep = .8 - min(0.8, time() - t0)
                    sleep(to_sleep)
                    connections = 0
                    t0 = time()
            else:
                coords = trace[crm][tad['end']][
                    'syntenic at' if synteny else 'mapped to']
            if type(coords) is dict and GOOD_CRM.match(coords['chr']):
                new_genome.setdefault(
                    coords['chr'], Chromosome(coords['chr'],
                                              species=crm.species,
                                              assembly=crm.assembly))
                if not exp.name in [e.name for e in
                                    new_genome[coords['chr']].experiments]:
                    new_genome[coords['chr']].add_experiment(
                        exp.name, cell_type=exp.cell_type,
                        identifier=exp.identifier, enzyme=exp.enzyme,
                        resolution=exp.resolution, no_warn=True)
                    new_genome[coords['chr']].experiments[exp.name]._tads = {
                        'end': [], 'start': [], 'score': []}
                tads = new_genome[coords['chr']].experiments[exp.name]._tads
                tads['start'].append(
                    (tads['end'][-1] + 1) if tads['end'] else 0.0)
                tads['end'  ].append(float(round(
                    float(coords['end']) / exp.resolution)))
                tads['score'].append(tad['score'])
                sys.stdout.write('.')
            else:
                if type(coords) is dict:
                    coords = 'Region %s:%s-%s not nice' % (
                        coords['chr'], coords['start'], coords['end'])
                if not coords in log:
                    log.append(coords)
                sys.stdout.write('E')
            if not (i + 1) % 50:
                sys.stdout.write (' %3s\n         ' % (i + 1))
            elif not (i+1) % 10:
                sys.stdout.write(' ')
                sys.stdout.flush()
            else:
                sys.stdout.flush()
        print ''
    if log:
        log.sort(key=lambda x: int(x.split(':')[1].split('-')[0]))
        sys.stdout.write('              '+'\n              '.join(log))
        sys.stdout.write('\n              ')
    log = []
    return new_crm


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
            ' '*10 + opts.genome)
    genome = load_genome(opts.genome, res=opts.res, verbose=opts.verbose
)
    # remap TADs
    if opts.verbose:
        print '\nCreating new TADbit chromosomes with new coordinates\n'
    new_genome = {}
    trace = {}
    global HTTP
    HTTP = httplib2.Http(".cache")

    # create directory for output
    rootdir = os.path.abspath(opts.out_path)
    if not os.path.exists(rootdir):
        os.mkdir(rootdir)
    
    # start iteration over chromosomes
    for crm in genome:
        print '\n   Chromosome:', crm
        new_crm = None
        # remap and syntenic region
        if opts.original_assembly:
            print '\n     remapping from %s to %s:' % (
                opts.original_assembly, opts.target_assembly)
            print '        and searching syntenic regions from %s to %s' %(
                opts.original_species.capitalize().replace('_', ' '),
                opts.target_species.capitalize().replace('_', ' '))
            convert_chromosome(genome[crm], new_genome, opts.original_species,
                               from_map=opts.original_assembly,
                               to_map=opts.target_assembly,
                               to_species=opts.target_species,
                               synteny=True, mapping=True, trace=trace)
        if opts.log:
            log = open
        for t in sorted(trace[crm]):
            try:
                print '%4s : %2s:%9s-%9s -> %2s:%9s-%9s -> %2s:%9s-%9s' %(
                    int(t),
                    trace[crm][t]['from']['crm']       , trace[crm][t]['from']['start']       , trace[crm][t]['from']['end'],
                    trace[crm][t]['mapped to']['chr']  , trace[crm][t]['mapped to']['start']  , trace[crm][t]['mapped to']['end'],
                    trace[crm][t]['syntenic at']['chr'], trace[crm][t]['syntenic at']['start'], trace[crm][t]['syntenic at']['end'],
                    )
            except KeyError:
                print '%4s : %2s:%9s-%9s -> %22s -> %22s' %(
                    int(t),
                    trace[crm][t]['from']['crm']       , trace[crm][t]['from']['start']       , trace[crm][t]['from']['end'],
                    'None',
                    'None',
                    )

    # save new chromosomes
    if opts.check:
        for crm in genome:
            new_crm = genome[crm]
            for exp in new_crm.experiments:
                tadcnt = 0
                new_tads = {}
                for tad in exp.tads:
                    cond = 'syntenic at' if opts.target_species else 'mapped to'
                    if trace[crm][exp.tads[tad]['end']][cond]['chr'] is None:
                        continue
                    new_tads[tadcnt] = exp.tads[tad]
                    tadcnt += 1
                exp.tads = new_tads
            crmdir = os.path.join(rootdir, crm)
            if not os.path.exists(crmdir):
                os.mkdir(crmdir)
            new_crm.save_chromosome(os.path.join(crmdir, crm + '.tdb'))
    else:
        for crm in new_genome:
            new_crm = new_genome[crm]
            for exp in new_crm.experiments:
                tads, norm = parse_tads(exp._tads)
                last = max(tads.keys())
                if not exp.size:
                    exp.size = tads[last]['end']
                exp.tads = tads
                exp.norm  = norm
            crmdir = os.path.join(rootdir, crm)
            if not os.path.exists(crmdir):
                os.mkdir(crmdir)
            new_crm.save_chromosome(os.path.join(crmdir, crm + '.tdb'))


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
