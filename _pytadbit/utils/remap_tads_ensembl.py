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
from pytadbit.utils.file_handling import check_pik
from pytadbit                     import load_chromosome, Chromosome
from time                         import sleep, time
from pytadbit.parsers.tad_parser  import parse_tads


GOOD_CRM = re.compile("^[A-Za-z]?[0-9]{0,3}[XVI]{0,3}(?:ito)?[A-Z-a-z]?$")


def decode_resolution(res):
    """
    Convertion from human to machine readable resolution
    """
    if 'Kb' in res:
        mult = 1000
    elif 'Mb' in res:
        mult = 1000000
    else:
        raise NotImplementedError('%s not know' % (res[-2:]))
    return int(res[:-2]) * mult


def load_genome(genome_path, res=None, verbose=False):
    """
    Search, at a given path, for chromosome folders containing TADbit chromosome
    objects saved as files.

    :param genome_path: Path where to search for TADbit chomosomes
    :param None res: Resolution at which should be saved chromosomes
    :param False verbose:

    :returns: a dictionary with all TADbit chromosomes found
    """
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
    Given a seed genomic segment find the coordinate of its syntenic segment on
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
    Given a seed genomic segment find its coordinates in a different assembly
    version.
    
    :param crm: chromosome name
    :param beg: region start position
    :param end: region end position
    :param species: name of the seed species (from which com the original
       coordinates). Scientific name with underscore instead of spaces
    :param None from_map: name of the original assembly version
    :param None to_map: name of the target assembly version
    :param 'http://beta.rest.ensembl.org' server: see http://beta.rest.ensembl.org/documentation/info/genomic_alignment_block_region
       for more explanations

    :returns: a dict with the coordinates of the remapped segment found (chr,
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
    Converts coordinates of 1 TAD border in a given chromosome. Convertion in
    terms of change in assembly version, or in terms of syntenic regions between
    species (or both).
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
    global HTTP
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
                HTTP = httplib2.Http(".cache")
    return coords


def convert_chromosome(crm, new_genome, from_species, synteny=True,
                       mapping=True, trace=None, **kwargs):
    """
    Converts coordinates of TAD borders in a given chromosome. Convertion in
    terms of change in assembly version, or in terms of syntenic regions between
    species (or both).
    """
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


def remap_genome(genome, original_assembly, target_assembly,
                 original_species, target_species, write_log=False):
    """
    Given a set of chromosomes from one species with a given assembly version.
    Remap the coordinates of TAD borders to a new assembly and/or to a new
    species.

    :para genome: a dict containing all chromosomes computed
    :param original_assembly: i.e.: 'NCBI36'
    :param target_assembly:  i.e.: 'GRCh37', if None, no remapping will be done
    :param original_species: i.e.: 'Homo_sapiens'
    :param target_species: i.e.: 'mus_musculus', if None, no search for syntenic
       regions will be done
    :param False write_log: path where to write a log file

    :returns: new genome dictionary with remapped TAD borders, and a trace
       dictionnary that allows to reconstruct the process of remapping and finding
       syntenic regions.
    """
    global HTTP
    HTTP = httplib2.Http(".cache")
    trace = {}

    new_genome = {}
    for crm in genome:
        print '\n   Chromosome:', crm
        # remap and syntenic region
        if original_assembly:
            print '\n     remapping from %s to %s:' % (
                original_assembly, target_assembly)
            print '        and searching syntenic regions from %s to %s' %(
                original_species.capitalize().replace('_', ' '),
                target_species.capitalize().replace('_', ' '))
            convert_chromosome(genome[crm], new_genome, original_species,
                               from_map=original_assembly,
                               to_map=target_assembly,
                               to_species=target_species,
                               synteny=True, mapping=True, trace=trace)
        if write_log:
            log = open(write_log, 'w')
        for t in sorted(trace[crm]):
            try:
                log.write('%4s : %2s:%9s-%9s -> %2s:%9s-%9s -> %2s:%9s-%9s\n' %(
                    int(t),
                    trace[crm][t]['from']['crm']       , trace[crm][t]['from']['start']       , trace[crm][t]['from']['end'],
                    trace[crm][t]['mapped to']['chr']  , trace[crm][t]['mapped to']['start']  , trace[crm][t]['mapped to']['end'],
                    trace[crm][t]['syntenic at']['chr'], trace[crm][t]['syntenic at']['start'], trace[crm][t]['syntenic at']['end'],
                    ))
            except KeyError:
                log.write('%4s : %2s:%9s-%9s -> %22s -> %22s\n' %(
                    int(t),
                    trace[crm][t]['from']['crm']       , trace[crm][t]['from']['start']       , trace[crm][t]['from']['end'],
                    'None',
                    'None',
                    ))
    return new_genome, trace


def save_new_genome(genome, trace, check=False, target_species=None, rootdir='./'):
    """
    Save new chromosomes with remapped or check TAD borders into a new folder.

    :para genome: a dict containing all chromosomes computed
    :param trace: dictionary containing a trace of all mapped TAD boundaries
    :param False check: if no remapping have to be done (only check TAD borders)
    :param None target_species: name of the target species, if None, it is
       assumed, that only a remapping has been done.
    :param './' rootdir: path where to write directories for remapped/checked
       chromosomes.
    """
    for crm in genome:
        new_crm = genome[crm]
        for exp in new_crm.experiments:
            if check:
                tadcnt = 0
                new_tads = {}
                for tad in exp.tads:
                    cond = 'syntenic at' if target_species else 'mapped to'
                    if trace[crm][exp.tads[tad]['end']][cond]['chr'] is None:
                        continue
                    new_tads[tadcnt] = exp.tads[tad]
                    tadcnt += 1
            else:
                tads, norm = parse_tads(exp._tads)
                last = max(tads.keys())
                if not exp.size:
                    exp.size = tads[last]['end']
                exp.norm  = norm
            exp.tads = new_tads
        crmdir = os.path.join(rootdir, crm)
        if not os.path.exists(crmdir):
            os.mkdir(crmdir)
        new_crm.save_chromosome(os.path.join(crmdir, crm + '.tdb'))
    
