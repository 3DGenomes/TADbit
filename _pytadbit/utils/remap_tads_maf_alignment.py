"""
12 Jan 2016

This script uses MAF alignment files to convert genomic coordinates of one
species to another species.
As borders from one chromosome in one species might be distribute along several
chromosome in an other species, this script needs to be used over whole genomes.
"""

import sys, re, os
from os                           import listdir
from os.path                      import isdir, isfile
from pytadbit.utils.file_handling import check_pik
from pytadbit                     import load_chromosome, Chromosome
from time                         import sleep, time
from pytadbit.parsers.tad_parser  import parse_tads
from cPickle import dump, load

GOOD_CRM = re.compile("^(?:chr)?[A-Za-z]?[0-9]{0,3}[XVI]{0,3}(?:ito)?[A-Z-a-z]?$")


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


def load_genome_from_tad_def(genome_path, res, verbose=False):
    """
    Search, at a given path, for chromosome folders containing TAD
    definitions in tsv files.

    :param genome_path: Path where to search for TADbit chromosomes
    :param res: Resolution at were saved chromosomes
    :param False verbose:

    :returns: a dictionary with all TADbit chromosomes found
    """
    ref_genome = {}
    for crm in listdir(genome_path):
        crm_path = os.path.join(genome_path, crm)
        if not isfile(crm_path):
            continue
        if crm in ref_genome:
            raise Exception('More than 1 TAD definition file found\n')
        crm = crm.replace('.tsv', '').replace('chr', '')
        if verbose:
            print '  Chromosome:', crm
        crmO = Chromosome(crm)
        crmO.add_experiment('sample', res)
        crmO.experiments[0].load_tad_def(crm_path)
        ref_genome[crm] = crmO
    return ref_genome


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
        crm_path = os.path.join(genome_path, crm)
        if not isdir(crm_path):
            continue
        for crm_fh in listdir(crm_path):
            crm_pik = os.path.join(crm_path, crm_fh)
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

def get_synteny_defitinion(fname):
    """
    From a tab separated file with this columns:
     1- Reference chromosome name
     2- Reference start position
     3- Reference end position
     4- Target chromosome name
     5- Target start position
     6- Target end position
     7- Target orientation

    :param fname: path to file with columns

    :returns: dictionary with synthenic block coordinates
    """
    synteny = {}
    for line in open(fname):
        if line.startswith('#'):
            continue
        at_crm, at_beg, at_end, to_crm, to_beg, to_end, to_std = line.split()
        at_beg, at_end, to_beg, to_end, to_std = map(
            int, [at_beg, at_end, to_beg, to_end, to_std])
        synteny.setdefault(at_crm, []).append({'from': {'crm': at_crm,
                                                        'beg': at_beg,
                                                        'end': at_end},
                                               'targ': {'crm': to_crm,
                                                        'beg': to_beg,
                                                        'end': to_end,
                                                        'std': to_std}})
    return synteny


def filter_non_syntenic_alignments(synteny_file):
    synteny = get_synteny_defitinion(synteny_file)
    goods = set()
    for at_crm in synteny:
        for blk in synteny[at_crm]:
            to_crm = blk['targ']['crm']
            # to_std = blk['targ']['std']
            to_pos = set(range(blk['targ']['beg']/100, blk['targ']['end']/100 + 1))
            at_pos = set(range(blk['from']['beg']/100, blk['from']['end']/100 + 1))
            for num, ali in enumerate(ALIGNMENTS):
                if ali['from']['crm'] != at_crm:
                    continue
                if not at_pos.intersection(set(range(ali['from']['beg']/100,
                                                     ali['from']['end']/100 + 1))):
                    continue
                # print 'ok: %3s %11s %11s [%5s] %2s %3s %11s %11s' % (
                #     ali['targ']['crm'], ali['targ']['beg'], ali['targ']['end'],
                #     ali['targ']['len'], ali['targ']['std'],
                #     ali['from']['crm'], ali['from']['beg'], ali['from']['end'])
                if ali['targ']['crm'] != to_crm:
                    continue
                # if ali['targ']['crm'] != to_crm or ali['targ']['std'] != to_std:
                #     print " " * 39, "'-> STRAND!!!"
                #     continue
                if not to_pos.intersection(set(range(ali['targ']['beg']/100,
                                                     ali['targ']['end']/100 + 1))):
                    print 'HA! %3s %11s %11s [%5s] %2s %3s %11s %11s' % (
                        blk['targ']['crm'], blk['targ']['beg'], blk['targ']['end'],
                        blk['targ']['end'] - blk['targ']['beg'], blk['targ']['std'],
                        blk['from']['crm'], blk['from']['beg'], blk['from']['end'])
                    print '=> POSITION??'
                    continue
                goods.add(num)
    for num in xrange(len(ALIGNMENTS) - 1, -1, -1):
        if not num in goods:
            ALIGNMENTS.pop(num)
                

def get_alignments(seed, targ, maf_file):
    """
    From a MAF file extract alignment of a pair of species. Tris to extend
    the alignment taking advantage of missing species.
    """
    species = [seed, targ]
    # get alignments
    pre_alignments = []
    pick_path = maf_file.replace('.maf', '%s-%s.pick' % (seed, targ))
    if isfile(pick_path):
        print '       Found pre-computed alignment at %s' % pick_path
        print '         -> loading it...'
        return load(open(pick_path))
    maf_handler = open(maf_file)
    while True:
        try:
            found = False
            while True:
                line = maf_handler.next()
                if not line.startswith('s '):
                    if line.startswith('a '):
                        ali = {}
                        break
                    if line.startswith('i '):
                        if not spe in species:
                            continue
                        if '_' in crm or not GOOD_CRM.match(crm):
                            continue
                        if spe in ali:
                            raise Exception()
                        found = True
                        ali[spe] = (crm, pos, slen, strand, seq,
                                    ' '.join(line.split()[-4:]))
                    continue
                _, spe_crm, pos, slen, strand, clen, seq = line.split()
                try:
                    spe, crm = spe_crm.split('.chr', 1)
                except ValueError:
                    spe, crm = spe_crm.split('.', 1)
                slen = int(slen)
                if strand == '+':
                    pos  = int(pos)
                    strand = 1
                else:
                    pos = int(clen) - int(pos)
                    strand = -1
                if spe == seed:
                    ali[spe] = crm, pos, slen, strand, seq, 'C 0 C 0'
            if found:
                pre_alignments.append(ali)
        except StopIteration:
            break

    # reduce alignments to our species, and merge chunks if possible
    alignments = []
    dico1 = {}
    for k, dico2 in enumerate(pre_alignments[1:]):
        if not dico1:
            dico1 = pre_alignments[k].copy()
        species1 = sorted(dico1.keys())
        species2 = sorted(dico2.keys())
        if species1 == species2:
            vals2  = [dico2[spe] for spe in species2]
            if all([val[-1].startswith('C 0') for val in vals2]):
                # sumup sequences and lengths
                for spe in dico1:
                    crm, pos, slen1, strand, seq1, _    = dico1[spe]
                    _  , _  , slen2, _     , seq2, comp = dico2[spe]
                    dico1[spe] = crm, pos, slen1 + slen2, strand, seq1 + seq2, comp
                continue
        alignments.append(dico1)
        dico1 = dico2
    alignments.append(dico1)

    bads = []
    for i in xrange(len(alignments)):
        try:
            alignments[i] = {
                'from': {'crm': alignments[i][seed][0],
                         'beg': alignments[i][seed][1],
                         'end': alignments[i][seed][1] + alignments[i][seed][2],
                         'len': alignments[i][seed][2],
                         'std': alignments[i][seed][3],
                         'seq': alignments[i][seed][4]},
                'targ': {'crm': alignments[i][targ][0],
                         'beg': alignments[i][targ][1],
                         'end': alignments[i][targ][1] + alignments[i][targ][2],
                         'len': alignments[i][targ][2],
                         'std': alignments[i][targ][3],
                         'seq': alignments[i][targ][4]}}
        except KeyError:
            bads.append(i)
    for i in bads[::-1]:
        alignments.pop(i)
    print '      Filtering alignment'
    # remove alignments that are not inside syntenic regions
    filter_non_syntenic_alignments(synteny_file)



    print 'Dumping alignment to file', pick_path
    out = open(pick_path, 'w')
    dump(alignments, out)
    out.close()
    return alignments


def syntenic_segment(crm, beg, end):
    ins = set(range(beg / 10, end / 10 + 1))
    matching = [(0, {})]
    for ali in ALIGNMENTS:
        if ali['from']['crm'] != crm:
            continue
        val = len(ins.intersection(range(ali['from']['beg'] / 10,
                                         ali['from']['end'] / 10 + 1)))
        if val > 1: # at least an overlap of 10 nucleotides
            matching.append((val, ali))
    return max(matching)[1]


def map_tad(tad, crm, resolution, trace=None):
    """
    Converts coordinates of 1 TAD border in a given chromosome. Convertion in
    terms of change in assembly version, or in terms of syntenic regions between
    species (or both).
    """
    beg = int(tad['end']       * resolution)
    end = int((tad['end'] + 1) * resolution)
    ## keep trace
    try:
        trace.setdefault(crm, {})
    except AttributeError:
        trace = {crm: {}}
    if not tad['end'] in trace[crm]:
        trace[crm][tad['end']] = {'from': {'crm': crm, 'beg': beg, 'end': end}}
    coords = {}
    try:
        coords = syntenic_segment(crm, beg, end)
        if not coords: # found nothing
            return coords
        diff_beg = beg - coords['from']['beg']
        diff_end = end - coords['from']['end']
        if coords['targ']['std'] == 1:
            coords['targ']['beg'] += coords['targ']['beg'] + diff_beg
            coords['targ']['end'] += diff_end
        else:
            coords['targ']['beg'] -= diff_end # use diff end it's all reverted
            coords['targ']['end'] -= diff_beg
    except Exception, e:
        print e
        raise Exception('ERROR: not able to find synteny %s:%s-%s\n' % (crm, beg, end))
    trace[crm][tad['end']]['syntenic at'] = coords['targ']
    return coords['targ']


def convert_chromosome(crm, new_genome, trace=None):
    """
    Converts coordinates of TAD borders in a given chromosome. Convertion in
    terms of change in assembly version, or in terms of syntenic regions between
    species (or both).
    """
    log = []
    crm_name = crm.name.replace('chr', '').replace('crm', '')
    for exp in crm.experiments:
        print '      Experiment %10s (%s TADs)' % (exp.name, len(exp.tads))
        sys.stdout.write('         ')
        for i, tad in enumerate(exp.tads.values()):
            trace.setdefault(crm, {})
            if not tad['end'] in trace[crm]:
                # MAP
                coords = map_tad(tad, crm_name, exp.resolution, trace=trace)
            else:
                print 'NEVER happening...'
                coords = trace[crm.name][tad['end']]
            if coords and GOOD_CRM.match(coords['crm']):
                # add new chromosome in the target genome, if not there
                new_genome.setdefault(
                    coords['crm'], Chromosome(coords['crm'],
                                              species=crm.species,
                                              assembly=crm.assembly))
                # add new experiment with empty dict of TADs, if not there
                if not exp.name in [e.name for e in
                                    new_genome[coords['crm']].experiments]:
                    new_genome[coords['crm']].add_experiment(
                        exp.name, cell_type=exp.cell_type,
                        identifier=exp.identifier, enzyme=exp.enzyme,
                        resolution=exp.resolution, no_warn=True)
                    new_genome[coords['crm']].experiments[exp.name]._tads = {
                        'end': [], 'start': [], 'score': []}
                tads = new_genome[coords['crm']].experiments[exp.name]._tads
                tads['start'].append(
                    (tads['end'][-1] + 1) if tads['end'] else 0.0)
                tads['end'  ].append(float(round(
                    float(coords['end']) / exp.resolution)))
                tads['score'].append(tad['score'])
                sys.stdout.write('.')
            else:
                if 'crm' in coords:
                    coords = 'Target region %s:%s-%s not nice' % (
                        coords['crm'], coords['start'], coords['end'])
                else:
                    coords = 'Seed region %s:%s-%s not found ' % (
                        crm_name, tad['end'] * exp.resolution,
                        (tad['end'] + 1) * exp.resolution)
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
        log.sort(key=lambda x: int(float(x.split(':')[1].split('-')[0])))
        sys.stdout.write('              '+'\n              '.join(log))
        sys.stdout.write('\n              ')


def remap_genome(seed_species, target_species, maf_file, genome, synteny_file=None,
                 write_log=False):
    """
    Given a set of chromosomes from one species with a given assembly version.
    Remap the coordinates of TAD borders to a new assembly and/or to a new
    species.

    :para genome: a dict containing all chromosomes computed of the target species
    :param original_species: i.e.: 'Homo_sapiens'
    :param target_species: i.e.: 'mus_musculus', if None, no search for syntenic
       regions will be done
    :param False write_log: path where to write a log file

    :returns: new genome dictionary with remapped TAD borders, and a trace
       dictionary that allows to reconstruct the process of remapping and finding
       syntenic regions.
    """
    trace = {}
    if write_log:
        log = open(write_log, 'w')
    else:
        log = sys.stdout

    new_genome = {}
    for crm in genome:
        print '\n   Chromosome:', crm
        print '      Loading MAF alignment'
        global ALIGNMENTS
        ALIGNMENTS = get_alignments(seed_species, target_species, maf_file % crm)
        # remap and syntenic region
        print '           Searching syntenic regions from %s to %s' %(
            seed_species.capitalize().replace('_', ' '),
            target_species.capitalize().replace('_', ' '))
        convert_chromosome(genome[crm], new_genome, trace=trace)

        crm = genome[crm]
        new_genome[crm.name] = Chromosome(crm.name, species=crm.species,
                                          assembly=crm.assembly)
        for t in sorted(trace[crm]):
            try:
                log.write('%4s : %2s:%9s-%9s -> %2s:%9s-%9s\n' %(
                    int(t)                             ,
                    trace[crm][t]['from']['crm']       ,
                    trace[crm][t]['from']['beg']       ,
                    trace[crm][t]['from']['end']       ,
                    trace[crm][t]['syntenic at']['crm'],
                    trace[crm][t]['syntenic at']['beg'],
                    trace[crm][t]['syntenic at']['end'],
                    ))
            except KeyError:
                log.write('%4s : %2s:%9s-%9s -> %22s -> %22s\n' %(
                    int(t),
                    trace[crm][t]['from']['crm'],
                    trace[crm][t]['from']['beg'],
                    trace[crm][t]['from']['end'],
                    'None',
                    'None',
                    ))
    if write_log:
        log.close()
        
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
                    try:
                        if trace[crm][exp.tads[tad]['end']][cond]['crm'] is None:
                            continue
                    except KeyError:
                        print ('Not found:', crm, exp.tads[tad]['end'],
                               trace[crm][exp.tads[tad]['end']])
                        continue
                    new_tads[tadcnt] = exp.tads[tad]
                    tadcnt += 1
                exp.tads = new_tads
            else:
                tads, norm = parse_tads(exp._tads)
                last = max(tads.keys())
                if not exp.size:
                    exp.size = tads[last]['end']
                exp.norm  = norm
                exp.tads = tads
        crmdir = os.path.join(rootdir, crm)
        if not os.path.exists(crmdir):
            os.mkdir(crmdir)
        new_crm.save_chromosome(os.path.join(crmdir, crm + '.tdb'))
    

