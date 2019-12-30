"""
12 Jan 2016

This script uses MAF alignment files to convert genomic coordinates of one
species to another species.
As borders from one chromosome in one species might be distribute along several
chromosome in an other species, this script needs to be used over whole genomes.
"""
from __future__ import print_function

import sys, re, os
from os                           import listdir
from os.path                      import isfile
from pytadbit                     import Chromosome
from pytadbit.parsers.tad_parser  import parse_tads
from cPickle import dump, load

GOOD_CRM = re.compile("^(?:chr)?[A-Za-z]?[0-9]{0,3}[XVI]{0,3}(?:ito)?[A-Z-a-z]?$")


def decode_resolution(res):
    """
    Conversion from human to machine readable resolution
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
        crm = crm.replace('.tsv', '').replace('chr', '').upper()
        if verbose:
            print('  Chromosome:', crm)
        crmO = Chromosome(crm)
        crmO.add_experiment('sample', res)
        crmO.experiments[0].load_tad_def(crm_path)
        ref_genome[crm] = crmO
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
        to_crm, to_beg, to_end, at_crm, at_beg, at_end, at_std = line.split()
        to_beg, to_end, at_beg, at_end, at_std = map(
            int, [to_beg, to_end, at_beg, at_end, at_std])
        synteny.setdefault(at_crm.upper(), []).append(
            {'targ': {'crm': to_crm.upper(),
                      'beg': to_beg,
                      'end': to_end},
             'from': {'crm': at_crm.upper(),
                      'beg': at_beg,
                      'end': at_end,
                      'std': at_std}})
    return synteny


def filter_non_syntenic_alignments(synteny_file, synteny_reso=0):
    synteny = get_synteny_defitinion(synteny_file)
    goods = set()
    for at_crm in synteny:
        for blk in synteny[at_crm]:
            to_crm = blk['targ']['crm']
            # to_std = blk['targ']['std']
            at_pos = set(range((blk['from']['beg'] - synteny_reso) / 100,
                               (blk['from']['end'] + synteny_reso) / 100 + 1))
            to_pos = set(range((blk['targ']['beg'] - synteny_reso) / 100,
                               (blk['targ']['end'] + synteny_reso) / 100 + 1))
            for num, ali in enumerate(ALIGNMENTS):
                if ali['targ']['crm'] != to_crm:
                    continue
                if not to_pos.intersection(set(range(ali['targ']['beg']/100,
                                                     ali['targ']['end']/100 + 1))):
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
                # check that the region in close from definition of any syntenic region
                if not at_pos.intersection(set(range(ali['from']['beg']/100,
                                                     ali['from']['end']/100 + 1))):
                    print('Found target alignment outside any syntenic block', end=' ')
                    print(' (while seed alignment is inside)')
                    print('SYNT: %3s %10s %10s %3s %10s %10s [%9s] %2s' % (
                        blk['targ']['crm'],  blk['targ']['beg'], blk['targ']['end'],
                        blk['from']['crm'],  blk['from']['beg'], blk['from']['end'],
                        blk['from']['end'] - blk['from']['beg'], blk['from']['std']))
                    print('ALI:  %3s %10s %10s %3s %10s %10s [%9s] %2s' % (
                        ali['targ']['crm'],  ali['targ']['beg'], ali['targ']['end'],
                        ali['from']['crm'],  ali['from']['beg'], ali['from']['end'],
                        ali['from']['end'] - ali['from']['beg'], ali['from']['std']))
                    continue
                goods.add(num)
    for num in xrange(len(ALIGNMENTS) - 1, -1, -1):
        if not num in goods:
            ALIGNMENTS.pop(num)


def get_alignments(seed, targ, maf_path, synteny_file, synteny_reso=0,
                   clean_all=False, breaker='i'):
    """
    From a MAF file extract alignment of a pair of species. Also extends
    the alignment taking advantage of missing species.
    :param 'i' breaker: first character of the lines thata marks the end of the
       sequence. Use '\n' for pairwise alignments
    """
    # get alignments
    pre_alignments = []
    pick_path = os.path.join(maf_path, 'alignments_%s-%s.pick' % (seed, targ))
    if isfile(pick_path) and not clean_all:
        print('       Found pre-computed alignment at %s' % pick_path)
        print('         -> loading it...')
        global ALIGNMENTS
        ALIGNMENTS = load(open(pick_path))
        return
    crm_num = 1
    for crm_file in os.listdir(maf_path):
        if (not os.path.isfile(os.path.join(maf_path, crm_file))
            or not crm_file.endswith('.maf')):
            continue
        print('     %2s- loading MAF file' % crm_num, crm_file)
        maf_handler = open(os.path.join(maf_path, crm_file))
        for line in maf_handler:
            if line.startswith('s'):
                _, spe_crm, pos, slen, strand, clen, _ = line.split()
                try:
                    spe, crm = spe_crm.split('.chr', 1)
                except ValueError:
                    spe, crm = spe_crm.split('.', 1)
                crm = crm.upper()
                slen = int(slen)
                if spe == seed:   # seed always on positive strand
                    ali = {}
                    ali[spe] = crm, int(pos), slen, 1, True
                elif spe != targ: # skip other species and their "i " line
                    _ = next(maf_handler)
                    continue
                elif strand == '+':
                    pos = int(pos)
                    strand = 1
                else:
                    pos = int(clen) - int(pos)
                    strand = -1
            elif line.startswith(breaker):
                if not GOOD_CRM.match(crm):
                    continue
                ali[spe] = (crm, pos, slen, strand, line.endswith('C 0 C 0\n'))
                # store this species and seed
                pre_alignments.append(ali)
                # skip other species
                for line in maf_handler:
                    if line.startswith('a'):
                        break
        crm_num += 1

    # reduce alignments to our species, and merge chunks if possible
    global ALIGNMENTS
    ALIGNMENTS = []
    dico1 = {}
    for k, dico2 in enumerate(pre_alignments[1:]):
        if not dico1:
            dico1 = pre_alignments[k].copy()
        species1 = sorted(dico1.keys())
        species2 = sorted(dico2.keys())
        if species1 == species2:
            vals2  = [dico2[spe] for spe in species2]
            if all([val[-1] for val in vals2]):
                # sumup sequences and lengths
                for spe in dico1:
                    # crm, pos, slen1, strand, seq1, _    = dico1[spe]
                    # _  , _  , slen2, _     , seq2, comp = dico2[spe]
                    # dico1[spe] = crm, pos, slen1 + slen2, strand, seq1 + seq2, comp
                    crm, pos, slen1, strand, _    = dico1[spe]
                    _  , _  , slen2, _     , comp = dico2[spe]
                    dico1[spe] = crm, pos, slen1 + slen2, strand, comp
                continue
        ALIGNMENTS.append(dico1)
        dico1 = dico2
    ALIGNMENTS.append(dico1)

    bads = []
    for i in xrange(len(ALIGNMENTS)):
        try:
            ALIGNMENTS[i] = {
                'targ': {'crm': ALIGNMENTS[i][seed][0],
                         'beg': ALIGNMENTS[i][seed][1],
                         'end': ALIGNMENTS[i][seed][1] + ALIGNMENTS[i][seed][2],
                         'len': ALIGNMENTS[i][seed][2],
                         'std': ALIGNMENTS[i][seed][3],
                         #'seq': ALIGNMENTS[i][seed][4]
                         },
                'from': {'crm': ALIGNMENTS[i][targ][0],
                         'beg': ALIGNMENTS[i][targ][1],
                         'end': ALIGNMENTS[i][targ][1] + ALIGNMENTS[i][targ][2],
                         'len': ALIGNMENTS[i][targ][2],
                         'std': ALIGNMENTS[i][targ][3],
                         #'seq': ALIGNMENTS[i][targ][4]
                         }}
        except KeyError:
            bads.append(i)
    for i in bads[::-1]:
        ALIGNMENTS.pop(i)
    print('      Filtering alignment')
    # remove alignments that are not inside syntenic regions
    if synteny_file:
        filter_non_syntenic_alignments(synteny_file, synteny_reso)

    print('Dumping alignment to file', pick_path)
    out = open(pick_path, 'w')
    dump(ALIGNMENTS, out)
    out.close()


def syntenic_segment(crm, beg, end):
    matching = [(0, {})]
    for ali in ALIGNMENTS:
        if ali['from']['crm'] != crm:
            continue
        if beg > ali['from']['end'] or end < ali['from']['beg']:
            continue
        val = min([end, ali['from']['end']]) - max([beg, ali['from']['beg']])
        if val > 10: # at least an overlap of 10 nucleotide
            matching.append((val, ali))
            if val == end - beg:
                break
    return max(matching)[1] # take alignment only, not val


def map_bin(crm, pos, resolution, trace=None):
    """
    Converts coordinates of 1 TAD border in a given chromosome. Convertion in
    terms of change in assembly version, or in terms of syntenic regions between
    species (or both).
    """
    beg = int(pos       * resolution)
    end = int((pos + 1) * resolution)
    ## keep trace
    try:
        trace.setdefault(crm, {})
    except AttributeError:
        trace = {crm: {}}
    if not end in trace[crm]:
        trace[crm][end] = {'from': {'crm': crm, 'beg': beg, 'end': end}}
    coords = {}
    try:
        coords = syntenic_segment(crm, beg, end)
        if not coords: # found nothing
            return coords
        diff_beg = beg - coords['from']['beg']
        diff_end = end - coords['from']['end']
        if coords['targ']['std'] == 1:
            coords['targ']['beg'] += diff_beg
            coords['targ']['end'] += diff_end
        else:
            coords['targ']['beg'] -= diff_end # use diff end it's all reverted
            coords['targ']['end'] -= diff_beg
    except Exception as e:
        print(e)
        raise Exception('ERROR: not able to find synteny %s:%s-%s\n' % (crm, beg, end))
    trace[crm][end]['syntenic at'] = coords['targ']
    return coords['targ']


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
            coords['targ']['beg'] += diff_beg
            coords['targ']['end'] += diff_end
        else:
            coords['targ']['beg'] -= diff_end # use diff end it's all reverted
            coords['targ']['end'] -= diff_beg
    except Exception as e:
        print(e)
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
        print('      Experiment %10s (%s TADs)' % (exp.name, len(exp.tads)))
        sys.stdout.write('         ')
        for i, tad in enumerate(exp.tads.values()):
            trace.setdefault(crm_name, {})
            if not tad['end'] in trace[crm_name]:
                # MAP
                coords = map_tad(tad, crm_name, exp.resolution, trace=trace)
            else: # already found in an other experiment
                coords = trace[crm_name][tad['end']]
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
            elif not (i + 1) % 10:
                sys.stdout.write(' ')
                sys.stdout.flush()
            else:
                sys.stdout.flush()
        print('')
    if log:
        log.sort(key=lambda x: int(float(x.split(':')[1].split('-')[0])))
        sys.stdout.write('              '+'\n              '.join(log))
        sys.stdout.write('\n              ')


def remap_genome(seed_species, target_species, maf_path, genome, synteny_file=None,
                 write_log=False, synteny_reso=0, clean_all=False):
    """
    Given a set of chromosomes from one species with a given assembly version.
    Remap the coordinates of TAD borders to a new assembly and/or to a new
    species.

    :para genome: a dict containing all chromosomes computed of the target species
    :param original_species: i.e.: 'Homo_sapiens'
    :param target_species: i.e.: 'Mus_musculus', if None, no search for syntenic
       regions will be done
    :param False write_log: path where to write a log file
    :param False clean_all: don't take into account pre-parsed MAF files

    :returns: new genome dictionary with remapped TAD borders, and a trace
       dictionary that allows to reconstruct the process of remapping and finding
       syntenic regions.
    """
    trace = {}
    if write_log:
        log = open(write_log, 'w')
    else:
        log = sys.stdout

    print('\n  Loading MAF alignment')
    get_alignments(seed_species, target_species, maf_path,
                   synteny_file, synteny_reso, clean_all)

    new_genome = {}
    for crm_num, crm in enumerate(genome, 1):
        print('\n   %2s- Chromosome:' % (crm_num), crm)
        # remap and syntenic region
        print('           Searching syntenic regions from %s to %s' %(
            seed_species.capitalize().replace('_', ' '),
            target_species.capitalize().replace('_', ' ')))
        convert_chromosome(genome[crm], new_genome, trace=trace)

        for t in sorted(trace[crm]):
            try:
                log.write('      %4s : %2s:%9s-%9s -> %2s:%9s-%9s\n' %(
                    int(t)                             ,
                    trace[crm][t]['from']['crm']       ,
                    trace[crm][t]['from']['beg']       ,
                    trace[crm][t]['from']['end']       ,
                    trace[crm][t]['syntenic at']['crm'],
                    trace[crm][t]['syntenic at']['beg'],
                    trace[crm][t]['syntenic at']['end'],
                    ))
            except KeyError:
                log.write('      %4s : %2s:%9s-%9s -> %22s\n' %(
                    int(t),
                    trace[crm][t]['from']['crm'],
                    trace[crm][t]['from']['beg'],
                    trace[crm][t]['from']['end'],
                    'None',
                    ))
    if write_log:
        log.close()
    return new_genome, trace


def save_new_genome(genome, trace, check=False, rootdir='./'):
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
            # reorder the TADs in increasing order of their end position
            end, _, score = zip(*sorted(zip(
                *[exp._tads[k] for k in ['end', 'start', 'score']])))
            exp._tads['start'] = [0.] + [v - 1 for v in end[:-1]]
            exp._tads['end'  ] = list(end)
            exp._tads['score'] = list(score)
            if check:
                # check TADs that have synteny
                tadcnt = 0
                new_tads = {}
                for tad in exp.tads:
                    try:
                        if trace[crm][
                            exp.tads[tad]['end']]['syntenic at']['crm'] is None:
                            continue
                    except KeyError:
                        print('Not found:', crm, exp.tads[tad]['end'],
                               trace[crm][exp.tads[tad]['end']])
                        continue
                    new_tads[tadcnt] = exp.tads[tad]
                    tadcnt += 1
                exp.tads = new_tads
            else:
                # create new genome on which are mapped the old coordinates
                tads, norm = parse_tads(exp._tads)
                last = max(tads.keys())
                if not exp.size:
                    exp.size = tads[last]['end']
                exp.norm  = norm
                exp.tads = tads
            if not os.path.exists(rootdir):
                os.mkdir(rootdir)
            exp.write_tad_borders(density=False,
                                  savedata=os.path.join(rootdir, new_crm.name + '.tsv'))
