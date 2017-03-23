"""
14 nov. 2014

Definition and mapping of restriction enymes
"""

from re import compile
from warnings import warn

def count_re_fragments(fnam):
    frag_count = {}
    fhandler = open(fnam)
    line = fhandler.next()
    while line.startswith('#'):
        line = fhandler.next()
    try:
        while True:
            _, cr1, _, _, _, rs1, _, cr2, _, _, _, rs2, _ = line.split('\t')
            try:
                frag_count[(cr1, rs1)] += 1
            except KeyError:
                frag_count[(cr1, rs1)] = 1
            try:
                frag_count[(cr2, rs2)] += 1
            except KeyError:
                frag_count[(cr2, rs2)] = 1
            line = fhandler.next()
    except StopIteration:
        pass
    return frag_count


def map_re_sites_nochunk(enzyme_name, genome_seq, verbose=False):
    """
    map all restriction enzyme (RE) sites of a given enzyme in a genome.
    Position of a RE site is defined as the genomic coordinate of the first
    nucleotide after the first cut (genomic coordinate starts at 1).
    In the case of HindIII the genomic coordinate is this one:

    123456 789...
           |
           v
    -----A|AGCT T--------------
    -----T TCGA|A--------------

    In this example the coordinate of the RE site would be 7.


    :param enzyme_name: name of the enzyme to map (upper/lower case are
       important)
    :param genome_seq: a dictionary containing the genomic sequence by
       chromosome
    """
    warn('WARNING: not reviewed since multiple-cut branch, and the use of regexpinstead of index')
    if isinstance(enzyme_name, str):
        enzyme_names = [enzyme_name]
    elif isinstance(enzyme_name, list):
        enzyme_names = enzyme_name
    enzymes = {}
    for name in enzyme_names:
        enzymes[name] = RESTRICTION_ENZYMES[name]

    # we match the full cut-site but report the position after the cut site
    # (third group of the regexp)
    restring = ('%s') % ('|'.join(['((%s)(%s))' % tuple(enzymes[n].split('|'))
                                   for n in enzymes]))
    # IUPAC conventions
    restring.replace('R', '[AG]')
    restring.replace('Y', '[CT]')
    restring.replace('M', '[AC]')
    restring.replace('K', '[GT]')
    restring.replace('S', '[CG]')
    restring.replace('W', '[AT]')
    restring.replace('H', '[ACT]')
    restring.replace('B', '[CGT]')
    restring.replace('V', '[ACG]')
    restring.replace('D', '[AGT]')
    restring.replace('N', '[ATGC]')

    enz_pattern = compile(restring)

    frags = {}
    count = 0
    for crm in genome_seq:
        seq = genome_seq[crm]
        frags[crm] = [1]
        for match in enz_pattern.finditer(seq):
            pos = match.start(3) + 1  # get 3rd group of regex (after the cut)
            frags[crm].append(pos)
            count += 1
        # at the end of last chunk we add the chromosome length
        frags[crm].append(len(seq))
    if verbose:
        print 'Found %d RE sites' % count
    return frags


def map_re_sites(enzyme_name, genome_seq, frag_chunk=100000, verbose=False):
    """
    map all restriction enzyme (RE) sites of a given enzyme in a genome.
    Position of a RE site is defined as the genomic coordinate of the first
    nucleotide after the first cut (genomic coordinate starts at 1).
    In the case of HindIII the genomic coordinate is this one:

    123456 789...
           |
           v
    -----A|AGCT T--------------
    -----T TCGA|A--------------

    In this example the coordinate of the RE site would be 7.


    :param enzyme_name: name of the enzyme to map (upper/lower case are
       important)
    :param genome_seq: a dictionary containing the genomic sequence by
       chromosome
    :param 100000 frag_chunk: in order to optimize the search for nearby RE
       sites, each chromosome is splitted into chunks.
    """
    if isinstance(enzyme_name, str):
        enzyme_names = [enzyme_name]
    elif isinstance(enzyme_name, list):
        enzyme_names = enzyme_name
    enzymes = {}
    for name in enzyme_names:
        enzymes[name] = RESTRICTION_ENZYMES[name]
    # we match the full cut-site but report the position after the cut site
    # (third group of the regexp)
    restring = ('%s') % ('|'.join(['((%s)(%s))' % tuple(enzymes[n].split('|'))
                                   for n in enzymes]))
    # IUPAC conventions
    restring.replace('R', '[AG]')
    restring.replace('Y', '[CT]')
    restring.replace('M', '[AC]')
    restring.replace('K', '[GT]')
    restring.replace('S', '[CG]')
    restring.replace('W', '[AT]')
    restring.replace('H', '[ACT]')
    restring.replace('B', '[CGT]')
    restring.replace('V', '[ACG]')
    restring.replace('D', '[AGT]')
    restring.replace('N', '[ATGC]')

    enz_pattern = compile(restring)
    
    frags = {}
    count = 0
    for crm in genome_seq:
        seq = genome_seq[crm]
        frags[crm] = dict([(i, []) for i in xrange(len(seq) / frag_chunk + 1)])
        frags[crm][0] = [1]
        for match in enz_pattern.finditer(seq):
            pos = match.start(3) + 1
            frags[crm][pos / frag_chunk].append(pos)
            count += 1
        # at the end of last chunk we add the chromosome length
        frags[crm][len(seq) / frag_chunk].append(len(seq))
        # now we need to assign as first RE site of a fragment the last RE site
        # of previsou fragment, and as last RE site, the first RE site of the
        # next fragment.
        for i in xrange(len(seq) / frag_chunk + 1):
            try:
                try:
                    frags[crm][i].insert(0, frags[crm][i - 1][-2])
                except IndexError:
                    # in case there was no RE site in previous fragment
                    frags[crm][i].insert(0, frags[crm][i - 1][-1])
            except KeyError:
                # it is the very first chunk
                pass
            plus = 1
            while True:
                try:
                    frags[crm][i].append(frags[crm][i + plus][0])
                    break
                except IndexError:
                    # no RE site in this fragment, get "next RE site" from next
                    plus += 1
                except KeyError:
                    # end of the chromosome
                    break
    if verbose:
        print 'Found %d RE sites' % count
    return frags


def complementary(seq):
    trs = dict([(nt1, nt2) for nt1, nt2 in zip('ATGCN', 'TACGN')])
    return ''.join([trs[s] for s in seq[::-1]])


def repaired(r_enz):
    """
    returns the resulting sequence after reparation of two digested and repaired
    ends, marking dangling ends.
    """
    site = RESTRICTION_ENZYMES[r_enz]
    beg, end = site.split('|')
    site = site.replace('|', '')
    return complementary(beg + site[min(len(beg), len(end)) :
                                    max(len(beg), len(end))])


def religateds(r_enzs):
    """
    returns the resulting list of all possible sequences after religation of two
    digested and repaired ends.
    """
    ligations = {}
    for r_enz1 in r_enzs:
        for r_enz2 in r_enzs:
            site1 = RESTRICTION_ENZYMES[r_enz1]
            site2 = RESTRICTION_ENZYMES[r_enz2]
            beg1, end1 = site1.split('|')
            _, end2 = site2.split('|')
            site1 = site1.replace('|', '')
            site2 = site2.replace('|', '')
            ligations[(r_enz1, r_enz2)] = beg1 + end1[:len(end1)-len(beg1)] + end2
    return ligations


class RE_dict(dict):
    def __getitem__(self, i):
        try:
            return super(RE_dict, self).__getitem__(i)
        except KeyError:
            for nam in self:
                if nam.lower() == i.lower():
                    return self[nam]
            raise KeyError('Restriction Enzyme %s not found\n' % (i))
    

RESTRICTION_ENZYMES = RE_dict([('AanI'       , 'TTA|TAA'                     ),
                               ('AarI'       , 'CACCTGC|'                    ),
                               ('AasI'       , 'GACNNNN|NNGTC'               ),
                               ('AatII'      , 'GACGT|C'                     ),
                               ('AbaSI'      , 'C|'                          ),
                               ('AbsI'       , 'CC|TCGAGG'                   ),
                               ('Acc16I'     , 'TGC|GCA'                     ),
                               ('Acc36I'     , 'ACCTGC|'                     ),
                               ('Acc65I'     , 'G|GTACC'                     ),
                               ('AccB1I'     , 'G|GYRCC'                     ),
                               ('AccB7I'     , 'CCANNNN|NTGG'                ),
                               ('AccBSI'     , 'CCG|CTC'                     ),
                               ('AccI'       , 'GT|MKAC'                     ),
                               ('AccII'      , 'CG|CG'                       ),
                               ('AccIII'     , 'T|CCGGA'                     ),
                               ('AceIII'     , 'CAGCTC|'                     ),
                               ('AciI'       , 'C|CGC'                       ),
                               ('AclI'       , 'AA|CGTT'                     ),
                               ('AclWI'      , 'GGATC|'                      ),
                               ('AcoI'       , 'Y|GGCCR'                     ),
                               ('AcsI'       , 'R|AATTY'                     ),
                               ('AcuI'       , 'CTGAAG|'                     ),
                               ('AcvI'       , 'CAC|GTG'                     ),
                               ('AcyI'       , 'GR|CGYC'                     ),
                               ('AdeI'       , 'CACNNN|GTG'                  ),
                               ('AfaI'       , 'GT|AC'                       ),
                               ('AfeI'       , 'AGC|GCT'                     ),
                               ('AfiI'       , 'CCNNNNN|NNGG'                ),
                               ('AflII'      , 'C|TTAAG'                     ),
                               ('AflIII'     , 'A|CRYGT'                     ),
                               ('AgeI'       , 'A|CCGGT'                     ),
                               ('AgsI'       , 'TTS|AA'                      ),
                               ('AhaIII'     , 'TTT|AAA'                     ),
                               ('AhdI'       , 'GACNNN|NNGTC'                ),
                               ('AhlI'       , 'A|CTAGT'                     ),
                               ('AjiI'       , 'CAC|GTC'                     ),
                               ('AjnI'       , '|CCWGG'                      ),
                               ('AjuI'       , 'GAANNNN|NNNTTGG'             ),
                               ('AleI'       , 'CACNN|NNGTG'                 ),
                               ('AlfI'       , 'GC|ANNNNNNTGC'               ),
                               ('AloI'       , 'GAACNN|NNNNTCC'              ),
                               ('AluBI'      , 'AG|CT'                       ),
                               ('AluI'       , 'AG|CT'                       ),
                               ('Alw21I'     , 'GWGCW|C'                     ),
                               ('Alw26I'     , 'GTCTC|'                      ),
                               ('Alw44I'     , 'G|TGCAC'                     ),
                               ('AlwFI'      , 'GAAAYNNNNNRTG|GAAAYNNNNNRTG' ),
                               ('AlwI'       , 'GGATC|'                      ),
                               ('AlwNI'      , 'CAGNNN|CTG'                  ),
                               ('Ama87I'     , 'C|YCGRG'                     ),
                               ('Aor13HI'    , 'T|CCGGA'                     ),
                               ('Aor51HI'    , 'AGC|GCT'                     ),
                               ('AoxI'       , '|GGCC'                       ),
                               ('ApaBI'      , 'GCANNNNN|TGC'                ),
                               ('ApaI'       , 'GGGCC|C'                     ),
                               ('ApaLI'      , 'G|TGCAC'                     ),
                               ('ApeKI'      , 'G|CWGC'                      ),
                               ('ApoI'       , 'R|AATTY'                     ),
                               ('ApyPI'      , 'ATCGAC|'                     ),
                               ('AquII'      , 'GCCGNAC|'                    ),
                               ('AquIII'     , 'GAGGAG|'                     ),
                               ('AquIV'      , 'GRGGAAG|'                    ),
                               ('ArsI'       , 'GACNN|NNNNTTYG'              ),
                               ('AscI'       , 'GG|CGCGCC'                   ),
                               ('AseI'       , 'AT|TAAT'                     ),
                               ('Asi256I'    , 'G|ATC'                       ),
                               ('AsiGI'      , 'A|CCGGT'                     ),
                               ('AsiSI'      , 'GCGAT|CGC'                   ),
                               ('Asp700I'    , 'GAANN|NNTTC'                 ),
                               ('Asp718I'    , 'G|GTACC'                     ),
                               ('AspA2I'     , 'C|CTAGG'                     ),
                               ('AspBHI'     , 'YSCNS|'                      ),
                               ('AspLEI'     , 'GCG|C'                       ),
                               ('AspS9I'     , 'G|GNCC'                      ),
                               ('AssI'       , 'AGT|ACT'                     ),
                               ('AsuC2I'     , 'CC|SGG'                      ),
                               ('AsuHPI'     , 'GGTGA|'                      ),
                               ('AsuI'       , 'G|GNCC'                      ),
                               ('AsuII'      , 'TT|CGAA'                     ),
                               ('AsuNHI'     , 'G|CTAGC'                     ),
                               ('AvaI'       , 'C|YCGRG'                     ),
                               ('AvaII'      , 'G|GWCC'                      ),
                               ('AvaIII'     , 'ATGCAT|ATGCAT'               ),
                               ('AvrII'      , 'C|CTAGG'                     ),
                               ('AxyI'       , 'CC|TNAGG'                    ),
                               ('BaeGI'      , 'GKGCM|C'                     ),
                               ('BaeI'       , 'A|CNNNNGTAYC'                ),
                               ('BalI'       , 'TGG|CCA'                     ),
                               ('BamHI'      , 'G|GATCC'                     ),
                               ('BanI'       , 'G|GYRCC'                     ),
                               ('BanII'      , 'GRGCY|C'                     ),
                               ('BarI'       , 'GAAGNN|NNNNTAC'              ),
                               ('BasI'       , 'CCANNNN|NTGG'                ),
                               ('BauI'       , 'C|ACGAG'                     ),
                               ('Bbr7I'      , 'GAAGAC|'                     ),
                               ('BbrPI'      , 'CAC|GTG'                     ),
                               ('BbsI'       , 'GAAGAC|'                     ),
                               ('Bbv12I'     , 'GWGCW|C'                     ),
                               ('BbvCI'      , 'CC|TCAGC'                    ),
                               ('BbvI'       , 'GCAGC|'                      ),
                               ('BbvII'      , 'GAAGAC|'                     ),
                               ('BccI'       , 'CCATC|'                      ),
                               ('Bce83I'     , 'CTTGAG|'                     ),
                               ('BceAI'      , 'ACGGC|'                      ),
                               ('BcefI'      , 'ACGGC|'                      ),
                               ('BcgI'       , 'CG|ANNNNNNTGC'               ),
                               ('BciT130I'   , 'CC|WGG'                      ),
                               ('BciVI'      , 'GTATCC|'                     ),
                               ('BclI'       , 'T|GATCA'                     ),
                               ('BcnI'       , 'CC|SGG'                      ),
                               ('BcoDI'      , 'GTCTC|'                      ),
                               ('BcuI'       , 'A|CTAGT'                     ),
                               ('BdaI'       , 'TG|ANNNNNNTCA'               ),
                               ('BetI'       , 'W|CCGGW'                     ),
                               ('BfaI'       , 'C|TAG'                       ),
                               ('BfiI'       , 'ACTGGG|'                     ),
                               ('BfmI'       , 'C|TRYAG'                     ),
                               ('BfoI'       , 'RGCGC|Y'                     ),
                               ('BfrI'       , 'C|TTAAG'                     ),
                               ('BfuAI'      , 'ACCTGC|'                     ),
                               ('BfuCI'      , '|GATC'                       ),
                               ('BfuI'       , 'GTATCC|'                     ),
                               ('BglI'       , 'GCCNNNN|NGGC'                ),
                               ('BglII'      , 'A|GATCT'                     ),
                               ('BinI'       , 'GGATC|'                      ),
                               ('BisI'       , 'GC|NGC'                      ),
                               ('BlnI'       , 'C|CTAGG'                     ),
                               ('BlpI'       , 'GC|TNAGC'                    ),
                               ('BlsI'       , 'GCN|GC'                      ),
                               ('BmcAI'      , 'AGT|ACT'                     ),
                               ('Bme1390I'   , 'CC|NGG'                      ),
                               ('Bme18I'     , 'G|GWCC'                      ),
                               ('BmeDI'      , 'C|'                          ),
                               ('BmeRI'      , 'GACNNN|NNGTC'                ),
                               ('BmeT110I'   , 'C|YCGRG'                     ),
                               ('BmgBI'      , 'CAC|GTC'                     ),
                               ('BmgI'       , 'GKGCCC|GKGCCC'               ),
                               ('BmgT120I'   , 'G|GNCC'                      ),
                               ('BmiI'       , 'GGN|NCC'                     ),
                               ('BmrFI'      , 'CC|NGG'                      ),
                               ('BmrI'       , 'ACTGGG|'                     ),
                               ('BmsI'       , 'GCATC|'                      ),
                               ('BmtI'       , 'GCTAG|C'                     ),
                               ('BmuI'       , 'ACTGGG|'                     ),
                               ('BoxI'       , 'GACNN|NNGTC'                 ),
                               ('BpiI'       , 'GAAGAC|'                     ),
                               ('BplI'       , 'GAG|NNNNNCTC'                ),
                               ('BpmI'       , 'CTGGAG|'                     ),
                               ('Bpu10I'     , 'CC|TNAGC'                    ),
                               ('Bpu1102I'   , 'GC|TNAGC'                    ),
                               ('Bpu14I'     , 'TT|CGAA'                     ),
                               ('BpuEI'      , 'CTTGAG|'                     ),
                               ('BpuMI'      , 'CC|SGG'                      ),
                               ('BpvUI'      , 'CGAT|CG'                     ),
                               ('Bsa29I'     , 'AT|CGAT'                     ),
                               ('BsaAI'      , 'YAC|GTR'                     ),
                               ('BsaBI'      , 'GATNN|NNATC'                 ),
                               ('BsaHI'      , 'GR|CGYC'                     ),
                               ('BsaI'       , 'GGTCTC|'                     ),
                               ('BsaJI'      , 'C|CNNGG'                     ),
                               ('BsaWI'      , 'W|CCGGW'                     ),
                               ('BsaXI'      , 'AC|NNNNNCTCC'                ),
                               ('BsbI'       , 'CAACAC|'                     ),
                               ('Bsc4I'      , 'CCNNNNN|NNGG'                ),
                               ('BscAI'      , 'GCATC|'                      ),
                               ('BscGI'      , 'CCCGT|CCCGT'                 ),
                               ('Bse118I'    , 'R|CCGGY'                     ),
                               ('Bse1I'      , 'ACTGG|'                      ),
                               ('Bse21I'     , 'CC|TNAGG'                    ),
                               ('Bse3DI'     , 'GCAATG|'                     ),
                               ('Bse8I'      , 'GATNN|NNATC'                 ),
                               ('BseAI'      , 'T|CCGGA'                     ),
                               ('BseBI'      , 'CC|WGG'                      ),
                               ('BseCI'      , 'AT|CGAT'                     ),
                               ('BseDI'      , 'C|CNNGG'                     ),
                               ('BseGI'      , 'GGATG|'                      ),
                               ('BseJI'      , 'GATNN|NNATC'                 ),
                               ('BseLI'      , 'CCNNNNN|NNGG'                ),
                               ('BseMI'      , 'GCAATG|'                     ),
                               ('BseMII'     , 'CTCAG|'                      ),
                               ('BseNI'      , 'ACTGG|'                      ),
                               ('BsePI'      , 'G|CGCGC'                     ),
                               ('BseRI'      , 'GAGGAG|'                     ),
                               ('BseSI'      , 'GKGCM|C'                     ),
                               ('BseX3I'     , 'C|GGCCG'                     ),
                               ('BseXI'      , 'GCAGC|'                      ),
                               ('BseYI'      , 'C|CCAGC'                     ),
                               ('BsgI'       , 'GTGCAG|'                     ),
                               ('Bsh1236I'   , 'CG|CG'                       ),
                               ('Bsh1285I'   , 'CGRY|CG'                     ),
                               ('BshFI'      , 'GG|CC'                       ),
                               ('BshNI'      , 'G|GYRCC'                     ),
                               ('BshTI'      , 'A|CCGGT'                     ),
                               ('BshVI'      , 'AT|CGAT'                     ),
                               ('BsiEI'      , 'CGRY|CG'                     ),
                               ('BsiHKAI'    , 'GWGCW|C'                     ),
                               ('BsiHKCI'    , 'C|YCGRG'                     ),
                               ('BsiI'       , 'C|ACGAG'                     ),
                               ('BsiSI'      , 'C|CGG'                       ),
                               ('BsiWI'      , 'C|GTACG'                     ),
                               ('BsiYI'      , 'CCNNNNN|NNGG'                ),
                               ('BslFI'      , 'GGGAC|'                      ),
                               ('BslI'       , 'CCNNNNN|NNGG'                ),
                               ('BsmAI'      , 'GTCTC|'                      ),
                               ('BsmBI'      , 'CGTCTC|'                     ),
                               ('BsmFI'      , 'GGGAC|'                      ),
                               ('BsmI'       , 'GAATGC|'                     ),
                               ('BsnI'       , 'GG|CC'                       ),
                               ('Bso31I'     , 'GGTCTC|'                     ),
                               ('BsoBI'      , 'C|YCGRG'                     ),
                               ('Bsp119I'    , 'TT|CGAA'                     ),
                               ('Bsp120I'    , 'G|GGCCC'                     ),
                               ('Bsp1286I'   , 'GDGCH|C'                     ),
                               ('Bsp13I'     , 'T|CCGGA'                     ),
                               ('Bsp1407I'   , 'T|GTACA'                     ),
                               ('Bsp143I'    , '|GATC'                       ),
                               ('Bsp1720I'   , 'GC|TNAGC'                    ),
                               ('Bsp19I'     , 'C|CATGG'                     ),
                               ('Bsp24I'     , 'GACN|NNNNNTGG'               ),
                               ('Bsp68I'     , 'TCG|CGA'                     ),
                               ('BspACI'     , 'C|CGC'                       ),
                               ('BspCNI'     , 'CTCAG|'                      ),
                               ('BspD6I'     , 'GACTC|'                      ),
                               ('BspDI'      , 'AT|CGAT'                     ),
                               ('BspEI'      , 'T|CCGGA'                     ),
                               ('BspFNI'     , 'CG|CG'                       ),
                               ('BspGI'      , 'CTGGAC|CTGGAC'               ),
                               ('BspHI'      , 'T|CATGA'                     ),
                               ('BspLI'      , 'GGN|NCC'                     ),
                               ('BspLU11I'   , 'A|CATGT'                     ),
                               ('BspMI'      , 'ACCTGC|'                     ),
                               ('BspMII'     , 'T|CCGGA'                     ),
                               ('BspNCI'     , 'CCAGA|CCAGA'                 ),
                               ('BspOI'      , 'GCTAG|C'                     ),
                               ('BspPI'      , 'GGATC|'                      ),
                               ('BspQI'      , 'GCTCTTC|'                    ),
                               ('BspT104I'   , 'TT|CGAA'                     ),
                               ('BspT107I'   , 'G|GYRCC'                     ),
                               ('BspTI'      , 'C|TTAAG'                     ),
                               ('BsrBI'      , 'CCG|CTC'                     ),
                               ('BsrDI'      , 'GCAATG|'                     ),
                               ('BsrFI'      , 'R|CCGGY'                     ),
                               ('BsrGI'      , 'T|GTACA'                     ),
                               ('BsrI'       , 'ACTGG|'                      ),
                               ('BsrSI'      , 'ACTGG|'                      ),
                               ('BssAI'      , 'R|CCGGY'                     ),
                               ('BssECI'     , 'C|CNNGG'                     ),
                               ('BssHII'     , 'G|CGCGC'                     ),
                               ('BssKI'      , '|CCNGG'                      ),
                               ('BssMI'      , '|GATC'                       ),
                               ('BssNAI'     , 'GTA|TAC'                     ),
                               ('BssNI'      , 'GR|CGYC'                     ),
                               ('BssSI'      , 'C|ACGAG'                     ),
                               ('BssT1I'     , 'C|CWWGG'                     ),
                               ('Bst1107I'   , 'GTA|TAC'                     ),
                               ('Bst2BI'     , 'C|ACGAG'                     ),
                               ('Bst2UI'     , 'CC|WGG'                      ),
                               ('Bst4CI'     , 'ACN|GT'                      ),
                               ('Bst6I'      , 'CTCTTC|'                     ),
                               ('BstACI'     , 'GR|CGYC'                     ),
                               ('BstAFI'     , 'C|TTAAG'                     ),
                               ('BstAPI'     , 'GCANNNN|NTGC'                ),
                               ('BstAUI'     , 'T|GTACA'                     ),
                               ('BstBAI'     , 'YAC|GTR'                     ),
                               ('BstBI'      , 'TT|CGAA'                     ),
                               ('BstC8I'     , 'GCN|NGC'                     ),
                               ('BstDEI'     , 'C|TNAG'                      ),
                               ('BstDSI'     , 'C|CRYGG'                     ),
                               ('BstEII'     , 'G|GTNACC'                    ),
                               ('BstENI'     , 'CCTNN|NNNAGG'                ),
                               ('BstF5I'     , 'GGATG|'                      ),
                               ('BstFNI'     , 'CG|CG'                       ),
                               ('BstH2I'     , 'RGCGC|Y'                     ),
                               ('BstHHI'     , 'GCG|C'                       ),
                               ('BstKTI'     , 'GAT|C'                       ),
                               ('BstMAI'     , 'GTCTC|'                      ),
                               ('BstMBI'     , '|GATC'                       ),
                               ('BstMCI'     , 'CGRY|CG'                     ),
                               ('BstMWI'     , 'GCNNNNN|NNGC'                ),
                               ('BstNI'      , 'CC|WGG'                      ),
                               ('BstNSI'     , 'RCATG|Y'                     ),
                               ('BstOI'      , 'CC|WGG'                      ),
                               ('BstPAI'     , 'GACNN|NNGTC'                 ),
                               ('BstPI'      , 'G|GTNACC'                    ),
                               ('BstSCI'     , '|CCNGG'                      ),
                               ('BstSFI'     , 'C|TRYAG'                     ),
                               ('BstSLI'     , 'GKGCM|C'                     ),
                               ('BstSNI'     , 'TAC|GTA'                     ),
                               ('BstUI'      , 'CG|CG'                       ),
                               ('BstV1I'     , 'GCAGC|'                      ),
                               ('BstV2I'     , 'GAAGAC|'                     ),
                               ('BstX2I'     , 'R|GATCY'                     ),
                               ('BstXI'      , 'CCANNNNN|NTGG'               ),
                               ('BstYI'      , 'R|GATCY'                     ),
                               ('BstZ17I'    , 'GTA|TAC'                     ),
                               ('BstZI'      , 'C|GGCCG'                     ),
                               ('Bsu15I'     , 'AT|CGAT'                     ),
                               ('Bsu36I'     , 'CC|TNAGG'                    ),
                               ('BsuI'       , 'GTATCC|'                     ),
                               ('BsuRI'      , 'GG|CC'                       ),
                               ('BtgI'       , 'C|CRYGG'                     ),
                               ('BtgZI'      , 'GCGATG|'                     ),
                               ('BthCI'      , 'GCNG|C'                      ),
                               ('BtrI'       , 'CAC|GTC'                     ),
                               ('BtsCI'      , 'GGATG|'                      ),
                               ('BtsI'       , 'GCAGTG|'                     ),
                               ('BtsIMutI'   , 'CAGTG|'                      ),
                               ('BtuMI'      , 'TCG|CGA'                     ),
                               ('BveI'       , 'ACCTGC|'                     ),
                               ('Cac8I'      , 'GCN|NGC'                     ),
                               ('CaiI'       , 'CAGNNN|CTG'                  ),
                               ('CauII'      , 'CC|SGG'                      ),
                               ('CchII'      , 'GGARGA|'                     ),
                               ('CchIII'     , 'CCCAAG|'                     ),
                               ('CciI'       , 'T|CATGA'                     ),
                               ('CciNI'      , 'GC|GGCCGC'                   ),
                               ('Cdi630V'    , 'CAAAAA|CAAAAA'               ),
                               ('CdiI'       , 'CATC|G'                      ),
                               ('CdpI'       , 'GCGGAG|'                     ),
                               ('CfoI'       , 'GCG|C'                       ),
                               ('Cfr10I'     , 'R|CCGGY'                     ),
                               ('Cfr13I'     , 'G|GNCC'                      ),
                               ('Cfr42I'     , 'CCGC|GG'                     ),
                               ('Cfr9I'      , 'C|CCGGG'                     ),
                               ('CfrI'       , 'Y|GGCCR'                     ),
                               ('Cgl13032I'  , 'GGCGCA|GGCGCA'               ),
                               ('Cgl13032II' , 'ACGABGG|ACGABGG'             ),
                               ('ChaI'       , 'GATC|'                       ),
                               ('CjeFIII'    , 'GCAAGG|GCAAGG'               ),
                               ('CjeFV'      , 'GGRCA|GGRCA'                 ),
                               ('CjeI'       , 'CCA|NNNNNNGT'                ),
                               ('CjeNII'     , 'GAGNNNNNGT|GAGNNNNNGT'       ),
                               ('CjeNIII'    , 'GKAAYG|'                     ),
                               ('CjeP659IV'  , 'CACNNNNNNNGAA|CACNNNNNNNGAA' ),
                               ('CjePI'      , 'CCANN|NNNNNTC'               ),
                               ('CjuI'       , 'CAYNNNNNRTG|CAYNNNNNRTG'     ),
                               ('CjuII'      , 'CAYNNNNNCTC|CAYNNNNNCTC'     ),
                               ('ClaI'       , 'AT|CGAT'                     ),
                               ('CpoI'       , 'CG|GWCCG'                    ),
                               ('CseI'       , 'GACGC|'                      ),
                               ('CsiI'       , 'A|CCWGGT'                    ),
                               ('Csp6I'      , 'G|TAC'                       ),
                               ('CspAI'      , 'A|CCGGT'                     ),
                               ('CspCI'      , 'C|AANNNNNGTGG'               ),
                               ('CspI'       , 'CG|GWCCG'                    ),
                               ('CstMI'      , 'AAGGAG|'                     ),
                               ('CviAII'     , 'C|ATG'                       ),
                               ('CviJI'      , 'RG|CY'                       ),
                               ('CviKI_1'    , 'RG|CY'                       ),
                               ('CviQI'      , 'G|TAC'                       ),
                               ('CviRI'      , 'TG|CA'                       ),
                               ('DdeI'       , 'C|TNAG'                      ),
                               ('DinI'       , 'GGC|GCC'                     ),
                               ('DpnI'       , 'GA|TC'                       ),
                               ('DpnII'      , '|GATC'                       ),
                               ('DraI'       , 'TTT|AAA'                     ),
                               ('DraII'      , 'RG|GNCCY'                    ),
                               ('DraIII'     , 'CACNNN|GTG'                  ),
                               ('DraRI'      , 'CAAGNAC|'                    ),
                               ('DrdI'       , 'GACNNNN|NNGTC'               ),
                               ('DrdII'      , 'GAACCA|GAACCA'               ),
                               ('DriI'       , 'GACNNN|NNGTC'                ),
                               ('DsaI'       , 'C|CRYGG'                     ),
                               ('DseDI'      , 'GACNNNN|NNGTC'               ),
                               ('EaeI'       , 'Y|GGCCR'                     ),
                               ('EagI'       , 'C|GGCCG'                     ),
                               ('Eam1104I'   , 'CTCTTC|'                     ),
                               ('Eam1105I'   , 'GACNNN|NNGTC'                ),
                               ('EarI'       , 'CTCTTC|'                     ),
                               ('EciI'       , 'GGCGGA|'                     ),
                               ('Ecl136II'   , 'GAG|CTC'                     ),
                               ('EclXI'      , 'C|GGCCG'                     ),
                               ('Eco105I'    , 'TAC|GTA'                     ),
                               ('Eco130I'    , 'C|CWWGG'                     ),
                               ('Eco147I'    , 'AGG|CCT'                     ),
                               ('Eco24I'     , 'GRGCY|C'                     ),
                               ('Eco31I'     , 'GGTCTC|'                     ),
                               ('Eco32I'     , 'GAT|ATC'                     ),
                               ('Eco47I'     , 'G|GWCC'                      ),
                               ('Eco47III'   , 'AGC|GCT'                     ),
                               ('Eco52I'     , 'C|GGCCG'                     ),
                               ('Eco53kI'    , 'GAG|CTC'                     ),
                               ('Eco57I'     , 'CTGAAG|'                     ),
                               ('Eco57MI'    , 'CTGRAG|'                     ),
                               ('Eco72I'     , 'CAC|GTG'                     ),
                               ('Eco81I'     , 'CC|TNAGG'                    ),
                               ('Eco88I'     , 'C|YCGRG'                     ),
                               ('Eco91I'     , 'G|GTNACC'                    ),
                               ('EcoHI'      , '|CCSGG'                      ),
                               ('EcoICRI'    , 'GAG|CTC'                     ),
                               ('EcoNI'      , 'CCTNN|NNNAGG'                ),
                               ('EcoO109I'   , 'RG|GNCCY'                    ),
                               ('EcoO65I'    , 'G|GTNACC'                    ),
                               ('EcoRI'      , 'G|AATTC'                     ),
                               ('EcoRII'     , '|CCWGG'                      ),
                               ('EcoRV'      , 'GAT|ATC'                     ),
                               ('EcoT14I'    , 'C|CWWGG'                     ),
                               ('EcoT22I'    , 'ATGCA|T'                     ),
                               ('EcoT38I'    , 'GRGCY|C'                     ),
                               ('EgeI'       , 'GGC|GCC'                     ),
                               ('EheI'       , 'GGC|GCC'                     ),
                               ('ErhI'       , 'C|CWWGG'                     ),
                               ('EsaBC3I'    , 'TC|GA'                       ),
                               ('EsaSSI'     , 'GACCAC|GACCAC'               ),
                               ('Esp3I'      , 'CGTCTC|'                     ),
                               ('EspI'       , 'GC|TNAGC'                    ),
                               ('FaeI'       , 'CATG|'                       ),
                               ('FaiI'       , 'YA|TR'                       ),
                               ('FalI'       , 'AAG|NNNNNCTT'                ),
                               ('FaqI'       , 'GGGAC|'                      ),
                               ('FatI'       , '|CATG'                       ),
                               ('FauI'       , 'CCCGC|'                      ),
                               ('FauNDI'     , 'CA|TATG'                     ),
                               ('FbaI'       , 'T|GATCA'                     ),
                               ('FblI'       , 'GT|MKAC'                     ),
                               ('FinI'       , 'GGGAC|GGGAC'                 ),
                               ('FmuI'       , 'GGNC|C'                      ),
                               ('Fnu4HI'     , 'GC|NGC'                      ),
                               ('FnuDII'     , 'CG|CG'                       ),
                               ('FokI'       , 'GGATG|'                      ),
                               ('FriOI'      , 'GRGCY|C'                     ),
                               ('FseI'       , 'GGCCGG|CC'                   ),
                               ('Fsp4HI'     , 'GC|NGC'                      ),
                               ('FspAI'      , 'RTGC|GCAY'                   ),
                               ('FspBI'      , 'C|TAG'                       ),
                               ('FspEI'      , 'CC|'                         ),
                               ('FspI'       , 'TGC|GCA'                     ),
                               ('GauT27I'    , 'CGCGCAGG|CGCGCAGG'           ),
                               ('GdiII'      , 'C|GGCCR'                     ),
                               ('GlaI'       , 'GC|GC'                       ),
                               ('GluI'       , 'GC|NGC'                      ),
                               ('GsaI'       , 'CCCAG|C'                     ),
                               ('GsuI'       , 'CTGGAG|'                     ),
                               ('HaeI'       , 'WGG|CCW'                     ),
                               ('HaeII'      , 'RGCGC|Y'                     ),
                               ('HaeIII'     , 'GG|CC'                       ),
                               ('HapII'      , 'C|CGG'                       ),
                               ('HauII'      , 'TGGCCA|'                     ),
                               ('HgaI'       , 'GACGC|'                      ),
                               ('HgiAI'      , 'GWGCW|C'                     ),
                               ('HgiCI'      , 'G|GYRCC'                     ),
                               ('HgiEII'     , 'ACCNNNNNNGGT|ACCNNNNNNGGT'   ),
                               ('HgiJII'     , 'GRGCY|C'                     ),
                               ('HhaI'       , 'GCG|C'                       ),
                               ('Hin1I'      , 'GR|CGYC'                     ),
                               ('Hin1II'     , 'CATG|'                       ),
                               ('Hin4I'      , 'GAY|NNNNNVTC'                ),
                               ('Hin4II'     , 'CCTTC|'                      ),
                               ('Hin6I'      , 'G|CGC'                       ),
                               ('HinP1I'     , 'G|CGC'                       ),
                               ('HincII'     , 'GTY|RAC'                     ),
                               ('HindII'     , 'GTY|RAC'                     ),
                               ('HindIII'    , 'A|AGCTT'                     ),
                               ('HinfI'      , 'G|ANTC'                      ),
                               ('HpaI'       , 'GTT|AAC'                     ),
                               ('HpaII'      , 'C|CGG'                       ),
                               ('HphI'       , 'GGTGA|'                      ),
                               ('Hpy166II'   , 'GTN|NAC'                     ),
                               ('Hpy178III'  , 'TC|NNGA'                     ),
                               ('Hpy188I'    , 'TCN|GA'                      ),
                               ('Hpy188III'  , 'TC|NNGA'                     ),
                               ('Hpy8I'      , 'GTN|NAC'                     ),
                               ('Hpy99I'     , 'CGWCG|'                      ),
                               ('Hpy99XIII'  , 'GCCTA|GCCTA'                 ),
                               ('Hpy99XIV'   , 'GGWTAA|GGWTAA'               ),
                               ('HpyAV'      , 'CCTTC|'                      ),
                               ('HpyCH4III'  , 'ACN|GT'                      ),
                               ('HpyCH4IV'   , 'A|CGT'                       ),
                               ('HpyCH4V'    , 'TG|CA'                       ),
                               ('HpyF10VI'   , 'GCNNNNN|NNGC'                ),
                               ('HpyF3I'     , 'C|TNAG'                      ),
                               ('HpySE526I'  , 'A|CGT'                       ),
                               ('Hsp92I'     , 'GR|CGYC'                     ),
                               ('Hsp92II'    , 'CATG|'                       ),
                               ('HspAI'      , 'G|CGC'                       ),
                               ('Jma19592I'  , 'GTATNAC|GTATNAC'             ),
                               ('KasI'       , 'G|GCGCC'                     ),
                               ('KflI'       , 'GG|GWCCC'                    ),
                               ('Kpn2I'      , 'T|CCGGA'                     ),
                               ('KpnI'       , 'GGTAC|C'                     ),
                               ('KroI'       , 'G|CCGGC'                     ),
                               ('Ksp22I'     , 'T|GATCA'                     ),
                               ('Ksp632I'    , 'CTCTTC|'                     ),
                               ('KspAI'      , 'GTT|AAC'                     ),
                               ('KspI'       , 'CCGC|GG'                     ),
                               ('Kzo9I'      , '|GATC'                       ),
                               ('LguI'       , 'GCTCTTC|'                    ),
                               ('LpnI'       , 'RGC|GCY'                     ),
                               ('LpnPI'      , 'CCDG|'                       ),
                               ('Lsp1109I'   , 'GCAGC|'                      ),
                               ('LweI'       , 'GCATC|'                      ),
                               ('MabI'       , 'A|CCWGGT'                    ),
                               ('MaeI'       , 'C|TAG'                       ),
                               ('MaeII'      , 'A|CGT'                       ),
                               ('MaeIII'     , '|GTNAC'                      ),
                               ('MalI'       , 'GA|TC'                       ),
                               ('MaqI'       , 'CRTTGAC|'                    ),
                               ('MauBI'      , 'CG|CGCGCG'                   ),
                               ('MbiI'       , 'CCG|CTC'                     ),
                               ('MboI'       , '|GATC'                       ),
                               ('MboII'      , 'GAAGA|'                      ),
                               ('McaTI'      , 'GCGC|GC'                     ),
                               ('McrI'       , 'CGRY|CG'                     ),
                               ('MfeI'       , 'C|AATTG'                     ),
                               ('MflI'       , 'R|GATCY'                     ),
                               ('MhlI'       , 'GDGCH|C'                     ),
                               ('MjaIV'      , 'GTNNAC|GTNNAC'               ),
                               ('MkaDII'     , 'GAGAYGT|GAGAYGT'             ),
                               ('MlsI'       , 'TGG|CCA'                     ),
                               ('MluCI'      , '|AATT'                       ),
                               ('MluI'       , 'A|CGCGT'                     ),
                               ('MluNI'      , 'TGG|CCA'                     ),
                               ('Mly113I'    , 'GG|CGCC'                     ),
                               ('MlyI'       , 'GAGTC|'                      ),
                               ('MmeI'       , 'TCCRAC|'                     ),
                               ('MnlI'       , 'CCTC|'                       ),
                               ('Mph1103I'   , 'ATGCA|T'                     ),
                               ('MreI'       , 'CG|CCGGCG'                   ),
                               ('MroI'       , 'T|CCGGA'                     ),
                               ('MroNI'      , 'G|CCGGC'                     ),
                               ('MroXI'      , 'GAANN|NNTTC'                 ),
                               ('MscI'       , 'TGG|CCA'                     ),
                               ('MseI'       , 'T|TAA'                       ),
                               ('MslI'       , 'CAYNN|NNRTG'                 ),
                               ('Msp20I'     , 'TGG|CCA'                     ),
                               ('MspA1I'     , 'CMG|CKG'                     ),
                               ('MspCI'      , 'C|TTAAG'                     ),
                               ('MspI'       , 'C|CGG'                       ),
                               ('MspJI'      , 'CNNR|'                       ),
                               ('MspR9I'     , 'CC|NGG'                      ),
                               ('MssI'       , 'GTTT|AAAC'                   ),
                               ('MstI'       , 'TGC|GCA'                     ),
                               ('MunI'       , 'C|AATTG'                     ),
                               ('Mva1269I'   , 'GAATGC|'                     ),
                               ('MvaI'       , 'CC|WGG'                      ),
                               ('MvnI'       , 'CG|CG'                       ),
                               ('MvrI'       , 'CGAT|CG'                     ),
                               ('MwoI'       , 'GCNNNNN|NNGC'                ),
                               ('NaeI'       , 'GCC|GGC'                     ),
                               ('NarI'       , 'GG|CGCC'                     ),
                               ('NciI'       , 'CC|SGG'                      ),
                               ('NcoI'       , 'C|CATGG'                     ),
                               ('NdeI'       , 'CA|TATG'                     ),
                               ('NdeII'      , '|GATC'                       ),
                               ('NgoAVIII'   , '|GACNNNNNTGA'                ),
                               ('NgoMIV'     , 'G|CCGGC'                     ),
                               ('NhaXI'      , 'CAAGRAG|CAAGRAG'             ),
                               ('NheI'       , 'G|CTAGC'                     ),
                               ('NlaCI'      , 'CATCAC|'                     ),
                               ('NlaIII'     , 'CATG|'                       ),
                               ('NlaIV'      , 'GGN|NCC'                     ),
                               ('Nli3877I'   , 'CYCGR|G'                     ),
                               ('NmeAIII'    , 'GCCGAG|'                     ),
                               ('NmeDI'      , '|RCCGGY'                     ),
                               ('NmuCI'      , '|GTSAC'                      ),
                               ('NotI'       , 'GC|GGCCGC'                   ),
                               ('NruI'       , 'TCG|CGA'                     ),
                               ('NsbI'       , 'TGC|GCA'                     ),
                               ('NsiI'       , 'ATGCA|T'                     ),
                               ('NspBII'     , 'CMG|CKG'                     ),
                               ('NspI'       , 'RCATG|Y'                     ),
                               ('NspV'       , 'TT|CGAA'                     ),
                               ('OliI'       , 'CACNN|NNGTG'                 ),
                               ('PabI'       , 'GTA|C'                       ),
                               ('PacI'       , 'TTAAT|TAA'                   ),
                               ('PaeI'       , 'GCATG|C'                     ),
                               ('PaeR7I'     , 'C|TCGAG'                     ),
                               ('PagI'       , 'T|CATGA'                     ),
                               ('PalAI'      , 'GG|CGCGCC'                   ),
                               ('PasI'       , 'CC|CWGGG'                    ),
                               ('PauI'       , 'G|CGCGC'                     ),
                               ('PceI'       , 'AGG|CCT'                     ),
                               ('PciI'       , 'A|CATGT'                     ),
                               ('PciSI'      , 'GCTCTTC|'                    ),
                               ('PcsI'       , 'WCGNNNN|NNNCGW'              ),
                               ('PctI'       , 'GAATGC|'                     ),
                               ('PdiI'       , 'GCC|GGC'                     ),
                               ('PdmI'       , 'GAANN|NNTTC'                 ),
                               ('PenI'       , 'GCAGT|GCAGT'                 ),
                               ('PfeI'       , 'G|AWTC'                      ),
                               ('Pfl1108I'   , 'TCGTAG|TCGTAG'               ),
                               ('Pfl23II'    , 'C|GTACG'                     ),
                               ('PflFI'      , 'GACN|NNGTC'                  ),
                               ('PflMI'      , 'CCANNNN|NTGG'                ),
                               ('PfoI'       , 'T|CCNGGA'                    ),
                               ('PinAI'      , 'A|CCGGT'                     ),
                               ('PlaDI'      , 'CATCAG|'                     ),
                               ('Ple19I'     , 'CGAT|CG'                     ),
                               ('PleI'       , 'GAGTC|'                      ),
                               ('PluTI'      , 'GGCGC|C'                     ),
                               ('PmaCI'      , 'CAC|GTG'                     ),
                               ('PmeI'       , 'GTTT|AAAC'                   ),
                               ('PmlI'       , 'CAC|GTG'                     ),
                               ('PpiI'       , 'GAACN|NNNNCTC'               ),
                               ('PpsI'       , 'GAGTC|'                      ),
                               ('Ppu10I'     , 'A|TGCAT'                     ),
                               ('Ppu21I'     , 'YAC|GTR'                     ),
                               ('PpuMI'      , 'RG|GWCCY'                    ),
                               ('PscI'       , 'A|CATGT'                     ),
                               ('PshAI'      , 'GACNN|NNGTC'                 ),
                               ('PshBI'      , 'AT|TAAT'                     ),
                               ('PsiI'       , 'TTA|TAA'                     ),
                               ('Psp03I'     , 'GGWC|C'                      ),
                               ('Psp124BI'   , 'GAGCT|C'                     ),
                               ('Psp1406I'   , 'AA|CGTT'                     ),
                               ('Psp5II'     , 'RG|GWCCY'                    ),
                               ('Psp6I'      , '|CCWGG'                      ),
                               ('PspCI'      , 'CAC|GTG'                     ),
                               ('PspEI'      , 'G|GTNACC'                    ),
                               ('PspGI'      , '|CCWGG'                      ),
                               ('PspLI'      , 'C|GTACG'                     ),
                               ('PspN4I'     , 'GGN|NCC'                     ),
                               ('PspOMI'     , 'G|GGCCC'                     ),
                               ('PspOMII'    , 'CGCCCAR|'                    ),
                               ('PspPI'      , 'G|GNCC'                      ),
                               ('PspPPI'     , 'RG|GWCCY'                    ),
                               ('PspPRI'     , 'CCYCAG|'                     ),
                               ('PspXI'      , 'VC|TCGAGB'                   ),
                               ('PsrI'       , 'GAACNN|NNNNTAC'              ),
                               ('PssI'       , 'RGGNC|CY'                    ),
                               ('PstI'       , 'CTGCA|G'                     ),
                               ('PstNI'      , 'CAGNNN|CTG'                  ),
                               ('PsuI'       , 'R|GATCY'                     ),
                               ('PsyI'       , 'GACN|NNGTC'                  ),
                               ('PteI'       , 'G|CGCGC'                     ),
                               ('PvuI'       , 'CGAT|CG'                     ),
                               ('PvuII'      , 'CAG|CTG'                     ),
                               ('R2_BceSIV'  , '|GCAGC'                      ),
                               ('RceI'       , 'CATCGAC|'                    ),
                               ('RdeGBI'     , 'CCGCAG|CCGCAG'               ),
                               ('RdeGBII'    , 'ACCCAG|'                     ),
                               ('RdeGBIII'   , '|TGRYCA'                     ),
                               ('RflFIII'    , 'CGCCAG|CGCCAG'               ),
                               ('RgaI'       , 'GCGAT|CGC'                   ),
                               ('RigI'       , 'GGCCGG|CC'                   ),
                               ('RlaI'       , 'VCW|VCW'                     ),
                               ('RleAI'      , 'CCCACA|'                     ),
                               ('RpaB5I'     , 'CGRGGAC|'                    ),
                               ('RpaBI'      , 'CCCGCAG|'                    ),
                               ('RpaI'       , 'GTYGGAG|'                    ),
                               ('RpaTI'      , 'GRTGGAG|GRTGGAG'             ),
                               ('RruI'       , 'TCG|CGA'                     ),
                               ('RsaI'       , 'GT|AC'                       ),
                               ('RsaNI'      , 'G|TAC'                       ),
                               ('RseI'       , 'CAYNN|NNRTG'                 ),
                               ('Rsr2I'      , 'CG|GWCCG'                    ),
                               ('RsrII'      , 'CG|GWCCG'                    ),
                               ('SacI'       , 'GAGCT|C'                     ),
                               ('SacII'      , 'CCGC|GG'                     ),
                               ('SalI'       , 'G|TCGAC'                     ),
                               ('SanDI'      , 'GG|GWCCC'                    ),
                               ('SapI'       , 'GCTCTTC|'                    ),
                               ('SaqAI'      , 'T|TAA'                       ),
                               ('SatI'       , 'GC|NGC'                      ),
                               ('Sau3AI'     , '|GATC'                       ),
                               ('Sau96I'     , 'G|GNCC'                      ),
                               ('SauI'       , 'CC|TNAGG'                    ),
                               ('SbfI'       , 'CCTGCA|GG'                   ),
                               ('ScaI'       , 'AGT|ACT'                     ),
                               ('SchI'       , 'GAGTC|'                      ),
                               ('SciI'       , 'CTC|GAG'                     ),
                               ('ScrFI'      , 'CC|NGG'                      ),
                               ('SdaI'       , 'CCTGCA|GG'                   ),
                               ('SdeAI'      , 'CAGRAG|'                     ),
                               ('SdeOSI'     , '|GACNNNNRTGA'                ),
                               ('SduI'       , 'GDGCH|C'                     ),
                               ('SecI'       , 'C|CNNGG'                     ),
                               ('SelI'       , '|CGCG'                       ),
                               ('SetI'       , 'ASST|'                       ),
                               ('SexAI'      , 'A|CCWGGT'                    ),
                               ('SfaAI'      , 'GCGAT|CGC'                   ),
                               ('SfaNI'      , 'GCATC|'                      ),
                               ('SfcI'       , 'C|TRYAG'                     ),
                               ('SfeI'       , 'C|TRYAG'                     ),
                               ('SfiI'       , 'GGCCNNNN|NGGCC'              ),
                               ('SfoI'       , 'GGC|GCC'                     ),
                               ('Sfr274I'    , 'C|TCGAG'                     ),
                               ('Sfr303I'    , 'CCGC|GG'                     ),
                               ('SfuI'       , 'TT|CGAA'                     ),
                               ('SgeI'       , 'CNNG|'                       ),
                               ('SgfI'       , 'GCGAT|CGC'                   ),
                               ('SgrAI'      , 'CR|CCGGYG'                   ),
                               ('SgrBI'      , 'CCGC|GG'                     ),
                               ('SgrDI'      , 'CG|TCGACG'                   ),
                               ('SgrTI'      , 'CCDS|'                       ),
                               ('SgsI'       , 'GG|CGCGCC'                   ),
                               ('SimI'       , 'GG|GTC'                      ),
                               ('SlaI'       , 'C|TCGAG'                     ),
                               ('SmaI'       , 'CCC|GGG'                     ),
                               ('SmiI'       , 'ATTT|AAAT'                   ),
                               ('SmiMI'      , 'CAYNN|NNRTG'                 ),
                               ('SmlI'       , 'C|TYRAG'                     ),
                               ('SmoI'       , 'C|TYRAG'                     ),
                               ('SnaBI'      , 'TAC|GTA'                     ),
                               ('SnaI'       , 'GTATAC|GTATAC'               ),
                               ('Sno506I'    , 'GGCCGAG|GGCCGAG'             ),
                               ('SpeI'       , 'A|CTAGT'                     ),
                               ('SphI'       , 'GCATG|C'                     ),
                               ('SplI'       , 'C|GTACG'                     ),
                               ('SpoDI'      , 'GCGGRAG|GCGGRAG'             ),
                               ('SrfI'       , 'GCCC|GGGC'                   ),
                               ('Sse232I'    , 'CG|CCGGCG'                   ),
                               ('Sse8387I'   , 'CCTGCA|GG'                   ),
                               ('Sse8647I'   , 'AG|GWCCT'                    ),
                               ('Sse9I'      , '|AATT'                       ),
                               ('SseBI'      , 'AGG|CCT'                     ),
                               ('SsiI'       , 'C|CGC'                       ),
                               ('SspD5I'     , 'GGTGA|'                      ),
                               ('SspDI'      , 'G|GCGCC'                     ),
                               ('SspI'       , 'AAT|ATT'                     ),
                               ('SstE37I'    , 'CGAAGAC|'                    ),
                               ('SstI'       , 'GAGCT|C'                     ),
                               ('Sth132I'    , 'CCCG|'                       ),
                               ('Sth302II'   , 'CC|GG'                       ),
                               ('StrI'       , 'C|TCGAG'                     ),
                               ('StsI'       , 'GGATG|'                      ),
                               ('StuI'       , 'AGG|CCT'                     ),
                               ('StyD4I'     , '|CCNGG'                      ),
                               ('StyI'       , 'C|CWWGG'                     ),
                               ('SwaI'       , 'ATTT|AAAT'                   ),
                               ('TaaI'       , 'ACN|GT'                      ),
                               ('TaiI'       , 'ACGT|'                       ),
                               ('TaqI'       , 'T|CGA'                       ),
                               ('TaqII'      , 'GACCGA|'                     ),
                               ('TasI'       , '|AATT'                       ),
                               ('TatI'       , 'W|GTACW'                     ),
                               ('TauI'       , 'GCSG|C'                      ),
                               ('TfiI'       , 'G|AWTC'                      ),
                               ('Tru1I'      , 'T|TAA'                       ),
                               ('Tru9I'      , 'T|TAA'                       ),
                               ('TscAI'      , 'CASTG|'                      ),
                               ('TseFI'      , '|GTSAC'                      ),
                               ('TseI'       , 'G|CWGC'                      ),
                               ('TsoI'       , 'TARCCA|'                     ),
                               ('Tsp45I'     , '|GTSAC'                      ),
                               ('Tsp4CI'     , 'ACN|GT'                      ),
                               ('TspDTI'     , 'ATGAA|'                      ),
                               ('TspEI'      , '|AATT'                       ),
                               ('TspGWI'     , 'ACGGA|'                      ),
                               ('TspMI'      , 'C|CCGGG'                     ),
                               ('TspRI'      , 'CASTG|'                      ),
                               ('TssI'       , 'GAGNNNCTC|GAGNNNCTC'         ),
                               ('TstI'       , 'CACN|NNNNNTCC'               ),
                               ('TsuI'       , 'GCGAC|GCGAC'                 ),
                               ('Tth111I'    , 'GACN|NNGTC'                  ),
                               ('Tth111II'   , 'CAARCA|'                     ),
                               ('UbaF11I'    , 'TCGTA|TCGTA'                 ),
                               ('UbaF12I'    , 'CTACNNNGTC|CTACNNNGTC'       ),
                               ('UbaF13I'    , 'GAGNNNNNNCTGG|GAGNNNNNNCTGG' ),
                               ('UbaF14I'    , 'CCANNNNNTCG|CCANNNNNTCG'     ),
                               ('UbaF9I'     , 'TACNNNNNRTGT|TACNNNNNRTGT'   ),
                               ('UbaPI'      , 'CGAACG|CGAACG'               ),
                               ('UcoMSI'     , '|GAGCTC'                     ),
                               ('UnbI'       , '|GGNCC'                      ),
                               ('Van91I'     , 'CCANNNN|NTGG'                ),
                               ('Vha464I'    , 'C|TTAAG'                     ),
                               ('VneI'       , 'G|TGCAC'                     ),
                               ('VpaK11AI'   , '|GGWCC'                      ),
                               ('VpaK11BI'   , 'G|GWCC'                      ),
                               ('VspI'       , 'AT|TAAT'                     ),
                               ('WviI'       , 'CACRAG|'                     ),
                               ('XagI'       , 'CCTNN|NNNAGG'                ),
                               ('XapI'       , 'R|AATTY'                     ),
                               ('XbaI'       , 'T|CTAGA'                     ),
                               ('XceI'       , 'RCATG|Y'                     ),
                               ('XcmI'       , 'CCANNNNN|NNNNTGG'            ),
                               ('XhoI'       , 'C|TCGAG'                     ),
                               ('XhoII'      , 'R|GATCY'                     ),
                               ('XmaI'       , 'C|CCGGG'                     ),
                               ('XmaIII'     , 'C|GGCCG'                     ),
                               ('XmaJI'      , 'C|CTAGG'                     ),
                               ('XmiI'       , 'GT|MKAC'                     ),
                               ('XmnI'       , 'GAANN|NNTTC'                 ),
                               ('XspI'       , 'C|TAG'                       ),
                               ('YkrI'       , 'C|'                          ),
                               ('ZraI'       , 'GAC|GTC'                     ),
                               ('ZrmI'       , 'AGT|ACT'                     ),
                               ('Zsp2I'      , 'ATGCA|T'                     )])
