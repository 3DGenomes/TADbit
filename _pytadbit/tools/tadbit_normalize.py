"""

information needed

 - path working directory with parsed reads

"""
from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
from argparse                             import HelpFormatter
from os                                   import path, remove, system
from sys                                  import exc_info, stdout
from string                               import ascii_letters
from random                               import random
from shutil                               import copyfile, rmtree
from collections                          import OrderedDict, defaultdict
from pickle                               import dump, load, HIGHEST_PROTOCOL
from traceback                            import print_exc
from multiprocessing                      import cpu_count
import multiprocessing  as mu
import sqlite3 as lite
import time

from pysam                                import AlignmentFile
from numpy                                import nanmean, isnan, nansum, nanpercentile, seterr
from matplotlib                           import pyplot as plt

from pytadbit                             import load_hic_data_from_bam
from pytadbit.utils.sqlite_utils          import already_run, digest_parameters
from pytadbit.utils.sqlite_utils          import add_path, get_jobid, print_db, retry
from pytadbit.utils.file_handling         import mkdir
from pytadbit.mapping.analyze             import plot_distance_vs_interactions
from pytadbit.mapping.filter              import MASKED
from pytadbit.utils                       import printime
from pytadbit.parsers.hic_bam_parser      import print_progress
from pytadbit.parsers.hic_bam_parser      import filters_to_bin
from pytadbit.parsers.bed_parser          import parse_mappability_bedGraph
from pytadbit.utils.extraviews            import nicer
# from pytadbit.utils.hic_filtering         import filter_by_local_ratio
from pytadbit.utils.hic_filtering         import plot_filtering
# from pytadbit.utils.hic_filtering         import filter_by_zero_count
from pytadbit.utils.normalize_hic         import oneD
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES
from pytadbit.parsers.genome_parser       import parse_fasta, get_gc_content
from functools import reduce

# removes annoying message when normalizing...
seterr(invalid='ignore')


DESC = 'normalize Hi-C data and generates array of biases'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()

    param_hash = digest_parameters(opts)
    if opts.bam:
        mreads = path.realpath(opts.bam)
    else:
        mreads = path.join(opts.workdir, load_parameters_fromdb(opts))

    filter_exclude = opts.filter

    outdir = path.join(opts.workdir, '04_normalization')
    mkdir(outdir)

    mappability = gc_content = n_rsites = None
    if opts.normalization == 'oneD':
        if not opts.fasta:
            raise Exception('ERROR: missing path to FASTA for oneD normalization')
        if not opts.renz:
            raise Exception('ERROR: missing restriction enzyme name for oneD normalization')
        if not opts.mappability:
            raise Exception('ERROR: missing path to mappability for oneD normalization')
        bamfile = AlignmentFile(mreads, 'rb')
        refs = bamfile.references
        bamfile.close()

        # get genome sequence ~1 min
        printime('  - parsing FASTA')
        genome = parse_fasta(opts.fasta, verbose=False)

        fas = set(genome.keys())
        bam = set(refs)
        if fas - bam:
            print('WARNING: %d extra chromosomes in FASTA (removing them)' % (len(fas - bam)))
            if len(fas - bam) <= 50:
                print('\n'.join([('  - ' + c) for c in (fas - bam)]))
        if bam - fas:
            txt = ('\n'.join([('  - ' + c) for c in (bam - fas)])
                   if len(bam - fas) <= 50 else '')
            raise Exception('ERROR: %d extra chromosomes in BAM (remove them):\n%s\n' % (
                len(bam - fas), txt))
        refs = [crm for crm in refs if crm in genome]
        if len(refs) == 0:
            raise Exception("ERROR: chromosomes in FASTA different the ones"
                            " in BAM")

        # get mappability ~2 min
        printime('  - Parsing mappability')
        mappability = parse_mappability_bedGraph(
            opts.mappability, opts.reso,
            wanted_chrom=refs[0] if len(refs)==1 else None)
        # resize chomosomes
        for c in refs:
            if not c in mappability:
                mappability[c] = [float('nan')] * (len(refs) // opts.reso + 1)
            if len(mappability[c]) < len(refs) // opts.reso + 1:
                mappability[c] += [float('nan')] * (
                    (len(refs) // opts.reso + 1) - len(mappability[c]))
        # concatenates
        mappability = reduce(lambda x, y: x + y,
                             (mappability.get(c, []) for c in refs))

        printime('  - Computing GC content per bin (removing Ns)')
        gc_content = get_gc_content(genome, opts.reso, chromosomes=refs,
                                    n_cpus=opts.cpus)
        # pad mappability at the end if the size is close to gc_content
        if len(mappability)<len(gc_content) and len(mappability)/len(gc_content) > 0.95:
            mappability += [float('nan')] * (len(gc_content)-len(mappability))

        # compute r_sites ~30 sec
        # TODO: read from DB
        printime('  - Computing number of RE sites per bin (+/- 200 bp)')
        n_rsites  = []
        re_site = RESTRICTION_ENZYMES[opts.renz].replace('|', '')
        for crm in refs:
            for pos in range(200, len(genome[crm]) + 200, opts.reso):
                seq = genome[crm][pos-200:pos + opts.reso + 200]
                n_rsites.append(seq.count(re_site))

        ## CHECK TO BE REMOVED
        # out = open('tmp_mappability.txt', 'w')
        # i = 0
        # for crm in refs:
        #     for pos in xrange(len(genome[crm]) / opts.reso + 1):
        #         out.write('%s\t%d\t%d\t%f\n' % (crm, pos * opts.reso, pos * opts.reso + opts.reso, mappability[i]))
        #         i += 1`
        # out.close()
        # compute GC content ~30 sec
        # TODO: read from DB
    biases, decay, badcol, raw_cisprc, norm_cisprc = read_bam(
        mreads, filter_exclude, opts.reso, min_count=opts.min_count, sigma=2,
        factor=1, outdir=outdir, extra_out=param_hash, ncpus=opts.cpus,
        normalization=opts.normalization, mappability=mappability,
        p_fit=opts.p_fit, cg_content=gc_content, n_rsites=n_rsites,
        seed=opts.seed,
        normalize_only=opts.normalize_only, max_njobs=opts.max_njobs,
        extra_bads=opts.badcols, biases_path=opts.biases_path, 
        cis_limit=opts.cis_limit, trans_limit=opts.trans_limit, 
        min_ratio=opts.ratio_limit, cistrans_filter=opts.cistrans_filter)

    inter_vs_gcoord = path.join(opts.workdir, '04_normalization',
                                'interactions_vs_genomic-coords.png_%s_%s.png' % (
                                    opts.reso, param_hash))

    # get and plot decay
    if not opts.normalize_only:
        printime('  - Computing interaction decay vs genomic distance')
        (_, _, _), (a2, _, _), (_, _, _) = plot_distance_vs_interactions(
            decay, max_diff=10000, resolution=opts.reso, normalized=not opts.filter_only,
            savefig=inter_vs_gcoord)

        print ('    -> Decay slope 0.7-10 Mb\t%s' % a2)
    else:
        a2 = 0.

    printime('  - Saving biases and badcol columns')
    # biases
    bias_file = path.join(outdir, 'biases_%s_%s.pickle' % (
        nicer(opts.reso).replace(' ', ''), param_hash))
    out = open(bias_file, 'wb')

    dump({'biases'    : biases,
          'decay'     : decay,
          'badcol'    : badcol,
          'resolution': opts.reso}, out, HIGHEST_PROTOCOL)
    out.close()

    finish_time = time.localtime()

    try:
        save_to_db(opts, bias_file, mreads, len(badcol),
                   len(biases), raw_cisprc, norm_cisprc,
                   inter_vs_gcoord, a2, opts.filter,
                   launch_time, finish_time)
    except:
        # release lock anyway
        print_exc()
        try:
            remove(path.join(opts.workdir, '__lock_db'))
        except OSError:
            pass
        exit(1)


@retry(lite.OperationalError, tries=20, delay=2)
def save_to_db(opts, bias_file, mreads,
               nbad_columns, ncolumns, raw_cisprc, norm_cisprc,
               inter_vs_gcoord, a2, bam_filter,
               launch_time, finish_time):
    if 'tmpdb' in opts and opts.tmpdb:
        # check lock
        while path.exists(path.join(opts.workdir, '__lock_db')):
            time.sleep(0.5)
        # close lock
        open(path.join(opts.workdir, '__lock_db'), 'a').close()
        # tmp file
        dbfile = opts.tmpdb
        try: # to copy in case read1 was already mapped for example
            copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
        except IOError:
            pass
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='JOBs'""")
        if not cur.fetchall():
            cur.execute("""
            create table PATHs
               (Id integer primary key,
                JOBid int, Path text, Type text,
                unique (Path))""")
            cur.execute("""
            create table JOBs
               (Id integer primary key,
                Parameters text,
                Launch_time text,
                Finish_time text,
                Type text,
                Parameters_md5 text,
                unique (Parameters_md5))""")
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='NORMALIZE_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
            create table NORMALIZE_OUTPUTs
               (Id integer primary key,
                JOBid int,
                Input int,
                N_columns int,
                N_filtered int,
                BAM_filter int,
                Cis_percentage_Raw real,
                Cis_percentage_Norm real,
                Slope_700kb_10Mb real,
                Resolution int,
                Normalization text,
                Factor int,
                unique (JOBid))""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            param_hash = digest_parameters(opts, get_md5=True )
            cur.execute("""
            insert into JOBs
            (Id  , Parameters, Launch_time, Finish_time, Type , Parameters_md5)
            values
            (NULL,       '%s',        '%s',        '%s', 'Normalize',           '%s')
            """ % (parameters,
                   time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
                   time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass
        jobid = get_jobid(cur)
        add_path(cur, bias_file       , 'BIASES'     , jobid, opts.workdir)
        add_path(cur, inter_vs_gcoord , 'FIGURE'     , jobid, opts.workdir)
        if opts.bam:
            add_path(cur, path.realpath(opts.bam), 'EXT_2D_BAM' , jobid, opts.workdir)
        if opts.mappability:
            add_path(cur, path.realpath(opts.mappability), 'EXT_MAPPABILITY' , jobid, opts.workdir)
        if opts.fasta:
            add_path(cur, path.realpath(opts.fasta), 'EXT_FASTA' , jobid, opts.workdir)
        # get pathid of input
        cur.execute("select id from paths where path = '%s'" % (path.relpath(mreads, opts.workdir)))
        input_bed = cur.fetchall()[0][0]

        a2 = 0 if isnan(a2) else a2
        try:
            cur.execute("""
            insert into NORMALIZE_OUTPUTs
            (Id  , JOBid,     Input, N_columns,   N_filtered, BAM_filter, Cis_percentage_Raw, Cis_percentage_Norm, Slope_700kb_10Mb,   Resolution,      Normalization,      Factor)
            values
            (NULL,    %d,        %d,        %d,           %d,         %d,                 %f,                  %f,               %f,           %d,               '%s',          %f)
            """ % (jobid, input_bed,  ncolumns, nbad_columns, bam_filter,   100 * raw_cisprc,   100 * norm_cisprc,               a2,    opts.reso, opts.normalization, opts.factor))
        except lite.OperationalError:
            try:
                cur.execute("""
                insert into NORMALIZE_OUTPUTs
                (Id  , JOBid,     Input, N_columns,   N_filtered, BAM_filter,      Cis_percentage_Raw, Cis_percentage_Norm, Slope_700kb_10Mb,   Resolution,     Normalization,       Factor)
                values
                (NULL,    %d,        %d,        %d,           %d,         %d,                      %f,                  %f,               %f,           %d,               '%s',          %f)
                """ % (jobid, input_bed,  ncolumns, nbad_columns, bam_filter,        100 * raw_cisprc,   100 * norm_cisprc,               a2,    opts.reso, opts.normalization, opts.factor))
            except lite.OperationalError:
                print('WANRING: Normalized table not written!!!')

        print_db(cur, 'PATHs')
        print_db(cur, 'JOBs')
        try:
            print_db(cur, 'FILTER_OUTPUTs')
            print_db(cur, 'INTERSECTION_OUTPUTs')
            print_db(cur, 'MAPPED_INPUTs')
            print_db(cur, 'MAPPED_OUTPUTs')
            print_db(cur, 'PARSED_OUTPUTs')
            print_db(cur, 'FILTER_OUTPUTs')
        except lite.OperationalError:
            pass
        print_db(cur, 'NORMALIZE_OUTPUTs')
    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
    # release lock
    try:
        remove(path.join(opts.workdir, '__lock_db'))
    except OSError:
        pass

def load_parameters_fromdb(opts, what='bam'):
    if 'tmpdb' in opts and opts.tmpdb:
        dbfile = opts.tmpdb
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        if not opts.jobid:
            # get the JOBid of the parsing job
            cur.execute("""
            select distinct Id from JOBs
            where Type = 'Filter' or Type = 'Merge' or Type = 'Import'
            """)
            jobids = cur.fetchall()
            if len(jobids) > 1:
                raise Exception('ERROR: more than one possible input found, use'
                                '"tadbit describe" and select corresponding '
                                'jobid with --jobid')
            try:
                parse_jobid = jobids[0][0]
            except IndexError:
                raise Exception('ERROR: no BAM file found... is it filtered?')
        else:
            parse_jobid = opts.jobid
        # fetch path to parsed BED files
        if what == 'bed':
            cur.execute("""
            select distinct path from paths
            inner join filter_outputs on filter_outputs.pathid = paths.id
            where filter_outputs.name = 'valid-pairs' and paths.jobid = %s
            """ % parse_jobid)
        elif what == 'bam':
            cur.execute("""
            select distinct path from paths
            where paths.type = 'HIC_BAM' and paths.jobid = %s
            """ % parse_jobid)
        bam = cur.fetchall()[0][0]
        return bam


def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: SmartFormatter(prog, width=95,
                                                       max_help_position=27)

    oblopt = parser.add_argument_group('Required options')
    glopts = parser.add_argument_group('General options')
    bfiltr = parser.add_argument_group('Bin filtering options')
    rfiltr = parser.add_argument_group('Read filtering options')
    normpt = parser.add_argument_group('Normalization options')

    oblopt.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    oblopt.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int, required=True,
                        help='''resolution at which to output matrices''')

    glopts.add_argument('--bam', dest='bam', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to a TADbit-generated BAM file with
                        all reads (other wise the tool will guess from the
                        working directory database)''')

    glopts.add_argument('-j', '--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')

    glopts.add_argument('--max_njobs', dest='max_njobs', metavar="INT",
                        action='store', default=100, type=int,
                        help='''[%(default)s] Define maximum number of jobs
                        for reading BAM file (set to higher numbers for large files
                        and low RAM memory).''')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument("-C", "--cpus", dest="cpus", type=int,
                        default=cpu_count(), help='''[%(default)s] Maximum
                        number of CPU cores  available in the execution host.
                        If higher than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    glopts.add_argument('--normalize_only', dest='normalize_only', action='store_true',
                        default=False,
                        help=
                        'skip calculation of Cis-percentage and decay')

    glopts.add_argument('--noX', action='store_true', help='no display server (X screen)')

    normpt.add_argument('--normalization', dest='normalization', metavar="STR",
                        action='store', default='Vanilla', type=str,
                        choices=['Vanilla', 'ICE', 'SQRT', 'oneD', 'custom'],
                        help='''[%(default)s] normalization(s) to apply.
                        Order matters. Choices: %(choices)s''')

    normpt.add_argument('--biases_path', dest='biases_path', type=str,
                        default=None, help='''biases file to compute decay.
                        REQUIRED with "custom" normalization. Format: single
                        column with header''')

    normpt.add_argument('--mappability', dest='mappability', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''R|Path to mappability bedGraph file, required for oneD normalization.
Mappability file can be generated with GEM (example from the genomic FASTA file hg38.fa):\n
     gem-indexer -i hg38.fa -o hg38
     gem-mappability -I hg38.gem -l 50 -o hg38.50mer -T 8
     gem-2-wig -I hg38.gem -i hg38.50mer.mappability -o hg38.50mer
     wigToBigWig hg38.50mer.wig hg38.50mer.sizes hg38.50mer.bw
     bigWigToBedGraph hg38.50mer.bw  hg38.50mer.bedGraph\n''')

    normpt.add_argument('--fasta', dest='fasta', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''Path to FASTA file with genome sequence, to compute
                        GC content and number of restriction sites per bin.
                        Required for oneD normalization''')

    normpt.add_argument('--renz', dest='renz', metavar="STR",
                        type=str, required=False,
                        help='''restriction enzyme name(s). Required for oneD
                        normalization''')

    normpt.add_argument('--factor', dest='factor', metavar="NUM",
                        action='store', default=1, type=float,
                        help='''[%(default)s] target mean value of a cell after
                        normalization (can be used to weight experiments before
                        merging)''')

    normpt.add_argument('--prop_data', dest='p_fit', metavar="FLOAT",
                        action='store', default=None, type=float,
                        help=('''[1] Only for oneD normalization: proportion of
                        data to be used in fitting (for very large datasets).
                        Number between 0 and 1.'''))

    normpt.add_argument('--seed', dest='seed', metavar="INT",
                        action='store', default=1, type=int,
                        help=('''[%(default)s] Only for oneD normalization: seed number for
                        the random picking of data when using the "prop_data"
                        parameter'''))

    bfiltr.add_argument('--min_count', dest='min_count', metavar="INT",
                        action='store', default=None, type=float,
                        help=('''[%(default)s] minimum number of reads mapped to
                        a bin (recommended value could be 2500). If set this
                        option overrides the perc_zero filtering... This option is
                        slightly slower.'''))

    bfiltr.add_argument('--cis_limit', dest='cis_limit', action='store',
                        default=None, type=int,
                        help='''Maximum distance in bins at which to consider an 
                        interaction cis for the filtering. By default it is the 
                        number of bins corresponding to 1Mb''')

    bfiltr.add_argument('--trans_limit', dest='trans_limit', action='store',
                        default=None, type=int,
                        help='''Maximum distance in bins at which to consider an 
                        interaction trans for the filtering. By default it is 
                        five times the cis_limit (if also default, it would 
                        correspond to the number of bins needed to reach 5Mb).''')

    bfiltr.add_argument('--ratio_limit', dest='ratio_limit', action='store',
                        default=1.0, type=float,
                        help='''[%(default)s] Minimum cis/trans (as defined with 
                        cis_limit and trans_limit parameters) to filter out bins.''')

    bfiltr.add_argument('--cistrans_filter', dest='cistrans_filter', action='store_true',
                        default=False,
                        help='''filter using cis-trans ratio.''')

    # bfiltr.add_argument('--min_perc', dest='min_perc', metavar="INT",
    #                     action='store', default=None, type=float,
    #                     help=('''[%(default)s] lower percentile from which
    #                     consider bins as good.'''))

    # bfiltr.add_argument('--max_perc', dest='max_perc', metavar="INT",
    #                     action='store', default=None, type=float,
    #                     help=('''[%(default)s] upper percentile until which
    #                     consider bins as good.'''))

    bfiltr.add_argument('--filter_only', dest='filter_only', action='store_true',
                        default=False,
                        help='skip normalization')

    bfiltr.add_argument('-B', '--badcols', dest='badcols', nargs='+',
                        type=str, default=None, metavar='CHR:POS1-POS2',
                        help=('extra regions to be added to bad-columns (in'
                              'genomic position). e.g.: --badcols'
                              ' 1:150000000-160000000 2:1200000-1300000'))

    rfiltr.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                        choices = list(range(1, 11)),
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))

    rfiltr.add_argument('--valid', dest='only_valid', action='store_true',
                        default=False,
                        help='input BAM file contains only valid pairs (already filtered).')


def check_options(opts):
    mkdir(opts.workdir)

    # transform filtering reads option
    opts.filter = filters_to_bin(opts.filter)

    # check custom normalization
    if opts.normalization=='custom':
        if not opts.biases_path:
            raise IOError('ERROR: biases file required for "custom" normalization.')
        elif not path.exists(opts.biases_path):
            raise IOError('ERROR: biases not found at path: %s' % opts.biases_path)

    # check filtering options
    if not opts.min_count and not opts.cistrans_filter:
        opts.cistrans_filter = True

    # check resume
    if not path.exists(opts.workdir):
        raise IOError('ERROR: workdir not found.')

    # for lustre file system....
    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)
        try:
            copyfile(path.join(opts.workdir, 'trace.db'), opts.tmpdb)
        except IOError:
            pass

    # number of cpus
    if opts.cpus == 0:
        opts.cpus = cpu_count()
    else:
        opts.cpus = min(opts.cpus, cpu_count())

    # check if job already run using md5 digestion of parameters
    try:
        if already_run(opts):
            if 'tmpdb' in opts and opts.tmpdb:
                remove(path.join(dbdir, dbfile))
            exit('WARNING: exact same job already computed, see JOBs table above')
    except IOError:  # new working directory
        pass


def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)


################################################################################
## TODO: This should be handled in the hic bam parser

def read_bam_frag_valid(inbam, filter_exclude, all_bins, sections,
                        resolution, outdir, extra_out,region, start, end,
                        next_position=1, last_position=None):
    
    if last_position is None:
        last_position = next_position * 5
    
    bamfile = AlignmentFile(inbam, 'rb')
    refs = bamfile.references
    try:
        dico = {}
        for r in bamfile.fetch(region=region,
                               start=start - (1 if start else 0), end=end,  # coords starts at 0
                               multiple_iterators=True):
            crm1 = r.reference_name
            pos1 = r.reference_start + 1
            crm2 = refs[r.mrnm]
            pos2 = r.mpos + 1
            try:
                pos1 = sections[(crm1, pos1 // resolution)]
                pos2 = sections[(crm2, pos2 // resolution)]
            except KeyError:
                continue  # not in the subset matrix we want
            try:
                dico[(pos1, pos2)] += 1
            except KeyError:
                dico[(pos1, pos2)] = 1
        cisprc = {}
        for (i, j), v in dico.items():
            if all_bins[i][0] == all_bins[j][0]:
                diff = abs(j - i)
                try:
                    cisprc[i][0] += v  # add to cis interactions
                    cisprc[i][1] += v  # add to total interactions
                    if diff <= next_position:
                        cisprc[i][2] += v  # add to total interactions
                    elif diff <= last_position:
                        cisprc[i][3] += v  # add to total interactions
                except KeyError:
                    if diff <= next_position:
                        cisprc[i] = [v, v, v, 0]
                    elif diff <= last_position:
                        cisprc[i] = [v, v, 0, v]
                    else:
                        cisprc[i] = [v, v, 0, 0]
            else:
                try:
                    cisprc[i][1] += v  # only add to total interactions
                except KeyError:
                    cisprc[i] = [0, v, 0, 0]
        out = open(path.join(outdir,
                             'tmp_%s:%d-%d_%s.pickle' % (region, start, end, extra_out)), 'wb')
        dump(dico, out, HIGHEST_PROTOCOL)
        out.close()
        out = open(path.join(outdir,
                             'tmp_bins_%s:%d-%d_%s.pickle' % (region, start, end, extra_out)), 'wb')
        dump(cisprc, out, HIGHEST_PROTOCOL)
        out.close()
    except Exception as e:
        exc_type, exc_obj, exc_tb = exc_info()
        fname = path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(e)
        print(exc_type, fname, exc_tb.tb_lineno)


def read_bam_frag_filter(inbam, filter_exclude, all_bins, sections,
                         resolution, outdir, extra_out,region, start, end,
                         next_position=1, last_position=None):

    if last_position is None:
        last_position = next_position * 5

    bamfile = AlignmentFile(inbam, 'rb')
    refs = bamfile.references
    try:
        dico = {}
        for r in bamfile.fetch(region=region,
                               start=start - (1 if start else 0), end=end,  # coords starts at 0
                               multiple_iterators=True):
            if r.flag & filter_exclude:
                continue
            crm1 = r.reference_name
            pos1 = r.reference_start + 1
            crm2 = refs[r.mrnm]
            pos2 = r.mpos + 1
            try:
                pos1 = sections[(crm1, pos1 // resolution)]
                pos2 = sections[(crm2, pos2 // resolution)]
            except KeyError:
                continue  # not in the subset matrix we want
            try:
                dico[(pos1, pos2)] += 1
            except KeyError:
                dico[(pos1, pos2)] = 1
        cisprc = {}
        for (i, j), v in dico.items():
            if all_bins[i][0] == all_bins[j][0]:  # same chromosome
                diff = abs(j - i)
                try:
                    cisprc[i][0] += v  # add to cis interactions
                    cisprc[i][1] += v  # add to total interactions
                    if diff <= next_position:
                        cisprc[i][2] += v  # add to total interactions
                    elif diff <= last_position:
                        cisprc[i][3] += v  # add to total interactions
                except KeyError:
                    if diff <= next_position:
                        cisprc[i] = [v, v, v, 0]
                    elif diff <= last_position:
                        cisprc[i] = [v, v, 0, v]
                    else:
                        cisprc[i] = [v, v, 0, 0]
            else:
                try:
                    cisprc[i][1] += v  # only add to total interactions
                except KeyError:
                    cisprc[i] = [0, v ,0, 0]
        out = open(path.join(outdir,
                             'tmp_%s:%d-%d_%s.pickle' % (region, start, end, extra_out)), 'wb')
        dump(dico, out, HIGHEST_PROTOCOL)
        out.close()
        out = open(path.join(outdir, 'tmp_bins_%s:%d-%d_%s.pickle' % (
            region, start, end, extra_out)), 'wb')
        dump(cisprc, out, HIGHEST_PROTOCOL)
        out.close()
    except Exception as e:
        exc_type, exc_obj, exc_tb = exc_info()
        fname = path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(e)
        print(exc_type, fname, exc_tb.tb_lineno)


def read_bam(inbam, filter_exclude, resolution, min_count=2500, biases_path='',
             normalization='Vanilla', mappability=None, n_rsites=None,
             cg_content=None, sigma=2, ncpus=8, factor=1, outdir='.', seed=1,
             extra_out='', only_valid=False, normalize_only=False, p_fit=None,
             max_njobs=100, extra_bads=None, 
             cis_limit=1, trans_limit=5, min_ratio=1.0, cistrans_filter=False):
    bamfile = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(list(zip(bamfile.references,
                               [x // resolution + 1 for x in bamfile.lengths])))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]
    bins = []
    for crm in sections:
        len_crm = sections[crm]
        bins.extend([(crm, i) for i in range(len_crm)])

    start_bin = 0
    end_bin   = len(bins)
    total     = len(bins)

    regs = []
    begs = []
    ends = []
    njobs = min(total, max_njobs) + 1
    nbins = total // njobs + 1
    for i in range(start_bin, end_bin, nbins):
        if i + nbins > end_bin:  # make sure that we stop
            nbins = end_bin - i
        try:
            (crm1, beg1), (crm2, end2) = bins[i], bins[i + nbins - 1]
        except IndexError:
            try:
                (crm1, beg1), (crm2, end2) = bins[i], bins[-1]
            except IndexError:
                break
        if crm1 != crm2:
            end1 = sections[crm1]
            beg2 = 0
            regs.append(crm1)
            regs.append(crm2)
            begs.append(beg1 * resolution)
            begs.append(beg2 * resolution)
            ends.append(end1 * resolution + resolution)  # last nt included
            ends.append(end2 * resolution + resolution - 1)  # last nt not included (overlap with next window)
        else:
            regs.append(crm1)
            begs.append(beg1 * resolution)
            ends.append(end2 * resolution + resolution - 1)
    ends[-1] += 1  # last nucleotide included

    # print '\n'.join(['%s %d %d' % (a, b, c) for a, b, c in zip(regs, begs, ends)])
    printime('  - Parsing BAM (%d chunks)' % (len(regs)))
    # define limits for cis and trans interactions if not given
    if cis_limit is None:
        cis_limit = int(1000000 / resolution)
    print('      -> cis interactions are defined as being bellow {}'.format(
        nicer(cis_limit * resolution)))
    if trans_limit is None:
        trans_limit = cis_limit * 5
    print('      -> trans interactions are defined as being bellow {}'.format(
        nicer(trans_limit * resolution)))

    bins_dict = dict([(j, i) for i, j in enumerate(bins)])
    pool = mu.Pool(ncpus)
    procs = []
    read_bam_frag = read_bam_frag_valid if only_valid else read_bam_frag_filter

    num_regs = len(regs)
    regs = map(lambda r: r.replace("|", "_"), regs)

    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        num_regs += 1
        procs.append(pool.apply_async(
            read_bam_frag, args=(inbam, filter_exclude, bins, bins_dict,
                                 resolution, outdir, extra_out,
                                 region, start, end, cis_limit, trans_limit)))
    pool.close()
    print_progress(procs)
    pool.join()
    ## COLLECT RESULTS
    cisprc = {}
    printime('  - Collecting cis and total interactions per bin (%d chunks)' % num_regs)
    stdout.write('     ')
    for countbin, (region, start, end) in enumerate(zip(regs, begs, ends)):
        if not countbin % 10 and countbin:
            stdout.write(' ')
        if not countbin % 50 and countbin:
            stdout.write(' %9s\n     ' % ('%s/%s' % (countbin , num_regs)))
        stdout.write('.')
        stdout.flush()

        fname = path.join(outdir,
                          'tmp_bins_%s:%d-%d_%s.pickle' % (region, start, end, extra_out))
        tmp_cisprc = load(open(fname,'rb'))
        system('rm -f %s' % fname)
        cisprc.update(tmp_cisprc)
    stdout.write('\n')

    # get cis/trans ratio
    for k in cisprc:
        try:
            cisprc[k][3] = cisprc[k][2] / cisprc[k][3]
        except ZeroDivisionError:
            cisprc[k][3] = 0

    # BIN FILTERINGS
    printime('  - Removing columns with too few or too much interactions')
    
    # define filter for minimum interactions per bin
    if cistrans_filter:
        if min_count is None:
            min_count = nanpercentile(
                [cisprc[k][2] for k in range(total) 
                if cisprc.get(k, [0, 0, 0, 0])[3] < min_ratio 
                and cisprc.get(k, [0, 0, 0, 0])[2] >= 1], 95)  # harcoded parameter we are filtering
                                                                # out bins with no interactions in cis

        print('      -> too few interactions defined as less than %9d '
               'interactions' % (min_count))
        badcol = dict((k, True) for k in range(total) 
                    if cisprc.get(k, [0, 0, 0, 0])[3] < min_ratio
                    or cisprc[k][2] < min_count)
        print('      -> removed %d columns of %d (%.1f%%)' % (
            len(badcol), total, float(len(badcol)) / total * 100))
    else:
    # if len(bamfile.references) == 1 and min_count is None:
    #     raise Exception("ERROR: only one chromosome can't filter by "
    #                     "cis-percentage, set min_count instead")
    # elif min_count is None and len(bamfile.references) > 1:
        # badcol = filter_by_cis_percentage(
        #     cisprc, sigma=sigma, verbose=True, min_perc=min_perc, max_perc=max_perc,
        #     size=total, savefig=None)

        print('      -> too few interactions defined as less than %9d '
               'interactions' % (min_count))
        badcol = {}
        countL = 0
        countZ = 0
        for c in range(total):
            if cisprc.get(c, [0, 0, 0, 0])[1] < min_count:
                badcol[c] = cisprc.get(c, [0, 0, 0, 0])[1]
                countL += 1
                if not c in cisprc:
                    countZ += 1
        print('      -> removed %d columns (%d/%d null/high counts) of %d (%.1f%%)' % (
            len(badcol), countZ, countL, total, float(len(badcol)) / total * 100))

    # Plot
    plot_filtering(dict((k, cisprc[k][2]) for k in cisprc), 
                   dict((k, cisprc[k][3]) for k in cisprc), total, min_count, min_ratio, 
                   path.join(outdir, 'filtering_summary_plot_{}_{}.png'.format(nicer(resolution, sep=''),
                    extra_out)),
                   base_position=0, next_position=cis_limit, last_position=trans_limit, resolution=resolution,
                   legend='Filtered {} of {} bins'.format(len(badcol), total))

    # no mappability will result in NaNs, better to filter out these columns
    if mappability:
        badcol.update((i, True) for i, m in enumerate(mappability) if not m)

    # add manually columns to bad columns
    if extra_bads:
        removed_manually = 0
        for ebc in extra_bads:
            c, ebc = ebc.split(':')
            b, e = list(map(int, ebc.split('-')))
            b = b // resolution + section_pos[c][0]
            e = e // resolution + section_pos[c][0]
            removed_manually += (e - b)
            badcol.update(dict((p, 'manual') for p in range(b, e)))
        printime('  - Removed %d columns manually.' % removed_manually)
    raw_cisprc = sum(float(cisprc[k][0]) / cisprc[k][1]
                     for k in cisprc if not k in badcol) / (len(cisprc) - len(badcol))

    printime('  - Rescaling sum of interactions per bins')
    size = len(bins)
    biases = [float('nan') if k in badcol else cisprc.get(k, [0, 1., 0, 0])[1]
              for k in range(size)]

    if normalization == 'ICE':
        printime('  - ICE normalization')
        hic_data = load_hic_data_from_bam(
            inbam, resolution, filter_exclude=filter_exclude,
            tmpdir=outdir, ncpus=ncpus, nchunks=max_njobs)
        hic_data.bads = badcol
        hic_data.normalize_hic(iterations=100, max_dev=0.000001)
        biases = hic_data.bias.copy()
        del(hic_data)
    elif normalization == 'Vanilla':
        printime('  - Vanilla normalization')
        mean_col = nanmean(biases)
        biases   = dict((k, b / mean_col * mean_col**0.5)
                        for k, b in enumerate(biases))
    elif normalization == 'SQRT':
        printime('  - Vanilla-SQRT normalization')
        biases = [b**0.5 for b in biases]
        mean_col = nanmean(biases)
        biases   = dict((k, b / mean_col * mean_col**0.5)
                        for k, b in enumerate(biases))
    elif normalization == 'oneD':
        printime('  - oneD normalization')
        if len(set([len(biases), len(mappability), len(n_rsites), len(cg_content)])) > 1:
            print("biases", "mappability", "n_rsites", "cg_content")
            print(len(biases), len(mappability), len(n_rsites), len(cg_content))
            raise Exception('Error: not all arrays have the same size')
        tmp_oneD = path.join(outdir,'tmp_oneD_%s' % (extra_out))
        mkdir(tmp_oneD)
        biases = oneD(tmp_dir=tmp_oneD, p_fit=p_fit, tot=biases, map=mappability,
                      res=n_rsites, cg=cg_content, seed=seed)
        biases = dict((k, b) for k, b in enumerate(biases))
        rmtree(tmp_oneD)
    elif normalization == 'custom':
        n_pos = 0
        biases = {}
        print('Using provided biases...')
        with open(biases_path, 'r') as r:
            next(r)
            for line in r:
                if line[0] == 'N':
                    #b = float('nan')
                    badcol[n_pos] = 0
                    biases[n_pos] = float('nan')
                else:
                    b = float(line)
                    if b == 0:
                        badcol[n_pos] = 0
                        biases[n_pos] = float('nan')
                    else:
                        biases[n_pos] = b
                n_pos += 1
        for add in range(max(biases.keys()), total + 1):
            biases[add] = float('nan')
    else:
        raise NotImplementedError('ERROR: method %s not implemented' %
                                  normalization)

    # collect subset-matrices and write genomic one
    # out = open(os.path.join(outdir,
    #                         'hicdata_%s.abc' % (nicer(resolution).replace(' ', ''))), 'w')
    printime('  - Getting sum of normalized bins')
    pool = mu.Pool(ncpus)
    procs = []
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        fname = path.join(outdir,
                          'tmp_%s:%d-%d_%s.pickle' % (region, start, end, extra_out))
        procs.append(pool.apply_async(sum_nrm_matrix,
                                      args=(fname, biases,)))
    pool.close()
    print_progress(procs)
    pool.join()

    # to correct biases
    sumnrm = sum(p.get() for p in procs)

    target = (sumnrm / float(size * size * factor))**0.5
    biases = dict([(b, biases[b] * target) for b in biases])

    if not normalize_only:
        printime('  - Computing Cis percentage')
        # Calculate Cis percentage

        pool = mu.Pool(ncpus)
        procs = []
        for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
            fname = path.join(outdir,
                              'tmp_%s:%d-%d_%s.pickle' % (region, start, end, extra_out))
            procs.append(pool.apply_async(get_cis_perc,
                                          args=(fname, biases, badcol, bins)))
        pool.close()
        print_progress(procs)
        pool.join()

        # collect results
        cis = total = 0
        for proc in procs:
            c, t = proc.get()
            cis += c
            total += t
        norm_cisprc = float(cis) / total
        print('    * Cis-percentage: %.1f%%' % (norm_cisprc * 100))
    else:
        norm_cisprc = 0.

    printime('  - Rescaling decay')
    # normalize decay by size of the diagonal, and by Vanilla correction
    # (all cells must still be equals to 1 in average)

    pool = mu.Pool(ncpus)
    procs = []
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        fname = path.join(outdir,
                          'tmp_%s:%d-%d_%s.pickle' % (region, start, end, extra_out))
        procs.append(pool.apply_async(sum_dec_matrix,
                                      args=(fname, biases, badcol, bins)))
    pool.close()
    print_progress(procs)
    pool.join()

    # collect results
    nrmdec = {}
    rawdec = {}
    for proc in procs:
        tmpnrm, tmpraw = proc.get()
        for c, d in tmpnrm.items():
            for k, v in d.items():
                try:
                    nrmdec[c][k] += v
                    rawdec[c][k] += tmpraw[c][k]
                except KeyError:
                    try:
                        nrmdec[c][k]  = v
                        rawdec[c][k] = tmpraw[c][k]
                    except KeyError:
                        nrmdec[c] = {k: v}
                        rawdec[c] = {k: tmpraw[c][k]}
    # count the number of cells per diagonal
    # TODO: parallelize
    # find largest chromosome
    len_crms = dict((c, section_pos[c][1] - section_pos[c][0]) for c in section_pos)
    # initialize dictionary
    ndiags = dict((c, dict((k, 0) for k in range(len_crms[c]))) for c in sections)
    for crm in section_pos:
        beg_chr, end_chr = section_pos[crm][0], section_pos[crm][1]
        chr_size = end_chr - beg_chr
        thesebads = [b for b in badcol if beg_chr <= b <= end_chr]
        for dist in range(1, chr_size):
            ndiags[crm][dist] += chr_size - dist
            # from this we remove bad columns
            # bad columns will only affect if they are at least as distant from
            # a border as the distance between the longest diagonal and the
            # current diagonal.
            bad_diag = set()  # 2 bad rows can point to the same bad cell in diagonal
            maxp = end_chr - dist
            minp = beg_chr + dist
            for b in thesebads:
                if b < maxp:  # not inclusive!!
                    bad_diag.add(b)
                if b >= minp:
                    bad_diag.add(b - dist)
            ndiags[crm][dist] -= len(bad_diag)
        # different behavior for longest diagonal:
        ndiags[crm][0] += chr_size - sum(beg_chr <= b < end_chr for b in thesebads)

    # normalize sum per diagonal by total number of cells in diagonal
    signal_to_noise = 0.05
    min_n = signal_to_noise ** -2. # equals 400 when default
    for crm in sections:
        if not crm in nrmdec:
            nrmdec[crm] = {}
            rawdec[crm] = {}
        tmpdec = 0  # store count by diagonal
        tmpsum = 0  # store count by diagonal
        ndiag  = 0
        val    = 0
        previous = [] # store diagonals to be summed in case not reaching the minimum
        for k in ndiags[crm]:
            tmpdec += nrmdec[crm].get(k, 0.)
            tmpsum += rawdec[crm].get(k, 0.)
            previous.append(k)
            if tmpsum > min_n:
                ndiag = sum(ndiags[crm][k] for k in previous)
                val = tmpdec  # backup of tmpdec kept for last ones outside the loop
                try:
                    ratio = val / ndiag
                    for l in previous:
                        nrmdec[crm][l] = ratio
                except ZeroDivisionError:  # all columns at this distance are "bad"
                    pass
                previous = []
                tmpdec = 0
                tmpsum = 0
        # last ones we average with previous result
        if  len(previous) == len(ndiags[crm]):
            nrmdec[crm] = {}
        elif tmpsum < min_n:
            ndiag += sum(ndiags[crm][k] for k in previous)
            val += tmpdec
            try:
                ratio = val / ndiag
                for k in previous:
                    nrmdec[crm][k] = ratio
            except ZeroDivisionError:  # all columns at this distance are "bad"
                pass
    return biases, nrmdec, badcol, raw_cisprc, norm_cisprc


def sum_dec_matrix(fname, biases, badcol, bins):
    dico = load(open(fname,'rb'))
    rawdec = {}
    nrmdec = {}
    for (i, j), v in dico.items():
        if i < j:
            continue
        # different chromosome
        c = bins[i][0]
        if c != bins[j][0]:
            continue
        if i in badcol or j in badcol:
            continue
        k = i - j
        val = v / biases[i] / biases[j]
        try:
            nrmdec[c][k] += val
            rawdec[c][k] += v
        except KeyError:
            try:
                nrmdec[c][k] = val
                rawdec[c][k] = v
            except KeyError:
                nrmdec[c] = {k: val}
                rawdec[c] = {k: v}
    system('rm -f %s' % (fname))
    return nrmdec, rawdec


def get_cis_perc(fname, biases, badcol, bins):
    dico = load(open(fname,'rb'))
    cis = total = 0
    for (i, j), v in dico.items():
        if i <= j:
            continue
        if i in badcol or j in badcol:
            continue
        val = v / biases[i] / biases[j]
        # same chromosome
        if bins[i][0] == bins[j][0]:
            cis += val
        total += val
    return cis, total


def sum_nrm_matrix(fname, biases):
    dico = load(open(fname,'rb'))
    sumnrm = nansum([v / biases[i] / biases[j]
                     for (i, j), v in dico.items()])
    return sumnrm


class SmartFormatter(HelpFormatter):
    """
    https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
    """
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return HelpFormatter._split_lines(self, text, width)
