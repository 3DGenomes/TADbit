"""

information needed

 - path working directory with parsed reads

"""
from argparse                             import HelpFormatter
from os                                   import path, remove, system
from sys                                  import exc_info, stdout
from string                               import ascii_letters
from random                               import random
from shutil                               import copyfile, rmtree
from collections                          import OrderedDict
from cPickle                              import dump, load
# from warnings                             import filterwarnings
from multiprocessing                      import cpu_count
import multiprocessing  as mu
import sqlite3 as lite
import time

from pysam                                import AlignmentFile
from numpy                                import nanmean, isnan, nansum, seterr

from pytadbit.utils.sqlite_utils          import already_run, digest_parameters
from pytadbit.utils.sqlite_utils          import add_path, get_jobid, print_db
from pytadbit.utils.file_handling         import mkdir
from pytadbit.mapping.analyze             import plot_distance_vs_interactions
from pytadbit.mapping.filter              import MASKED
from pytadbit.parsers.hic_bam_parser      import printime, print_progress
from pytadbit.parsers.hic_bam_parser      import filters_to_bin
from pytadbit.utils.extraviews            import nicer
from pytadbit.utils.hic_filtering         import filter_by_cis_percentage
from pytadbit.utils.normalize_hic         import oneD
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES
from pytadbit.parsers.genome_parser       import parse_fasta, get_gc_content

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
            print 'WARNING: %d extra chromosomes in FASTA (removing them)' % (len(fas - bam))
            if len(fas - bam) <= 50:
                print '\n'.join([('  - ' + c) for c in (fas - bam)])
        if bam - fas:
            txt = ('\n'.join([('  - ' + c) for c in (bam - fas)])
                   if len(bam - fas) <= 50 else '')
            raise Exception('ERROR: %d extra chromosomes in BAM (remove them):\n%s\n' % (
                len(bam - fas), txt))
        refs = [crm for crm in refs if crm in genome]
        if len(refs) == 0:
            raise Exception("ERROR: chromosomes in FASTA different the ones in BAM")

        # get mappability ~2 min
        printime('  - Parsing mappability')
        fh = open(opts.mappability)
        mappability = []
        ordered_crm = True
        line = fh.next()
        crmM, begM, endM, val = line.split()
        for crm in refs:
            try:
                while line and crm != crmM:
                    ordered_crm = False
                    line = fh.next()
                    crmM, begM, endM, val = line.split()
            except EOFError:
                break
            for begB in xrange(0, len(genome[crm]), opts.reso):
                endB = begB + opts.reso
                tmp = 0
                try:
                    while True:
                        crmM, begM, endM, val = line.split()
                        begM = int(begM)
                        endM = int(endM)
                        val = float(val)
                        weight = min(endM, endB) - max(begM, begB)
                        if weight < 0 or crm != crmM:
                            break
                        if endM > endB:
                            tmp += weight * val
                            break
                        tmp += weight * val
                        line = fh.next()
                except EOFError:
                    continue
                except StopIteration:
                    pass
                mappability.append(tmp / opts.reso)
            if not ordered_crm:
                fh.seek(0, 0)

        printime('  - Computing GC content per bin (removing Ns)')
        gc_content = get_gc_content(genome, opts.reso, chromosomes=refs,
                                    n_cpus=opts.cpus)
        # compute r_sites ~30 sec
        # TODO: read from DB
        printime('  - Computing number of RE sites per bin (+/- 200 bp)')
        n_rsites  = []
        re_site = RESTRICTION_ENZYMES[opts.renz].replace('|', '')
        for crm in refs:
            for pos in xrange(200, len(genome[crm]) + 200, opts.reso):
                seq = genome[crm][pos-200:pos + opts.reso + 200]
                n_rsites.append(seq.count(re_site))

        ## CHECK TO BE REMOVED
        # out = open('tmp_mappability.txt', 'w')
        # i = 0
        # for crm in refs:
        #     for pos in xrange(len(genome[crm]) / opts.reso + 1):
        #         out.write('%s\t%d\t%d\t%f\n' % (crm, pos * opts.reso, pos * opts.reso + opts.reso, mappability[i]))
        #         i += 1
        # out.close()
        # compute GC content ~30 sec
        # TODO: read from DB
    biases, decay, badcol, raw_cisprc, norm_cisprc = read_bam(
        mreads, filter_exclude, opts.reso, min_count=opts.min_count, sigma=2,
        factor=1, outdir=outdir, extra_out=param_hash, ncpus=opts.cpus,
        normalization=opts.normalization, mappability=mappability,
        cg_content=gc_content, n_rsites=n_rsites, min_perc=opts.min_perc, max_perc=opts.max_perc,
        normalize_only=opts.normalize_only, max_njobs=opts.max_njobs)

    bad_col_image = path.join(outdir, 'filtered_bins_%s_%s.png' % (
        nicer(opts.reso).replace(' ', ''), param_hash))

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
    out = open(bias_file, 'w')

    dump({'biases'    : biases,
          'decay'     : decay,
          'badcol'    : badcol,
          'resolution': opts.reso}, out)
    out.close()

    finish_time = time.localtime()

    save_to_db(opts, bias_file, mreads, bad_col_image,
               len(badcol), len(biases), raw_cisprc, norm_cisprc,
               inter_vs_gcoord, a2, opts.filter,
               launch_time, finish_time)


def save_to_db(opts, bias_file, mreads, bad_col_image,
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
        add_path(cur, bad_col_image   , 'FIGURE'     , jobid, opts.workdir)
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
                print 'WANRING: Normalized table not written!!!'

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
            where Type = 'Filter' or Type = 'Merge'
            """)
            jobids = cur.fetchall()
            if len(jobids) > 1:
                raise Exception('ERROR: more than one possible input found, use'
                                '"tadbit describe" and select corresponding '
                                'jobid with --jobid')
            parse_jobid = jobids[0][0]
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

    glopts.add_argument('--force', dest='force', action='store_true',
                        default=False,
                        help='overwrite previously run job')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument("-C", "--cpus", dest="cpus", type=int,
                        default=0, help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    glopts.add_argument('--normalize_only', dest='normalize_only', action='store_true',
                        default=False,
                        help='skip calculation of Cis-percentage and decay')

    normpt.add_argument('--normalization', dest='normalization', metavar="STR",
                        action='store', default='Vanilla', type=str,
                        choices=['Vanilla', 'oneD'],
                        help='''[%(default)s] normalization(s) to apply.
                        Order matters. Choices: [%(choices)s]''')

    normpt.add_argument('--mappability', dest='mappability', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''R|Path to mappability bedGraph file, required for oneD normalization.
Mappability file can be generated with GEM (example from the genomic fasta file hg38.fa):\n
     gem-indexer -i hg38.fa -o hg38
     gem-mappability -I hg38.gem -l 50 -o hg38.50mer -T 8
     gem-2-wig -I hg38.gem -i hg38.50mer.mappability -o hg38.50mer
     wigToBigWig hg38.50mer.wig hg38.50mer.sizes hg38.50mer.bw
     bigWigToBedGraph hg38.50mer.bw  hg38.50mer.bedGraph\n''')

    normpt.add_argument('--fasta', dest='fasta', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''Path to fasta file with genome sequence, to compute
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

    bfiltr.add_argument('--perc_zeros', dest='perc_zeros', metavar="FLOAT",
                        action='store', default=95, type=float,
                        help=('[%(default)s%%] maximum percentage of zeroes '
                              'allowed per column.'))

    bfiltr.add_argument('--min_count', dest='min_count', metavar="INT",
                        action='store', default=None, type=float,
                        help=('''[%(default)s] minimum number of reads mapped to
                        a bin (recommended value could be 2500). If set this
                        option overrides the perc_zero filtering... This option is
                        slightly slower.'''))

    bfiltr.add_argument('--min_perc', dest='min_perc', metavar="INT",
                        action='store', default=None, type=float,
                        help=('''[%(default)s] lower percentile from which
                        consider bins as good.'''))

    bfiltr.add_argument('--max_perc', dest='max_perc', metavar="INT",
                        action='store', default=None, type=float,
                        help=('''[%(default)s] upper percentile until which
                        consider bins as good.'''))

    bfiltr.add_argument('--filter_only', dest='filter_only', action='store_true',
                        default=False,
                        help='skip normalization')

    bfiltr.add_argument('--fast_filter', dest='fast_filter', action='store_true',
                        default=False,
                        help='''only filter according to the percentage of zero
                        count or minimum count of reads''')

    rfiltr.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                        choices = range(1, 11),
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
    if already_run(opts):
        if 'tmpdb' in opts and opts.tmpdb:
            remove(path.join(dbdir, dbfile))
        exit('WARNING: exact same job already computed, see JOBs table above')

def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)


################################################################################
## TODO: This should be handled in the hic bam parser

def read_bam_frag_valid(inbam, filter_exclude, all_bins, sections,
                  resolution, outdir, extra_out,region, start, end):
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
                pos1 = sections[(crm1, pos1 / resolution)]
                pos2 = sections[(crm2, pos2 / resolution)]
            except KeyError:
                continue  # not in the subset matrix we want
            try:
                dico[(pos1, pos2)] += 1
            except KeyError:
                dico[(pos1, pos2)] = 1
        cisprc = {}
        for (i, j), v in dico.iteritems():
            # out.write('%d\t%d\t%d\n' % (i, j, v))
            try:
                if all_bins[i][0] == all_bins[j][0]:
                    cisprc[i][0] += v
                cisprc[i][1] += v
            except KeyError:
                if all_bins[i][0] == all_bins[j][0]:
                    cisprc[i] = [v, v]
                else:
                    cisprc[i] = [0, v]
        out = open(path.join(outdir,
                             'tmp_%s:%d-%d_%s.pickle' % (region, start, end, extra_out)), 'w')
        dump(dico, out)
        out.close()
        out = open(path.join(outdir,
                             'tmp_bins_%s:%d-%d_%s.pickle' % (region, start, end, extra_out)), 'w')
        dump(cisprc, out)
        out.close()
    except Exception, e:
        exc_type, exc_obj, exc_tb = exc_info()
        fname = path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print e
        print(exc_type, fname, exc_tb.tb_lineno)


def read_bam_frag_filter(inbam, filter_exclude, all_bins, sections,
                         resolution, outdir, extra_out,region, start, end):
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
                pos1 = sections[(crm1, pos1 / resolution)]
                pos2 = sections[(crm2, pos2 / resolution)]
            except KeyError:
                continue  # not in the subset matrix we want
            try:
                dico[(pos1, pos2)] += 1
            except KeyError:
                dico[(pos1, pos2)] = 1
        cisprc = {}
        for (i, j), v in dico.iteritems():
            # out.write('%d\t%d\t%d\n' % (i, j, v))
            try:
                if all_bins[i][0] == all_bins[j][0]:
                    cisprc[i][0] += v
                cisprc[i][1] += v
            except KeyError:
                if all_bins[i][0] == all_bins[j][0]:
                    cisprc[i] = [v, v]
                else:
                    cisprc[i] = [0, v]
        out = open(path.join(outdir,
                             'tmp_%s:%d-%d_%s.pickle' % (region, start, end, extra_out)), 'w')
        dump(dico, out)
        out.close()
        out = open(path.join(outdir,
                             'tmp_bins_%s:%d-%d_%s.pickle' % (region, start, end, extra_out)), 'w')
        dump(cisprc, out)
        out.close()
    except Exception, e:
        exc_type, exc_obj, exc_tb = exc_info()
        fname = path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print e
        print(exc_type, fname, exc_tb.tb_lineno)


def read_bam(inbam, filter_exclude, resolution, min_count=2500,
             normalization='Vanilla', mappability=None, n_rsites=None,
             cg_content=None, sigma=2, ncpus=8, factor=1, outdir='.',
             extra_out='', only_valid=False, normalize_only=False,
             max_njobs=100, min_perc=None, max_perc=None):
    bamfile = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references,
                               [x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]
    bins = []
    for crm in sections:
        len_crm = sections[crm]
        bins.extend([(crm, i) for i in xrange(len_crm)])

    start_bin = 0
    end_bin   = len(bins)
    total = len(bins)

    total = end_bin - start_bin
    regs = []
    begs = []
    ends = []
    njobs = min(total, max_njobs) + 1
    nbins = total / njobs + 1
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
    bins_dict = dict([(j, i) for i, j in enumerate(bins)])
    pool = mu.Pool(ncpus)
    procs = []
    read_bam_frag = read_bam_frag_valid if only_valid else read_bam_frag_filter
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        procs.append(pool.apply_async(
            read_bam_frag, args=(inbam, filter_exclude, bins, bins_dict,
                                 resolution, outdir, extra_out,
                                 region, start, end,)))
    pool.close()
    print_progress(procs)
    pool.join()
    ## COLLECT RESULTS
    cisprc = {}
    printime('  - Collecting cis and total interactions per bin (%d chunks)' % (len(regs)))
    stdout.write('     ')
    for countbin, (region, start, end) in enumerate(zip(regs, begs, ends)):
        if not countbin % 10 and countbin:
            stdout.write(' ')
        if not countbin % 50 and countbin:
            stdout.write(' %9s\n     ' % ('%s/%s' % (countbin , len(regs))))
        stdout.write('.')
        stdout.flush()

        fname = path.join(outdir,
                          'tmp_bins_%s:%d-%d_%s.pickle' % (region, start, end, extra_out))
        tmp_cisprc = load(open(fname))
        system('rm -f %s' % fname)
        cisprc.update(tmp_cisprc)
    print '%s %9s\n' % (' ' * (54 - (countbin % 50) - (countbin % 50) / 10),
                        '%s/%s' % (len(regs),len(regs)))

    printime('  - Removing columns with too few or too much interactions')
    if len(bamfile.references) == 1 and min_count is None:
        raise Exception("ERROR: only one chromosome can't filter by "
                        "cis-percentage, set min_count instead")
    elif min_count is None and len(bamfile.references) > 1:
        badcol = filter_by_cis_percentage(
            cisprc, sigma=sigma, verbose=True, min_perc=min_perc, max_perc=max_perc,
            savefig=path.join(outdir, 'filtered_bins_%s_%s.png' % (
                nicer(resolution).replace(' ', ''), extra_out)))
    else:
        print '      -> too few interactions defined as less than %9d interactions' % (
            min_count)
        badcol = {}
        countL = 0
        countZ = 0
        for c in xrange(total):
            if cisprc.get(c, [0, 0])[1] < min_count:
                badcol[c] = cisprc.get(c, [0, 0])[1]
                countL += 1
                if not c in cisprc:
                    countZ += 1
        print '      -> removed %d columns (%d/%d null/high counts) of %d (%.1f%%)' % (
            len(badcol), countZ, countL, total, float(len(badcol)) / total * 100)

    raw_cisprc = sum(float(cisprc[k][0]) / cisprc[k][1]
                     for k in cisprc if not k in badcol) / (len(cisprc) - len(badcol))

    printime('  - Rescaling sum of interactions per bins')
    size = len(bins)
    biases = [float('nan') if k in badcol else cisprc.get(k, [0, 1.])[1]
              for k in xrange(size)]

    if normalization=='Vanilla':
        printime('  - Vanilla normalization')
        mean_col = nanmean(biases)
        biases   = dict((k, b / mean_col * mean_col**0.5)
                        for k, b in enumerate(biases))
    elif normalization=='oneD':
        printime('  - oneD normalization')
        if len(set([len(biases), len(mappability), len(n_rsites), len(cg_content)])) > 1:
            print "biases", "mappability", "n_rsites", "cg_content"
            print len(biases), len(mappability), len(n_rsites), len(cg_content)
            raise Exception('Error: not all arrays have the same size')
        tmp_oneD = path.join(outdir,'tmp_oneD_%s' % (extra_out))
        mkdir(tmp_oneD)
        biases = oneD(tmp_dir=tmp_oneD, tot=biases, map=mappability, res=n_rsites, cg=cg_content)
        biases = dict((k, b) for k, b in enumerate(biases))
        rmtree(tmp_oneD)
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
        print '    * Cis-percentage: %.1f%%' % (norm_cisprc * 100)
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
        for c, d in tmpnrm.iteritems():
            for k, v in d.iteritems():
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
    ndiags = dict((c, dict((k, 0) for k in xrange(len_crms[c]))) for c in sections)
    for crm in section_pos:
        beg_chr, end_chr = section_pos[crm][0], section_pos[crm][1]
        chr_size = end_chr - beg_chr
        thesebads = [b for b in badcol if beg_chr <= b <= end_chr]
        for dist in xrange(1, chr_size):
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
                    for k in previous:
                        nrmdec[crm][k] = ratio
                except ZeroDivisionError:  # all columns at this distance are "bad"
                    pass
                previous = []
                tmpdec = 0
                tmpsum = 0
        # last ones we average with previous result
        if tmpsum < min_n:
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
    dico = load(open(fname))
    rawdec = {}
    nrmdec = {}
    for (i, j), v in dico.iteritems():
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
    dico = load(open(fname))
    cis = total = 0
    for (i, j), v in dico.iteritems():
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
    dico = load(open(fname))
    sumnrm = nansum([v / biases[i] / biases[j]
                     for (i, j), v in dico.iteritems()])
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
