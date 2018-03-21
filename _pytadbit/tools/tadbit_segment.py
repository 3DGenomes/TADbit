"""

information needed

 - path working directory with parsed reads

"""
from argparse                       import HelpFormatter
from os                             import path, remove
from shutil                         import copyfile
from string                         import ascii_letters
from random                         import random
from warnings                       import warn
from cPickle                        import load
from multiprocessing                import cpu_count
from traceback                      import print_exc
import sqlite3 as lite
import time

from pytadbit                       import load_hic_data_from_bam
from pytadbit                       import tadbit
from pytadbit.utils.sqlite_utils    import already_run, digest_parameters
from pytadbit.utils.sqlite_utils    import add_path, get_jobid, print_db
from pytadbit.utils.file_handling   import mkdir
from pytadbit.parsers.tad_parser    import parse_tads
from pytadbit.parsers.genome_parser import parse_fasta, get_gc_content
from pytadbit.mapping.filter        import MASKED


DESC = 'Finds TAD or compartment segmentation in Hi-C data.'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()
    param_hash = digest_parameters(opts, get_md5=True)

    if opts.nosql:
        biases = opts.biases
        mreads = opts.mreads
        inputs = []
    else:
        biases, mreads, biases_id, mreads_id = load_parameters_fromdb(opts)
        inputs = [biases_id, mreads_id]
        # store path ids to be saved in database
        mreads = path.join(opts.workdir, mreads)
        biases = path.join(opts.workdir, biases)

    reso   = opts.reso

    mkdir(path.join(opts.workdir, '06_segmentation'))

    print 'loading %s \n    at resolution %s' % (mreads, nice(reso))
    region = None
    if opts.crms and len(opts.crms) == 1:
        region = opts.crms[0]
    hic_data = load_hic_data_from_bam(mreads, reso, ncpus=opts.cpus,
                                      region=region,
                                      biases=None if opts.all_bins else biases,
                                      filter_exclude=opts.filter)

    # compartments
    cmp_result = {}
    richA_stats = {}
    firsts = {}
    if not opts.only_tads:
        print 'Searching compartments'
        cmprt_dir = path.join(opts.workdir, '06_segmentation',
                              'compartments_%s' % (nice(reso)))
        mkdir(cmprt_dir)
        if opts.fasta:
            print '  - Computing GC content to label compartments'
            rich_in_A = get_gc_content(parse_fasta(opts.fasta, chr_filter=opts.crms), reso,
                                       chromosomes=opts.crms,
                                       by_chrom=True, n_cpus=opts.cpus)
        elif opts.rich_in_A:
            rich_in_A = opts.rich_in_A
        else:
            rich_in_A = None
        n_evs = opts.n_evs if opts.n_evs > 0 else 3
        firsts, richA_stats = hic_data.find_compartments(
            crms=opts.crms, savefig=cmprt_dir, verbose=True, suffix=param_hash,
            rich_in_A=rich_in_A, show_compartment_labels=rich_in_A is not None,
            savecorr=cmprt_dir if opts.savecorr else None,
            max_ev=n_evs,
            ev_index=opts.ev_index,
            vmin=None if opts.fix_corr_scale else 'auto',
            vmax=None if opts.fix_corr_scale else 'auto')

        for ncrm, crm in enumerate(opts.crms or hic_data.chromosomes):
            if not crm in firsts:
                continue
            ev_file = open(path.join(
                cmprt_dir, '%s_EigVect%d_%s.tsv' % (
                    crm, opts.ev_index[ncrm] if opts.ev_index else 1,
                    param_hash)), 'w')
            ev_file.write('# %s\n' % ('\t'.join(
                'EV_%d (%.4f)' % (i, v)
                for i, v in enumerate(firsts[crm][0], 1))))
            ev_file.write('\n'.join(['\t'.join([str(v) for v in vs])
                                     for vs in zip(*firsts[crm][1])]))
            ev_file.close()

        for ncrm, crm in enumerate(opts.crms or hic_data.chromosomes):
            cmprt_file1 = path.join(cmprt_dir, '%s_%s.tsv' % (crm, param_hash))
            cmprt_file2 = path.join(cmprt_dir, '%s_EigVect%d_%s.tsv' % (
                crm, opts.ev_index[ncrm] if opts.ev_index else 1, param_hash))
            cmprt_image = path.join(cmprt_dir, '%s_EV%d_%s.%s' % (
                crm, opts.ev_index[ncrm] if opts.ev_index else 1,
                param_hash, opts.format))
            if opts.savecorr:
                cormat_file = path.join(cmprt_dir, '%s_corr-matrix%s.tsv' %
                                       (crm, param_hash))
            else:
                cormat_file = None
            hic_data.write_compartments(cmprt_file1, chroms=[crm])
            cmp_result[crm] = {'path_cmprt1': cmprt_file1,
                               'path_cmprt2': cmprt_file2,
                               'path_cormat': cormat_file,
                               'image_cmprt': cmprt_image,
                               'num' : len(hic_data.compartments[crm])}

    # TADs
    tad_result = {}
    if not opts.only_compartments:
        print 'Searching TADs'
        tad_dir = path.join(opts.workdir, '06_segmentation',
                             'tads_%s' % (nice(reso)))
        mkdir(tad_dir)
        for crm in hic_data.chromosomes:
            if opts.crms and not crm in opts.crms:
                continue
            print '  - %s' % crm
            matrix = hic_data.get_matrix(focus=crm)
            beg, end = hic_data.section_pos[crm]
            size = len(matrix)
            if size < 10:
                print "     Chromosome too short (%d bins), skipping..." % size
                continue
            # transform bad column in chromosome referential
            if hic_data.bads:
                to_rm = tuple([1 if i in hic_data.bads else 0 for i in xrange(beg, end)])
            else:
                to_rm = None
            # maximum size of a TAD
            max_tad_size = (size - 1) if opts.max_tad_size is None else opts.max_tad_size
            result = tadbit([matrix], remove=to_rm,
                            n_cpus=opts.cpus, verbose=opts.verbose,
                            max_tad_size=max_tad_size,
                            no_heuristic=False)

            # use normalization to compute height on TADs called
            if opts.all_bins:
                if opts.nosql:
                    biases = load(open(biases))
                else:
                    biases = load(open(path.join(opts.workdir, biases)))
                hic_data.bads = biases['badcol']
                hic_data.bias = biases['biases']
            tads = load_tad_height(result, size, beg, end, hic_data)
            table = ''
            table += '%s\t%s\t%s\t%s\t%s\n' % ('#', 'start', 'end', 'score', 'density')
            for tad in tads:
                table += '%s\t%s\t%s\t%s%s\n' % (
                    tad, int(tads[tad]['start'] + 1), int(tads[tad]['end'] + 1),
                    abs(tads[tad]['score']), '\t%s' % (round(
                        float(tads[tad]['height']), 3)))
            out_tad = path.join(tad_dir, '%s_%s.tsv' % (crm, param_hash))
            out = open(out_tad, 'w')
            out.write(table)
            out.close()
            tad_result[crm] = {'path' : out_tad,
                               'num': len(tads)}

    finish_time = time.localtime()

    if not opts.nosql:
        try:
            save_to_db(opts, cmp_result, tad_result, reso, inputs,
                       richA_stats, firsts, param_hash,
                       launch_time, finish_time)
        except:
            # release lock anyway
            print_exc()
            try:
                remove(path.join(opts.workdir, '__lock_db'))
            except OSError:
                pass
            exit(1)


def save_to_db(opts, cmp_result, tad_result, reso, inputs,
               richA_stats, firsts, param_hash,
               launch_time, finish_time):
    if 'tmpdb' in opts and opts.tmpdb:
        # check lock
        while path.exists(path.join(opts.workdir, '__lock_db')):
            time.sleep(0.5)
        # close lock
        open(path.join(opts.workdir, '__lock_db'), 'a').close()
        # tmp file
        dbfile = opts.tmpdb
        copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='SEGMENT_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
            create table SEGMENT_OUTPUTs
               (Id integer primary key,
                JOBid int,
                Inputs text,
                TADs int,
                Compartments int,
                richA_corr real,
                EV_index int,
                EValue real,
                Chromosome text,
                Resolution int)""")
        try:
            parameters = digest_parameters(opts, get_md5=False, extra=['fasta'])
            cur.execute("""
            insert into JOBs
            (Id  , Parameters, Launch_time, Finish_time, Type , Parameters_md5)
            values
            (NULL,       '%s',        '%s',        '%s', 'Segment',       '%s')
            """ % (parameters,
                   time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
                   time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass
        jobid = get_jobid(cur)
        for ncrm, crm in enumerate(max(cmp_result.keys(), tad_result.keys(), key=len)):
            if crm in cmp_result:
                add_path(cur, cmp_result[crm]['path_cmprt1'], 'COMPARTMENT',
                         jobid, opts.workdir)
                add_path(cur, cmp_result[crm]['path_cmprt2'], 'COMPARTMENT',
                         jobid, opts.workdir)
                add_path(cur, cmp_result[crm]['image_cmprt'], 'FIGURE',
                         jobid, opts.workdir)
                if opts.savecorr:
                    add_path(cur, cmp_result[crm]['path_cormat'],
                             'CROSS_CORR_MAT', jobid, opts.workdir)
            if crm in tad_result:
                add_path(cur, tad_result[crm]['path'], 'TAD', jobid, opts.workdir)
            if opts.rich_in_A:
                add_path(cur, opts.rich_in_A, 'BED', jobid, opts.workdir)

            if crm in firsts:
                evalue = firsts[crm][0][(opts.ev_index[ncrm] - 1) if opts.ev_index else 0]
                eindex = opts.ev_index[ncrm] if opts.ev_index else 1
            else:
                evalue = 'NULL'
                eindex = 'NULL'
            try:
                cur.execute("""
                insert into SEGMENT_OUTPUTs
                (Id  , JOBid, Inputs, TADs, Compartments, richA_corr, EV_index, EValue, Chromosome, Resolution)
                values
                (NULL,    %d,   '%s',   %s,           %s,         %s,       %s,     %s,       '%s',         %d)
                """ % (jobid,
                       ','.join([str(i) for i in inputs]),
                       tad_result[crm]['num'] if crm in tad_result else 'NULL',
                       cmp_result[crm]['num'] if crm in cmp_result else 'NULL',
                       (richA_stats[crm] if crm in richA_stats
                        and richA_stats[crm] is not None else 'NULL'),
                       eindex, evalue, crm, reso))
            except lite.OperationalError:  # TODO: remove this
                print_exc()
                try:
                    cur.execute("alter table SEGMENT_OUTPUTs add column 'richA_corr' 'real'")
                except:
                    pass
                try:
                    cur.execute("alter table SEGMENT_OUTPUTs add column 'EValue' 'real'")
                except:
                    pass
                try:
                    cur.execute("alter table SEGMENT_OUTPUTs add column 'EV_index', 'int'")
                except:
                    pass
                cur.execute("""
                insert into SEGMENT_OUTPUTs
                (Id  , JOBid, Inputs, TADs, Compartments, richA_corr, EV_index, EValue, Chromosome, Resolution)
                values
                (NULL,    %d,   '%s',   %d,           %d,         %s,       %s,     %s,       '%s',         %d)
                """ % (jobid,
                       ','.join([str(i) for i in inputs]),
                       tad_result[crm]['num'] if crm in tad_result else 0,
                       cmp_result[crm]['num'] if crm in cmp_result else 0,
                       (richA_stats[crm] if crm in richA_stats
                        and richA_stats[crm] is not None else 'NULL'),
                       eindex, evalue, crm, reso))
        print_db(cur, 'PATHs')
        print_db(cur, 'JOBs')
        print_db(cur, 'SEGMENT_OUTPUTs')
    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
        # release lock
        remove(path.join(opts.workdir, '__lock_db'))


def load_parameters_fromdb(opts):
    if 'tmpdb' in opts and opts.tmpdb:
        dbfile = opts.tmpdb
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        if not opts.jobid:
            # get the JOBid of the parsing job
            try:
                cur.execute("""
                select distinct Id from JOBs
                where Type = 'Normalize'
                """)
                jobids = cur.fetchall()
                parse_jobid = jobids[0][0]
            except IndexError:
                cur.execute("""
                select distinct Id from JOBs
                where Type = '%s'
                """ % ('Filter'))
                jobids = cur.fetchall()
                try:
                    parse_jobid = jobids[0][0]
                except IndexError:
                    parse_jobid = 1
            if len(jobids) > 1:
                found = False
                cur.execute("""
                select distinct JOBid from NORMALIZE_OUTPUTs
                where Resolution = %d
                """ % (opts.reso))
                jobs = cur.fetchall()
                try:
                    parse_jobid = jobs[0][0]
                    found = True
                except IndexError:
                    found = False
                if len(jobs ) > 1:
                    found = False
                if not found:
                    raise Exception('ERROR: more than one possible input found, use'
                                    '"tadbit describe" and select corresponding '
                                    'jobid with --jobid')

        # fetch path to BAM files
        # try:
        biases = mreads = reso = None
        try:
            cur.execute("""
            select distinct Path, PATHs.id from PATHs
            where paths.jobid = %s and paths.Type = 'BIASES'
            """ % parse_jobid)
            biases, biases_id = cur.fetchall()[0]

            cur.execute("""
            select distinct Path, PATHs.id from PATHs
            inner join NORMALIZE_OUTPUTs on PATHs.Id = NORMALIZE_OUTPUTs.Input
            where NORMALIZE_OUTPUTs.JOBid = %d;
            """ % parse_jobid)
            mreads, mreads_id = cur.fetchall()[0]

            cur.execute("""
            select distinct Resolution from NORMALIZE_OUTPUTs
            where NORMALIZE_OUTPUTs.JOBid = %d;
            """ % parse_jobid)
            reso = int(cur.fetchall()[0][0])
            if reso != opts.reso:
                warn('WARNING: input resolution does not match '
                     'the one of the precomputed normalization')
        except IndexError:
            raise Exception('ERROR: normalization not found')
        return biases, mreads, biases_id, mreads_id


def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')
    cmopts = parser.add_argument_group('Compartment calling options')
    tdopts = parser.add_argument_group('TAD calling options')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument('--nosql', dest='nosql',
                        action='store_true', default=False,
                        help='do not load/store data from/in sqlite database')

    glopts.add_argument('--all_bins', dest='all_bins',
                        action='store_true', default=False,
                        help='skip the filtering of bins for detection of TADs')

    glopts.add_argument('--mreads', dest='mreads', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path valid-pairs file (TADbit BAM format)''')

    glopts.add_argument('--biases',   dest='biases', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to file with precalculated biases by
                        columns''')

    glopts.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int, required=True,
                        help='''resolution at which to output matrices''')

    glopts.add_argument('--norm_matrix', dest='norm_matrix', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to a matrix file with normalized read
                        counts''')

    glopts.add_argument('--raw_matrix', dest='raw_matrix', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to a matrix file with raw read
                        counts''')

    glopts.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                        choices = range(1, 11),
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))

    cmopts.add_argument('--rich_in_A', dest='rich_in_A', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to a BAD or bedGraph file with list of
                        protein coding gene or other active epigenetic mark,
                        to be used to label compartments instead of using
                        the relative interaction count.''')

    cmopts.add_argument('--fasta', dest='fasta', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''Path to fasta file with genome sequence, to compute
                        GC content and use it to label compartments''')

    cmopts.add_argument('--savecorr', dest='savecorr', action='store_true', default=False,
                        help='''Save correlation matrix used to predict compartments.''')

    cmopts.add_argument('--fix_corr_scale', dest='fix_corr_scale', action='store_true',
                        default=False,
                        help='''Correlation matrix plot scaled between correlation 1 and -1
                        instead of maximum observed values.''')

    cmopts.add_argument('--format', dest='format', action='store', default='png',
                        type=str, help='[%(default)s] file format for figures')

    cmopts.add_argument('--n_evs', dest='n_evs', metavar="INT",
                        action='store', default=3, type=int,
                        help='''[%(default)s] Number of eigenvectors to
                        store. if "-1" all eigenvectors will be calculated''')

    cmopts.add_argument('--ev_index', dest='ev_index', nargs='+',
                        type=int, metavar='INT', default=None,
                        help="""list of indexes of eigenvectors capturing
                        compartments signal (one index per chromosome, in the
                        same order as chromosomes in fasta file). Example
                        picking the first eigenvector for all chromosomes but
                        for chromosome 3:
                        '--ev_index 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1""")

    glopts.add_argument('--only_compartments', dest='only_compartments',
                        action='store_true', default=False,
                        help='''search A/B compartments using first eigen vector
                        of the correlation matrix''')

    glopts.add_argument('--only_tads', dest='only_tads',
                        action='store_true', default=False,
                        help='''search TAD boundaries break-point detection
                        algorithm''')

    glopts.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='''print more messages''')

    glopts.add_argument('-j', '--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')

    glopts.add_argument('-c', '--chromosomes', dest='crms', metavar="STR",
                        action='store', default=None, type=str, nargs='+',
                        help='''Name of the chromosomes on which to search
                        for TADs or compartments.''')

    tdopts.add_argument('--max_tad_size', dest='max_tad_size', metavar="INT",
                        action='store', default=None, type=int,
                        help='''an integer defining the maximum size of TAD. Default
                        defines it as the number of rows/columns''')

    glopts.add_argument("-C", "--cpu", dest="cpus", type=int,
                        default=cpu_count(), help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    glopts.add_argument('--force', dest='force', action='store_true',
                        default=False,
                        help='overwrite previously run job')

    parser.add_argument_group(glopts)


def check_options(opts):

    # check resume
    if not path.exists(opts.workdir):
        raise IOError('ERROR: %s does not exists' % opts.workdir)

    # number of cpus
    if opts.cpus == 0:
        opts.cpus = cpu_count()
    else:
        opts.cpus = min(opts.cpus, cpu_count())

    if opts.rich_in_A and opts.fasta:
        raise Exception('ERROR: should choose one of FASTA or rich_in_A')

    # rich_in_A
    if opts.fasta:
        if opts.rich_in_A:
            raise Exception(('ERROR: if you input a FASTA file, GC content will'
                             'will be used as "rich in A" metric to infer '
                             'compartments.'))
        opts.rich_in_A = opts.fasta

    # N EVs
    if opts.ev_index:
        if max(opts.ev_index) > opts.n_evs:
            warn('WARNING: increasing number of calculated eigenvectors to %d, '
                 'to match the u=input eigenvectors indices' % max(opts.ev_index))
            opts.n_evs = max(opts.ev_index)

    # tmp folder
    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)
        copyfile(path.join(opts.workdir, 'trace.db'), opts.tmpdb)

    if already_run(opts) and not opts.force:
        if 'tmpdb' in opts and opts.tmpdb:
            remove(path.join(dbdir, dbfile))
        exit('WARNING: exact same job already computed, see JOBs table above')


def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)


def load_tad_height(tad_def, size, beg, end, hic_data):
    bias, zeros = hic_data.bias, hic_data.bads
    if bias:
        norm = lambda i, j: bias[i] * bias[j]
    else:
        norm = lambda i, j: 1 # Non-normalized height, keep in mind!
    tads, _ = parse_tads(tad_def)
    diags = []
    for k in xrange(1, size):
        try:
            diags.append(sum(
                hic_data[i, i + k] / norm(i + k, i)
                for i in xrange(beg, end - k) if not i in zeros and not i + k in zeros
                ) / float(sum(1 for i in range(beg, end - k)
                              if not i in zeros and not i + k in zeros)))
        except ZeroDivisionError:
            diags.append(0.)
    for tad in tads:
        start, final = (int(tads[tad]['start']) + 1,
                        int(tads[tad]['end']) + 1)
        matrix = sum(
            hic_data[i, j] / norm(j, i)
            for i in xrange(beg + start - 1, beg + final - 1) if not i in zeros
            for j in xrange(i + 1          , beg + final - 1) if not j in zeros)
        try:
            height = float(matrix) / sum(
                [diags[i - 1] * (final - start - i)
                 for i in xrange(1, final - start)])
        except ZeroDivisionError:
            height = 0.
        tads[tad]['height'] = height
    return tads
