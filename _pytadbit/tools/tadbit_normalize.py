"""

information needed

 - path working directory with parsed reads

"""
from argparse                     import HelpFormatter
from pytadbit                     import load_hic_data_from_reads
from pytadbit.utils.sqlite_utils  import already_run, digest_parameters
from pytadbit.utils.sqlite_utils  import add_path, get_jobid, print_db
from pytadbit.utils.file_handling import mkdir
from pytadbit.mapping.analyze     import plot_distance_vs_interactions, hic_map
from os                           import path
import sqlite3 as lite
import time

DESC = 'normalize Hi-C data and write results to file as matrices'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()


    param_hash = digest_parameters(opts)
    if opts.bed:
        mreads = path.realpath(opts.bed)
    else:
        mreads = path.join(opts.workdir, load_parameters_fromdb(opts))

    print 'loading', mreads
    hic_data = load_hic_data_from_reads(mreads, opts.reso)

    print 'Get poor bins...'
    try:
        hic_data.filter_columns(perc_zero=opts.perc_zeros, draw_hist=True,
                                by_mean=not opts.fast_filter, savefig=path.join(
                                    opts.workdir, '04_normalization',
                                    'bad_columns_%d_.pdf' % opts.perc_zeros) if
                                not opts.fast_filter else None)
    except ValueError:
        hic_data.filter_columns(perc_zero=100, draw_hist=True,
                                by_mean=not opts.fast_filter, savefig=path.join(
                                    opts.workdir, '04_normalization',
                                    'bad_columns_%d_.pdf' % opts.perc_zeros) if
                                not opts.fast_filter else None)

    mkdir(path.join(opts.workdir, '04_normalization'))

    # bad columns
    bad_columns_file = path.join(opts.workdir, '04_normalization',
                                 'bad_columns_%s.tsv' %param_hash)
    out_bad = open(bad_columns_file, 'w')
    out_bad.write('\n'.join([str(i) for i in hic_data.bads.keys()]))
    out_bad.close()

    # Identify biases
    print 'Get biases using ICE...'
    hic_data.normalize_hic(silent=False, max_dev=0.1, iterations=0,
                           factor=opts.factor)

    print 'Getting cis/trans...'
    cis_trans = hic_data.cis_trans_ratio(normalized=True)
        
    print 'Cis/Trans ratio of normalized matrix including the diagonal', cis_trans

    # Plot genomic distance vs interactions
    print 'Plot genomic distance vs interactions...'
    inter_vs_gcoord = path.join(opts.workdir, '04_normalization',
                                'interactions_vs_genomic-coords.pdf_%s.pdf' % (
                                    param_hash))
    (_, _, _), (a2, _, _), (_, _, _) = plot_distance_vs_interactions(
        hic_data, max_diff=10000, resolution=opts.reso, normalized=True,
        savefig=inter_vs_gcoord)
    
    print 'Decay slope 0.7-10 Mb\t%s' % a2

    inter_vs_gcoord = intra_dir_fig = intra_dir_txt = inter_dir_fig = inter_dir_txt = genome_map_fig = genome_map_txt = None
    
    if "intra" in opts.keep:
        print "  Saving normalized intra chromosomal matrix..."
        if opts.only_txt:
            intra_dir_fig = path.join(opts.workdir, '04_normalization',
                                      'intra_chromosome_map_images_%s' % (param_hash))
        else:
            intra_dir_fig = None
        intra_dir_txt = path.join(opts.workdir, '04_normalization',
                                  'intra_chromosome_map_matrices_%s' % (param_hash))
        hic_map(hic_data, normalized=True, by_chrom='intra', cmap='jet',
                name = path.split(opts.workdir)[-1],
                savefig =intra_dir_fig, savedata=intra_dir_txt)

    if "inter" in opts.keep:
        print "  Saving normalized inter chromosomal matrix..."
        if opts.only_txt:
            inter_dir_fig = path.join(opts.workdir, '04_normalization',
                                      'inter_chromosome_map_images_%s' % (param_hash))
        else:
            inter_dir_fig = None
        inter_dir_txt = path.join(opts.workdir, '04_normalization',
                                  'inter_chromosome_map_matrices_%s' % (param_hash))
        hic_map(hic_data, normalized=True, by_chrom='inter', cmap='jet',
                name = path.split(opts.workdir)[-1],
                savefig =inter_dir_fig, savedata=intra_dir_txt)

    if "genome" in opts.keep:
        print "  Saving normalized genomic matrix..."
        if opts.only_txt:
            genome_map_fig = path.join(opts.workdir, '04_normalization',
                                       'genomic_maps_%s.pdf' % (param_hash))
        else:
            genome_map_fig = None
        genome_map_txt = path.join(opts.workdir, '04_normalization',
                                   'genomic_maps_%s.tsv' % (param_hash))
        hic_map(hic_data, normalized=True, cmap='jet',
                name = path.split(opts.workdir)[-1],
                savefig =genome_map_fig, savedata=genome_map_txt)

    finish_time = time.localtime()

    save_to_db (opts, cis_trans, a2, bad_columns_file,
                inter_vs_gcoord, intra_dir_fig, intra_dir_txt, inter_dir_fig,
                inter_dir_txt, genome_map_fig, genome_map_txt,
                launch_time, finish_time)

def save_to_db(opts, cis_trans, a2, bad_columns_file,
               inter_vs_gcoord, intra_dir_fig, intra_dir_txt, inter_dir_fig,
               inter_dir_txt, genome_map_fig, genome_map_txt,
               launch_time, finish_time):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='NORMALIZE_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
            create table NORMALIZE_OUTPUTs
               (Id integer primary key,
                JOBid int,
                CisTrans real,
                Slope_700kb_10Mb real,
                Resolution int,
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
        add_path(cur, bad_columns_file, 'TSV', jobid)
        add_path(cur, inter_vs_gcoord, 'FIGURE', jobid)
        add_path(cur, intra_dir_fig, 'FIGURE_DIR', jobid)
        add_path(cur, intra_dir_txt, 'MATRiX_DIR', jobid)
        add_path(cur, inter_dir_fig, 'FIGURE_DIR', jobid)
        add_path(cur, inter_dir_txt, 'MATRIX_DIR', jobid)
        add_path(cur, genome_map_fig, 'FIGURE', jobid)
        add_path(cur, genome_map_txt, 'MATRIX', jobid)

        cur.execute("""
        insert into NORMALIZE_OUTPUTs
        (Id  , JOBid, CisTrans, Slope_700kb_10Mb, Resolution,      Factor)
        values
        (NULL,    %d,        %f,             %f,          %d,          %f)
        """ % (jobid, cis_trans,             a2,   opts.reso, opts.factor))
        print_db(cur, 'MAPPED_INPUTs')
        print_db(cur, 'PATHs')
        print_db(cur, 'MAPPED_OUTPUTs')
        print_db(cur, 'PARSED_OUTPUTs')
        print_db(cur, 'JOBs')
        print_db(cur, 'INTERSECTION_OUTPUTs')        
        print_db(cur, 'FILTER_OUTPUTs')
        print_db(cur, 'NORMALIZE_OUTPUTs')

def load_parameters_fromdb(opts):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        if not opts.jobid:
            # get the JOBid of the parsing job
            cur.execute("""
            select distinct Id from JOBs
            where Type = 'Filter'
            """)
            jobids = cur.fetchall()
            if len(jobids) > 1:
                raise Exception('ERROR: more than one possible input found, use'
                                '"tadbit describe" and select corresponding '
                                'jobid with --jobid')
        parse_jobid = jobids[0][0]
        # fetch path to parsed BED files
        cur.execute("""
        select distinct path from paths
        inner join filter_outputs on filter_outputs.pathid = paths.id
        where filter_outputs.name = 'valid-pairs' and paths.jobid = %s
        """ % parse_jobid)
        return cur.fetchall()[0][0]

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    glopts.add_argument('--bed', dest='bed', metavar="PATH",
                        action='store', default=None, type=str, 
                        help='''path to a TADbit-generated BED file with
                        filtered reads (other wise the tool will guess from the
                        working directory database)''')

    glopts.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int, required=True,
                        help='''resolution at which to output matrices''')

    glopts.add_argument('--perc_zeros', dest='perc_zeros', metavar="FLOAT",
                        action='store', default=99, type=float, 
                        help='maximum percentage of zeroes allowed per column')

    glopts.add_argument('--normalization', dest='resolution', metavar="STR",
                        action='store', default='ICE', nargs='+', type=str,
                        choices=['ICE', 'EXP'],
                        help='''[%(default)s] normalization(s) to apply. Order matters.''')

    glopts.add_argument('--factor', dest='factor', metavar="NUM",
                        action='store', default=1, type=float,
                        help='''[%(default)s] target mean value of a cell after
                        normalization (can be used to weight experiments before
                        merging)''')

    glopts.add_argument('--save', dest='save', metavar="STR",
                        action='store', default='genome', nargs='+', type=str,
                        choices=['genome', 'chromosomes'],
                        help='''[%(default)s] save genomic or chromosomic matrix.''')

    glopts.add_argument('--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')    

    glopts.add_argument('--keep', dest='keep', action='store',
                        default=['intra', 'genome'], nargs='+',
                        choices = ['intra', 'genome'],
                        help='''%(default)s Matrices to save, choices are
                        "intra" to keep intra-chromosomal matrices, "inter" to
                        keep inter-chromosomal matrices and "genome", to keep
                        genomic matrices.''')

    glopts.add_argument('--only_txt', dest='only_txt', action='store_true',
                      default=False,
                      help='Save only text file for matrices, not images')

    glopts.add_argument('--fast_filter', dest='fast_filter', action='store_true',
                      default=False,
                      help='only filter according to the percentage of zero count')

    glopts.add_argument('--force', dest='force', action='store_true',
                      default=False,
                      help='overwrite previously run job')

    parser.add_argument_group(glopts)

def check_options(opts):

    # check resume
    if not path.exists(opts.workdir) and opts.resume:
        print ('WARNING: can use output files, found, not resuming...')
        opts.resume = False

    if already_run(opts) and not opts.force:
        exit('WARNING: exact same job already computed, see JOBs table above')
