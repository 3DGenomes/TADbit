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
from os                           import path, remove
from string                       import ascii_letters
from random                       import random
from shutil                       import copyfile
from cPickle                      import dump
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

    mkdir(path.join(opts.workdir, '04_normalization'))

    print 'Get poor bins...'
    try:
        hic_data.filter_columns(perc_zero=opts.perc_zeros, draw_hist=True,
                                by_mean=not opts.fast_filter, savefig=path.join(
                                    opts.workdir, '04_normalization',
                                    'bad_columns_%s_%d_%s.pdf' % (
                                        opts.reso, opts.perc_zeros, param_hash)) if
                                not opts.fast_filter else None)
    except ValueError:
        hic_data.filter_columns(perc_zero=100, draw_hist=True,
                                by_mean=not opts.fast_filter, savefig=path.join(
                                    opts.workdir, '04_normalization',
                                    'bad_columns_%s_%d_%s.pdf' % (
                                        opts.reso, opts.perc_zeros, param_hash)) if
                                not opts.fast_filter else None)

    # bad columns
    bad_columns_file = path.join(opts.workdir, '04_normalization',
                                 'bad_columns_%s_%d_%s.tsv' % (
                                     opts.reso, opts.perc_zeros, param_hash))
    out_bad = open(bad_columns_file, 'w')
    out_bad.write('\n'.join([str(i) for i in hic_data.bads.keys()]))
    out_bad.close()

    # Identify biases
    if not opts.filter_only:
        print 'Get biases using ICE...'
        hic_data.normalize_hic(silent=False, max_dev=0.1, iterations=0,
                               factor=opts.factor)

    print 'Getting cis/trans...'
    cis_trans_N_D = cis_trans_N_d = float('nan')
    if not opts.filter_only:
        cis_trans_N_D = hic_data.cis_trans_ratio(normalized=True , diagonal=True )
        cis_trans_N_d = hic_data.cis_trans_ratio(normalized=True , diagonal=False)
    cis_trans_n_D = hic_data.cis_trans_ratio(normalized=False, diagonal=True )
    cis_trans_n_d = hic_data.cis_trans_ratio(normalized=False, diagonal=False)
        
    if not opts.filter_only:
        print 'Cis/Trans ratio of normalized matrix including the diagonal', cis_trans_N_D
        print 'Cis/Trans ratio of normalized matrix excluding the diagonal', cis_trans_N_d
    print 'Cis/Trans ratio of raw matrix including the diagonal', cis_trans_n_D
    print 'Cis/Trans ratio of raw matrix excluding the diagonal', cis_trans_n_d

    # Plot genomic distance vs interactions
    print 'Plot genomic distance vs interactions...'
    inter_vs_gcoord = path.join(opts.workdir, '04_normalization',
                                'interactions_vs_genomic-coords.pdf_%s_%s.pdf' % (
                                    opts.reso, param_hash))
    (_, _, _), (a2, _, _), (_, _, _) = plot_distance_vs_interactions(
        hic_data, max_diff=10000, resolution=opts.reso, normalized=not opts.filter_only,
        savefig=inter_vs_gcoord)
    
    print 'Decay slope 0.7-10 Mb\t%s' % a2

    # write biases
    bias_file = path.join(opts.workdir, '04_normalization',
                          'bias_%s_%s.tsv' % (opts.reso, param_hash))
    out_bias = 'NA'
    if not opts.filter_only:
        out_bias = open(bias_file, 'w')
        out_bias.write('\n'.join(['%d\t%f' % (i, hic_data.bias[i])
                                  for i in hic_data.bias])
                       + '\n')
        out_bias.close()


    # pickle the HiC-data object
    print 'Saving genomic matrix'
    pickle_path = path.join(opts.workdir, '04_normalization',
                            'hic-data_%s_%s.pickle' % (nice(opts.reso), param_hash))
    out = open(pickle_path, 'w')
    dump(hic_data, out)
    out.close()

    # to feed the save_to_db funciton
    intra_dir_nrm_fig = intra_dir_nrm_txt = None
    inter_dir_nrm_fig = inter_dir_nrm_txt = None
    genom_map_nrm_fig = genom_map_nrm_txt = None
    intra_dir_raw_fig = intra_dir_raw_txt = None
    inter_dir_raw_fig = inter_dir_raw_txt = None
    genom_map_raw_fig = genom_map_raw_txt = None

    if "intra" in opts.keep:
        print "  Saving intra chromosomal raw and normalized matrices..."
        if opts.only_txt:
            intra_dir_nrm_fig = None
            intra_dir_raw_fig = None
        else:
            intra_dir_nrm_fig = path.join(opts.workdir, '04_normalization',
                                          'intra_chromosome_nrm_images_%s_%s' % (opts.reso, param_hash))
            intra_dir_raw_fig = path.join(opts.workdir, '04_normalization',
                                          'intra_chromosome_raw_images_%s_%s' % (opts.reso, param_hash))
        intra_dir_nrm_txt = path.join(opts.workdir, '04_normalization',
                                      'intra_chromosome_nrm_matrices_%s_%s' % (opts.reso, param_hash))
        intra_dir_raw_txt = path.join(opts.workdir, '04_normalization',
                                      'intra_chromosome_raw_matrices_%s_%s' % (opts.reso, param_hash))
        if not opts.filter_only:
            hic_map(hic_data, normalized=True, by_chrom='intra', cmap='jet',
                    name=path.split(opts.workdir)[-1],
                    savefig=intra_dir_nrm_fig, savedata=intra_dir_nrm_txt)
        hic_map(hic_data, normalized=False, by_chrom='intra', cmap='jet',
                name=path.split(opts.workdir)[-1],
                savefig=intra_dir_raw_fig, savedata=intra_dir_raw_txt)

    if "inter" in opts.keep:
        print "  Saving inter chromosomal raw and normalized matrices..."
        if opts.only_txt:
            inter_dir_nrm_fig = None
            inter_dir_raw_fig = None
        else:
            if not opts.filter_only:
                inter_dir_nrm_fig = path.join(opts.workdir, '04_normalization',
                                              'inter_chromosome_nrm_images_%s_%s' % (opts.reso, param_hash))
            inter_dir_raw_fig = path.join(opts.workdir, '04_normalization',
                                      'inter_chromosome_raw_images_%s_%s' % (opts.reso, param_hash))
        if not opts.filter_only:
            inter_dir_nrm_txt = path.join(opts.workdir, '04_normalization',
                                          'inter_chromosome_nrm_matrices_%s_%s' % (opts.reso, param_hash))
        inter_dir_raw_txt = path.join(opts.workdir, '04_normalization',
                                  'inter_chromosome_raw_matrices_%s_%s' % (opts.reso, param_hash))
        if not opts.filter_only:
            hic_map(hic_data, normalized=True, by_chrom='inter', cmap='jet',
                    name=path.split(opts.workdir)[-1],
                    savefig=inter_dir_nrm_fig, savedata=inter_dir_nrm_txt)
        hic_map(hic_data, normalized=False, by_chrom='inter', cmap='jet',
                name=path.split(opts.workdir)[-1],
                savefig=inter_dir_raw_fig, savedata=inter_dir_raw_txt)

    if "genome" in opts.keep:
        print "  Saving normalized genomic matrix..."
        if opts.only_txt:
            if not opts.filter_only:
                genom_map_nrm_fig = path.join(opts.workdir, '04_normalization',
                                              'genomic_maps_nrm_%s_%s.pdf' % (opts.reso, param_hash))
            genom_map_raw_fig = path.join(opts.workdir, '04_normalization',
                                          'genomic_maps_raw_%s_%s.pdf' % (opts.reso, param_hash))
        else:
            genom_map_nrm_fig = None
            genom_map_raw_fig = None
        if not opts.filter_only:
            genom_map_nrm_txt = path.join(opts.workdir, '04_normalization',
                                          'genomic_nrm_%s_%s.tsv' % (opts.reso, param_hash))
        genom_map_raw_txt = path.join(opts.workdir, '04_normalization',
                                      'genomic_raw_%s_%s.tsv' % (opts.reso, param_hash))
        if not opts.filter_only:
            hic_map(hic_data, normalized=True, cmap='jet',
                    name=path.split(opts.workdir)[-1],
                savefig=genom_map_nrm_fig, savedata=genom_map_nrm_txt)
        hic_map(hic_data, normalized=False, cmap='jet',
                name=path.split(opts.workdir)[-1],
                savefig=genom_map_raw_fig, savedata=genom_map_raw_txt)

    finish_time = time.localtime()

    save_to_db (opts, cis_trans_N_D, cis_trans_N_d, cis_trans_n_D, cis_trans_n_d,
                a2, bad_columns_file, bias_file, inter_vs_gcoord, mreads,
                len(hic_data.bads.keys()), len(hic_data),
                intra_dir_nrm_fig, intra_dir_nrm_txt,
                inter_dir_nrm_fig, inter_dir_nrm_txt,
                genom_map_nrm_fig, genom_map_nrm_txt,
                intra_dir_raw_fig, intra_dir_raw_txt,
                inter_dir_raw_fig, inter_dir_raw_txt,
                genom_map_raw_fig, genom_map_raw_txt,
                pickle_path, launch_time, finish_time)

def save_to_db(opts, cis_trans_N_D, cis_trans_N_d, cis_trans_n_D, cis_trans_n_d,
               a2, bad_columns_file, bias_file, inter_vs_gcoord, mreads,
               nbad_columns, ncolumns,
               intra_dir_nrm_fig, intra_dir_nrm_txt,
               inter_dir_nrm_fig, inter_dir_nrm_txt,
               genom_map_nrm_fig, genom_map_nrm_txt,
               intra_dir_raw_fig, intra_dir_raw_txt,
               inter_dir_raw_fig, inter_dir_raw_txt,
               genom_map_raw_fig, genom_map_raw_txt,
               pickle_path, launch_time, finish_time):
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
                       type='table' AND name='NORMALIZE_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
            create table NORMALIZE_OUTPUTs
               (Id integer primary key,
                JOBid int,
                Input int,
                N_columns int,
                N_filtered int,
                CisTrans_nrm_all real,
                CisTrans_nrm_out real,
                CisTrans_raw_all real,
                CisTrans_raw_out real,
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
        add_path(cur, pickle_path     , 'PICKLE'     , jobid, opts.workdir)
        add_path(cur, bad_columns_file, 'BAD_COLUMNS', jobid, opts.workdir)
        add_path(cur, bias_file       , 'BIASES'     , jobid, opts.workdir)
        add_path(cur, inter_vs_gcoord , 'FIGURE'     , jobid, opts.workdir)
        add_path(cur, mreads          , '2D_BED'     , jobid, opts.workdir)
        # get pathid of input
        cur.execute("select id from paths where path = '%s'" % (path.relpath(mreads, opts.workdir)))
        input_bed = cur.fetchall()[0][0]
        if intra_dir_nrm_fig:
            add_path(cur, intra_dir_nrm_fig, 'FIGURES', jobid, opts.workdir)
        if intra_dir_nrm_fig:
            add_path(cur, intra_dir_nrm_txt, 'NRM_MATRICES', jobid, opts.workdir)
        if inter_dir_nrm_fig:
            add_path(cur, inter_dir_nrm_fig, 'FIGURES', jobid, opts.workdir)
        if inter_dir_nrm_fig:
            add_path(cur, inter_dir_nrm_txt, 'NRM_MATRICES', jobid, opts.workdir)
        if genom_map_nrm_fig:
            add_path(cur, genom_map_nrm_fig, 'FIGURE', jobid, opts.workdir)
        if genom_map_nrm_txt:
            add_path(cur, genom_map_nrm_txt, 'NRM_MATRIX', jobid, opts.workdir)
        if intra_dir_raw_fig:
            add_path(cur, intra_dir_raw_fig, 'FIGURES', jobid, opts.workdir)
        if intra_dir_raw_fig:
            add_path(cur, intra_dir_raw_txt, 'RAW_MATRICES', jobid, opts.workdir)
        if inter_dir_raw_fig:
            add_path(cur, inter_dir_raw_fig, 'FIGURES', jobid, opts.workdir)
        if inter_dir_raw_fig:
            add_path(cur, inter_dir_raw_txt, 'RAW_MATRICES', jobid, opts.workdir)
        if genom_map_raw_fig:
            add_path(cur, genom_map_raw_fig, 'FIGURE', jobid, opts.workdir)
        if genom_map_raw_txt:
            add_path(cur, genom_map_raw_txt, 'RAW_MATRIX', jobid, opts.workdir)

        try:
            cur.execute("""
            insert into NORMALIZE_OUTPUTs
            (Id  , JOBid,     Input, N_columns,   N_filtered, CisTrans_nrm_all,   CisTrans_nrm_out,   CisTrans_raw_all,   CisTrans_raw_out, Slope_700kb_10Mb,   Resolution,      Factor)
            values
            (NULL,    %d,        %d,        %d,           %d,               %f,                 %f,                 %f,                 %f,               %f,           %d,          %f)
            """ % (jobid, input_bed,  ncolumns, nbad_columns,    cis_trans_N_D,      cis_trans_N_d,      cis_trans_n_D,      cis_trans_n_d,               a2,    opts.reso, opts.factor))
        except lite.OperationalError:
            try:
                cur.execute("""
                insert into NORMALIZE_OUTPUTs
                (Id  , JOBid,     Input, N_columns,   N_filtered,  CisTrans_raw_all,   CisTrans_raw_out, Slope_700kb_10Mb,   Resolution,      Factor)
                values
                (NULL,    %d,        %d,        %d,           %d,                %f,                 %f,               %f,           %d,          %f)
                """ % (jobid, input_bed,  ncolumns, nbad_columns,     cis_trans_n_D,      cis_trans_n_d,               a2,    opts.reso, opts.factor))
            except lite.OperationalError:
                print 'WANRING: Normalized table not written!!!'
            
        print_db(cur, 'MAPPED_INPUTs')
        print_db(cur, 'PATHs')
        print_db(cur, 'MAPPED_OUTPUTs')
        print_db(cur, 'PARSED_OUTPUTs')
        print_db(cur, 'JOBs')
        print_db(cur, 'INTERSECTION_OUTPUTs')        
        print_db(cur, 'FILTER_OUTPUTs')
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
        else:
            parse_jobid = opts.jobid
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
                        action='store', default=95, type=float, 
                        help=('[%(default)s%%] maximum percentage of zeroes '
                              'allowed per column.'))

    glopts.add_argument('--normalization', dest='resolution', metavar="STR",
                        action='store', default='ICE', nargs='+', type=str,
                        choices=['ICE', 'EXP'],
                        help='''[%(default)s] normalization(s) to apply.
                        Order matters.''')

    glopts.add_argument('--factor', dest='factor', metavar="NUM",
                        action='store', default=1, type=float,
                        help='''[%(default)s] target mean value of a cell after
                        normalization (can be used to weight experiments before
                        merging)''')

    glopts.add_argument('-j', '--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')    

    glopts.add_argument('--keep', dest='keep', action='store',
                        default=['intra', 'genome'], nargs='+',
                        choices = ['intra', 'inter', 'genome', 'none'],
                        help='''%(default)s Matrices to save, choices are
                        "intra" to keep intra-chromosomal matrices, "inter" to
                        keep inter-chromosomal matrices and "genome", to keep
                        genomic matrices.''')

    glopts.add_argument('--only_txt', dest='only_txt', action='store_true',
                      default=False,
                      help='Save only text file for matrices, not images')

    glopts.add_argument('--filter_only', dest='filter_only', action='store_true',
                      default=False,
                      help='skip normalization')

    glopts.add_argument('--fast_filter', dest='fast_filter', action='store_true',
                      default=False,
                      help='only filter according to the percentage of zero count')

    glopts.add_argument('--force', dest='force', action='store_true',
                      default=False,
                      help='overwrite previously run job')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    parser.add_argument_group(glopts)

def check_options(opts):

    # check resume
    if not path.exists(opts.workdir):
        raise IOError('ERROR: wordir not found.')

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

    # check if job already run using md5 digestion of parameters
    if already_run(opts):
        if 'tmpdb' in opts and opts.tmpdb:
            remove(path.join(dbdir, dbfile))
        exit('WARNING: exact same job already computed, see JOBs table above')

def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)
