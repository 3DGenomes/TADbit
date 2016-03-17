"""

information needed

 - path working directory with parsed reads

"""
from argparse                     import HelpFormatter
from pytadbit                     import load_hic_data_from_reads
from pytadbit                     import Chromosome
from pytadbit.utils.sqlite_utils  import already_run, digest_parameters
from pytadbit.utils.sqlite_utils  import add_path, get_jobid, print_db
from pytadbit.utils.file_handling import mkdir
from pytadbit.mapping.analyze     import plot_distance_vs_interactions, hic_map
from os                           import path
import sqlite3 as lite
import time

DESC = 'Finds TAD or compartment segmentation in Hi-C data.'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()
    param_hash = digest_parameters(opts)

    raw_matrices, raw_matrix, nrm_matrices, nrm_matrix = path.join(
        opts.workdir, load_parameters_fromdb(opts))

    for crm_name in raw_matrices:
        crm = Chromosome('crm_name')
        crm.add_experiment(
            xnam, exp_type='Hi-C', enzyme=opts.enzyme,
            cell_type=opts.cell,
            identifier=opts.identifier, # general descriptive fields
            project=opts.project, # user descriptions
            resolution=opts.res,
            hic_data=xpath,
            norm_data=xnorm)

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
            where Type = 'Normalize'
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
        select distinct path from paths where paths.jobid = %s and paths.Type = 'BAD_COLUMNS'
        """ % parse_jobid)
        bad_columns = cur.fetchall()[0][0]
        cur.execute("""
        select distinct path from paths where paths.jobid = %s and paths.Type = 'BIASES'
        """ % parse_jobid)
        biases = cur.fetchall()[0][0]
        return bad_columns, biases

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

    glopts.add_argument('--norm_matrix', dest='norm_matrix', metavar="PATH",
                        action='store', default=None, type=str, 
                        help='''path to a matrix file with normalized read
                        counts''')
    
    glopts.add_argument('--raw_matrix', dest='raw_matrix', metavar="PATH",
                        action='store', default=None, type=str, 
                        help='''path to a matrix file with raw read
                        counts''')

    glopts.add_argument('--perc_zeros', dest='perc_zeros', metavar="FLOAT",
                        action='store', default=95, type=float, 
                        help='maximum percentage of zeroes allowed per column')

    glopts.add_argument('-j', '--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')    

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
