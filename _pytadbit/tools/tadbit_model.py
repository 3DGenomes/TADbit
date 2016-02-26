"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""

from argparse                       import HelpFormatter
from os                             import path
from pytadbit.imp.imp_modelling import generate_3d_models
from pytadbit import Chromosome
from pytadbit.utils.file_handling   import mkdir
from pytadbit.utils.sqlite_utils    import get_path_id, add_path, print_db, get_jobid
from pytadbit.utils.sqlite_utils    import already_run, digest_parameters
from itertools import product
import time
import logging
import fcntl
import sqlite3 as lite
from warnings import warn

DESC = ("Generates 3D models given an input interaction matrix and a set of "
        "input parameters")


def run(opts):
    check_options(opts)

    launch_time = time.localtime()

    if opts.xname:
        xnames = opts.xname
    elif opts.data[0]:
        xnames = [path.split(d)[-1] for d in opts.data]
    else:
        xnames = [path.split(d)[-1] for d in opts.norm]

    # load data
    print 'loading data'
    crm = load_hic_data(opts, xnames)
    exp = crm.experiments[0]
    

    beg, end = opts.beg or 1, opts.end or exp.size
    zscores, values, zeros = exp._sub_experiment_zscore(beg, end)
    nloci = end - beg + 1
    coords = {"crm"  : opts.crm,
              "start": opts.beg,
              "end"  : opts.end}
    
    print 'Start modeling'
    for m, u, l, d, s in product(opts.maxdist, opts.upfreq, opts.lowfreq,
                                 opts.dcutoff, opts.scale):

        optpar = {'maxdist': m,
                  'upfreq' : u,
                  'lowfreq': l,
                  'dcutoff': d,
                  'scale'  : s,
                  'kforce' : 5}
        
        # compute models
        models =  generate_3d_models(zscores, opts.res, nloci,
                                     values=values, n_models=opts.nmodels_mod,
                                     n_keep=opts.nkeep_mod,
                                     n_cpus=opts.ncpus, keep_all=True,
                                     first=opts.first, container=opts.container,
                                     config=optpar, coords=coords, zeros=zeros)

        # Save models
        models.save_models(
            path.join(opts.outdir, "%s", "%s" + ".models"))
    


    finish_time = time.localtime()

    # save all job information to sqlite DB
    # save_to_db(opts, counts, multis, f_names1, f_names2, out_file1, out_file2,
    #            launch_time, finish_time)


def save_to_db(opts, counts, multis, f_names1, f_names2, out_file1, out_file2,
               launch_time, finish_time):
    con = lite.connect(path.join(opts.workdir, 'trace.db'))
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='PARSED_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
        create table MAPPED_OUTPUTs
           (Id integer primary key,
            PATHid int,
            BEDid int,
            Uniquely_mapped int,
            unique (PATHid, BEDid))""")
            cur.execute("""
        create table PARSED_OUTPUTs
           (Id integer primary key,
            PATHid int,
            Total_interactions int,
            Multiples int,
            unique (PATHid))""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            param_hash = digest_parameters(opts, get_md5=True )
            cur.execute("""
    insert into JOBs
     (Id  , Parameters, Launch_time, Finish_time,    Type, Parameters_md5)
    values
     (NULL,       '%s',        '%s',        '%s', 'Parse',           '%s')
     """ % (parameters,
            time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
            time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass
        jobid = get_jobid(cur)
        add_path(cur, out_file1, 'BED', jobid, opts.workdir)
        for genome in opts.genome:
            add_path(cur, genome, 'FASTA', jobid, opts.workdir)
        if out_file2:
            add_path(cur, out_file2, 'BED', jobid, opts.workdir)
        fnames = f_names1, f_names2
        outfiles = out_file1, out_file2
        for count in counts:
            try:
                sum_reads = 0
                for i, item in enumerate(counts[count]):
                    cur.execute("""
                    insert into MAPPED_OUTPUTs
                    (Id  , PATHid, BEDid, Uniquely_mapped)
                    values
                    (NULL,    %d,     %d,      %d)
                    """ % (get_path_id(cur, fnames[count][i], opts.workdir),
                           get_path_id(cur, outfiles[count], opts.workdir),
                           counts[count][item]))
                    sum_reads += counts[count][item]
            except lite.IntegrityError:
                print 'WARNING: already parsed (MAPPED_OUTPUTs)'
            try:
                cur.execute("""
                insert into PARSED_OUTPUTs
                (Id  , PATHid, Total_interactions, Multiples)
                values
                (NULL,     %d,      %d,        %d)
                """ % (get_path_id(cur, outfiles[count], opts.workdir),
                       sum_reads, multis[count]))
            except lite.IntegrityError:
                print 'WARNING: already parsed (PARSED_OUTPUTs)'
        print_db(cur, 'MAPPED_INPUTs')
        print_db(cur, 'PATHs')
        print_db(cur, 'MAPPED_OUTPUTs')
        print_db(cur, 'PARSED_OUTPUTs')
        print_db(cur, 'JOBs')

def load_parameters_fromdb(workdir, reads=None, jobids=None):
    con = lite.connect(path.join(workdir, 'trace.db'))
    fnames = {1: [], 2: []}
    ids = []
    with con:
        cur = con.cursor()
        # fetch file names to parse
        if not jobids:
            jobids = {}
            for read in reads:
                cur.execute("""
                select distinct JOBs.Id from JOBs
                   inner join PATHs on (JOBs.Id = PATHs.JOBid)
                   inner join MAPPED_INPUTs on (PATHs.Id = MAPPED_INPUTs.MAPPED_OUTPUTid)
                 where MAPPED_INPUTs.Read = %d
                """ % read)
                jobids[read] = [j[0] for j in cur.fetchall()]
                if len(jobids[read]) > 1:
                    warn(('WARNING: more than one possible input found for read %d '
                          '(jobids: %s), use "tadbit describe" and select corresponding '
                          'jobid with --jobids option') % (
                             read, ', '.join([str(j) for j in jobids[read]])))
        for read in reads:
            for jobid in jobids[read]:
                cur.execute("""
                select distinct PATHs.Id,PATHs.Path from PATHs
                inner join MAPPED_INPUTs on PATHs.Id = MAPPED_INPUTs.MAPPED_OUTPUTid
                where MAPPED_INPUTs.Read = %d and PATHs.JOBid = %d
                """ % (read, jobid))
                for fname in cur.fetchall():
                    ids.append(fname[0])
                    fnames[read].append(path.join(workdir, fname[1]))
        # GET enzyme name
        enzymes = []
        for fid in ids:
            cur.execute("""
            select distinct MAPPED_INPUTs.Enzyme from MAPPED_INPUTs
            where MAPPED_INPUTs.MAPPED_OUTPUTid=%d
            """ % fid)
            enzymes.extend(cur.fetchall())
        if len(set(reduce(lambda x, y: x+ y, enzymes))) != 1:
            raise Exception(
                'ERROR: different enzymes used to generate these files')
        renz = enzymes[0][0]
    return fnames[1], fnames[2], renz
        

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    # glopts.add_argument('--qc_plot', dest='quality_plot', action='store_true',
    #                   default=False,
    #                   help='generate a quality plot of FASTQ and exits')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    glopts.add_argument('--optimize', dest='optimize', 
                        default=False, action="store_true",
                        help='''optimization run, store less info about models''')

    glopts.add_argument('--rand', dest='rand', metavar="INT",
                        type=int, default=1, 
                        help='''[%(default)s] random initial number. NOTE:
                        when running single model at the time, should be
                        different for each run''')

    glopts.add_argument('--skip', dest='skip', action='store_true',
                      default=False,
                      help='[DEBUG] in case already mapped.')

    glopts.add_argument('--crm', dest='crm', metavar="NAME",
                        help='chromosome name')

    glopts.add_argument('--beg', dest='beg', metavar="INT", type=float,
                        default=None,
                        help='genomic coordinate from which to start modeling')

    glopts.add_argument('--end', dest='end', metavar="INT", type=float,
                        help='genomic coordinate where to end modeling')

    glopts.add_argument('--res', dest='res', metavar="INT", type=int,
                        help='resolution of the Hi-C experiment')

    glopts.add_argument('--input_matrix', dest='matrix', metavar="PATH",
                        type=str,
                        help='''In case input was not generated with the TADbit
                        tools''')

    glopts.add_argument('--nmodels_run', dest='nmodels_run', metavar="INT",
                        default=None, type=int,
                        help='[ALL] number of models to run with this call')

    glopts.add_argument('--nmodels_mod', dest='nmodels_mod', metavar="INT",
                        default='5000', type=int,
                        help=('[%(default)s] number of models to generate for' +
                              ' modeling'))

    glopts.add_argument('--nkeep_mod', dest='nkeep_mod', metavar="INT",
                        default='1000', type=int,
                        help=('[%(default)s] number of models to keep for ' +
                        'modeling'))

    glopts.add_argument('--maxdist', action='store', metavar="LIST",
                        default='400', dest='maxdist',
                        help='range of numbers for maxdist' +
                        ', i.e. 400:1000:100 -- or just a number')
    glopts.add_argument('--upfreq', dest='upfreq', metavar="LIST",
                        default='0',
                        help='range of numbers for upfreq' +
                        ', i.e. 0:1.2:0.3 --  or just a number')
    glopts.add_argument('--lowfreq', dest='lowfreq', metavar="LIST",
                        default='0',
                        help='range of numbers for lowfreq' +
                        ', i.e. -1.2:0:0.3 -- or just a number')
    glopts.add_argument('--scale', dest='scale', metavar="LIST",
                        default="0.01",
                        help='[%(default)s] range of numbers to be test as ' +
                        'optimal scale value, i.e. 0.005:0.01:0.001 -- Can ' +
                        'also pass only one number')
    glopts.add_argument('--dcutoff', dest='dcutoff', metavar="LIST",
                        default="2",
                        help='[%(default)s] range of numbers to be test as ' +
                        'optimal distance cutoff parameter (distance, in ' +
                        'number of beads, from which to consider 2 beads as ' +
                        'being close), i.e. 1:5:0.5 -- Can also pass only one' +
                        ' number')
    glopts.add_argument('--nmodels_opt', dest='nmodels_opt', metavar="INT",
                        default='500', type=int,
                        help='[%(default)s] number of models to generate for ' +
                        'optimization')
    glopts.add_argument('--nkeep_opt', dest='nkeep_opt', metavar="INT",
                        default='100', type=int,
                        help='[%(default)s] number of models to keep for ' +
                        'optimization')

    parser.add_argument_group(glopts)

def check_options(opts):
    # do the division to bins
    if not opts.tad_only:
        try:
            opts.beg = int(float(opts.beg) / opts.res)
            opts.end = int(float(opts.end) / opts.res)
            if opts.end - opts.beg <= 2:
                raise Exception('"beg" and "end" parameter should be given in ' +
                                'genomic coordinates, not bin')
        except TypeError:
            pass

    # turn options into lists
    opts.scale   = (tuple([float(i) for i in opts.scale.split(':')  ])
                    if ':' in opts.scale   else [float(opts.scale  )])
    
    opts.maxdist = (tuple([int  (i) for i in opts.maxdist.split(':')])
                    if ':' in opts.maxdist else [int  (opts.maxdist)])

    opts.upfreq  = (tuple([float(i) for i in opts.upfreq.split(':') ])
                    if ':' in opts.upfreq  else [float(opts.upfreq )])

    opts.lowfreq = (tuple([float(i) for i in opts.lowfreq.split(':')])
                    if ':' in opts.lowfreq else [float(opts.lowfreq)])

    opts.dcutoff = (tuple([float(i) for i in opts.dcutoff.split(':')])
                    if ':' in opts.dcutoff else [float(opts.dcutoff)])


def load_hic_data(opts, xnames):
    """
    Load Hi-C data
    """
    # Start reading the data
    crm = Chromosome(opts.crm, species=(
        opts.species.split('_')[0].capitalize() + opts.species.split('_')[1]
                          if '_' in opts.species else opts.species),
                          centromere_search=opts.centromere,
                          assembly=opts.assembly) # Create chromosome object

    logging.info("\tReading input data...")
    for xnam, xpath, xnorm in zip(xnames, opts.data, opts.norm):
        crm.add_experiment(
            xnam, exp_type='Hi-C', enzyme=opts.enzyme,
            cell_type=opts.cell,
            identifier=opts.identifier, # general descriptive fields
            project=opts.project, # user descriptions
            resolution=opts.res,
            hic_data=xpath,
            norm_data=xnorm)
        if not xnorm:
            crm.experiments[xnam].filter_columns(diagonal=not opts.nodiag)
            logging.info("\tNormalizing HiC data of %s...", xnam)
            crm.experiments[xnam].normalize_hic(iterations=10, max_dev=0.1)
    if opts.beg > crm.experiments[-1].size:
        raise Exception('ERROR: beg parameter is larger than chromosome size.')
    if opts.end > crm.experiments[-1].size:
        logging.info('WARNING: end parameter is larger than chromosome ' +
                     'size. Setting end to %s.\n' % (crm.experiments[-1].size *
                                                     opts.res))
        opts.end = crm.experiments[-1].size
    return crm
