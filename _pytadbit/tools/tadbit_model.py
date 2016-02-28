"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""

from argparse                       import HelpFormatter
from os                             import path, listdir
from pytadbit.imp.imp_modelling import generate_3d_models, TADbitModelingOutOfBound
from pytadbit import load_structuralmodels
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
from numpy import arange
from cPickle import load

DESC = ("Generates 3D models given an input interaction matrix and a set of "
        "input parameters")


def run(opts):
    check_options(opts)

    launch_time = time.localtime()

    # load data
    print 'loading data'
    crm = load_hic_data(opts)
    exp = crm.experiments[0]
    
    mkdir(path.join(opts.workdir, '06_model'))

    beg, end = opts.beg or 1, opts.end or exp.size

    outdir = path.join(opts.workdir, '06_model',
		       'chr%s_%s-%s' % (opts.crm, beg, end))
    mkdir(outdir)

    models = compile_models(outdir)
    print models

    zscores, values, zeros = exp._sub_experiment_zscore(beg, end)
    nloci = end - beg + 1
    coords = {"crm"  : opts.crm,
              "start": opts.beg,
              "end"  : opts.end}
    
    print 'Start modeling'
    print ('# %3s %6s %7s %7s %6s %7s\n' % (
        "num", "upfrq", "lowfrq", "maxdist",
        "scale", "cutoff"))

    if opts.job_list:
        job_file_handler = open(path.join(outdir, 'job_list.q'), 'w')
    for m, u, l, d, s in product(opts.maxdist, opts.upfreq, opts.lowfreq,
                                 opts.dcutoff, opts.scale):
        
        print('%5s %6s %7s %7s %6s %7s  ' % ('x', u, l, m, s, d))
        mkdir(path.join(outdir, 'cfg_%s_%s_%s_%s_%s' % (m, u, l, d, s)))

        optpar = {'maxdist': m,
                  'upfreq' : u,
                  'lowfreq': l,
                  'dcutoff': d,
                  'scale'  : s,
                  'kforce' : 5}

	# write list of jobs to be run separately
        if opts.job_list:
            for rand in xrange(1, opts.nmodels + 1, opts.nmodels_run):
                job_file_handler.write(('tadbit model --input_matrix %s '
                                        '--maxdist %s --upfreq %s --lowfreq=%s '
                                        '--dcutoff %s --scale %s --rand %s '
                                        '--nmodels_run %s\n') % (
                                           opts.matrix, m, u, l, d, s, rand,
					   opts.nmodels_run))
            continue

        # compute models
        try:
            models =  generate_3d_models(zscores, opts.reso, nloci,
                                         values=values, n_models=opts.nmodels,
                                         n_keep=opts.nkeep,
                                         n_cpus=opts.cpus, keep_all=True,
                                         first=opts.rand, container=None,
                                         config=optpar, coords=coords,
					 zeros=zeros)

        except TADbitModelingOutOfBound:
            warn('WARNING: scale (here %s) x resolution (here %d) should be '
                 'lower than maxdist (here %d instead of at least: %d)' % (
                     s, opts.reso, m, s * opts.reso))
            continue

        # Save models
        models.save_models(
            path.join(outdir, 'cfg_%s_%s_%s_%s_%s' % (m, u, l, d, s),
                      ('models_%s-%s.pick' % (opts.rand, opts.rand + opts.nmodels))
                      if opts.nmodels > 1 else 
                      ('model_%s.pick' % (opts.rand))))

    if opts.optimize:
	finish_time = time.localtime()
	return



    finish_time = time.localtime()

    # save all job information to sqlite DB
    # save_to_db(opts, counts, multis, f_names1, f_names2, out_file1, out_file2,
    #            launch_time, finish_time)

def compile_models(outdir):
    models = {}
    for cfg_dir in listdir(outdir):
        if not cfg_dir.startswith('cfg_'):
            continue
        _, m, u, l, d, s = cfg_dir.split('_')
        m, u, l, d, s = int(m), float(u), float(l), float(d), float(s)
        for fmodel in listdir(path.join(outdir, cfg_dir)):
            if not fmodel.startswith('models_'):
                continue
            if not (m, u, l, d, s) in models:
                models[(m, u, l, d, s)] = load_structuralmodels(path.join(
                    outdir, cfg_dir, fmodel))
            else:
                sm = load(path.join(outdir, cfg_dir, fmodel))
                if models[(m, u, l, d, s)]._config != sm['config']:
                    raise Exception('ERROR: clean directory, hetergoneous data')
                models[(m, u, l, d, s)]._extend_models(sm['models'])
                models[(m, u, l, d, s)]._bad_models.extend(sm['bad_models'])
    return models

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

    glopts.add_argument('--job_list', dest='job_list', action='store_true',
                      default=False,
                      help='generate a a file with a list of jobs to be run in '
                        'a cluster')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
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
    glopts.add_argument('-r', '--reso', dest='reso', metavar="INT", type=int,
                        help='resolution of the Hi-C experiment')
    glopts.add_argument('--input_matrix', dest='matrix', metavar="PATH",
                        type=str,
                        help='''In case input was not generated with the TADbit
                        tools''')

    glopts.add_argument('--nmodels_run', dest='nmodels_run', metavar="INT",
                        default=None, type=int,
                        help='[ALL] number of models to run with this call')

    glopts.add_argument('--nmodels', dest='nmodels', metavar="INT",
                        default=5000, type=int,
                        help=('[%(default)s] number of models to generate for' +
                              ' modeling'))

    glopts.add_argument('--nkeep', dest='nkeep', metavar="INT",
                        default=1000, type=int,
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
    glopts.add_argument("-C", "--cpu", dest="cpus", type=int,
                        default=1, help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    parser.add_argument_group(glopts)

def check_options(opts):
    # do the division to bins
    try:
        opts.beg = int(float(opts.beg) / opts.reso)
        opts.end = int(float(opts.end) / opts.reso)
        if opts.end - opts.beg <= 2:
            raise Exception('"beg" and "end" parameter should be given in ' +
                            'genomic coordinates, not bin')
    except TypeError:
        pass

    # turn options into lists
    opts.scale   = (tuple(arange(*[float(s) for s in opts.scale.split(':')  ]))
                    if ':' in opts.scale   else [float(opts.scale  )])
    
    opts.maxdist = (tuple(range (*[int  (i) for i in opts.maxdist.split(':')]))
                    if ':' in opts.maxdist else [int  (opts.maxdist)])

    opts.upfreq  = (tuple(arange(*[float(i) for i in opts.upfreq.split(':') ]))
                    if ':' in opts.upfreq  else [float(opts.upfreq )])

    opts.lowfreq = (tuple(arange(*[float(i) for i in opts.lowfreq.split(':')]))
                    if ':' in opts.lowfreq else [float(opts.lowfreq)])

    opts.dcutoff = (tuple(arange(*[float(i) for i in opts.dcutoff.split(':')]))
                    if ':' in opts.dcutoff else [float(opts.dcutoff)])

    opts.nmodels_run = opts.nmodels_run or opts.nmodels

    mkdir(opts.workdir)

def load_hic_data(opts):
    """
    Load Hi-C data
    """
    # Start reading the data
    crm = Chromosome(opts.crm) # Create chromosome object

    crm.add_experiment('test', exp_type='Hi-C', resolution=opts.reso,
                       norm_data=opts.matrix)
    if opts.beg > crm.experiments[-1].size:
        raise Exception('ERROR: beg parameter is larger than chromosome size.')
    if opts.end > crm.experiments[-1].size:
        print ('WARNING: end parameter is larger than chromosome ' +
               'size. Setting end to %s.\n' % (crm.experiments[-1].size *
                                               opts.reso))
        opts.end = crm.experiments[-1].size
    return crm
