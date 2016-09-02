"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""

from argparse                     import HelpFormatter
from os                           import path, listdir, remove
from string                       import ascii_letters
from random                       import random
from shutil                       import copyfile
from pytadbit.imp.imp_modelling   import generate_3d_models, TADbitModelingOutOfBound
from pytadbit                     import load_structuralmodels
from pytadbit                     import Chromosome
from pytadbit.utils.file_handling import mkdir
from pytadbit.utils.extraviews    import nicer
from pytadbit.utils.sqlite_utils  import get_path_id, add_path, print_db, get_jobid
from pytadbit.utils.sqlite_utils  import digest_parameters
from pytadbit                     import load_hic_data_from_reads
from pytadbit                     import get_dependencies_version
from pytadbit.parsers.hic_parser  import optimal_reader
from itertools                    import product
from warnings                     import warn
from numpy                        import arange
from cPickle                      import load
from hashlib                      import md5
import sqlite3 as lite
import time



DESC = ("Generates 3D models given an input interaction matrix and a set of "
        "input parameters")


def write_one_job(opts, rand, m, u, l, s, job_file_handler, nmodels_run=None):
    if nmodels_run is None:
        nmodels_run = opts.nmodels_run
    (m, u, l, s) = tuple(map(my_round, (m, u, l, s)))
    job_file_handler.write(
        'tadbit model -w %s -r %d -C %d --input_matrix %s '
        '--maxdist %s --upfreq %s --lowfreq=%s --scale %s '
        '--rand %s --nmodels %s --nkeep %s '
        '--crm %s --beg %d --end %d --perc_zero %s %s\n' % (
            opts.workdir, opts.reso,
            min(opts.cpus, opts.nmodels_run), opts.matrix,
            m, u, l, s, rand,
            nmodels_run,  # equal to nmodels if not defined
            nmodels_run,  # keep, equal to nmodel_run if defined
            opts.crm, opts.beg * opts.reso, opts.end * opts.reso,
            opts.perc_zero,
            '--optimize ' if opts.optimize else ''))


def run_batch_job(exp, opts, m, u, l, s, outdir):
    zscores, values, zeros = exp._sub_experiment_zscore(opts.beg, opts.end)
    zeros = tuple([i not in zeros for i in xrange(opts.end - opts.beg + 1)])
    nloci = opts.end - opts.beg + 1
    coords = {"crm"  : opts.crm,
              "start": opts.beg,
              "end"  : opts.end}

    optpar = {'maxdist': float(m),
              'upfreq' : float(u),
              'lowfreq': float(l),
              'scale'  : float(s),
              'kforce' : 5}

    models = generate_3d_models(zscores, opts.reso, nloci,
                                values=values, n_models=opts.nmodels,
                                n_keep=opts.nkeep,
                                n_cpus=opts.cpus, keep_all=True,
                                start=int(opts.rand), container=None,
                                config=optpar, coords=coords,
                                zeros=zeros)
    # Save models
    muls = tuple(map(my_round, (m, u, l, s)))
    dirname = 'cfg_%s_%s_%s_%s' % muls
    runned = [int(mod['rand_init']) for mod in models]
    if not len(runned):
        raise Exception(("\n\n\nNothing to be done.\n\n"
                         "   All models asked for are already run.\n"
                         "    - ask for more models\n"
                         "    - use higher random initial number\n"
                         "    - go ahead with the analysis!"))
    models.save_models(
        path.join(outdir, dirname,
                  ('models_%s-%s.pick' % (min(runned), max(runned)))
                  if len(runned) > 1 else
                  ('model_%s.pick' % (runned[0]))))


def my_round(num, val=4):
    num = round(float(num), val)
    return str(int(num) if num == int(num) else num)


def optimization(exp, opts, job_file_handler, outdir):
    models = compile_models(opts, outdir)
    print('\nOptimizing parameters...')
    print('# %3s %6s %7s %7s %6s\n' % (
        "num", "upfrq", "lowfrq", "maxdist",
        "scale"))
    for m, u, l, s in product(opts.maxdist, opts.upfreq, opts.lowfreq, opts.scale):
        muls = tuple(map(my_round, (m, u, l, s)))
        if muls in models:
            print('%5s %6s %7s %7s %6s  ' % ('x', u, l, m, s))
            continue
        elif opts.job_list:
            print('%5s %6s %7s %7s %6s  ' % ('o', u, l, m, s))
        else:
            print('%5s %6s %7s %7s %6s  ' % ('-', u, l, m, s))
        mkdir(path.join(outdir, 'cfg_%s_%s_%s_%s' % muls))

        # write list of jobs to be run separately
        if opts.job_list:
            for rand in xrange(1, opts.nmodels + 1, opts.nmodels_run):
                write_one_job(opts, rand, m, u, l, s, job_file_handler)
            continue

        # compute models
        try:
            run_batch_job(exp, opts, m, u, l, s, outdir)
        except TADbitModelingOutOfBound:
            warn('WARNING: scale (here %s) x resolution (here %d) should be '
                 'lower than maxdist (here %d instead of at least: %d)' % (
                     s, opts.reso, m, s * opts.reso))
            continue

    if opts.job_list:
        job_file_handler.close()


def correlate_models(opts, outdir, exp, corr='spearman', off_diag=1,
                     verbose=True):
    models = compile_models(opts, outdir, exp=exp)
    results = {}
    num = 1
    if verbose:
        print('\n\n# %13s %6s %7s %7s %6s %7s %7s\n' % (
            "Optimization", "UpFreq", "LowFreq", "MaxDist",
            "scale", "cutoff", "| Correlation"))
    for m, u, l, s in models:
        muls = tuple(map(my_round, (m, u, l, s)))
        result = 0
        d = float('nan')
        for d in opts.dcutoff:
            try:
                result = models[muls].correlate_with_real_data(
                    cutoff=(int(d * opts.reso * float(s))),
                    corr=corr, off_diag=off_diag)[0]
            except Exception, e:
                print '  SKIPPING: %s' % e
                result = 0
            name = tuple(map(my_round, (m, u, l, d, s)))
            if verbose:
                print(('%8s/%-6s %6s %7s %7s %6s %7s | %.4f' %
                       (num, len(models), u, l, m, s, d, result)))
            num += 1
            results[name] = {'corr'   : result,
                             'nmodels': (len(models[muls]) +
                                         len(models[muls]._bad_models)),
                             'kept'   : len(models[muls])}
    # get the best combination
    best = ({'corr': None}, [None, None, None, None, None])
    for m, u, l, d, s in results:
        if results[(m, u, l, d, s)]['corr'] > best[0]['corr']:
            best = results[(m, u, l, d, s)], [u, l, m, s, d]
    if verbose:
        print '\nBest combination:'
        print('  %5s     %6s %7s %7s %6s %6s %.4f\n' % tuple(
            ['=>'] + best[1] + [best[0]['corr']]))

    u, l, m, s, d = best[1]
    optpar = {'maxdist': m,
              'upfreq' : u,
              'lowfreq': l,
              'scale'  : s,
              'kforce' : 5}
    return optpar, results


def big_run(exp, opts, job_file_handler, outdir, optpar):
    m, u, l, s = (optpar['maxdist'], 
                  optpar['upfreq' ],
                  optpar['lowfreq'],
                  optpar['scale'  ])
    muls = tuple(map(my_round, (m, u, l, s)))
    # this to take advantage of previously runned models
    models = compile_models(opts, outdir, wanted=muls)[muls]
    models.define_best_models(len(models) + len(models._bad_models))
    
    if opts.rand == '1':
        # In this case we can use models already_run runned
        runned = [int(mod['rand_init']) for mod in models]
        start = str(max(runned) + 1)
        # reduce number of models to run
        opts.nmodels -= len(runned)
        print 'Using %s pre-calculated models with m:%s u:%s l:%s s:%s ' % (
            start, m, u, l, s)
        opts.rand = str(int(start) + 1)
    else:
        start = 1
    if opts.rand != '1' and int(opts.rand) < int(start):
        raise Exception('ERROR: found %s pre-computed models, use a higher '
                        'rand. init. number or delete the files' % (start))

    print 'Computing %s models' % opts.nmodels
    if opts.job_list:
        # write jobs
        nmodels_run = opts.nmodels_run
        print 'START', int(start), int(start) + opts.nmodels
        for rand in xrange(int(start), int(start) + opts.nmodels,
                           opts.nmodels_run):
            nmodels_run = min(opts.nmodels_run,
                              int(start) + opts.nmodels - rand)
            write_one_job(opts, rand, m, u, l, s, job_file_handler,
                          nmodels_run)
        job_file_handler.close()
        return

    # compute models
    try:
        run_batch_job(exp, opts, m, u, l, s, outdir)
    except TADbitModelingOutOfBound:
        warn('WARNING: scale (here %s) x resolution (here %d) should be '
             'lower than maxdist (here %d instead of at least: %d)' % (
                 s, opts.reso, m, s * opts.reso))


def run(opts):
    check_options(opts)

    launch_time = time.localtime()

    print('''
%s%s

  - Region: Chromosome %s from %d to %d at resolution %s (%d particles)


    ''' % ('Preparing ' if opts.job_list else '',
           ('Optimization\n' + '*' * (21 if opts.job_list else 11)) if opts.optimize else
           ('Modeling\n' + '*' * (18 if opts.job_list else 8)),
           opts.crm, opts.ori_beg, opts.ori_end, nicer(opts.reso),
           opts.end - opts.beg))

    # load data
    if opts.matrix:
        crm = load_hic_data(opts)
    else:
        # FIXME: copied from somewhere else
        (bad_co, bad_co_id, biases, biases_id,
         mreads, mreads_id, reso) = load_parameters_fromdb(opts)
        hic_data = load_hic_data_from_reads(mreads, reso)
        hic_data.bads = dict((int(l.strip()), True) for l in open(bad_co))
        hic_data.bias = dict((int(l.split()[0]), float(l.split()[1]))
                             for l in open(biases))

    exp = crm.experiments[0]
    opts.beg, opts.end = opts.beg or 1, opts.end or exp.size

    # prepare output folders
    batch_job_hash = digest_parameters(opts, get_md5=True , extra=[
        'maxdist', 'upfreq', 'lowfreq', 'scale', 'dcutoff',
        'nmodels_run', 'job_list', 'rand', 'nmodels', 'nkeep', 'optimize',
        'optimization_id', 'cpus', 'workdir', 'matrix', 'ori_beg', 'ori_end'])

    mkdir(path.join(opts.workdir, '06_model'))
    outdir = path.join(opts.workdir, '06_model',
                       '%s_chr%s_%s-%s' % (batch_job_hash,
                                           opts.crm, opts.beg, opts.end))
    mkdir(outdir)

    # in case we are not going to run
    if opts.job_list:
        job_file_handler = open(path.join(
            outdir, 'job_list_%s.q' % ('optimization' if
                                       opts.optimize else 'modeling')), 'w')
    else:
        job_file_handler = None

    ###############
    # Optimization
    if opts.optimize:
        optimization(exp, opts, job_file_handler, outdir)
        finish_time = time.localtime()
        print('\n optimization done')
        # correlate all optimization and get best set of parameters

    if not (opts.optimize and opts.job_list):
        optpar, results = correlate_models(opts, outdir, exp)
    else:
        results = []

    ###########
    # Modeling
    if not opts.optimize:
        big_run(exp, opts, job_file_handler, outdir, optpar)

    finish_time = time.localtime()

    # save all job information to sqlite DB
    save_to_db(opts, outdir, results, batch_job_hash,
               launch_time, finish_time)


def save_to_db(opts, outdir, results, batch_job_hash, 
               launch_time, finish_time):
    if 'tmpdb' in opts and opts.tmpdb:
        # check lock
        while path.exists(path.join(opts.workdir, '__lock_db')):
            time.sleep(0.5)
        # close lock
        open(path.join(opts.workdir, '__lock_db'), 'a').close()
        # tmp file
        dbfile = opts.tmpdb
        try:  # to copy in case read1 was already mapped for example
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
                       type='table' AND name='MODELED_REGIONs'""")
        if not cur.fetchall():
            cur.execute("""
        create table MODELED_REGIONs
           (Id integer primary key,
            PATHid int,
            PARAM_md5 text,
            RESO int,
            BEG int,
            END int,
            unique (PARAM_md5))""")
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='MODELs'""")
        if not cur.fetchall():
            cur.execute("""
        create table MODELs
           (Id integer primary key,
            REGIONid int,
            JOBid int,
            OPTPAR_md5 text,
            MaxDist int,
            UpFreq int,
            LowFreq int,
            Scale int,
            Cutoff int,
            Nmodels int,
            Kept int,
            Correlation int)""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            # In case optimization or modeling  is split in different computers
            param_hash = digest_parameters(opts, get_md5=True)
            cur.execute("""
    insert into JOBs
     (Id  , Parameters, Launch_time, Finish_time,    Type, Parameters_md5)
    values
     (NULL,       '%s',        '%s',        '%s',    '%s',           '%s')
     """ % ((parameters, time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
             time.strftime("%d/%m/%Y %H:%M:%S", finish_time),
             (('PRE_' if opts.job_list else '') +
              ('OPTIM' if opts.optimize else 'MODEL')), param_hash)))
        except lite.IntegrityError:
            pass
        ##### STORE OPTIMIZATION RESULT
        jobid = get_jobid(cur)
        add_path(cur, outdir, 'DIR', jobid, opts.workdir)
        pathid = get_path_id(cur, outdir, opts.workdir)
        # models = compile_models(opts, outdir, exp=exp, ngood=opts.nkeep)
        ### STORE GENERAL OPTIMIZATION INFO
        try:
            cur.execute("""
            insert into MODELED_REGIONs
            (Id  , PATHid, PARAM_md5, RESO, BEG, END)
            values
            (NULL,     %d,      "%s",   %d,  %d,  %d)
            """ % (pathid, batch_job_hash, opts.reso,
                   opts.beg, opts.end))
        except lite.IntegrityError:
            pass
        ### STORE EACH OPTIMIZATION
        cur.execute("SELECT Id from MODELED_REGIONs where PARAM_md5='%s'" % (
            batch_job_hash))
        optimid = cur.fetchall()[0][0]
        for m, u, l, d, s in results:
            optpar_md5 = md5('%s%s%s%s%s' %
                             (m, u, l, d, s)).hexdigest()[:12]
            cur.execute(("SELECT Id from MODELs where "
                         "OPTPAR_md5='%s' and REGIONid='%s'") % (
                             optpar_md5, optimid))
            if not cur.fetchall():
                cur.execute("""
                insert into MODELs
                (Id  , REGIONid, JOBid, OPTPAR_md5, MaxDist, UpFreq, LowFreq, Cutoff, Scale, Nmodels, Kept, Correlation)
                values
                (NULL,             %d,    %d,      '%s',      %s,     %s,      %s,     %s,    %s,      %d,   %d,          %f)
                """ % ((optimid, jobid, optpar_md5, m, u, l, d, s,
                        results[(m, u, l, d, s)]['nmodels'],
                        results[(m, u, l, d, s)]['kept'],
                        results[(m, u, l, d, s)]['corr'])))
            else:
                cur.execute(("update MODELs "
                             "set Nmodels = %d, Kept = %d, Correlation = %f "
                             "where "
                             "OPTPAR_md5='%s' and REGIONid='%s'") % (
                                 results[(m, u, l, d, s)]['nmodels'],
                                 results[(m, u, l, d, s)]['kept'],
                                 results[(m, u, l, d, s)]['corr'],
                                 optpar_md5, optimid))

        ### MODELING
        if not opts.optimization_id:
            cur.execute("SELECT Id from MODELED_REGIONs")
            optimid = cur.fetchall()[0]
            if len(optimid) > 1:
                raise IndexError("ERROR: more than 1 optimization in folder "
                                 "choose with 'tadbit describe' and "
                                 "--optimization_id")
            optimid = optimid[0]
        else:
            cur.execute("SELECT Id from MODELED_REGIONs where Id=%d" % (
                opts.optimization_id))
            optimid = cur.fetchall()[0][0]
        
    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
    # release lock
    try:
        remove(path.join(opts.workdir, '__lock_db'))
    except OSError:
        pass


def compile_models(opts, outdir, exp=None, ngood=None, wanted=None):
    if exp:
        zscores, _, zeros = exp._sub_experiment_zscore(opts.beg, opts.end)
        zeros = tuple([i not in zeros for i in xrange(opts.end - opts.beg + 1)])
    models = {}
    for cfg_dir in listdir(outdir):
        if not cfg_dir.startswith('cfg_'):
            continue
        _, m, u, l, s = cfg_dir.split('_')
        muls = m, u, l, s
        m, u, l, s = int(m), float(u), float(l), float(s)
        if wanted and muls != wanted:
            continue
        for fmodel in listdir(path.join(outdir, cfg_dir)):
            if not fmodel.startswith('models_'):
                continue
            if muls not in models:
                models[muls] = load_structuralmodels(path.join(
                    outdir, cfg_dir, fmodel))
            else:
                sm = load(open(path.join(outdir, cfg_dir, fmodel)))
                for k in sm['config']:
                    if not isinstance(sm['config'][k], float):
                        continue
                    sm['config'][k] = float(my_round(sm['config'][k]))
                if models[muls]._config != sm['config']:
                    print 'Different configuration found in this directory:'
                    print '=' * 80
                    print sm['config']
                    print '-' * 80
                    print models[muls]._config
                    print '=' * 80
                    raise Exception('ERROR: clean directory, '
                                    'heterogeneous data')
                models[muls]._extend_models(sm['models'])
                models[muls]._extend_models(sm['bad_models'])
            if exp:
                models[muls].experiment = exp
                models[muls]._zscores   = zscores
                models[muls]._zeros     = zeros
            if ngood:
                models[muls].define_best_models(ngood)
    return models


def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class = lambda prog: HelpFormatter(prog, width=95,
                                                        max_help_position=27)

    glopts = parser.add_argument_group('General options')
    reopts = parser.add_argument_group('Modeling preparation')
    opopts = parser.add_argument_group('Parameter optimization')
    anopts = parser.add_argument_group('Analysis')
    ruopts = parser.add_argument_group('Computation')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool TADbit mapper)''')
    glopts.add_argument('--input_matrix', dest='matrix', metavar="PATH",
                        type=str,
                        help='''In case input was not generated with the TADbit
                        tools''')
    glopts.add_argument('--rand', dest='rand', metavar="INT",
                        type=str, default='1', 
                        help='''[%(default)s] random initial number. NOTE:
                        when running single model at the time, should be
                        different for each run''')
    glopts.add_argument('--nmodels', dest='nmodels', metavar="INT",
                        default=5000, type=int,
                        help=('[%(default)s] number of models to generate for' +
                              ' modeling'))
    glopts.add_argument('--nkeep', dest='nkeep', metavar="INT",
                        default=1000, type=int,
                        help=('[%(default)s] number of models to keep for ' +
                              'modeling'))
    glopts.add_argument('--optimization_id', dest='optimization_id', metavar="INT",
                        type=float, default=None,
                        help="[%(default)s] ID of a pre-run optimization batch job")

    reopts.add_argument('--crm', dest='crm', metavar="NAME",
                        help='chromosome name')
    reopts.add_argument('--beg', dest='beg', metavar="INT", type=float,
                        required=True,
                        help='genomic coordinate from which to start modeling')
    reopts.add_argument('--end', dest='end', metavar="INT", type=float,
                        required=True,
                        help='genomic coordinate where to end modeling')
    reopts.add_argument('-r', '--reso', dest='reso', metavar="INT", type=int,
                        help='resolution of the Hi-C experiment')
    reopts.add_argument('--perc_zero', dest='perc_zero', metavar="FLOAT",
                        type=float, default=90.0)

    opopts.add_argument('--optimize', dest='optimize', 
                        default=False, action="store_true",
                        help='''optimization run, store less info about models''')
    opopts.add_argument('--maxdist', action='store', metavar="LIST",
                        default='400', dest='maxdist',
                        help='range of numbers for maxdist' +
                        ', i.e. 400:1000:100 -- or just a number')
    opopts.add_argument('--upfreq', dest='upfreq', metavar="LIST",
                        default='0',
                        help='range of numbers for upfreq' +
                        ', i.e. 0:1.2:0.3 --  or just a number')
    opopts.add_argument('--lowfreq', dest='lowfreq', metavar="LIST",
                        default='0',
                        help='range of numbers for lowfreq' +
                        ', i.e. -1.2:0:0.3 -- or just a number')
    opopts.add_argument('--scale', dest='scale', metavar="LIST",
                        default="0.01",
                        help='[%(default)s] range of numbers to be test as ' +
                        'optimal scale value, i.e. 0.005:0.01:0.001 -- Can ' +
                        'also pass only one number')
    opopts.add_argument('--dcutoff', dest='dcutoff', metavar="LIST",
                        default="2",
                        help='[%(default)s] range of numbers to be test as ' +
                        'optimal distance cutoff parameter (distance, in ' +
                        'number of beads, from which to consider 2 beads as ' +
                        'being close), i.e. 1:5:0.5 -- Can also pass only one' +
                        ' number')

    anopts.add_argument('--analyze', dest='analyze', 
                        default=False, action="store_true",
                        help='''analyze models.''')
    
    ruopts.add_argument('--nmodels_run', dest='nmodels_run', metavar="INT",
                        default=None, type=int,
                        help='[ALL] number of models to run with this call')
    ruopts.add_argument("-C", "--cpu", dest="cpus", type=int,
                        default=1, help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')
    ruopts.add_argument('--job_list', dest='job_list', action='store_true',
                        default=False,
                        help=('generate a list of commands stored in a file '
                              'named joblist_HASH.q (where HASH is replaced by '
                              'a string specific to the parameters used)'))
    ruopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    parser.add_argument_group(glopts)


def check_options(opts):
    # check resume
    if not path.exists(opts.workdir):
        warn('ERROR: workdir not found, creating it')
        mkdir(opts.workdir)
        # write version log
        vlog_path = path.join(opts.workdir, 'TADbit_and_dependencies_versions.log')
        dependencies = get_dependencies_version()
        if not path.exists(vlog_path) or open(vlog_path).readlines() != dependencies:
            print('Writing versions of TADbit and dependencies')
            vlog = open(vlog_path, 'w')
            vlog.write(dependencies)
            vlog.close()
    # do the division to bins
    try:
        opts.ori_beg = opts.beg
        opts.ori_end = opts.end
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

    if opts.matrix:
        opts.matrix  = path.abspath(opts.matrix)
    opts.workdir = path.abspath(opts.workdir)

    mkdir(opts.workdir)
    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)


def load_parameters_fromdb(opts):
    """
    TODO: should load optimization specific parameters like nkeep nmodels etc.. 
          to ensure they are always the same.
    """
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
        select distinct Path, PATHs.Id from PATHs
        where paths.jobid = %s and paths.Type = 'BAD_COLUMNS'
        """ % parse_jobid)
        bad_co, bad_co_id  = cur.fetchall()[0]
        cur.execute("""
        select distinct Path, PATHs.Id from PATHs
        where paths.jobid = %s and paths.Type = 'BIASES'
        """ % parse_jobid)
        biases, biases_id = cur.fetchall()[0]
        cur.execute("""
        select distinct Path, PATHs.Id from PATHs
        inner join NORMALIZE_OUTPUTs on PATHs.Id = NORMALIZE_OUTPUTs.Input
        where NORMALIZE_OUTPUTs.JOBid = %d;
        """ % parse_jobid)
        mreads, mreads_id = cur.fetchall()[0]
        cur.execute("""
        select distinct Resolution from NORMALIZE_OUTPUTs
        where NORMALIZE_OUTPUTs.JOBid = %d;
        """ % parse_jobid)
        reso = int(cur.fetchall()[0][0])
        return (bad_co, bad_co_id, biases, biases_id,
                mreads, mreads_id, reso)


def load_hic_data(opts):
    """
    Load Hi-C data
    """
    # Start reading the data
    crm = Chromosome(opts. crm)  # Create chromosome object
    try:
        global HIC_DATA
        HIC_DATA = False # we are reading a normalized matrix
        hic = optimal_reader(open(opts.matrix), resolution=opts.reso)
        crm.add_experiment('test', exp_type='Hi-C', resolution=opts.reso,
                           norm_data=hic)
    except:
         crm.add_experiment('test', exp_type='Hi-C', resolution=opts.reso,
                            norm_data=opts.matrix)
    # TODO: if not bad columns:...
    crm.experiments[-1].filter_columns(perc_zero=opts.perc_zero)
    if opts.beg > crm.experiments[-1].size:
        raise Exception('ERROR: beg parameter is larger than chromosome size.')
    if opts.end > crm.experiments[-1].size:
        print ('WARNING: end parameter is larger than chromosome ' +
               'size. Setting end to %s.\n' % (crm.experiments[-1].size *
                                               opts.reso))
        opts.end = crm.experiments[-1].size
    return crm
