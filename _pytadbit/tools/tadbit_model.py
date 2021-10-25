"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""
from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
from argparse                         import HelpFormatter
from os                               import path, remove, system, rename, makedirs
from string                           import ascii_letters
from math                             import ceil
from random                           import random
from shutil                           import copyfile
from itertools                        import product
from warnings                         import warn
from pickle                           import dump, HIGHEST_PROTOCOL
from hashlib                          import md5
from functools                        import partial
from multiprocessing                  import cpu_count, TimeoutError, Pool
from multiprocessing.dummy            import Pool as ThreadPool
import sqlite3 as lite
import subprocess
import time
import sys
import logging
import collections

from numpy                            import arange

from pytadbit                         import load_structuralmodels
from pytadbit.modelling.impoptimizer  import IMPoptimizer
from pytadbit                         import Chromosome
from pytadbit.utils.file_handling     import mkdir
from pytadbit.utils.extraviews        import nicer
from pytadbit.utils.sqlite_utils      import get_path_id, add_path, get_jobid
from pytadbit.utils.sqlite_utils      import digest_parameters, retry
from pytadbit                         import get_dependencies_version
from pytadbit.parsers.hic_parser      import read_matrix

try:
    basestring
except NameError:
    basestring = str

DESC = ("Generates 3D models given an input interaction matrix and a set of "
        "input parameters")

## Define analysis actions:
actions = {0  : "do nothing",
           1  : "optimization plot",
           2  : "correlation real/models",
           3  : "z-score plot",
           4  : "constraints",
           5  : "objective function",
           6  : "centroid",
           7  : "consistency",
           8  : "density",
           9  : "contact map",
           10 : "walking angle",
           11 : "persistence length",
           12 : "accessibility",
           13 : "interaction"}

def abortable_worker(func, *args, **kwargs):
    timeout = kwargs.get('timeout', None)
    p = ThreadPool(1)
    res = p.apply_async(func, args=args)
    try:
        out = res.get(timeout)  # Wait timeout seconds for func to complete.
        return out
    except TimeoutError:
        print("Model took more than %s seconds to complete ... canceling" % str(timeout))
        p.terminate()
        raise
    except:
        print("Unknown error with process")
        p.terminate()
        raise

def convert_from_unicode(data):
    if isinstance(data, basestring):
        return str(data)
    if isinstance(data, collections.Mapping):
        return dict(list(map(convert_from_unicode, iter(data.items()))))
    if isinstance(data, collections.Iterable):
        return type(data)(list(map(convert_from_unicode, data)))
    return data

def prepare_distributed_jobs(exp, opts, m, u, l, s, outdir, batch_job_hash):
    zscores, values, zeros = exp._sub_experiment_zscore(opts.beg - opts.offset + 1,
                                                        opts.end - opts.offset)
    zeros = tuple([i not in zeros for i in range(opts.end - opts.beg)])
    nloci = opts.end - opts.beg
    if exp.norm and exp.norm[0].chromosomes:
        coords = []
        tot = 0
        chrs = []
        chrom_offset_start = opts.beg
        chrom_offset_end = 0
        for k, v in exp.norm[0].chromosomes.items():
            tot += v
            if opts.beg > tot:
                chrom_offset_start = opts.beg - tot
            if opts.end <= tot:
                chrom_offset_end = tot - opts.end
                chrs.append(k)
                break
            if opts.beg < tot and opts.end >= tot:
                chrs.append(k)
        for k in chrs:
            coords.append({'crm'  : k,
                  'start': 1,
                  'end'  : exp.norm[0].chromosomes[k]})
        coords[0]['start'] = chrom_offset_start + 1
        coords[-1]['end'] -= chrom_offset_end
    else:
        coords = {"crm"  : opts.crm,
                  "start": opts.beg + 1,
                  "end"  : opts.end}

    optpar = {'maxdist': float(m),
              'upfreq' : float(u),
              'lowfreq': float(l),
              'scale'  : float(s),
              'kforce' : 5}

    muls = tuple(map(my_round, (m, u, l, s)))
    dirname = path.join(outdir, 'cfg_%s_%s_%s_%s' % muls)

    n_jobs = int(ceil(opts.nmodels/opts.nmodels_per_job))
    n_last = n_jobs*opts.nmodels_per_job - opts.nmodels
    paramsfile = path.join(dirname,'_tmp_common_params.pickle')
    tmp_params = open(paramsfile, 'wb')
    dump(exp, tmp_params)
    dump(zscores, tmp_params)
    dump(zeros, tmp_params)
    dump(values, tmp_params)
    dump(opts, tmp_params)
    dump(optpar, tmp_params)
    dump(coords, tmp_params)
    tmp_params.close()
    for n_job in range(n_jobs):
        nmodels_per_job = opts.nmodels_per_job
        if n_job == n_jobs - 1:
            nmodels_per_job -= n_last
        job_dir = path.join(dirname,'_tmp_results_%s_%s_%s' % (n_job, opts.rand, batch_job_hash))
        mkdir(job_dir)
        scriptname = path.join(job_dir,'_tmp_optim.py')
        tmp = open(scriptname, 'w')
        tmp.write('''
from pickle import load, dump
from os                               import path
from pytadbit.modelling.imp_modelling import generate_3d_models

params_file = open("%s","rb")
exp = load(params_file)
zscores = load(params_file)
zeros = load(params_file)
values = load(params_file)
opts = load(params_file)
optpar = load(params_file)
coords = load(params_file)
params_file.close()

try:
    models = generate_3d_models(zscores, opts.reso, %s,
                                        values=values, n_models=%s,
                                        n_keep=%s,
                                        n_cpus=opts.cpus_per_job, keep_all=True,
                                        start=int(opts.rand)+%s, container=None,
                                        config=optpar, coords=coords, experiment=exp,
                                        zeros=zeros)

    models.save_models(path.join("%s",'results.models'),minimal=%s)
except Exception as e:
    print(e)
    open(path.join("%s",'failed.flag'), 'a').close()
    ''' % (paramsfile, nloci, nmodels_per_job, nmodels_per_job,
           n_job*opts.nmodels_per_job, job_dir,
           '()' if n_job==0 else '["restraints", "zscores", "original_data"]',
           job_dir))

        tmp.close()

def run_distributed_job(job_dir, script_cmd , script_args):

    scriptname = path.join(job_dir,'_tmp_optim.py')
    logname = path.join(job_dir,'_tmp_log.log')
    with open(logname, 'a') as f:
        f.write('Log %s\n' % job_dir)
        f.flush()
        subprocess.Popen([script_cmd, script_args, scriptname], stdout = f, stderr = f,
                         universal_newlines=True)
    f.close()
    results_file = path.join(job_dir,'results.models')
    failed_flag = path.join(job_dir,'failed.flag')
    while not (path.exists(results_file) or path.exists(failed_flag)):
        time.sleep(1)

def run_distributed_jobs(opts, m, u, l, s, outdir, batch_job_hash, job_file_handler = None,
                         exp = None, script_cmd = 'python',
                         script_args = '', verbose = True):

    muls = tuple(map(my_round, (m, u, l, s)))
    dirname = path.join(outdir, 'cfg_%s_%s_%s_%s' % muls)
    modelsfile = path.join(outdir, dirname, 'models_%s_%s_%s_%s.models' % muls)

    if path.exists(modelsfile):
        models = load_structuralmodels(modelsfile)
    else:
        n_jobs = int(ceil(opts.nmodels/opts.nmodels_per_job))
        pool = Pool(processes=opts.cpus, maxtasksperchild=opts.concurrent_jobs)
        jobs = {}
        for n_job in range(n_jobs):
            job_dir = path.join(dirname, '_tmp_results_%s_%s_%s' % (n_job, opts.rand, batch_job_hash))
            if job_file_handler:
                scriptname = path.join(job_dir,'_tmp_optim.py')
                job_file_handler.write('%s %s %s\n'%(script_cmd, script_args, scriptname))
            else:
                jobs[n_job] = partial(abortable_worker, run_distributed_job)
                pool.apply_async(jobs[n_job], args=(job_dir, script_cmd , script_args))
        pool.close()
        pool.join()

        if job_file_handler:
            return None, None

        models = None
        for n_job in range(n_jobs):
            try:
                job_dir = path.join(dirname, '_tmp_results_%s_%s_%s' % (n_job, opts.rand, batch_job_hash))
                results_file = path.join(job_dir,'results.models')
                if path.isfile(results_file):
                    results = load_structuralmodels(results_file)
                    if models:
                        models._extend_models(results, nbest=len(models)+len(results))
                    else:
                        models = results
                failed_flag = path.join(job_dir,'failed.flag')
                logname = path.join(job_dir,'_tmp_log.log')
                if path.isfile(failed_flag):
                    f = open(logname, 'r')
                    logging.error(f.read())
                    f.close()
                system('rm -rf %s' % (job_dir))
            except TimeoutError:
                logging.info("Model took more than %s seconds to complete ... canceling"
                             % str(opts.timeout_job))
                jobs[n_job].cancel()
            except Exception as error:
                logging.info("Function raised %s" % error)
                jobs[n_job].cancel()
        paramsfile = path.join(dirname,'_tmp_common_params.pickle')
        system('rm %s' % (paramsfile))

        models.define_best_models(opts.nkeep)
        if isinstance(models.description['start'],list):
            models.description['start'] = [(st + opts.matrix_beg * opts.reso) 
                                        for st in models.description['start']]
            models.description['end'] = [(st + opts.matrix_beg * opts.reso) 
                                        for st in models.description['end']]
        else:
            models.description['start'] += opts.matrix_beg * opts.reso
            models.description['end'] += opts.matrix_beg * opts.reso  
        for model in models:
            model['description']['start'] = models.description['start']
            model['description']['end'] = models.description['end']
        models.save_models(modelsfile)

    num=1
    results_corr = {}
    cuts = dict([(d, int(d * opts.reso * float(s))) for d in opts.dcutoff])
    for d, cut in sorted(cuts.items()):
        try:
            result = models.correlate_with_real_data(
                cutoff=cut, corr=opts.corr)[0]
        except Exception as e:
            logging.info('  SKIPPING correlation: %s' % e)
            result = 0
        name = tuple(map(my_round, (m, u, l, d, s)))
        if verbose:
            logging.info(('%8s/%-6s %6s %7s %7s %6s %7s | %.4f' %
                   (num, len(cuts), u, l, m, s, d, result)))
        num += 1
        results_corr[name] = {'corr'   : result,
                              'nmodels': (len(models) +
                                          len(models._bad_models)),
                              'kept'   : len(models)}
    if exp: #Store more data for the models
        out = open(path.join(outdir, dirname, 'constraints.txt'),
                   'w')
        out.write('# Harmonic\tpart1\tpart2\tdist\tkforce\n')
        out.write('\n'.join(['%s\t%s\t%s\t%.1f\t%.3f' % (
            harm, p1, p2, dist, kforce)
                             for (p1, p2), (harm, dist, kforce)
                             in models._restraints.items()]) + '\n')
        out.close()

    return results_corr, modelsfile


def my_round(num, val=4):
    num = round(float(num), val)
    return str(int(num) if num == int(num) else num)


def optimization_distributed(exp, opts, outdir, batch_job_hash, job_file_handler = None,
                             script_cmd = 'python', script_args = '', verbose=True):
    logging.info('\nOptimizing parameters...')
    if verbose:
        logging.info('\n\n# %13s %6s %7s %7s %6s %7s %7s\n' % (
            "Optimization", "UpFreq", "LowFreq", "MaxDist",
            "scale", "cutoff", "| Correlation"))
    for m, u, l, s in product(opts.maxdist, opts.upfreq, opts.lowfreq, opts.scale):
        m, u, l, s = list(map(my_round, (m, u, l, s)))
        muls = tuple((m, u, l, s))
        cfgfolder = path.join(outdir, 'cfg_%s_%s_%s_%s' % muls)
        modelsfile = path.join(cfgfolder,'models_%s_%s_%s_%s.models' % muls)
        if not path.exists(modelsfile):
            mkdir(cfgfolder)
            prepare_distributed_jobs(exp, opts, m, u, l, s, outdir, batch_job_hash)

    # get the best combination
    best = ({'corr': 0}, [0, 0, 0, 0, 0])
    results = {}
    for m, u, l, s in product(opts.maxdist, opts.upfreq, opts.lowfreq, opts.scale):
        m, u, l, s = list(map(my_round, (m, u, l, s)))
        muls_results, _ = run_distributed_jobs(opts, m, u, l, s, outdir, batch_job_hash, job_file_handler, 
                            script_cmd = script_cmd, script_args = script_args, verbose=verbose)
        if muls_results:
            results.update(muls_results)
    if not job_file_handler:
        for m, u, l, d, s in results:
            if results[(m, u, l, d, s)]['corr'] > best[0]['corr']:
                best = results[(m, u, l, d, s)], [u, l, m, s, d]
    else:
        return None, None
    if verbose:
        logging.info( '\nBest combination:')
        logging.info('  %5s     %6s %7s %7s %6s %6s %.4f\n' % tuple(
            ['=>'] + best[1] + [best[0]['corr']]))

    u, l, m, s, d = best[1]
    optpar = {'maxdist': m,
              'upfreq' : u,
              'lowfreq': l,
              'scale'  : s,
              'kforce' : 5}

    return optpar, results

def run_distributed(exp, batch_job_hash, opts, outdir, optpar,
                    job_file_handler = None,
                    script_cmd = 'python', script_args = ''):
    m, u, l, s = (optpar['maxdist'],
                  optpar['upfreq' ],
                  optpar['lowfreq'],
                  optpar['scale'  ])
    muls = tuple(map(my_round, (m, u, l, s)))
    cfgfolder = path.join(outdir, 'cfg_%s_%s_%s_%s' % muls)
    # if path.exists(cfgfolder):
    #     logging.info( '\nJob already run. Please use tadbit clean if you want to redo it.')
    #     return []
    mkdir(cfgfolder)
    prepare_distributed_jobs(exp, opts, m, u, l, s, outdir, batch_job_hash)
    results, modelsfile = run_distributed_jobs(opts, m, u, l, s, outdir, batch_job_hash,
                                               job_file_handler=job_file_handler,
                                               exp=exp, script_cmd=script_cmd,
                                               script_args=script_args)
    if not job_file_handler:
        rename(modelsfile, path.join(outdir, batch_job_hash + '%s_%s.models' % (batch_job_hash, opts.rand)))
    return results

def run(opts):
    check_options(opts)

    launch_time = time.localtime()

    # prepare output folders
    batch_job_hash = digest_parameters(opts, get_md5=True , extra=[
        'maxdist', 'upfreq', 'lowfreq', 'scale', 'dcutoff',
        'job_list', 'rand', 'optimize',
        'optimization_id', 'cpus', 'workdir', 'matrix'])

    # write log
    if opts.optimize:
        log_format = '[OPTIMIZATION {}_{}_{}_{}_{}]   %(message)s'.format(
            opts.maxdist, opts.upfreq, opts.lowfreq, opts.scale, opts.dcutoff)
    elif opts.analyze:
        log_format = '[ANALYZE]   %(message)s'
    else:
        log_format = '[DEFAULT]   %(message)s'
    try:
        logging.basicConfig(filename=path.join(opts.workdir, batch_job_hash + '.log'),
                            level=logging.INFO, format=log_format)
    except IOError:
        logging.basicConfig(filename=path.join(opts.workdir, batch_job_hash + '.log2'),
                            level=logging.INFO, format=log_format)
    logging.getLogger().addHandler(logging.StreamHandler())

    if opts.optimize or opts.model:
        # load data
        if opts.matrix:
            opts.matrix = convert_from_unicode(path.realpath(opts.matrix))
        else:
            opts.matrix, _ = load_matrix_path_fromdb(opts)
            opts.matrix = convert_from_unicode(path.join(opts.workdir, opts.matrix))

        crm = load_hic_data(opts)
        exp = crm.experiments[0]
        opts.beg, opts.end = opts.beg or 0, opts.end or exp.size

        name = '{0}_{1}_{2}'.format(opts.crm if opts.crm else 'all',
                                    int(opts.beg + 1), int(opts.end))

        mkdir(path.join(opts.workdir, '07_model'))
        outdir = path.join(opts.workdir, '07_model',
                           '%s_%s_%s-%s' % (batch_job_hash,
                                opts.crm if opts.crm else 'all',
                                (opts.beg + 1) if opts.beg is not None else '',
                                opts.end if opts.end else ''))
        if path.exists(outdir):
            logging.info ('     o WARNING: reusing existing folder, \
                                    please check if you need to remove')
        mkdir(outdir)

        if opts.crm:
            logging.info('''
        %s%s

          - Region: Chromosome %s from %d to %d at resolution %s (%d particles)\n''' % (
              'Preparing ' if opts.job_list else '',
                   ('Optimization\n        ' + '*' * (21 if opts.job_list else 11))
                   if opts.optimize else
                   ('Modeling\n' + '*' * (18 if opts.job_list else 8)),
                   opts.crm, opts.ori_beg, opts.ori_end, nicer(opts.reso),
                   opts.end - opts.beg))
        else:
            logging.info('''
        %s%s

          - Region: full genome at resolution %s (%d particles)\n''' % (
              'Preparing ' if opts.job_list else '',
                   ('Optimization\n        ' + '*' * (21 if opts.job_list else 11))
                   if opts.optimize else
                   ('Modeling\n' + '*' * (18 if opts.job_list else 8)),
                   nicer(opts.reso),
                   opts.end - opts.beg))
        # in case we are not going to run
        if opts.job_list:
            job_file_handler = open(path.join(
                outdir, 'job_list_%s.q' % ('optimization' if
                                           opts.optimize else 'modeling')), 'w')
        else:
            job_file_handler = None

        optpar = None
        results = []
        ###############
        # Optimization
        if opts.optimize:
            logging.info ('     o Optimizing parameters')
            optpar, results = optimization_distributed(exp, opts, outdir, batch_job_hash,
                                                       job_file_handler = job_file_handler,
                                                       script_cmd = opts.script_cmd,
                                                       script_args = opts.script_args)
            if not opts.job_list and "optimization plot" in opts.analyze_list:
                if optpar:
                    optimizer = IMPoptimizer(exp, opts.beg - opts.offset + 1, 
                                             opts.end - opts.offset)
                    optimizer.scale_range    = [i for i in opts.scale]
                    optimizer.kbending_range = [0.0]
                    optimizer.maxdist_range  = [i for i in opts.maxdist]
                    optimizer.lowfreq_range  = [i for i in opts.lowfreq]
                    optimizer.upfreq_range   = [i for i in opts.upfreq]
                    optimizer.dcutoff_range  = [i for i in opts.dcutoff]
                    optimizer.results = dict((
                        (float(s),0.0,float(m),float(l), float(u), float(d)),
                        results[(m,u,l,d,s)]['corr'])
                                             for m,u,l,d,s in results)
                    optimizer.plot_2d(show_best=20,
                                savefig="%s/optimal_params_%s.%s" % (
                                    outdir,''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]),
                                    opts.fig_format))
                logging.info('\n optimization done')

        if opts.model:
            if not opts.force:
                logging.info( '     o Loading optimized parameters')
                if not optpar:
                    optpar = load_optpar_fromdb(opts)
            else:
                optpar = {'maxdist': opts.maxdist[0],
                          'upfreq' : opts.upfreq[0],
                          'lowfreq': opts.lowfreq[0],
                          'scale'  : opts.scale[0],
                          'kforce' : 5}

            ###########
            # Modeling
            results = run_distributed(exp, batch_job_hash, opts, outdir, optpar,
                            job_file_handler = job_file_handler,
                            script_cmd = opts.script_cmd, script_args = opts.script_args)

        finish_time = time.localtime()
        # save all job information to sqlite DB
        save_to_db(opts, outdir, results, batch_job_hash,
               launch_time, finish_time)

    if opts.analyze and not opts.job_list:
        outdir, _ = load_models_path_fromdb(opts)
        batch_job_hash, opts.crm, beg_end = list(map(str,(outdir.split('/')[-1]).split('_')))
        opts.beg,opts.end = list(map(int,beg_end.split('-')))
        name = '{0}_{1}_{2}'.format(opts.crm if opts.crm else 'all',
                                    int(opts.beg), int(opts.end))
        outdir = path.join(opts.workdir,outdir)
        models = load_structuralmodels(path.join(outdir, '%s_%s.models' % (batch_job_hash, opts.rand)))
        opts.reso = models.description['resolution']

        logging.info('''
        %s

          - Region: Chromosome %s from %d to %d at resolution %s (%d particles)
            ''' % ('Analysis',
                   opts.crm, opts.beg, opts.end, nicer(opts.reso),
                   opts.end - opts.beg + 1))

        dcutoff = int((opts.dcutoff[0] if opts.dcutoff else 2) *
                  models._config['scale']   *
                  models.resolution)
        if "correlation real/models" in opts.analyze_list:
            # Calculate the correlation coefficient between a set of kept models and
            # the original HiC matrix
            logging.info("\tCorrelation with data...")
            rho, pval = models.correlate_with_real_data(
                cutoff=dcutoff, corr=opts.corr,
                savefig=path.join(outdir, batch_job_hash + '_corre_real.' + opts.fig_format),
                plot=True)
            logging.info("\t Correlation coefficient: %s [p-value: %s]", rho, pval)

        if "z-score plot" in opts.analyze_list:
            # zscore plots
            logging.info("\tZ-score plot...")
            models.zscore_plot(
                savefig=path.join(outdir, batch_job_hash + '_zscores.' + opts.fig_format))

        # Cluster models based on structural similarity
        logging.info("\tClustering all models into sets of structurally similar" +
                     " models...")
        ffact    = 0.95 # Fraction of particles that are within the dcutoff value
        clcutoff = dcutoff - 50 # RMSD cut-off to consider two models equivalent(nm)
        for ffact in [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5]:
            logging.info('   fact = ' + str(ffact))
            for clcutoff in [dcutoff / 2 , dcutoff, dcutoff * 1.5]:
                try:
                    logging.info('      cutoff = ' + str(clcutoff))
                    models.cluster_models(fact=ffact, dcutoff=clcutoff,
                                          n_cpus=int(opts.cpus))
                    break
                except:
                    continue
            else:
                continue
            break
        logging.info("\tSaving again the models this time with clusters...")
        models.save_models(path.join(outdir, '%s_%s.models' % (batch_job_hash, opts.rand)))
        # Plot the clustering
        try:
            models.cluster_analysis_dendrogram(
                color=True, savefig=path.join(
                    outdir, batch_job_hash + '_clusters.' + opts.fig_format))
        except:
            logging.info("\t\tWARNING: plot for clusters could not be made...")

        if not opts.not_write_json:
            models.write_json(path.join(outdir, batch_job_hash + '.json'),
                              title = opts.project+' '+name if opts.project else name,
                              infer_unrestrained=True)

        if not (opts.not_write_xyz and opts.not_write_cmm):
            # Save the clustered models into directories for easy visualization with
            # Chimera (http://www.cgl.ucsf.edu/chimera/)
            # Move into the cluster directory and run in the prompt
            # "chimera cl_1_superimpose.cmd"
            logging.info("\t\tWriting models, list and chimera files...")
            for cluster in models.clusters:
                logging.info("\t\tCluster #{0} has {1} models {2}".format(
                    cluster, len(models.clusters[cluster]),
                    models.clusters[cluster]))
                if not path.exists(path.join(
                    outdir, 'models', 'cl_' + str(cluster))):
                    makedirs(path.join(
                        outdir, 'models', 'cl_' + str(cluster)))
                if not opts.not_write_xyz:
                    models.write_xyz(directory=path.join(
                        outdir, 'models', 'cl_' + str(cluster)),
                                     cluster=cluster)
                if not opts.not_write_cmm:
                    models.write_cmm(directory=path.join(
                        outdir, 'models', 'cl_' + str(cluster)),
                                     cluster=cluster)
                # Write list file
                clslstfile = path.join(
                    outdir,
                    'models', 'cl_{}.lst'.format(str(cluster)))
                out = open(clslstfile,'w')
                for model_n in models.clusters[cluster]:
                    out.write("model.{0}\n".format(model_n))
                out.close()
                if not opts.not_write_cmm:
                    # Write chimera file
                    clschmfile = path.join(
                        outdir, 'models',
                        'cl_{}_superimpose.cmd'.format(str(cluster)))
                    out = open(clschmfile, 'w')
                    out.write("open " + " ".join(["cl_{0}/model.{1}.cmm".format(
                        cluster, model_n) for model_n in models.clusters[cluster]]))
                    out.write("\nlabel; represent wire; ~bondcolor\n")
                    for i in range(1, len(models.clusters[cluster]) + 1):
                        out.write("match #{0} #0\n".format(i-1))
                    out.close()
            # same with singletons
            singletons = [m['rand_init'] for m in models if m['cluster']=='Singleton']
            logging.info("\t\tSingletons has %s models %s", len(singletons),
                         singletons)
            if not path.exists(path.join(
                outdir, 'models', 'Singletons')):
                makedirs(path.join(
                    outdir, 'models', 'Singletons'))
            if not opts.not_write_xyz:
                models.write_xyz(directory=path.join(
                    outdir, 'models', 'Singletons'),
                                 models=singletons)
            if not opts.not_write_cmm:
                models.write_cmm(directory=path.join(
                    outdir, 'models', 'Singletons'),
                                 models=singletons)
            # Write best model and centroid model
            models[models.centroid_model()].write_cmm(
                directory=path.join(outdir, 'models'),
                filename='centroid.cmm')
            models[models.centroid_model()].write_xyz(
                directory=path.join(outdir, 'models'),
                filename='centroid.xyz')
            models[0].write_cmm(
                directory=path.join(outdir, 'models'),
                filename='best.cmm')
            models[0].write_xyz(
                directory=path.join(outdir, 'models'),
                filename='best.xyz')
            # Write list file
            clslstfile = path.join(
                outdir, 'models', 'Singletons.lst')
            out = open(clslstfile,'w')
            for model_n in singletons:
                out.write("model.{0}\n".format(model_n))
            out.close()
            if not opts.not_write_cmm:
                # Write chimera file
                clschmfile = path.join(
                    outdir, 'models', 'Singletons_superimpose.cmd')
                out = open(clschmfile, 'w')
                out.write("open " + " ".join(["Singletons/model.{0}.cmm".format(
                    model_n) for model_n in singletons]))
                out.write("\nlabel; represent wire; ~bondcolor\n")
                for i in range(1, len(singletons) + 1):
                    out.write("match #{0} #0\n".format(i-1))
                out.close()

        if "objective function" in opts.analyze_list:
            logging.info("\tPlotting objective function decay for vbest model...")
            models.objective_function_model(
                0, log=True, smooth=False,
                savefig=path.join(outdir, batch_job_hash + '_obj-func.' + opts.fig_format))

        if "centroid" in opts.analyze_list:
            # Get the centroid model of cluster #1
            logging.info("\tGetting centroid...")
            centroid = models.centroid_model(cluster=1)
            logging.info("\t\tThe model centroid (closest to the average) " +
                         "for cluster 1 is: {}".format(centroid))

        if "consistency" in opts.analyze_list:
            # Calculate a consistency plot for all models in cluster #1
            logging.info("\tGetting consistency data...")
            models.model_consistency(
                cluster=1, cutoffs=list(range(50, dcutoff + 50, 50)),
                savefig =path.join(outdir, batch_job_hash + '_consistency.' + opts.fig_format),
                savedata=path.join(outdir, batch_job_hash + '_consistency.dat'))

        if "density" in opts.analyze_list:
            # Calculate a DNA density plot
            logging.info("\tGetting density data...")
            models.density_plot(
                error=True, steps=(1,3,5,7),
                savefig =path.join(outdir, batch_job_hash + '_density.' + opts.fig_format),
                savedata=path.join(outdir, batch_job_hash + '_density.dat'))

        if "contact map" in opts.analyze_list:
            # Get a contact map at cut-off of 150nm for cluster #1
            logging.info("\tGetting a contact map...")
            models.contact_map(
                cluster=1, cutoff=dcutoff,
                savedata=path.join(outdir, batch_job_hash + '_contact.dat'))

        if "walking angle" in opts.analyze_list:
            # Get Dihedral angle plot for cluster #1
            logging.info("\tGetting angle data...")
            models.walking_angle(
                cluster=1, steps=(1,5),
                savefig = path.join(outdir, batch_job_hash + '_wang.' + opts.fig_format),
                savedata= path.join(outdir, batch_job_hash + '_wang.dat'))

        if "persistence length" in opts.analyze_list:
            # Get persistence length of all models
            logging.info("\tGetting persistence length data...")
            pltfile = path.join(outdir, batch_job_hash + '_pL.dat')
            f = open(pltfile,'w')
            f.write('#Model_Number\tpL\n')
            for model in models:
                try:
                    f.write('%s\t%.2f\n' % (model["rand_init"],
                                            model.persistence_length()))
                except:
                    sys.stderr.write('WARNING: failed to compute persistence ' +
                         'length for model %s' % model["rand_init"])

        if "accessibility" in opts.analyze_list:
            # Calculate a DNA density plot
            logging.info("\tGetting accessibility data...")
            radius = 75   # Radius of an object to calculate accessibility
            nump   = 30   # number of particles (resolution)
            logging.info("\tGetting accessibility data (this can take long)...")
            models.accessibility(radius, nump=nump,
                error=True,
                savefig =path.join(outdir, batch_job_hash + '_accessibility.' + opts.fig_format),
                savedata=path.join(outdir, batch_job_hash + '_accessibility.dat'))

        if "interaction" in opts.analyze_list:
            # Get interaction data of all models at 200 nm cut-off
            logging.info("\tGetting interaction data...")
            models.interactions(
                cutoff=dcutoff, steps=(1,3,5),
                savefig =path.join(outdir, batch_job_hash + '_interactions.' + opts.fig_format),
                savedata=path.join(outdir, batch_job_hash + '_interactions.dat'),
                error=True)


@retry(lite.OperationalError, tries=20, delay=2)
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
            JOBid int,
            Type text,
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
            Correlation text)""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            # In case optimization or modeling  is split in different computers
            cur.execute("""
    insert into JOBs
     (Id  , Parameters, Launch_time, Finish_time,    Type, Parameters_md5)
    values
     (NULL,       '%s',        '%s',        '%s',    '%s',           '%s')
     """ % ((parameters, time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
             time.strftime("%d/%m/%Y %H:%M:%S", finish_time),
             (('PRE_' if opts.job_list else '') +
              ('OPTIM' if opts.optimize else 'MODEL')), batch_job_hash)))
        except lite.IntegrityError:
            pass
        if not opts.job_list:
            ##### STORE OPTIMIZATION RESULT
            jobid = get_jobid(cur)
            add_path(cur, outdir, 'EXT_OPTIM_FOLDER' if opts.optimize else 'DIR',
                     jobid, opts.workdir)
            pathid = get_path_id(cur, outdir, opts.workdir)
            # models = compile_models(opts, outdir, exp=exp, ngood=opts.nkeep)
            ### STORE GENERAL OPTIMIZATION INFO
            try:
                cur.execute("""
                insert into MODELED_REGIONs
                (Id  , JOBid , Type, PATHid, PARAM_md5, RESO, BEG, END)
                values
                (NULL,   %d,     "%s",   %d,      "%s",   %d,  %d,  %d)
                """ % (jobid, 'OPTIM' if opts.optimize else 'MODEL', pathid,
                    batch_job_hash, opts.reso, opts.beg, opts.end))
            except lite.IntegrityError:
                pass
            ### STORE EACH OPTIMIZATION
            cur.execute("SELECT Id from MODELED_REGIONs where PARAM_md5='%s'" % (
                batch_job_hash))
            optimid = cur.fetchall()[0][0]
            for m, u, l, d, s in results:
                optpar_md5 = md5(('%s%s%s%s%s' %
                                 (m, u, l, d, s)).encode('utf-8')).hexdigest()[:12]
                cur.execute(("SELECT Id from MODELs where "
                             "OPTPAR_md5='%s' and REGIONid='%s'") % (
                                 optpar_md5, optimid))
                if not cur.fetchall():
                    cur.execute("""
                    insert into MODELs
                    (Id  , REGIONid, JOBid, OPTPAR_md5, MaxDist, UpFreq, LowFreq, Cutoff, Scale, Nmodels, Kept, Correlation)
                    values
                    (NULL,             %d,    %d,      '%s',      %s,     %s,      %s,     %s,    %s,      %d,   %d,   '%f')
                    """ % ((optimid, jobid, optpar_md5, m, u, l, d, s,
                            results[(m, u, l, d, s)]['nmodels'],
                            results[(m, u, l, d, s)]['kept'],
                            results[(m, u, l, d, s)]['corr'])))
                    muls = tuple(map(my_round, (m, u, l, s)))
                    dirname = path.join(outdir, 'cfg_%s_%s_%s_%s' % muls)
                    add_path(cur, dirname, 'MODELS', jobid, opts.workdir)
                else:
                    cur.execute(("update MODELs "
                                 "set Nmodels = %d, Kept = %d, Correlation = %f "
                                 "where "
                                 "OPTPAR_md5='%s' and REGIONid='%s'") % (
                                     results[(m, u, l, d, s)]['nmodels'],
                                     results[(m, u, l, d, s)]['kept'],
                                     results[(m, u, l, d, s)]['corr'],
                                     optpar_md5, optimid))

    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
    # release lock
    try:
        remove(path.join(opts.workdir, '__lock_db'))
    except OSError:
        pass

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class = lambda prog: HelpFormatter(prog, width=95,
                                                        max_help_position=27)

    glopts = parser.add_argument_group('General options')
    descro = parser.add_argument_group('Descriptive, optional arguments')
    reopts = parser.add_argument_group('Modeling preparation')
    opopts = parser.add_argument_group('Parameter optimization')
    analyz = parser.add_argument_group('Analysis')
    ruopts = parser.add_argument_group('Running jobs')

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
    glopts.add_argument('-j', '--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')
    glopts.add_argument('--optimization_id', dest='optimization_id', metavar="INT",
                        type=float, default=None,
                        help="[%(default)s] ID of a pre-run optimization batch job")
    glopts.add_argument('--fig_format', dest='fig_format', metavar="STR",
                        default="pdf",
                        help='''file format and extension for figures and plots
                        (can be any supported by matplotlib, png, eps...)''')
    glopts.add_argument('--noX', action='store_true', help='no display server (X screen)')
    glopts.add_argument('--corr', dest='corr', metavar="STR",
                        default="spearman",
                        help='''correlation method to compare contact maps and original matrix
                        (options are speraman, pearson, kendall, logpearson, chi2, scc )''')

    #########################################
    # DESCRIPTION
    descro.add_argument('--species', dest='species', metavar="STRING",
                        default='UNKNOWN',
                        help='species name, with no spaces, i.e.: homo_sapiens')
    descro.add_argument('--assembly', dest='assembly', metavar="STRING",
                        default=None,
                        help='''NCBI ID of the original assembly
                        (i.e.: NCBI36 for human)''')
    descro.add_argument('--cell', dest='cell', metavar="STRING",
                        help='cell type name')
    descro.add_argument('--exp_type', dest='exp_type', metavar="STRING",
                        help='experiment type name (i.e.: Hi-C)')
    descro.add_argument('--project', dest='project', metavar="STRING",
                        default=None,
                        help='''project name''')

    reopts.add_argument('--crm', dest='crm', metavar="NAME",default=None,
                        help='chromosome name')
    reopts.add_argument('--beg', dest='beg', metavar="INT", type=float,
                        default=None,
                        help='genomic coordinate from which to start modeling')
    reopts.add_argument('--end', dest='end', metavar="INT", type=float,
                        default=None,
                        help='genomic coordinate where to end modeling')
    reopts.add_argument('--matrix_beg', dest='matrix_beg', metavar="INT", type=int,
                        default=None,
                        help='genomic coordinate of the first row/column ' +
                        'of the input matrix. This has to be specified if ' +
                        'the input matrix is not the TADbit tools generated abc format')
    reopts.add_argument('-r', '--reso', dest='reso', metavar="INT", type=int,
                        help='resolution of the Hi-C experiment')
    reopts.add_argument('--perc_zero', dest='perc_zero', metavar="FLOAT",
                        type=float, default=90.0)

    opopts.add_argument('--optimize', dest='optimize',
                        default=False, action="store_true",
                        help='''optimization run, store less info about models''')
    opopts.add_argument('--model', dest='model',
                        default=False, action="store_true",
                        help='''modelling run''')
    opopts.add_argument('--force', dest='force',
                        default=False, action="store_true",
                        help='''use input parameters, and skip any precalculated
                        optimization''')
    opopts.add_argument('--maxdist', action='store', metavar="LIST",
                        default=[400], dest='maxdist', nargs='+',
                        help='range of numbers for maxdist' +
                        ', i.e. 400:1000:100 -- or just a number -- or a ' +
                        'list of numbers')
    opopts.add_argument('--upfreq', dest='upfreq', metavar="LIST",
                        default=[0], nargs='+',
                        help='range of numbers for upfreq' +
                        ', i.e. 0:1.2:0.3 -- or just a number -- or a ' +
                        'list of numbers')
    opopts.add_argument('--lowfreq', dest='lowfreq', metavar="LIST",
                        default=[0], nargs='+',
                        help='range of numbers for lowfreq' +
                        ', i.e. -1.2:0:0.3 -- or just a number -- or a ' +
                        'list of numbers')
    opopts.add_argument('--scale', dest='scale', metavar="LIST",
                        default=[0.01], nargs='+',
                        help='%(default)s range of numbers to be test as ' +
                        'optimal scale value, i.e. 0.005:0.01:0.001 -- Can ' +
                        'also pass only one number -- or a ' +
                        'list of numbers')
    opopts.add_argument('--dcutoff', dest='dcutoff', metavar="LIST",
                        default=[2], nargs='+',
                        help='%(default)s range of numbers to be test as ' +
                        'optimal distance cutoff parameter (distance, in ' +
                        'number of beads, from which to consider 2 beads as ' +
                        'being close), i.e. 1:1.5:0.5 -- Can also pass only one' +
                        ' number -- or a list of numbers')

    opopts.add_argument('--analyze', dest='analyze',
                        default=False, action="store_true",
                        help='''analyze models.''')

    ruopts.add_argument("-C", "--cpu", dest="cpus", type=int,
                        default=cpu_count(), help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')
    ruopts.add_argument('--job_list', dest='job_list', action='store_true',
                        default=False,
                        help=('generate a list of commands stored in a '
                              'file named joblist_HASH.q (where HASH is '
                              'replaced by a string specific to the parameters '
                              'used). note that dcutoff will never be split as '
                              'it does not require to re-run models.'))
    ruopts.add_argument('--nmodels_per_job', dest='nmodels_per_job', metavar="INT",
                        default=1, type=int,
                        help=('Number of models per distributed job.'))
    ruopts.add_argument('--cpus_per_job', dest='cpus_per_job',
                        metavar="INT", default=1, type=int,
                        help=('Number of cpu nodes per distributed job.'))
    ruopts.add_argument('--concurrent_jobs', dest='concurrent_jobs',
                        metavar="INT", default=cpu_count(), type=int,
                        help=('Number of concurrent jobs in distributed mode.'))
    ruopts.add_argument('--timeout_job', dest='timeout_job', metavar="INT", type=int,
                        default=5000,
                        help=('Time to wait for a concurrent jobs to finish before '
                              'canceling it in distributed mode.'))
    ruopts.add_argument('--script_cmd', dest='script_cmd', metavar="STR", type=str,
                        default='python',
                        help=('Command to call the jobs '
                              'in distributed mode.'))
    ruopts.add_argument('--script_args', dest='script_args', metavar="STR", type=str,
                        default='-u',
                        help=('Argumnets to script_cmd to call the jobs '
                              'in distributed mode.'))
    # ruopts.add_argument('--job_list', dest='job_list', metavar='LIST/nothing', nargs='*',
    #                     choices=['maxdist', 'upfreq', 'lowfreq', 'scale', 'dcutoff'],
    #                     default=None,
    #                     help=('[False] generate a list of commands stored in a '
    #                           'file named joblist_HASH.q (where HASH is '
    #                           'replaced by a string specific to the parameters '
    #                           'used). With no extra argument any combination of'
    #                           ' optimized parameters will generate a new job. '
    #                           'Otherwise the parameters to be splitted can be '
    #                           'specified (e.g. "--job_list maxdist lowfreq '
    #                           'upfreq scale" will run in a single job all '
    #                           'combinations of dcutoff -- which makes sense as '
    #                           'dcutoff does not require to re-run models). '
    #                           'Choices are: %(choices)s'))
    ruopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    #########################################
    # OUTPUT
    analyz.add_argument('--analyze_list', dest='analyze_list', nargs='+',
                        choices=list(range(len(actions))), type=int,
                        default=list(range(1, len(actions))), metavar='INT',
                        help=('''[%s] list of numbers representing the
                        analysis to be done. Choose between:
                        %s''' % (' '.join([str(i) for i in range(
                                  2, len(actions))]),
                                 '\n'.join(['%s) %s' % (k, actions[k])
                                            for k in actions]))))
    analyz.add_argument('--not_write_cmm', dest='not_write_cmm',
                        default=False, action='store_true',
                        help='''[%(default)s] do not generate cmm files for each
                        model (Chimera input)''')
    analyz.add_argument('--not_write_xyz', dest='not_write_xyz',
                        default=False, action='store_true',
                        help='''[%(default)s] do not generate xyz files for each
                        model (3D coordinates)''')
    analyz.add_argument('--not_write_json', dest='not_write_json',
                        default=False, action='store_true',
                        help='''[%(default)s] do not generate json file.''')

def check_options(opts):
    # check resume
    if not path.exists(opts.workdir):
        warn('ERROR: workdir not found, creating it')
        mkdir(opts.workdir)
        # write version log
        vlog_path = path.join(opts.workdir, 'TADbit_and_dependencies_versions.log')
        dependencies = get_dependencies_version()
        if not path.exists(vlog_path) or open(vlog_path).readlines() != dependencies:
            logging.info('Writing versions of TADbit and dependencies')
            vlog = open(vlog_path, 'w')
            vlog.write(dependencies)
            vlog.close()
    # do the division to bins
    try:
        opts.ori_beg = opts.beg
        opts.ori_end = opts.end
        if opts.matrix_beg:
            opts.matrix_beg = int(float(opts.matrix_beg) / opts.reso)
        opts.beg = int(float(opts.beg) / opts.reso)
        opts.end = int(float(opts.end) / opts.reso)
        if opts.end - opts.beg <= 2:
            raise Exception('"beg" and "end" parameter should be given in ' +
                            'genomic coordinates, not bin')
    except TypeError:
        pass

    # turn options into lists
    def _load_range(range_str, num=float, decs=2):
        try:
            beg, end, step = list(map(num, range_str[0].split(':')))
            return tuple([round(x,decs) for x in arange(beg, end + step / 2, step)])
        except (AttributeError, ValueError):
            return tuple([round(num(v),decs) for v in range_str])

    opts.scale   = _load_range(opts.scale, decs=4)
    opts.maxdist = _load_range(opts.maxdist, num=int)
    opts.upfreq  = _load_range(opts.upfreq)
    opts.lowfreq = _load_range(opts.lowfreq)
    opts.dcutoff = _load_range(opts.dcutoff)

    if opts.matrix:
        opts.matrix  = path.abspath(opts.matrix)
    opts.workdir = path.abspath(opts.workdir)

        # rename analysis actions
    for i, j in enumerate(opts.analyze_list):
        opts.analyze_list[i] = actions[int(j)]

    mkdir(opts.workdir)
    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)


def load_optpar_fromdb(opts):
    if 'tmpdb' in opts and opts.tmpdb:
        dbfile = opts.tmpdb
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        if not opts.optimization_id:
            cur.execute("SELECT Id from MODELED_REGIONs where RESO=%d and BEG=%d and END=%d and Type='OPTIM'" % (
            opts.reso, opts.beg, opts.end))
        else:
            cur.execute("SELECT Id from MODELED_REGIONs where JOBid=%d and RESO=%d and BEG=%d and END=%d and Type='OPTIM'" % (
                opts.optimization_id,opts.reso, opts.beg, opts.end))
        try:
            optimid = cur.fetchall()
        except IndexError:
            raise IndexError("ERROR: no optimization found. Run 'tadbit model' "
                                 "with --optimize or use --force to indicate optimal "
                                 "parameters ")
        if len(optimid) > 1:
            raise IndexError("ERROR: more than 1 optimization in folder "
                             "choose with 'tadbit describe' and "
                             "--optimization_id")
        optimid = optimid[0]
        cur.execute("""
                select Id , JOBid, OPTPAR_md5, MaxDist, UpFreq, LowFreq, 
                    Cutoff, Scale, Nmodels, Kept, Correlation from MODELs
                    where REGIONid=%d ORDER BY Correlation DESC
                """ % (optimid))
        optpar = cur.fetchall()[0]
        logging.info(('Loaded UpFreq:%6s LowFreq:%7s MaxDist:%7s scale:%6s cutoff:%7s Correlation:%.4s' %
                (optpar[4], optpar[5], optpar[3], optpar[7], optpar[6], optpar[10])))
        optpar = {'maxdist': optpar[3],
                  'upfreq' : optpar[4],
                  'lowfreq': optpar[5],
                  'scale'  : optpar[7],
                  'kforce' : 5}
    return optpar

def load_models_path_fromdb(opts):
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
            where Type = 'MODEL'
            """)
            jobids = cur.fetchall()
            if len(jobids) > 1:
                raise Exception('ERROR: more than one possible input found, use'
                                '"tadbit describe" and select corresponding '
                                'jobid with --jobid')
            parse_jobid = jobids[0][0]
        else:
            parse_jobid = opts.jobid
        # fetch path to normalized matrix
        cur.execute("""
        select distinct Path, PATHs.Id from PATHs
        where paths.jobid = %s and paths.Type = 'DIR'
        """ % parse_jobid)
        jobids = cur.fetchall()
        if len(jobids) < 1:
            raise Exception('ERROR: no folders found in job.')
        fold, fold_id = jobids[0]
        return (fold, fold_id)

def load_matrix_path_fromdb(opts):
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
            where Type = 'Bin'
            """)
            jobids = cur.fetchall()
            if len(jobids) > 1:
                raise Exception('ERROR: more than one possible input found, use'
                                '"tadbit describe" and select corresponding '
                                'jobid with --jobid')
            parse_jobid = jobids[0][0]
        else:
            parse_jobid = opts.jobid
        # fetch path to normalized matrix
        cur.execute("""
        select distinct Path, PATHs.Id from PATHs
        where paths.jobid = %s and paths.Type = 'NRM_MATRIX'
        """ % parse_jobid)
        jobids = cur.fetchall()
        if len(jobids) < 1:
            raise Exception('ERROR: no normalized matrix found in job.')
        mat, mat_id = jobids[0]
        return (mat, mat_id)

def load_hic_data(opts):
    """
    Load Hi-C data
    """
    # Start reading the data
    crm = Chromosome(opts.crm, species=(
        opts.species.split('_')[0].capitalize() + opts.species.split('_')[1]
                          if '_' in opts.species else opts.species),
                          assembly=opts.assembly) # Create chromosome object
    logging.info( '     o Loading Hi-C matrix')
    opts.offset = 0
    try:
        hic = read_matrix(opts.matrix, hic=False,
                          resolution=opts.reso)
        hic_bads = hic.bads
        if opts.crm: # just to avoid loading the full chromosome to model a small region
            if opts.crm not in hic.chromosomes:
                raise Exception('ERROR: chromosome %s not in input matrix(%s).' % (opts.crm,
                                                    ','.join([h for h in hic.chromosomes])))
            hic_bads = {k: v for k, v in list(hic.bads.items()) if k >= opts.beg and k < opts.end}
            hic = hic.get_matrix(focus=(opts.beg+1,opts.end))
            opts.offset = opts.beg
        else:
            if len(hic.chromosomes) == 1: # we assume full chromosome
                logging.info('WARNING: assuming full chromosome %s.' % next(iter(hic.chromosomes)))
                opts.crm = next(iter(hic.chromosomes))
                opts.beg = 1
                opts.end = hic.chromosomes[opts.crm]
        opts.matrix_beg = 0
        crm.add_experiment('test',
            cell_type=opts.cell,
            project=opts.project, # user descriptions
            norm_data=[hic], resolution=opts.reso)
        crm.experiments[0].norm[0].bads = hic_bads
        crm.experiments[0]._zeros = hic_bads
        crm.experiments[0]._filtered_cols = True
    except Exception as e:
        #logging.info( str(e))
        warn('WARNING: failed to load data as TADbit standardized matrix\n')
        if opts.matrix_beg is None:
            raise Exception('"matrix_beg" parameter should be given if ' +
                            'input matrix is not in TADbit abc format')
        if opts.beg is None:
            opts.beg = 1
        if opts.end is None:
            opts.end = sum(1 for _ in open(opts.matrix))
        opts.beg -= opts.matrix_beg
        opts.end -= opts.matrix_beg
        crm.add_experiment('test',
            cell_type=opts.cell,
            project=opts.project, # user descriptions
            resolution=opts.reso,
            norm_data=opts.matrix)

    if opts.beg is not None:
        if opts.beg - opts.offset + 1 > crm.experiments[-1].size:
            raise Exception('ERROR: beg parameter is larger than chromosome size.')
        if opts.end - opts.offset > crm.experiments[-1].size:
            logging.info ('WARNING: end parameter is larger than chromosome ' +
                   'size. Setting end to %s.\n' % (crm.experiments[-1].size *
                                                   opts.reso))
            opts.end = crm.experiments[-1].size + opts.offset
    return crm
