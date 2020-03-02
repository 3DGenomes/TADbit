"""
16 Mar 2019


"""
from string import uppercase as uc, lowercase as lc
from os.path import exists
from random import randint, seed, random, sample, shuffle
from cPickle import load, dump
#from pebble import ProcessPool
#from concurrent.futures import TimeoutError
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
from math import atan2
from itertools import combinations, product
from shutil import copyfile

import sys
import copy
import os
import shutil
import multiprocessing

from numpy import sin, cos, arccos, sqrt, fabs, pi
import numpy as np


try:
    from pytadbit.modelling.imp_modelling import generate_3d_models
except ImportError:
    pass

try:
    from lammps import lammps
except ImportError:
    pass

from pytadbit.modelling import LAMMPS_CONFIG as CONFIG
from pytadbit.modelling.lammpsmodel import LAMMPSmodel
from pytadbit.modelling.structuralmodels import StructuralModels
from pytadbit.modelling.restraints import HiCBasedRestraints

class InitalConformationError(Exception):
    """
    Exception to handle failed initial conformation.
    """
    pass

def abortable_worker(func, *args, **kwargs):
    timeout = kwargs.get('timeout', None)
    failedSeedLog = kwargs.get('failedSeedLog', None)
    p = ThreadPool(1)
    res = p.apply_async(func, args=args)
    try:
        out = res.get(timeout)  # Wait timeout seconds for func to complete.
        return out
    except multiprocessing.TimeoutError:
        print "Model took more than %s seconds to complete ... canceling" % str(timeout)
        p.terminate()
        raise
    except:
        print "Unknown error with process"
        if failedSeedLog != None:
            failedSeedLog, k = failedSeedLog
            with open(failedSeedLog, 'a') as f:
                f.write('%s\t%s\n' %(k, 'Failed'))
        p.terminate()
        raise


def generate_lammps_models(zscores, resolution, nloci, start=1, n_models=5000,
                       n_keep=1000, close_bins=1, n_cpus=1,
                       verbose=0, outfile=None, config=None,
                       values=None, experiment=None, coords=None, zeros=None,
                       first=None, container=None,tmp_folder=None,timeout_job=10800,
                       initial_conformation=None, connectivity="FENE",
                       timesteps_per_k=10000,keep_restart_out_dir=None,
                       kfactor=1, adaptation_step=False, cleanup=False,
                       hide_log=True, remove_rstrn=[], initial_seed=0,
                       restart_path=False, store_n_steps=10,
                       useColvars=False):
    """
    This function generates three-dimensional models starting from Hi-C data.
    The final analysis will be performed on the n_keep top models.

    :param zscores: the dictionary of the Z-score values calculated from the
       Hi-C pairwise interactions
    :param resolution:  number of nucleotides per Hi-C bin. This will be the
       number of nucleotides in each model's particle
    :param nloci: number of particles to model (may not all be present in
       zscores)
    :param None experiment: experiment from which to do the modelling (used only
       for descriptive purpose)
    :param None coords: a dictionary like:
       ::

         {'crm'  : '19',
          'start': 14637,
          'end'  : 15689}

    :param 5000 n_models: number of models to generate
    :param 1000 n_keep: number of models used in the final analysis (usually
       the top 20% of the generated models). The models are ranked according to
       their objective function value (the lower the better)
    :param 1 close_bins: number of particles away (i.e. the bin number
       difference) a particle pair must be in order to be considered as
       neighbors (e.g. 1 means consecutive particles)
    :param n_cpus: number of CPUs to use
    :param False verbose: if set to True, information about the distance, force
       and Z-score between particles will be printed. If verbose is 0.5 than
       constraints will be printed only for the first model calculated.
    :param None values: the normalized Hi-C data in a list of lists (equivalent
       to a square matrix)
    :param None config: a dictionary containing the standard
       parameters used to generate the models. The dictionary should contain
       the keys kforce, lowrdist, maxdist, upfreq and lowfreq. Examples can be
       seen by doing:

       ::

         from pytadbit.modelling.CONFIG import CONFIG

         where CONFIG is a dictionary of dictionaries to be passed to this function:

       ::

         CONFIG = {
          'dmel_01': {
              # Paramaters for the Hi-C dataset from:
              'reference' : 'victor corces dataset 2013',

              # Force applied to the restraints inferred to neighbor particles
              'kforce'    : 5,

              # Space occupied by a nucleotide (nm)
              'scale'     : 0.005

              # Strength of the bending interaction
              'kbending'     : 0.0, # OPTIMIZATION:

              # Maximum experimental contact distance
              'maxdist'   : 600, # OPTIMIZATION: 500-1200

              # Minimum thresholds used to decide which experimental values have to be
              # included in the computation of restraints. Z-score values bigger than upfreq
              # and less that lowfreq will be include, whereas all the others will be rejected
              'lowfreq'   : -0.7 # OPTIMIZATION: min/max Z-score

              # Maximum threshold used to decide which experimental values have to be
              # included in the computation of restraints. Z-score values greater than upfreq
              # and less than lowfreq will be included, while all the others will be rejected
              'upfreq'    : 0.3 # OPTIMIZATION: min/max Z-score

              }
          }
    :param None first: particle number at which model should start (0 should be
       used inside TADbit)
    :param None container: restrains particle to be within a given object. Can
       only be a 'cylinder', which is, in fact a cylinder of a given height to
       which are added hemispherical ends. This cylinder is defined by a radius,
       its height (with a height of 0 the cylinder becomes a sphere) and the
       force applied to the restraint. E.g. for modeling E. coli genome (2
       micrometers length and 0.5 micrometer of width), these values could be
       used: ['cylinder', 250, 1500, 50], and for a typical mammalian nuclei
       (6 micrometers diameter): ['cylinder', 3000, 0, 50]
    :param None tmp_folder: path to a temporary file created during
        the clustering computation. Default will be created in /tmp/ folder
    :param 10800 timeout_job: maximum seconds a job can run in the multiprocessing
        of lammps before is killed
    :param initial_conformation: lammps input data file with the particles initial conformation.
    :param True hide_log: do not generate lammps log information
    :param FENE connectivity: use FENE for a fene bond or harmonic for harmonic
        potential for neighbours
    :param None keep_restart_out_dir: path to write files to restore LAMMPs
                session (binary)
    :param True cleanup: delete lammps folder after completion
    :param [] remove_rstrn: list of particles which must not have restrains
    :param 0 initial_seed: Initial random seed for modelling.
    :param False restart_path: path to files to restore LAMMPs session (binary)
    :param 10 store_n_steps: Integer with number of steps to be saved if 
        restart_file != False
    :param False useColvars: True if you want the restrains to be loaded by colvars

    :returns: a StructuralModels object
    """

    if not tmp_folder:
        tmp_folder = '/tmp/tadbit_tmp_%s/' % (
            ''.join([(uc + lc)[int(random() * 52)] for _ in xrange(4)]))
    else:
        if tmp_folder[-1] != '/':
            tmp_folder += '/'
        randk = ''.join([(uc + lc)[int(random() * 52)] for _ in xrange(4)])
        tmp_folder = '%s%s/' %(tmp_folder, randk)
    while os.path.exists(tmp_folder):
        randk = ''.join([(uc + lc)[int(random() * 52)] for _ in xrange(4)])
        tmp_folder = '%s%s/' %(tmp_folder[:-1], randk)
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)

    # Setup CONFIG
    if isinstance(config, dict):
        CONFIG.HiC.update(config)
    elif config:
        raise Exception('ERROR: "config" must be a dictionary')

    global RADIUS

    #RADIUS = float(resolution * CONFIG['scale']) / 2
    RADIUS = 0.5
    CONFIG.HiC['resolution'] = resolution
    CONFIG.HiC['maxdist'] = CONFIG.HiC['maxdist'] / (float(resolution * CONFIG.HiC['scale']))

    global LOCI
    # if z-scores are generated outside TADbit they may not start at zero
    if first is None:
        first = min([int(j) for i in zscores[0] for j in zscores[0][i]] +
                    [int(i) for i in zscores[0]])
    LOCI  = range(first, nloci + first)
    
    # random inital number
    global START
    START = start
    # verbose
    global VERBOSE
    VERBOSE = verbose
    #VERBOSE = 3

    HiCRestraints = [HiCBasedRestraints(nloci,RADIUS,CONFIG.HiC,resolution, zs,
                 chromosomes=coords, close_bins=close_bins,first=first,
                 remove_rstrn=remove_rstrn) for zs in zscores]

    run_time = 1000
    
    colvars = 'colvars.dat'

    steering_pairs = None
    time_dependent_steering_pairs = None
    if len(HiCRestraints) > 1:
        time_dependent_steering_pairs = {
            'colvar_input'              : HiCRestraints,
            'colvar_output'             : colvars,
            'chrlength'                 : nloci,
            'binsize'                   : resolution,
            'timesteps_per_k_change'    : [float(timesteps_per_k)]*6,
            'k_factor'                  : kfactor,
            'perc_enfor_contacts'       : 100.,
            'colvar_dump_freq'          : int(timesteps_per_k/100),
            'adaptation_step'           : adaptation_step,
        }
        if not initial_conformation:
            initial_conformation = 'tadbit'
    else:
        steering_pairs = {
            'colvar_input': HiCRestraints[0],
            'colvar_output': colvars,
            'binsize': resolution,
            'timesteps_per_k'           : timesteps_per_k,
            'k_factor'                  : kfactor,
            'colvar_dump_freq'          : int(timesteps_per_k/100),
            'timesteps_relaxation'      : int(timesteps_per_k*10)
        }
        if not initial_conformation:
            initial_conformation = 'random'

    if not container:
        container = ['cube',1000.0] # http://lammps.sandia.gov/threads/msg48683.html

    ini_model = ini_sm_model = None
    if initial_conformation != 'random':
        if isinstance(initial_conformation, dict):
            sm = [initial_conformation]
            sm[0]['x'] = sm[0]['x'][0:nloci]
            sm[0]['y'] = sm[0]['y'][0:nloci]
            sm[0]['z'] = sm[0]['z'][0:nloci]
        elif initial_conformation == 'tadbit':
            sm = generate_3d_models(zscores[0], resolution, nloci,
                  values=values[0], n_models=n_models, n_keep=n_keep, n_cpus=n_cpus,
                  verbose=verbose, first=first, close_bins=close_bins,
                  config=config, container=container,
                  coords=coords, zeros=zeros[0])
            print "Succesfully generated tadbit initial conformation \n"
        sm_diameter = float(resolution * CONFIG.HiC['scale'])
        for single_m in sm:
            for i in xrange(len(single_m['x'])):
                single_m['x'][i] /= sm_diameter
                single_m['y'][i] /= sm_diameter
                single_m['z'][i] /= sm_diameter
            cm0 = single_m.center_of_mass()
            for i in xrange(len(single_m['x'])):
                single_m['x'][i] -= cm0['x']
                single_m['y'][i] -= cm0['y']
                single_m['z'][i] -= cm0['z']
        ini_sm_model = [[single_sm.copy()] for single_sm in sm]
        ini_model = [single_sm.copy() for single_sm in sm]

    models = lammps_simulate(lammps_folder=tmp_folder, run_time=run_time,
                             initial_conformation=ini_sm_model,
                             connectivity=connectivity,
                             steering_pairs=steering_pairs,
                             time_dependent_steering_pairs=time_dependent_steering_pairs,
                             initial_seed=initial_seed,
                             n_models=n_keep, n_keep=n_keep, n_cpus=n_cpus,
                             keep_restart_out_dir=keep_restart_out_dir,
                             confining_environment=container, timeout_job=timeout_job,
                             cleanup=cleanup, to_dump=int(timesteps_per_k/100.),
                             hide_log=hide_log, restart_path=restart_path,
                             store_n_steps=store_n_steps,
                             useColvars=useColvars)

    try:
        xpr = experiment
        crm = xpr.crm
        description = {'identifier'        : xpr.identifier,
                       'chromosome'        : coords['crm'] if isinstance(coords,dict) \
                                             else [c['crm'] for c in coords],
                       'start'             : xpr.resolution * coords['start'] \
                                             if isinstance(coords,dict) \
                                             else [xpr.resolution*c['start'] for c in coords],
                       'end'               : xpr.resolution * coords['end'] \
                                             if isinstance(coords,dict) \
                                             else [xpr.resolution*c['end'] for c in coords],
                       'species'           : crm.species,
                       'restriction enzyme': xpr.enzyme,
                       'cell type'         : xpr.cell_type,
                       'experiment type'   : xpr.exp_type,
                       'resolution'        : xpr.resolution,
                       'assembly'          : crm.assembly}
        for desc in xpr.description:
            description[desc] = xpr.description[desc]
        for desc in crm.description:
            description[desc] = xpr.description[desc]
        for i, m in enumerate(models.values()):
            m['index'] = i
            m['description'] = description
    except AttributeError: # case we are doing optimization
        description = None
        for i, m in enumerate(models.values()):
            m['index'] = i
    if outfile:
        if exists(outfile):
            old_models, _ = load(open(outfile))
        else:
            old_models, _ = {}, {}
        models.update(old_models)
        out = open(outfile, 'w')
        dump((models), out)
        out.close()
    else:
        stages = {}
        timepoints = None
        allzeros = tuple([all([zero_stg[x] for zero_stg in zeros]) for x in xrange(len(zeros[0]))])
        if len(HiCRestraints)>1:
            #for timepoint in xrange(len(zeros)-1):
            #    allzeros = tuple([sum(x) for x in zip(allzeros, zeros[timepoint])])
            timepoints = time_dependent_steering_pairs['colvar_dump_freq']
            nbr_produced_models = len(models)/(timepoints*(len(HiCRestraints)-1))
            stages[0] = [i for i in xrange(nbr_produced_models)]

            for sm_id, single_m in enumerate(ini_model):
                for i in xrange(len(single_m['x'])):
                    single_m['x'][i] *= sm_diameter
                    single_m['y'][i] *= sm_diameter
                    single_m['z'][i] *= sm_diameter
                
                lammps_model = LAMMPSmodel({ 'x'          : single_m['x'],
                                              'y'          : single_m['y'],
                                              'z'          : single_m['z'],
                                              'cluster'    : 'Singleton',
                                              'objfun'     : single_m['objfun'],
                                              'log_objfun' : single_m['log_objfun'],
                                              'radius'     : float(CONFIG.HiC['resolution'] * \
                                                                   CONFIG.HiC['scale'])/2,
                                              'rand_init'  : str(sm_id+1+initial_seed)})

                models[sm_id] = lammps_model
            for timepoint in xrange((len(HiCRestraints)-1)*timepoints):
                stages[timepoint+1] = [(t+nbr_produced_models+timepoint*nbr_produced_models)
                                       for t in xrange(nbr_produced_models)]

        return StructuralModels(
            len(LOCI), models, {}, resolution, original_data=values \
                if len(HiCRestraints)>1 else values[0],
            zscores=zscores, config=CONFIG.HiC, experiment=experiment, zeros=allzeros,
            restraints=HiCRestraints[0]._get_restraints(),
            description=description, stages=stages, models_per_step=timepoints)
# Initialize the lammps simulation with standard polymer physics based
# interactions: chain connectivity (FENE) ; excluded volume (WLC) ; and
# bending rigidity (KP)
def init_lammps_run(lmp, initial_conformation,
                    neighbor=CONFIG.neighbor,
                    hide_log=True,
                    connectivity="FENE",
                    restart_file=False):

    """
    Initialise the parameters for the computation in lammps job

    :param lmp: lammps instance object.
    :param initial_conformation: lammps input data file with the particles initial conformation.
    :param CONFIG.neighbor neighbor: see LAMMPS_CONFIG.py.
    :param True hide_log: do not generate lammps log information
    :param FENE connectivity: use FENE for a fene bond or harmonic for harmonic 
        potential for neighbours
    :param False restart_file: path to file to restore LAMMPs session (binary)

    """

    if hide_log:
        lmp.command("log none")
    #os.remove("log.lammps")

    #######################################################
    # Box and units  (use LJ units and period boundaries) #
    #######################################################
    lmp.command("units %s" % CONFIG.units)
    lmp.command("atom_style %s" % CONFIG.atom_style) #with stiffness
    lmp.command("boundary %s" % CONFIG.boundary)
    """
    try:
        lmp.command("communicate multi")
    except:
        pass
    """

    ##########################
    # READ "start" data file #
    ##########################
    if restart_file == False :
        lmp.command("read_data %s" % initial_conformation)
    else:
        restart_time = int(restart_file.split('/')[-1].split('_')[4][:-8])
        print 'Previous unfinished LAMMPS steps found'
        print 'Loaded %s file' %restart_file
        lmp.command("read_restart %s" % restart_file)
        lmp.command("reset_timestep %i" % restart_time)
        
    lmp.command("mass %s" % CONFIG.mass)

    ##################################################################
    # Pair interactions require lists of neighbours to be calculated #
    ##################################################################
    lmp.command("neighbor %s" % neighbor)
    lmp.command("neigh_modify %s" % CONFIG.neigh_modify)
    
    ##############################################################
    # Sample thermodynamic info  (temperature, energy, pressure) #
    ##############################################################
    lmp.command("thermo %i" % CONFIG.thermo)
    
    ###############################
    # Stiffness term              #
    # E = K * (1+cos(theta)), K>0 #
    ###############################
    lmp.command("angle_style %s" % CONFIG.angle_style) # Write function for kinks 
    lmp.command("angle_coeff * %f" % CONFIG.persistence_length)
    
    ###################################################################
    # Pair interaction between non-bonded atoms                       #
    #                                                                 #
    #  Lennard-Jones 12-6 potential with cutoff:                      #
    #  potential E=4epsilon[ (sigma/r)^12 - (sigma/r)^6]  for r<r_cut #
    #  r_cut =1.12246 = 2^(1/6) is the minimum of the potential       #
    ###################################################################
    lmp.command("pair_style hybrid/overlay lj/cut %f morse 0.0" % CONFIG.PurelyRepulsiveLJcutoff)
    #lmp.command("pair_style lj/cut %f" % CONFIG.PurelyRepulsiveLJcutoff)
    
    ################################################################
    #  pair_modify shift yes adds a constant to the potential such #
    #  that E(r_cut)=0. Forces remains unchanged.                  #
    ################################################################
    lmp.command("pair_modify shift yes")
    
    ######################################
    #  pair_coeff for lj/cut, specify 4: #
    #    * atom type interacting with    #
    #    * atom type                     #
    #    * epsilon (energy units)        #
    #    * sigma (distance units)        #
    ######################################
    lmp.command("pair_coeff * * lj/cut %f %f %f" % (CONFIG.PurelyRepulsiveLJepsilon, 
                                                    CONFIG.PurelyRepulsiveLJsigma, 
                                                    CONFIG.PurelyRepulsiveLJcutoff/2.0))
 
    lmp.command("pair_coeff * * morse  %f %f %f" % (0.0, 0.0, 0.0))
    
    if connectivity == "FENE":
        #########################################################
        # Pair interaction between bonded atoms                 #
        #                                                       #
        # Fene potential + Lennard Jones 12-6:                  #
        #  E= - 0.5 K R0^2 ln[ 1- (r/R0)^2]                     #
        #     + 4epsilon[ (sigma/r)^12 - (sigma/r)^6] + epsilon #
        #########################################################
        lmp.command("bond_style fene")
    
        ########################################
        # For style fene, specify:             #
        #   * bond type                        #
        #   * K (energy/distance^2)            #
        #   * R0 (distance)                    #
        #   * epsilon (energy)  (LJ component) #
        #   * sigma (distance)  (LJ component) #
        ########################################
        lmp.command("bond_coeff * %f %f %f %f" % (CONFIG.FENEK, CONFIG.FENER0, CONFIG.FENEepsilon, CONFIG.FENEsigma))
        lmp.command("special_bonds fene") #<=== I M P O R T A N T (new command)
    if connectivity == "harmonic":
        lmp.command("bond_style harmonic")
        lmp.command("bond_coeff * 100.0 1.2")

    ##############################
    # set timestep of integrator #
    ##############################
    lmp.command("timestep %f" % CONFIG.timestep)

# This splits the lammps calculations on different processors:
def lammps_simulate(lammps_folder, run_time,
                    initial_conformation=None,
                    connectivity="FENE",
                    initial_seed=0, n_models=500, n_keep=100,
                    neighbor=CONFIG.neighbor, tethering=True,
                    minimize=True, compress_with_pbc=False,
                    compress_without_pbc=False,
                    keep_restart_out_dir=None, outfile=None, n_cpus=1,
                    confining_environment=['cube',100.],
                    steering_pairs=None,
                    time_dependent_steering_pairs=None,
                    loop_extrusion_dynamics=None, cleanup = True,
                    to_dump=100000, pbc=False, timeout_job=3600,
                    hide_log=True,
                    restart_path=False,
                    store_n_steps=10,
                    useColvars=False):

    """
    This function launches jobs to generate three-dimensional models in lammps
    
    :param initial_conformation: structural _models object with the particles initial conformation. 
            http://lammps.sandia.gov/doc/2001/data_format.html
    :param FENE connectivity: use FENE for a fene bond or harmonic for harmonic potential
        for neighbours (see init_lammps for details)
    :param run_time: # of timesteps.
    :param None steering_pairs: dictionary with all the info to perform
            steered molecular dynamics.
            steering_pairs = { 'colvar_input'              : "ENST00000540866.2chr7_clean_enMatch.txt",
                               'colvar_output'             : "colvar_list.txt",
                               'kappa_vs_genomic_distance' : "kappa_vs_genomic_distance.txt",
                               'chrlength'                 : 3182,
                               'copies'                    : ['A'],
                               'binsize'                   : 50000,
                               'number_of_kincrease'       : 1000,
                               'timesteps_per_k'           : 1000,
                               'timesteps_relaxation'      : 100000,
                               'perc_enfor_contacts'       : 10
                             }
            Should at least contain Chromosome, loci1, loci2 as 1st, 2nd and 3rd column

    :param None loop_extrusion_dynamics: dictionary with all the info to perform loop
            extrusion dynamics.
            loop_extrusion_dynamics = { 'target_loops_input'          : "target_loops.txt",
                                        'loop_extrusion_steps_output' : "loop_extrusion_steps.txt",
                                        'attraction_strength'         : 4.0,
                                        'equilibrium_distance'        : 1.0,
                                        'chrlength'                   : 3182,
                                        'copies'                      : ['A'],
                                        'timesteps_per_loop_extrusion_step' : 1000,
                                        'timesteps_relaxation'        : 100000,
                                        'perc_enfor_loops'            : 10
                             }

            Should at least contain Chromosome, loci1, loci2 as 1st, 2nd and 3rd column 

    :param 0 initial_seed: Initial random seed for modelling
    :param 500 n_models: number of models to generate.
    :param CONFIG.neighbor neighbor: see LAMMPS_CONFIG.py.
    :param True minimize: whether to apply minimize command or not. 
    :param None keep_restart_out_dir: path to write files to restore LAMMPs
                session (binary)
    :param None outfile: store result in outfile
    :param 1 n_cpus: number of CPUs to use.
    :param False restart_path: path to files to restore LAMMPs session (binary)
    :param 10 store_n_steps: Integer with number of steps to be saved if 
        restart_file != False
    :param False useColvars: True if you want the restrains to be loaded by colvars

    :returns: a StructuralModels object

    """
    
    if confining_environment[0] != 'cube' and pbc == True:
        print "ERROR: It is not possible to implement the pbc"
        print "for simulations inside a %s" % (confining_environment[0])    

    if initial_seed:
        seed(initial_seed)

    #pool = mu.Pool(n_cpus)
    timepoints = 1
    if time_dependent_steering_pairs:
        timepoints = (len(time_dependent_steering_pairs['colvar_input'])-1) * \
            time_dependent_steering_pairs['colvar_dump_freq']

    chromosome_particle_numbers = [int(x) for x in [len(LOCI)]]
    chromosome_particle_numbers.sort(key=int,reverse=True)

    kseeds = []
    for k in xrange(n_models):
        kseeds.append(k+1+initial_seed)
    #while len(kseeds) < n_models:
    #    rnd = randint(1,100000000)
    #    if all([(abs(ks - rnd) > timepoints) for ks in kseeds]):
    #        kseeds.append(rnd)

    #pool = ProcessPool(max_workers=n_cpus, max_tasks=n_cpus)
    pool = multiprocessing.Pool(processes=n_cpus, maxtasksperchild=n_cpus)

    results = []
    def collect_result(result):
        results.append((result[0], result[1]))

    jobs = {}
    for k_id, k in enumerate(kseeds):
        k_folder = lammps_folder + 'lammps_' + str(k) + '/'
        failedSeedLog = None
        # First we check if the modelling fails with this seed
        if restart_path != False:
            restart_file = restart_path + 'lammps_' + str(k) + '/'
            failedSeedLog = restart_file + 'runLog.txt'
            if os.path.exists(failedSeedLog):
                with open(failedSeedLog, 'r') as f:
                    for line in f:
                        prevRun = line.split()
                # add number of models done so dont repeat same seed
                if prevRun[1] == 'Failed':
                    k = int(prevRun[0]) + n_models
                    k_folder = lammps_folder + 'lammps_' + str(k) + '/'

        #print "#RandomSeed: %s" % k
        keep_restart_out_dir2 = None
        if keep_restart_out_dir != None:
            keep_restart_out_dir2 = keep_restart_out_dir + 'lammps_' + str(k) + '/'
            if not os.path.exists(keep_restart_out_dir2):
                os.makedirs(keep_restart_out_dir2)
        model_path = False
        if restart_path != False:
            # check presence of previously finished jobs
            model_path = restart_path + 'lammps_' + str(k) + '/finishedModel_%s.pickle' %k
        # define restart file by checking for finished jobs or last step
        if model_path != False and os.path.exists(model_path):
            with open(model_path, "rb") as input_file:
                m = load(input_file)
            results.append((m[0], m[1]))
        else:
            if restart_path != False:
                restart_file = restart_path + 'lammps_' + str(k) + '/'
                dirfiles = os.listdir(restart_file)
                # check for last k and step
                maxi = (0, 0, '')
                for f in dirfiles:
                    if f.startswith('restart_kincrease_'):
                        kincrease = int(f.split('_')[2])
                        step = int(f.split('_')[-1][:-8])
                        if kincrease > maxi[0]:
                            maxi = (kincrease, step, f)
                        elif kincrease == maxi[0] and step > maxi[1]:
                            maxi = (kincrease, step, f)
                # In case there is no restart file at all
                if maxi[2] == '':
                    print 'Could not find a LAMMPS restart file'
                    # will check later if we have a path or a file
                    getIniConf = True
                    #restart_file = False
                else:
                    restart_file = restart_file + maxi[2]
                    getIniConf = False
            else:
                restart_file = False
                getIniConf = True

            ini_conf = None
            if not os.path.exists(k_folder):
                os.makedirs(k_folder)
                if initial_conformation and getIniConf == True:
                    ini_conf = '%sinitial_conformation.dat' % k_folder
                    write_initial_conformation_file(initial_conformation[k_id],
                                                    chromosome_particle_numbers,
                                                    confining_environment,
                                                    out_file=ini_conf)
    #         jobs[k] = run_lammps(k, k_folder, run_time,
    #                                               initial_conformation, connectivity,
    #                                               neighbor,
    #                                               tethering, minimize,
    #                                               compress_with_pbc, compress_without_pbc,
    #                                               confining_environment,
    #                                               steering_pairs,
    #                                               time_dependent_steering_pairs,
    #                                               loop_extrusion_dynamics,
    #                                               to_dump, pbc,)
    #       jobs[k] = pool.schedule(run_lammps,
            jobs[k] = partial(abortable_worker, run_lammps, timeout=timeout_job,
                                failedSeedLog=[failedSeedLog, k])
            pool.apply_async(jobs[k],
                            args=(k, k_folder, run_time,
                                ini_conf, connectivity,
                                neighbor,
                                tethering, minimize,
                                compress_with_pbc, compress_without_pbc,
                                confining_environment,
                                steering_pairs,
                                time_dependent_steering_pairs,
                                loop_extrusion_dynamics,
                                to_dump, pbc, hide_log,
                                keep_restart_out_dir2,
                                restart_file,
                                model_path,
                                store_n_steps,
                                useColvars,), callback=collect_result)
    #                         , timeout=timeout_job)

    pool.close()
    pool.join()

#     for k in kseeds:
#         try:
#             #results.append((k, jobs[k]))
#             results.append((k, jobs[k].result()))
#         except TimeoutError:
#             print "Model took more than %s seconds to complete ... canceling" % str(timeout_job)
#             jobs[k].cancel()
#         except Exception as error:
#             print "Function raised %s" % error
#             jobs[k].cancel()

    models = {}
    if timepoints > 1:
        for t in xrange(timepoints):
            time_models = []
            for res in results:
                (k,restarr) = res
                time_models.append(restarr[t])
            for i, m in enumerate(time_models[:n_keep]):
                models[i+t*len(time_models[:n_keep])+n_keep] = m
            #for i, (_, m) in enumerate(
            #    sorted(time_models.items(), key=lambda x: x[1]['objfun'])[:n_keep]):
            #    models[i+t+1] = m

    else:
        for i, (_, m) in enumerate(
            sorted(results, key=lambda x: x[1][0]['objfun'])[:n_keep]):
            models[i] = m[0]

    if cleanup:
        for k in kseeds:
            k_folder = lammps_folder + '/lammps_' + str(k) + '/'
            if os.path.exists(k_folder):
                shutil.rmtree(k_folder)

    return models

    
    
# This performs the dynamics: I should add here: The steered dynamics (Irene and Hi-C based) ; 
# The loop extrusion dynamics ; the binders based dynamics (Marenduzzo and Nicodemi)...etc...
def run_lammps(kseed, lammps_folder, run_time,
               initial_conformation=None, connectivity="FENE",
               neighbor=CONFIG.neighbor,
               tethering=False, minimize=True,
               compress_with_pbc=None, compress_without_pbc=None,
               confining_environment=None,
               steering_pairs=None,
               time_dependent_steering_pairs=None,
               loop_extrusion_dynamics=None,
               to_dump=10000, pbc=False,
               hide_log=True,
               keep_restart_out_dir2=None,
               restart_file=False,
               model_path=False, 
               store_n_steps=10,
               useColvars=False):
    """
    Generates one lammps model
    
    :param kseed: Random number to identify the model.
    :param initial_conformation_folder: folder where to store lammps input 
        data file with the particles initial conformation. 
        http://lammps.sandia.gov/doc/2001/data_format.html
    :param FENE connectivity: use FENE for a fene bond or harmonic for harmonic
        potential for neighbours (see init_lammps_run) 
    :param run_time: # of timesteps.
    :param None initial_conformation: path to initial conformation file or None 
        for random walk initial start.
    :param CONFIG.neighbor neighbor: see LAMMPS_CONFIG.py.
    :param False tethering: whether to apply tethering command or not.
    :param True minimize: whether to apply minimize command or not. 
    :param None compress_with_pbc: whether to apply the compression dynamics in case of a
      system with cubic confinement and pbc. This compression step is usually apply 
      to obtain a system with the desired particle density. The input have to be a list 
      of three elements:
      0 - XXX;
      1 - XXX;
      2 - The compression simulation time span (in timesteps).
      e.g. compress_with_pbc=[0.01, 0.01, 100000]
    :param None compress_without_pbc: whether to apply the compression dynamics in case of a
      system with spherical confinement. This compression step is usually apply to obtain a 
      system with the desired particle density. The simulation shrinks/expands the initial 
      sphere to a sphere of the desired radius using many short runs. In each short run the
      radius is reduced by 0.1 box units. The input have to be a list of three elements:
      0 - Initial radius;
      1 - Final desired radius;
      2 - The time span (in timesteps) of each short compression run.
      e.g. compress_without_pbc=[300, 100, 100]
    :param None steering_pairs: particles contacts file from colvars fix 
      http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf. 
      steering_pairs = { 'colvar_input'              : "ENST00000540866.2chr7_clean_enMatch.txt",
                         'colvar_output'             : "colvar_list.txt",
                         'kappa_vs_genomic_distance' : "kappa_vs_genomic_distance.txt",
                         'chrlength'                 : 3182,
                         'copies'                    : ['A'],
                         'binsize'                   : 50000,
                         'number_of_kincrease'       : 1000,
                         'timesteps_per_k'           : 1000,
                         'timesteps_relaxation'      : 100000,
                         'perc_enfor_contacts'       : 10
                       }

    :param None loop_extrusion_dynamics: dictionary with all the info to perform loop 
            extrusion dynamics.
            loop_extrusion_dynamics = { 'target_loops_input'          : "target_loops.txt",
                                        'loop_extrusion_steps_output' : "loop_extrusion_steps.txt",
                                        'attraction_strength'         : 4.0,
                                        'equilibrium_distance'        : 1.0,
                                        'chrlength'                   : 3182,
                                        'copies'                      : ['A'],
                                        'timesteps_per_loop_extrusion_step' : 1000,
                                        'timesteps_relaxation'        : 100000,
                                        'perc_enfor_loops'            : 10
                             }

            Should at least contain Chromosome, loci1, loci2 as 1st, 2nd and 3rd column 
    :param None keep_restart_out_dir2: path to write files to restore LAMMPs
                session (binary)
    :param False restart_file: path to file to restore LAMMPs session (binary)
    :param False model_path: path to/for pickle with finished model (name included)
    :param 10 store_n_steps: Integer with number of steps to be saved if 
        restart_file != False
    :param False useColvars: True if you want the restrains to be loaded by colvars
    :returns: a LAMMPSModel object

    """

    lmp = lammps(cmdargs=['-screen','none','-log',lammps_folder+'log.lammps','-nocite'])
    # check if we have a restart file or a path to which restart
    if restart_file == False:
        doRestart = False
        saveRestart = False
    elif os.path.isdir(restart_file):
        doRestart = False
        saveRestart = True
    else:
        doRestart = True
        saveRestart = True
    if not initial_conformation and doRestart == False:    
        initial_conformation = lammps_folder+'initial_conformation.dat'
        generate_chromosome_random_walks_conformation ([len(LOCI)],
                                                       outfile=initial_conformation,
                                                       seed_of_the_random_number_generator=kseed,
                                                       confining_environment=confining_environment)
    
    # Just prepared the steps recovery for steering pairs
    if steering_pairs and doRestart == True:
        init_lammps_run(lmp, initial_conformation,
                neighbor=neighbor,
                hide_log=hide_log,
                connectivity=connectivity,
                restart_file=restart_file)
    else:
        init_lammps_run(lmp, initial_conformation,
                    neighbor=neighbor,
                    hide_log=hide_log,
                    connectivity=connectivity)

    lmp.command("dump    1       all    custom    %i   %slangevin_dynamics_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
    #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")

    # ##########################################################
    # # Generate RESTART file, SPECIAL format, not a .txt file #
    # # Useful if simulation crashes             
    # Prepared an optimisation for steering pairs, but not for the rest#
    # ##########################################################
    # create lammps restart files every x steps. 1000 is ok
    # There was the doubt of using text format session info (which allows use in other computers)
    # but since the binary can be converted later and this: "Because a data file is in text format, 
    # if you use a data file written out by this command to restart a simulation, the initial state 
    # of the new run will be slightly different than the final state of the old run (when the file 
    # was written) which was represented internally by LAMMPS in binary format. A new simulation 
    # which reads the data file will thus typically diverge from a simulation that continued 
    # in the original input script." will continue with binary. To convert use restart2data
    #if keep_restart_out_dir2:
    #    lmp.command("restart %i %s/relaxation_%i_*.restart" % (keep_restart_step, keep_restart_out_dir2, kseed))


    #######################################################
    # Set up fixes                                        #
    # use NVE ensemble                                    #
    # Langevin integrator Tstart Tstop 1/friction rndseed #
    # => sampling NVT ensamble                            #
    #######################################################
    # Define the langevin dynamics integrator
    lmp.command("fix 1 all nve")
    lmp.command("fix 2 all langevin 1.0  1.0  2.0 %i" % kseed)
    # Define the tethering to the center of the confining environment
    if tethering:
        lmp.command("fix 3 all spring tether 50.0 0.0 0.0 0.0 0.0")

    # Do a minimization step to prevent particles
    # clashes in the initial conformation
    if minimize:

        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %sminimization_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")
        
        print "Performing minimization run..."
        lmp.command("minimize 1.0e-4 1.0e-6 100000 100000")
        
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %slangevin_dynamics_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")        

    if compress_with_pbc:
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %scompress_with_pbc_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")

        # Re-setting the initial timestep to 0
        lmp.command("reset_timestep 0")

        lmp.command("unfix 1")
        lmp.command("unfix 2")

        # default as in PLoS Comp Biol Di Stefano et al. 2013 compress_with_pbc = [0.01, 0.01, 100000]
        lmp.command("fix 1 all   nph   iso   %s %s   2.0" % (compress_with_pbc[0], 
                                                             compress_with_pbc[1]))
        lmp.command("fix 2 all langevin 1.0  1.0  2.0 %i" % kseed)
        print "run %i" % compress_with_pbc[2]
        lmp.command("run %i" % compress_with_pbc[2])

        lmp.command("unfix 1")
        lmp.command("unfix 2")

        lmp.command("fix 1 all nve")
        lmp.command("fix 2 all langevin 1.0  1.0  2.0 %i" % kseed)        

        # Here We have to re-define the confining environment
        print "# Previous particle density (nparticles/volume)", lmp.get_natoms()/(confining_environment[1]**3)
        confining_environment[1] = lmp.extract_global("boxxhi",1) - lmp.extract_global("boxxlo",1)
        print ""
        print "# New cubic box dimensions after isotropic compression"
        print lmp.extract_global("boxxlo",1), lmp.extract_global("boxxhi",1)
        print lmp.extract_global("boxylo",1), lmp.extract_global("boxyhi",1)
        print lmp.extract_global("boxzlo",1), lmp.extract_global("boxzhi",1)
        print "# New confining environment", confining_environment
        print "# New particle density (nparticles/volume)", lmp.get_natoms()/(confining_environment[1]**3)
        print ""

        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %slangevin_dynamics_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")        

    if compress_without_pbc:
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %scompress_without_pbc_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")

        # Re-setting the initial timestep to 0
        lmp.command("reset_timestep 0")

        # default as in Sci Rep Di Stefano et al. 2016 
        # compress_without_pbc = [initial_radius, final_radius, timesteps_per_minirun] 
        # = [350, 161.74, 100]
        radius = compress_without_pbc[0]
        while radius > compress_without_pbc[1]:

            print "New radius %f" % radius
            if radius != compress_without_pbc[0]:
                lmp.command("region sphere delete")
            
            lmp.command("region sphere sphere 0.0 0.0 0.0 %f units box side in" % radius)

            # Performing the simulation
            lmp.command("fix 5 all  wall/region sphere lj126 1.0 1.0 1.12246152962189")
            lmp.command("run %i" % compress_without_pbc[2])

            radius -= 0.1

        # Here we have to re-define the confining environment
        volume = 4.*np.pi/3.0*(compress_without_pbc[0]**3)
        print "# Previous particle density (nparticles/volume)", lmp.get_natoms()/volume
        confining_environment[1] = compress_without_pbc[1]
        print ""
        volume = 4.*np.pi/3.0*(compress_without_pbc[1]**3)
        print "# New particle density (nparticles/volume)", lmp.get_natoms()/volume
        print ""

    timepoints = 1
    xc = []
    # Setup the pairs to co-localize using the COLVARS plug-in
    if steering_pairs:
        
        if doRestart == False:
            # Start relaxation step
            lmp.command("reset_timestep 0")   # cambiar para punto ionicial
            lmp.command("run %i" % steering_pairs['timesteps_relaxation'])
            lmp.command("reset_timestep %i" % 0)
        
            # Start Steered Langevin dynamics
            if to_dump:
                lmp.command("undump 1")
                lmp.command("dump    1       all    custom    %i   %ssteered_MD_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
                #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id")

        if 'number_of_kincrease' in steering_pairs:
            nbr_kincr = steering_pairs['number_of_kincrease']
        else:
            nbr_kincr = 1
        
        if doRestart == True:
            restart_k_increase = int(restart_file.split('/')[-1].split('_')[2])
            restart_time       = int(restart_file.split('/')[-1].split('_')[4][:-8])

        #steering_pairs['colvar_output'] = os.path.dirname(os.path.abspath(steering_pairs['colvar_output'])) + '/' + str(kseed) + '_'+ os.path.basename(steering_pairs['colvar_output'])    
        steering_pairs['colvar_output'] = lammps_folder+os.path.basename(steering_pairs['colvar_output'])
        for kincrease in xrange(nbr_kincr):
            # Write the file containing the pairs to constraint
            # steering_pairs should be a dictionary with:
            # Avoid to repeat calculations in case of restart
            if (doRestart == True) and (kincrease < restart_k_increase):
                continue

            if useColvars == True:
                
                generate_colvars_list(steering_pairs, kincrease+1)

                # Adding the colvar option
                #print "fix 4 all colvars %s output %s" % (steering_pairs['colvar_output'],lammps_folder)
                lmp.command("fix 4 all colvars %s output %sout" % (steering_pairs['colvar_output'],lammps_folder))

                if to_dump:
                    # lmp.command("thermo_style   custom   step temp epair emol")
                    lmp.command("thermo_style   custom   step temp epair emol pe ke etotal f_4")
                    lmp.command("thermo_modify norm no flush yes")
                    lmp.command("variable step equal step")
                    lmp.command("variable objfun equal f_4")
                    lmp.command('''fix 5 all print %s "${step} ${objfun}" file "%sobj_fun_from_time_point_%s_to_time_point_%s.txt" screen "no" title "#Timestep Objective_Function"''' % (steering_pairs['colvar_dump_freq'],lammps_folder,str(0), str(1)))

            # will load the bonds directly into LAMMPS
            else:
                bond_list = generate_bond_list(steering_pairs)
                for bond in bond_list:
                    lmp.command(bond)

                if to_dump:
                    lmp.command("thermo_style   custom   step temp etotal")
                    lmp.command("thermo_modify norm no flush yes")
                    lmp.command("variable step equal step")
                    lmp.command("variable objfun equal etotal")
                    lmp.command('''fix 5 all print %s "${step} ${objfun}" file "%sobj_fun_from_time_point_%s_to_time_point_%s.txt" screen "no" title "#Timestep Objective_Function"''' % (steering_pairs['colvar_dump_freq'],lammps_folder,str(0), str(1)))



            #lmp.command("reset_timestep %i" % 0)
            resettime = 0
            runtime   = steering_pairs['timesteps_per_k']
            if (doRestart == True) and (kincrease == restart_k_increase):
                resettime = restart_time 
                runtime   = steering_pairs['timesteps_per_k'] - restart_time

            # Create 10 restarts with name restart_kincrease_%s_time_%s.restart every
            if saveRestart == True:
                if os.path.isdir(restart_file):
                    restart_file_new = restart_file + 'restart_kincrease_%s_time_*.restart' %(kincrease)
                else:
                    restart_file_new = '/'.join(restart_file.split('/')[:-1]) + '/restart_kincrease_%s_time_*.restart' %(kincrease)
                print restart_file_new
                lmp.command("restart %i %s" %(int(steering_pairs['timesteps_per_k']/store_n_steps), restart_file_new))

            #lmp.command("reset_timestep %i" % resettime)
            lmp.command("run %i" % runtime)

    # Setup the pairs to co-localize using the COLVARS plug-in
    if time_dependent_steering_pairs:
        timepoints = time_dependent_steering_pairs['colvar_dump_freq']

        #if exists("objective_function_profile.txt"):
        #    os.remove("objective_function_profile.txt")

        #print "# Getting the time dependent steering pairs!"
        time_dependent_restraints = get_time_dependent_colvars_list(time_dependent_steering_pairs)
        time_points = sorted(time_dependent_restraints.keys())
        print "#Time_points",time_points        
        sys.stdout.flush()            

        time_dependent_steering_pairs['colvar_output'] = lammps_folder+os.path.basename(time_dependent_steering_pairs['colvar_output'])
        # Performing the adaptation step from TADbit to TADdyn excluded volume
        if time_dependent_steering_pairs['adaptation_step']:
            restraints = {}
            for time_point in time_points[0:1]:
                lmp.command("reset_timestep %i" % 0)    
                # Change to_dump with some way to load the conformations you want to store in TADbit
                # This Adaptation could be discarded in the trajectory files.
                if to_dump:
                    lmp.command("undump 1")
                    lmp.command("dump    1       all    custom    %i  %sadapting_MD_from_TADbit_to_TADdyn_at_time_point_%s.XYZ  id  xu yu zu" % (to_dump, lammps_folder, time_point))
                    lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")

                restraints[time_point] = {}
                print "# Step %s - %s" % (time_point, time_point)
                sys.stdout.flush()            
                for pair in time_dependent_restraints[time_point].keys():
                    # Strategy changing gradually the spring constant and the equilibrium distance
                    # Case 1: The restraint is present at time point 0 and time point 1:
                    if pair in time_dependent_restraints[time_point]:
                        # Case A: The restrainttype is the same at time point 0 and time point 1 ->
                        # The spring force changes, and the equilibrium distance is the one at time_point+1
                        restraints[time_point][pair] = [
                            # Restraint type
                            [time_dependent_restraints[time_point][pair][0]], 
                            # Initial spring constant 
                            [time_dependent_restraints[time_point][pair][1]*time_dependent_steering_pairs['k_factor']], 
                            # Final spring constant 
                            [time_dependent_restraints[time_point][pair][1]*time_dependent_steering_pairs['k_factor']], 
                            # Initial equilibrium distance
                            [time_dependent_restraints[time_point][pair][2]], 
                            # Final equilibrium distance
                            [time_dependent_restraints[time_point][pair][2]], 
                            # Number of timesteps for the gradual change
                            [int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point]*0.1)]]

                generate_time_dependent_colvars_list(restraints[time_point], time_dependent_steering_pairs['colvar_output'], time_dependent_steering_pairs['colvar_dump_freq'])
                copyfile(time_dependent_steering_pairs['colvar_output'], 
                         "colvar_list_from_time_point_%s_to_time_point_%s.txt" % 
                         (str(time_point), str(time_point)))

                lmp.command("velocity all create 1.0 %s" % kseed)
                # Adding the colvar option and perfoming the steering
                if time_point != time_points[0]:
                    lmp.command("unfix 4")
                print "#fix 4 all colvars %s" % time_dependent_steering_pairs['colvar_output']
                sys.stdout.flush()
                lmp.command("fix 4 all colvars %s tstat 2 output %sout" % (time_dependent_steering_pairs['colvar_output'],lammps_folder))
                lmp.command("run %i" % int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point]*0.1))

        # Time dependent steering
        restraints = {}
        #for i in xrange(time_points[0],time_points[-1]):
        for time_point in time_points[0:-1]:
            lmp.command("reset_timestep %i" % 0)    
            # Change to_dump with some way to load the conformations you want to store in TADbit
            if to_dump:
                lmp.command("undump 1")
                lmp.command("dump    1       all    custom    %i   %ssteered_MD_from_time_point_%s_to_time_point_%s.XYZ  id  xu yu zu" % (to_dump, lammps_folder, time_point, time_point+1))
                lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")

            restraints[time_point] = {}
            print "# Step %s - %s" % (time_point, time_point+1)
            sys.stdout.flush()            
            # Compute the current distance between any two particles
            xc_tmp = np.array(lmp.gather_atoms("x",1,3))
            current_distances = compute_particles_distance(xc_tmp)

            for pair in set(time_dependent_restraints[time_point].keys()+time_dependent_restraints[time_point+1].keys()):                
                r = 0
                
                # Strategy changing gradually the spring constant
                # Case 1: The restraint is present at time point 0 and time point 1:
                if pair     in time_dependent_restraints[time_point] and pair     in time_dependent_restraints[time_point+1]:
                    # Case A: The restrainttype is the same at time point 0 and time point 1 ->
                    # The spring force changes, and the equilibrium distance is the one at time_point+1
                    if time_dependent_restraints[time_point][pair][0]   == time_dependent_restraints[time_point+1][pair][0]:
                        r += 1
                        restraints[time_point][pair] = [
                            # Restraint type
                            [time_dependent_restraints[time_point+1][pair][0]], 
                            # Initial spring constant 
                            [time_dependent_restraints[time_point][pair][1]  *time_dependent_steering_pairs['k_factor']], 
                            # Final spring constant 
                            [time_dependent_restraints[time_point+1][pair][1]*time_dependent_steering_pairs['k_factor']], 
                            # Initial equilibrium distance
                            [time_dependent_restraints[time_point][pair][2]], 
                            # Final equilibrium distance
                            [time_dependent_restraints[time_point+1][pair][2]], 
                            # Number of timesteps for the gradual change
                            [int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point])]]
                    # Case B: The restrainttype is different between time point 0 and time point 1
                    if time_dependent_restraints[time_point][pair][0]   != time_dependent_restraints[time_point+1][pair][0]:
                        # Case a: The restrainttype is "Harmonic" at time point 0 
                        # and "LowerBoundHarmonic" at time point 1                        
                        if time_dependent_restraints[time_point][pair][0] == "Harmonic":
                            r += 1
                            restraints[time_point][pair] = [
                                # Restraint type
                                [time_dependent_restraints[time_point][pair][0], time_dependent_restraints[time_point+1][pair][0]], 
                                # Initial spring constant 
                                [time_dependent_restraints[time_point][pair][1]*time_dependent_steering_pairs['k_factor'], 0.0],
                                # Final spring constant 
                                [0.0, time_dependent_restraints[time_point+1][pair][1]*time_dependent_steering_pairs['k_factor']],
                                # Initial equilibrium distance
                                [time_dependent_restraints[time_point][pair][2], time_dependent_restraints[time_point][pair][2]],
                                # Final equilibrium distance
                                [time_dependent_restraints[time_point+1][pair][2], time_dependent_restraints[time_point+1][pair][2]],
                                # Number of timesteps for the gradual change
                                #[int(time_dependent_steering_pairs['timesteps_per_k_change']), int(time_dependent_steering_pairs['timesteps_per_k_change'])]]
                                [int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point]), int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point])]]
                        # Case b: The restrainttype is "LowerBoundHarmonic" at time point 0 
                        # and "Harmonic" at time point 1
                        if time_dependent_restraints[time_point][pair][0] == "HarmonicLowerBound":
                            r += 1
                            restraints[time_point][pair] = [
                                # Restraint type
                                [time_dependent_restraints[time_point][pair][0], time_dependent_restraints[time_point+1][pair][0]], 
                                # Initial spring constant 
                                [time_dependent_restraints[time_point][pair][1]*time_dependent_steering_pairs['k_factor'], 0.0],
                                # Final spring constant 
                                [0.0, time_dependent_restraints[time_point+1][pair][1]*time_dependent_steering_pairs['k_factor']],
                                # Initial equilibrium distance
                                [time_dependent_restraints[time_point][pair][2], time_dependent_restraints[time_point][pair][2]],
                                # Final equilibrium distance
                                [time_dependent_restraints[time_point+1][pair][2], time_dependent_restraints[time_point+1][pair][2]],
                                # Number of timesteps for the gradual change
                                #[int(time_dependent_steering_pairs['timesteps_per_k_change']), int(time_dependent_steering_pairs['timesteps_per_k_change'])]]
                                [int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point]), int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point])]]

                # Case 2: The restraint is not present at time point 0, but it is at time point 1:                            
                elif pair not in time_dependent_restraints[time_point] and pair     in time_dependent_restraints[time_point+1]:
                    # List content: restraint_type,kforce,distance
                    r += 1
                    restraints[time_point][pair] = [
                        # Restraint type -> Is the one at time point time_point+1
                        [time_dependent_restraints[time_point+1][pair][0]],
                        # Initial spring constant 
                        [0.0],
                        # Final spring constant 
                        [time_dependent_restraints[time_point+1][pair][1]*time_dependent_steering_pairs['k_factor']], 
                        # Initial equilibrium distance 
                        [time_dependent_restraints[time_point+1][pair][2]], 
                        # Final equilibrium distance 
                        [time_dependent_restraints[time_point+1][pair][2]], 
                        # Number of timesteps for the gradual change
                        [int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point])]] 

                # Case 3: The restraint is     present at time point 0, but it is not at time point 1:                            
                elif pair     in time_dependent_restraints[time_point] and pair not in time_dependent_restraints[time_point+1]:
                    # List content: restraint_type,kforce,distance
                    r += 1
                    restraints[time_point][pair] = [
                        # Restraint type -> Is the one at time point time_point
                        [time_dependent_restraints[time_point][pair][0]], 
                        # Initial spring constant 
                        [time_dependent_restraints[time_point][pair][1]*time_dependent_steering_pairs['k_factor']],                         
                        # Final spring constant 
                        [0.0],
                        # Initial equilibrium distance 
                        [time_dependent_restraints[time_point][pair][2]],                         
                        # Final equilibrium distance 
                        [time_dependent_restraints[time_point][pair][2]], 
                        # Number of timesteps for the gradual change
                        [int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point])]]
                
                    #current_distances[pair],                          
                else:
                    print "#ERROR None of the previous conditions is matched!"
                    if pair     in time_dependent_restraints[time_point]:
                        print "# Pair %s at timepoint %s %s  " % (pair, time_point, time_dependent_restraints[time_point][pair])
                    if pair     in time_dependent_restraints[time_point+1]:
                        print "# Pair %s at timepoint %s %s  " % (pair, time_point+1, time_dependent_restraints[time_point+1][pair])
                    continue

                if r > 1:
                    print "#ERROR Two of the previous conditions are matched!"

                #if pair     in time_dependent_restraints[time_point]:
                #    print "# Pair %s at timepoint %s %s  " % (pair, time_point, time_dependent_restraints[time_point][pair])
                #else:
                #    print "# Pair %s at timepoint %s None" % (pair, time_point)

                #if pair     in time_dependent_restraints[time_point+1]:
                #    print "# Pair %s at timepoint %s %s  " % (pair, time_point+1, time_dependent_restraints[time_point+1][pair])
                #else:
                #    print "# Pair %s at timepoint %s None" % (pair, time_point+1)
                #print restraints[pair]
                #print ""

            generate_time_dependent_colvars_list(restraints[time_point], time_dependent_steering_pairs['colvar_output'], time_dependent_steering_pairs['colvar_dump_freq'])
            copyfile(time_dependent_steering_pairs['colvar_output'], 
                     "%scolvar_list_from_time_point_%s_to_time_point_%s.txt" % 
                     (lammps_folder, str(time_point), str(time_point+1)))

            lmp.command("velocity all create 1.0 %s" % kseed)
            # Adding the colvar option and perfoming the steering
            if time_point != time_points[0]:
                lmp.command("unfix 4")
            print "#fix 4 all colvars %s" % time_dependent_steering_pairs['colvar_output']
            sys.stdout.flush()
            lmp.command("fix 4 all colvars %s tstat 2 output %sout" % (time_dependent_steering_pairs['colvar_output'],lammps_folder))
            if to_dump:
                lmp.command("thermo_style   custom   step temp epair emol pe ke etotal f_4")
                lmp.command("thermo_modify norm no flush yes")
                lmp.command("variable step equal step")
                lmp.command("variable objfun equal f_4")
                lmp.command('''fix 5 all print %s "${step} ${objfun}" file "%sobj_fun_from_time_point_%s_to_time_point_%s.txt" screen "no" title "#Timestep Objective_Function"''' % (time_dependent_steering_pairs['colvar_dump_freq'],lammps_folder,str(time_point), str(time_point+1)))
            
            lmp.command("run %i" % int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point]))
            
            if time_point > 0:
                    
                if exists("%sout.colvars.traj.BAK" % lammps_folder):

                    copyfile("%sout.colvars.traj.BAK" % lammps_folder, "%srestrained_pairs_equilibrium_distance_vs_timestep_from_time_point_%s_to_time_point_%s.txt" % (lammps_folder, str(time_point-1), str(time_point)))
            
                    os.remove("%sout.colvars.traj.BAK" % lammps_folder)

    # Setup the pairs to co-localize using the COLVARS plug-in
    if loop_extrusion_dynamics:

        # Start relaxation step
        lmp.command("reset_timestep 0")
        lmp.command("run %i" % loop_extrusion_dynamics['timesteps_relaxation'])

        lmp.command("reset_timestep 0")
        # Start Loop extrusion dynamics
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %sloop_extrusion_MD_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id")

        # List of target loops of the form [(loop1_start,loop1_stop),...,(loopN_start,loopN_stop)]
        target_loops  = read_target_loops_input(loop_extrusion_dynamics['target_loops_input'],
                                                loop_extrusion_dynamics['chrlength'],
                                                loop_extrusion_dynamics['perc_enfor_loops'])

        # Randomly extract starting point of the extrusion dynamics between start and stop
        initial_loops = draw_loop_extrusion_starting_points(target_loops,
                                                            loop_extrusion_dynamics['chrlength'])

        # Maximum number of particles to be extruded during the dynamics
        maximum_number_of_extruded_particles = get_maximum_number_of_extruded_particles(target_loops, initial_loops)
        print "Number of LES",maximum_number_of_extruded_particles

        # Loop extrusion steps
        for LES in xrange(1,maximum_number_of_extruded_particles):

            # Loop extrusion steps

            # Update the Lennard-Jones coefficients between extruded particles
            loop_number = 1
            for particle1,particle2 in initial_loops:

                print "# fix LE%i all restrain bond %i  %i %f %f %f" % (loop_number,
                                                                        particle1,
                                                                        particle2,
                                                                        loop_extrusion_dynamics['attraction_strength'],
                                                                        loop_extrusion_dynamics['attraction_strength'],
                                                                        loop_extrusion_dynamics['equilibrium_distance'])
                
                lmp.command("fix LE%i all restrain bond %i  %i %f %f %f" % (loop_number,
                                                                            particle1,
                                                                            particle2,
                                                                            loop_extrusion_dynamics['attraction_strength'],
                                                                            loop_extrusion_dynamics['attraction_strength'],
                                                                            loop_extrusion_dynamics['equilibrium_distance']))

                loop_number += 1

            # Doing the LES
            lmp.command("run %i" % loop_extrusion_dynamics['timesteps_per_loop_extrusion_step'])

            # Remove the loop extrusion restraint!
            loop_number = 1
            for particle1,particle2 in initial_loops:

                print "# unfix LE%i" % (loop_number)
                lmp.command("unfix LE%i" % (loop_number))

                loop_number += 1

            # Update the particles involved in the loop extrusion interaction:
            # decrease the intial_start by one until you get to start
            # increase the intial_stop by one until you get to stop
            for initial_loop, target_loop in zip(initial_loops,target_loops):

                if initial_loop[0] > target_loop[0]:
                    initial_loop[0] -= 1
                if initial_loop[1] < target_loop[1]:
                    initial_loop[1] += 1



    #if to_dump:
    #    lmp.command("undump 1")
    #    lmp.command("dump    1       all    custom    %i   langevin_dynamics_*.XYZ  id  xu yu zu" % to_dump)
    #    lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")

    # Post-processing analysis
    if time_dependent_steering_pairs:
        
        copyfile("%sout.colvars.traj" % lammps_folder, "%srestrained_pairs_equilibrium_distance_vs_timestep_from_time_point_%s_to_time_point_%s.txt" % (lammps_folder, str(time_point), str(time_point+1)))


        os.remove("%sout.colvars.traj" % lammps_folder)
        os.remove(time_dependent_steering_pairs['colvar_output'])
        for time_point in time_points[0:-1]:
            # Compute energy associated to the restraints: something like the IMP objective function
            #compute_the_objective_function("%srestrained_pairs_equilibrium_distance_vs_timestep_from_time_point_%s_to_time_point_%s.txt" % (lammps_folder, str(time_point), str(time_point+1)),
            #                               "%sobjective_function_profile_from_time_point_%s_to_time_point_%s.txt" % (lammps_folder, str(time_point), str(time_point+1)),
            #                               time_point,
            #                               time_dependent_steering_pairs['timesteps_per_k_change'][time_point])
        
            # Compute the % of satysfied constraints between 2. sigma = 2./sqrt(k)
            compute_the_percentage_of_satysfied_restraints("%srestrained_pairs_equilibrium_distance_vs_timestep_from_time_point_%s_to_time_point_%s.txt" % (lammps_folder, str(time_point), str(time_point+1)),
                                                           restraints[time_point],
                                                           "%spercentage_of_established_restraints_from_time_point_%s_to_time_point_%s.txt" % (lammps_folder, str(time_point), str(time_point+1)),
                                                           time_point,
                                                           time_dependent_steering_pairs['timesteps_per_k_change'][time_point])
            
        for time_point in time_points[0:-1]:        
            xc.append(np.array(read_trajectory_file("%ssteered_MD_from_time_point_%s_to_time_point_%s.XYZ" % (lammps_folder, time_point, time_point+1))))
    
    else:    
        # Managing the final model
        xc.append(np.array(lmp.gather_atoms("x",1,3)))
            
    lmp.close()    
        
    result = []
    for stg in xrange(len(xc)):
        log_objfun = read_objective_function("%sobj_fun_from_time_point_%s_to_time_point_%s.txt" % (lammps_folder, str(stg), str(stg+1)))
        for timepoint in range(1,timepoints+1):
            lammps_model = LAMMPSmodel({'x'          : [],
                                  'y'          : [],
                                  'z'          : [],
                                  'cluster'    : 'Singleton',
                                  'log_objfun' : log_objfun,
                                  'objfun'     : log_objfun[-1],
                                  'radius'     : float(CONFIG.HiC['resolution'] * CONFIG.HiC['scale'])/2,
                                  'rand_init'  : str(kseed+timepoint)})
        
            if pbc:
                store_conformation_with_pbc(xc[stg], lammps_model, confining_environment)    
            else:
                skip_first = 0
                if time_dependent_steering_pairs:
                    skip_first = 1
                for i in range((timepoint-1+skip_first)*len(LOCI)*3,(timepoint+skip_first)*len(LOCI)*3,3):
                    lammps_model['x'].append(xc[stg][i]*float(CONFIG.HiC['resolution'] * CONFIG.HiC['scale']))
                    lammps_model['y'].append(xc[stg][i+1]*float(CONFIG.HiC['resolution'] * CONFIG.HiC['scale']))
                    lammps_model['z'].append(xc[stg][i+2]*float(CONFIG.HiC['resolution'] * CONFIG.HiC['scale']))
            result.append(lammps_model)

    #os.remove("%slog.cite" % lammps_folder)
    # safe finished model
    if model_path != False:
        with open(model_path, "wb") as output_file:
            dump((kseed,result), output_file)
    ################### Special case for clusters with disk quota
    # Remove the saved steps
    if saveRestart == True:
        if os.path.isdir(restart_file):
            restart_path = restart_file
        else:
            restart_path = '/'.join(restart_file.split('/')[:-1]) + '/'
        for pathfile in os.listdir(restart_path):
            if pathfile.startswith('restart'):
                os.remove(restart_path + pathfile)    
    ##################################################################

    return (kseed,result)

def read_trajectory_file(fname):

    coords=[]
    fhandler = open(fname)
    line = fhandler.next()
    try:
        while True:
            if line.startswith('ITEM: TIMESTEP'):
                while not line.startswith('ITEM: ATOMS'):
                    line = fhandler.next()
                if line.startswith('ITEM: ATOMS'):
                    line = fhandler.next()
            line = line.strip()
            if len(line) == 0:
                continue
            line_vals = line.split()
            coords += [float(line_vals[1]),float(line_vals[2]),float(line_vals[3])]
            line = fhandler.next()
    except StopIteration:
        pass
    fhandler.close()        
            
    return coords    

########## Part to perform the restrained dynamics ##########
# I should add here: The steered dynamics (Irene's and Hi-C based models) ; 
# The loop extrusion dynamics ; the binders based dynamics (Marenduzzo and Nicodemi)...etc...

def linecount(filename):
    """
    Count valid lines of input colvars file
    
    :param filename: input colvars file.
    
    :returns: number of valid contact lines

    """
    
    k = 0
    tfp = open(filename)
    for i, line in enumerate(tfp):   

        if line.startswith('#'):
            continue
        cols_vals = line.split()
        if cols_vals[1] == cols_vals[2]:
            continue
        k += 1
        
    return k

##########

def generate_colvars_list(steering_pairs,
                          kincrease=0,
                          colvars_header='# collective variable: monitor distances\n\ncolvarsTrajFrequency 1000 # output every 1000 steps\ncolvarsRestartFrequency 10000000\n',
                          colvars_template='''

colvar {
  name %s
  # %s %s %i
  distance {
      group1 {
        atomNumbers %i
      }
      group2 {
        atomNumbers %i
      }
  }
}''',
                            colvars_tail = '''

harmonic {
  name h_pot_%s
  colvars %s
  centers %s
  forceConstant %f 
}\n''',                     colvars_harmonic_lower_bound_tail = '''

harmonicWalls {
  name hlb_pot_%s
  colvars %s
  lowerWalls %s
  forceConstant %f 
  lowerWallConstant 1.0
}\n'''
                            ):
                            
    """
    Generates lammps colvars file http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf
    
    :param dict steering_pairs: dictionary containing all the information to write down the
      the input file for the colvars implementation
    :param exisiting_template colvars_header: header template for colvars file.
    :param exisiting_template colvars_template: contact template for colvars file.
    :param exisiting_template colvars_tail: tail template for colvars file.

    """

    # Getting the input
    # XXXThe target_pairs could be also a list as the one in output of get_HiCbased_restraintsXXX
    target_pairs                 = steering_pairs['colvar_input'] 
    outfile                      = steering_pairs['colvar_output'] 
    if 'kappa_vs_genomic_distance' in steering_pairs:
        kappa_vs_genomic_distance    = steering_pairs['kappa_vs_genomic_distance']
    if 'chrlength' in steering_pairs:
        chrlength                    = steering_pairs['chrlength']
    else:
        chrlength                    = 0
    if 'copies' in steering_pairs:
        copies                       = steering_pairs['copies']
    else:
        copies                       = ['A']
    kbin                         = 10000000
    binsize                      = steering_pairs['binsize']
    if 'percentage_enforced_contacts' in steering_pairs:
        percentage_enforced_contacts = steering_pairs['perc_enfor_contacts']
    else:
        percentage_enforced_contacts = 100

    # Here we extract from all the restraints only 
    # a random sub-sample of percentage_enforced_contacts/100*totcolvars
    rand_lines = []
    i=0
    j=0
    if isinstance(target_pairs, str):    
        totcolvars = linecount(target_pairs)
        ncolvars = int(totcolvars*(float(percentage_enforced_contacts)/100))
        
        #print "Number of enforced contacts = %i over %i" % (ncolvars,totcolvars)
        rand_positions = sample(list(range(totcolvars)), ncolvars)
        rand_positions = sorted(rand_positions)
    
        tfp = open(target_pairs)
        with open(target_pairs) as f:
            for line in f:
                line = line.strip()
                if j >= ncolvars:
                    break
                if line.startswith('#'):
                    continue
             
                cols_vals = line.split()
                # Avoid to enforce restraints between the same bin
                if cols_vals[1] == cols_vals[2]:
                    continue
            
                if i == rand_positions[j]:
                    rand_lines.append(line)
                    j += 1
                i += 1
        tfp.close()
    elif isinstance(target_pairs, HiCBasedRestraints):
        
        rand_lines = target_pairs.get_hicbased_restraints()
        totcolvars = len(rand_lines)
        ncolvars = int(totcolvars*(float(percentage_enforced_contacts)/100))
        
        #print "Number of enforced contacts = %i over %i" % (ncolvars,totcolvars)
        rand_positions = sample(list(range(totcolvars)), ncolvars)
        rand_positions = sorted(rand_positions)
        
        
    else:
        print "Unknown target_pairs"
        return    
    
        
    
    #print rand_lines

    seqdists = {}
    poffset=0
    outf = open(outfile,'w')
    outf.write(colvars_header)
    for copy_nbr in copies:
        i = 1
        for line in rand_lines:
            if isinstance(target_pairs, str):   
                cols_vals = line.split()
            else:
                cols_vals = ['chr'] + line
                
            #print cols_vals
            
            if isinstance(target_pairs, HiCBasedRestraints) and cols_vals[3] != "Harmonic" and cols_vals[3] != "HarmonicLowerBound":
                continue
            
            part1_start = int(cols_vals[1])*binsize
            part1_end = (int(cols_vals[1])+1)*binsize
            #print part1_start, part1_end

            part2_start = int(cols_vals[2])*binsize
            part2_end = (int(cols_vals[2])+1)*binsize
            #print part2_start, part2_end

            name = str(i)+copy_nbr  
            seqdist = abs(part1_start-part2_start)
            #print seqdist

            region1 = cols_vals[0] + '_' + str(part1_start) + '_' + str(part1_end)
            region2 = cols_vals[0] + '_' + str(part2_start) + '_' + str(part2_end)

            particle1 = int(cols_vals[1]) + 1 + poffset
            particle2 = int(cols_vals[2]) + 1 + poffset

            seqdists[name] = seqdist

            outf.write(colvars_template % (name,region1,region2,seqdist,particle1,particle2))
            
            if isinstance(target_pairs, HiCBasedRestraints):
                # If the spring constant is zero we avoid to add the restraint!
                if cols_vals[4] == 0.0:
                    continue
                 
                centre                 = cols_vals[5]
                kappa                  = cols_vals[4]*steering_pairs['k_factor']
                 
                if cols_vals[3] == "Harmonic":
                    outf.write(colvars_tail % (name,name,centre,kappa))
         
                if cols_vals[3] == "HarmonicLowerBound":
                    outf.write(colvars_harmonic_lower_bound_tail % (name,name,centre,kappa)) 
            
            i += 1
        poffset += chrlength
            
    outf.flush()
    
    #===========================================================================
    # if isinstance(target_pairs, HiCBasedRestraints):
    #     for copy_nbr in copies:
    #         i = 1
    #         for line in rand_lines:
    #             cols_vals = line
    #                             
    #             if cols_vals[3] == 0.0:
    #                 continue
    #             
    #             name = str(i)+copy_nbr 
    #             
    #             centre                 = cols_vals[4]
    #             kappa                  = cols_vals[3]
    #             
    #             if cols_vals[2] == "Harmonic":
    #                 outf.write(colvars_tail % (name,name,centre,kappa))
    #     
    #             elif cols_vals[2] == "HarmonicLowerBound":
    #                 outf.write(colvars_harmonic_lower_bound_tail % (name,name,centre,kappa))
    #             
    #                  
    #             
    #             i += 1
    #         poffset += chrlength
    #             
    #     outf.flush()
    #===========================================================================
    
    if 'kappa_vs_genomic_distance' in steering_pairs:   
            
        kappa_values = {}
        with open(kappa_vs_genomic_distance) as kgd:
            for line in kgd:
                line_vals = line.split()
                kappa_values[int(line_vals[0])] = float(line_vals[1])
            
        for seqd in set(seqdists.values()):
            kappa = 0
            if seqd in kappa_values:
                kappa = kappa_values[seqd]*kincrease
            else:
                for kappa_key in sorted(kappa_values, key=int):
                    if int(kappa_key) > seqd:
                        break
                    kappa = kappa_values[kappa_key]*kincrease
            centres=''
            names=''
            for seq_name in seqdists:
                if seqdists[seq_name] == seqd:
                    centres += ' 1.0'
                    names += ' '+seq_name
          
            outf.write(colvars_tail % (str(seqd),names,centres,kappa))
                    
    outf.flush()
    
    outf.close()
        
    
def generate_bond_list(steering_pairs):
                            
    """
    Generates lammps bond commands
    
    :param dict steering_pairs: dictionary containing all the information to write down the
      the input file for the colvars implementation
   

    """

    # Getting the input
    # The target_pairs could be also a list as the one in output of get_HiCbased_restraintsXXX
    target_pairs                 = steering_pairs['colvar_input'] 
    if 'kappa_vs_genomic_distance' in steering_pairs:
        kappa_vs_genomic_distance    = steering_pairs['kappa_vs_genomic_distance']
    if 'chrlength' in steering_pairs:
        chrlength                    = steering_pairs['chrlength']
    else:
        chrlength                    = 0
    if 'copies' in steering_pairs:
        copies                       = steering_pairs['copies']
    else:
        copies                       = ['A']
    kbin                         = 10000000
    binsize                      = steering_pairs['binsize']
    if 'percentage_enforced_contacts' in steering_pairs:
        percentage_enforced_contacts = steering_pairs['perc_enfor_contacts']
    else:
        percentage_enforced_contacts = 100

    # Here we extract from all the restraints only 
    # a random sub-sample of percentage_enforced_contacts/100*totcolvars
    rand_lines = []
    i=0
    j=0
    if isinstance(target_pairs, str):    
        totcolvars = linecount(target_pairs)
        ncolvars = int(totcolvars*(float(percentage_enforced_contacts)/100))
        
        #print "Number of enforced contacts = %i over %i" % (ncolvars,totcolvars)
        rand_positions = sample(list(range(totcolvars)), ncolvars)
        rand_positions = sorted(rand_positions)
    
        tfp = open(target_pairs)
        with open(target_pairs) as f:
            for line in f:
                line = line.strip()
                if j >= ncolvars:
                    break
                if line.startswith('#'):
                    continue
             
                cols_vals = line.split()
                # Avoid to enforce restraints between the same bin
                if cols_vals[1] == cols_vals[2]:
                    continue
            
                if i == rand_positions[j]:
                    rand_lines.append(line)
                    j += 1
                i += 1
        tfp.close()
    elif isinstance(target_pairs, HiCBasedRestraints):
        
        rand_lines = target_pairs.get_hicbased_restraints()
        totcolvars = len(rand_lines)
        ncolvars = int(totcolvars*(float(percentage_enforced_contacts)/100))
        
        #print "Number of enforced contacts = %i over %i" % (ncolvars,totcolvars)
        rand_positions = sample(list(range(totcolvars)), ncolvars)
        rand_positions = sorted(rand_positions)
        
        
    else:
        print "Unknown target_pairs"
        return    
    
        
    
    #print rand_lines

    seqdists = {}
    poffset=0
    outf = []  #### a list
    for copy_nbr in copies:
        i = 1
        for line in rand_lines:
            if isinstance(target_pairs, str):   
                cols_vals = line.split()
            else:
                cols_vals = ['chr'] + line
                
            #print cols_vals
            
            if isinstance(target_pairs, HiCBasedRestraints) and cols_vals[3] != "Harmonic" and cols_vals[3] != "HarmonicLowerBound":
                continue
            
            part1_start = int(cols_vals[1])*binsize
            part1_end = (int(cols_vals[1])+1)*binsize
            #print part1_start, part1_end

            part2_start = int(cols_vals[2])*binsize
            part2_end = (int(cols_vals[2])+1)*binsize
            #print part2_start, part2_end

            name = str(i)+copy_nbr  
            seqdist = abs(part1_start-part2_start)
            #print seqdist

            region1 = cols_vals[0] + '_' + str(part1_start) + '_' + str(part1_end)
            region2 = cols_vals[0] + '_' + str(part2_start) + '_' + str(part2_end)

            particle1 = int(cols_vals[1]) + 1 + poffset
            particle2 = int(cols_vals[2]) + 1 + poffset

            seqdists[name] = seqdist

            
            if isinstance(target_pairs, HiCBasedRestraints):
                # If the spring constant is zero we avoid to add the restraint!
                if cols_vals[4] == 0.0:
                    continue
                 
                centre                 = cols_vals[5]
                kappa                  = cols_vals[4]*steering_pairs['k_factor']
                 
                bonType = None
                if cols_vals[3] == "Harmonic":
                    bonType = 'bond'
                elif cols_vals[3] == "HarmonicLowerBound":
                    bonType = 'lbond'

                if bonType:
                    outf.append('fix %s all restrain %s %d %d %f %f %f %f' %(
                        name, bonType, particle1, particle2, 0, kappa, 
                        centre, centre))

            
            i += 1
        poffset += chrlength

    return outf

##########

def generate_time_dependent_colvars_list(steering_pairs,
                                         outfile,
                                         colvar_dump_freq,
                                         colvars_header='# collective variable: monitor distances\n\ncolvarsTrajFrequency %i # output every %i steps\ncolvarsRestartFrequency 1000000\n',
                                         colvars_template='''

colvar {
  name %s
  # %s %s %i
  width 1.0
  distance {
      group1 {
        atomNumbers %i
      }
      group2 {
        atomNumbers %i
      }
  }
}''',
                                         colvars_harmonic_tail = '''

harmonic {
  name h_pot_%s
  colvars %s
  forceConstant %f       
  targetForceConstant %f 
  centers %s             
  targetCenters %s       
  targetNumSteps %s
  outputEnergy   yes
}\n''',
                                         colvars_harmonic_lower_bound_tail = '''
harmonicBound {
  name hlb_pot_%s
  colvars %s
  forceConstant %f
  targetForceConstant %f
  centers %f
  targetCenters %f
  targetNumSteps %s
  outputEnergy   yes
}\n'''
                            ):


    """
    harmonicWalls {
    name hlb_pot_%s
    colvars %s
    forceConstant %f       # This is the force constant at time_point
    targetForceConstant %f # This is the force constant at time_point+1
    centers %f             # This is the equilibrium distance at time_point+1
    targetCenters %f      # This is the equilibrium distance at time_point+1
    targetNumSteps %d      # This is the number of timesteps between time_point and time_point+1
    outputEnergy   yes
    }\n''',


    colvars_harmonic_lower_bound_tail = '''
    
    harmonicBound {
    name hlb_pot_%s
    colvars %s
    forceConstant %f       # This is the force constant at time_point
    targetForceConstant %f # This is the force constant at time_point+1
    centers %f             # This is the equilibrium distance at time_point+1
    targetCenters %f      # This is the equilibrium distance at time_point+1
    targetNumSteps %d      # This is the number of timesteps between time_point and time_point+1
    outputEnergy   yes
    }\n''',
    
    Generates lammps colvars file http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf
    
    :param dict steering_pairs: dictionary containing all the information to write down the
      the input file for the colvars implementation
    :param exisiting_template colvars_header: header template for colvars file.
    :param exisiting_template colvars_template: contact template for colvars file.
    :param exisiting_template colvars_tail: tail template for colvars file.

    """

    #restraints[pair] = [time_dependent_restraints[time_point+1][pair][0],     # Restraint type -> Is the one at time point time_point+1
    #time_dependent_restraints[time_point][pair][1]*10.,                       # Initial spring constant 
    #time_dependent_restraints[time_point+1][pair][1]*10.,                     # Final spring constant 
    #time_dependent_restraints[time_point][pair][2],                           # Initial equilibrium distance 
    #time_dependent_restraints[time_point+1][pair][2],                         # Final equilibrium distance 
    #int(time_dependent_steering_pairs['timesteps_per_k_change'][time_point])] # Number of timesteps for the gradual change

    outf = open(outfile,'w')
    #tfreq=10000
    #for pair in steering_pairs:
    #    if len(steering_pairs[pair][5]) >= 1:
    #        tfreq = int(steering_pairs[pair][5][0]/100)
    #        break
    
    tfreq = colvar_dump_freq
    outf.write(colvars_header % (tfreq, tfreq))
    # Defining the particle pairs
    for pair in steering_pairs:

        #print steering_pairs[pair]
        sys.stdout.flush()
        for i in xrange(len(steering_pairs[pair][0])):
            name    = "%s_%s_%s" % (i, int(pair[0])+1, int(pair[1])+1)
            seqdist = abs(int(pair[1])-int(pair[0])) 
            region1 = "particle_%s" % (int(pair[0])+1)
            region2 = "particle_%s" % (int(pair[1])+1)
            
            outf.write(colvars_template % (name,region1,region2,seqdist,int(pair[0])+1,int(pair[1])+1))

            restraint_type         = steering_pairs[pair][0][i]
            kappa_start            = steering_pairs[pair][1][i]
            kappa_stop             = steering_pairs[pair][2][i]
            centre_start           = steering_pairs[pair][3][i]
            centre_stop            = steering_pairs[pair][4][i]
            timesteps_per_k_change = steering_pairs[pair][5][i]     

            if restraint_type == "Harmonic":
                outf.write(colvars_harmonic_tail             % (name,name,kappa_start,kappa_stop,centre_start,centre_stop,timesteps_per_k_change))
                
            if restraint_type == "HarmonicLowerBound":
                outf.write(colvars_harmonic_lower_bound_tail % (name,name,kappa_start,kappa_stop,centre_start,centre_stop,timesteps_per_k_change))



                    
    outf.flush()
    
    outf.close()  

##########

def get_time_dependent_colvars_list(time_dependent_steering_pairs):
                            
    """
    Generates lammps colvars file http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf
    
    :param dict time_dependent_steering_pairs: dictionary containing all the information to write down the
      the input file for the colvars implementation
    """

    # Getting the input
    # XXXThe target_pairs_file could be also a list as the one in output of get_HiCbased_restraintsXXX
    target_pairs            = time_dependent_steering_pairs['colvar_input'] 
    outfile                      = time_dependent_steering_pairs['colvar_output']
    if 'chrlength' in time_dependent_steering_pairs:
        chrlength                    = time_dependent_steering_pairs['chrlength']
    binsize                      = time_dependent_steering_pairs['binsize']
    if 'percentage_enforced_contacts' in time_dependent_steering_pairs:
        percentage_enforced_contacts = time_dependent_steering_pairs['perc_enfor_contacts']
    else:
        percentage_enforced_contacts = 100

    # HiCbasedRestraints is a list of restraints returned by this function.
    # Each entry of the list is a list of 5 elements describing the details of the restraint:
    # 0 - particle_i
    # 1 - particle_j
    # 2 - type_of_restraint = Harmonic or HarmonicLowerBound or HarmonicUpperBound
    # 3 - the kforce of the restraint
    # 4 - the equilibrium (or maximum or minimum respectively) distance associated to the restraint
   
    # Here we extract from all the restraints only a random sub-sample
    # of percentage_enforced_contacts/100*totcolvars
    rand_lines = []
    i=0
    j=0
    if isinstance(target_pairs, str):    
        time_dependent_restraints = {}
        totcolvars = linecount(target_pairs)
        ncolvars = int(totcolvars*(float(percentage_enforced_contacts)/100))
        
        #print "Number of enforced contacts = %i over %i" % (ncolvars,totcolvars)
        rand_positions = sample(list(range(totcolvars)), ncolvars)
        rand_positions = sorted(rand_positions)
        
        with open(target_pairs) as f:
            for line in f:
                line = line.strip()
                if j >= ncolvars:
                    break
                if line.startswith('#') or line == "":
                    continue
                
                # Line format: timepoint,particle1,particle2,restraint_type,kforce,distance
                cols_vals = line.split()
                
                if int(cols_vals[1]) < int(cols_vals[2]):                
                    pair = (int(cols_vals[1]), int(cols_vals[2]))
                else:
                    pair = (int(cols_vals[2]), int(cols_vals[1]))
            
                try:
                    if pair    in time_dependent_restraints[int(cols_vals[0])]:
                        print "WARNING: Check your restraint list! pair %s is repeated in time point %s!" % (pair, int(cols_vals[0]))
                    # List content: restraint_type,kforce,distance
                    time_dependent_restraints[int(cols_vals[0])][pair] = [cols_vals[3],
                                                                          float(cols_vals[4]),
                                                                          float(cols_vals[5])]
                except:
                    time_dependent_restraints[int(cols_vals[0])] = {}
                    # List content: restraint_type,kforce,distance
                    time_dependent_restraints[int(cols_vals[0])][pair] = [cols_vals[3],
                                                                          float(cols_vals[4]),
                                                                          float(cols_vals[5])]
                if i == rand_positions[j]:
                    rand_lines.append(line)
                    j += 1
                i += 1
    elif isinstance(target_pairs, list):
        time_dependent_restraints = dict((i,{}) for i in xrange(len(target_pairs)))
        for time_point, HiCR in enumerate(target_pairs):
            rand_lines = HiCR.get_hicbased_restraints()
            totcolvars = len(rand_lines)
            ncolvars = int(totcolvars*(float(percentage_enforced_contacts)/100))
            
            #print "Number of enforced contacts = %i over %i" % (ncolvars,totcolvars)
            rand_positions = sample(list(range(totcolvars)), ncolvars)
            rand_positions = sorted(rand_positions)
            
            for cols_vals in rand_lines:
                
                if cols_vals[2] != "Harmonic" and cols_vals[2] != "HarmonicLowerBound":
                    continue
                if int(cols_vals[0]) < int(cols_vals[1]):                
                    pair = (int(cols_vals[0]), int(cols_vals[1]))
                else:
                    pair = (int(cols_vals[1]), int(cols_vals[0]))
            
                if pair in time_dependent_restraints[time_point]:
                    print "WARNING: Check your restraint list! pair %s is repeated in time point %s!" % (pair, time_point)
                # List content: restraint_type,kforce,distance
                time_dependent_restraints[time_point][pair] = [cols_vals[2],
                                                                      float(cols_vals[3]),
                                                                      float(cols_vals[4])]
        
    else:
        print "Unknown target_pairs"
        return 

#     for time_point in sorted(time_dependent_restraints.keys()):
#         for pair in time_dependent_restraints[time_point]:
#             print "#Time_dependent_restraints", time_point,pair, time_dependent_restraints[time_point][pair]
    return time_dependent_restraints

### TODO Add the option to add also spheres of different radii (e.g. to simulate nucleoli)
########## Part to generate the initial conformation ##########
def generate_chromosome_random_walks_conformation ( chromosome_particle_numbers ,
                                                    confining_environment=['sphere',100.] ,
                                                    particle_radius=0.5 ,
                                                    seed_of_the_random_number_generator=1 ,
                                                    number_of_conformations=1,
                                                    outfile="Initial_random_walk_conformation.dat",
                                                    pbc=False,
                                                    center=True):
    """
    Generates lammps initial conformation file by random walks
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param ['sphere',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
            ['sphere',radius]
            ['ellipsoid',x-semiaxes, y-semiaxes, z-semiaxes]
            ['cylinder', basal radius, height]
    :param 0.5 particle_radius: Radius of each particle.
    :param 1 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)

    for cnt in xrange(number_of_conformations):

        final_random_walks = generate_random_walks(chromosome_particle_numbers,
                                                   particle_radius,
                                                   confining_environment,
                                                   pbc=pbc,
                                                   center=center)

        # Writing the final_random_walks conformation
        #print "Succesfully generated conformation number %d\n" % (cnt+1)
        write_initial_conformation_file(final_random_walks,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        out_file=outfile)

##########
        
def generate_chromosome_rosettes_conformation ( chromosome_particle_numbers ,
                                                fractional_radial_positions=None,
                                                confining_environment=['sphere',100.] ,
                                                rosette_radius=12.0 , particle_radius=0.5 ,
                                                seed_of_the_random_number_generator=1 ,
                                                number_of_conformations=1,
                                                outfile = "Initial_rosette_conformation.dat"):
    """
    Generates lammps initial conformation file by rosettes conformation
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param None fractional_radial_positions: list with fractional radial positions for all the chromosomes.
    :param ['sphere',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
            ['sphere',radius]
            ['ellipsoid',x-semiaxes, y-semiaxes, z-semiaxes]
            ['cylinder', basal radius, height]
    :param 0.5 particle_radius: Radius of each particle.
    :param 1 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)    

    initial_rosettes , rosettes_lengths = generate_rosettes(chromosome_particle_numbers,
                                                            rosette_radius,
                                                            particle_radius)
    print rosettes_lengths

    
    # Constructing the rosettes conformations
    for cnt in xrange(number_of_conformations):

        particle_inside   = 0 # 0 means a particle is outside
        particles_overlap = 0 # 0 means two particles are overlapping
        while particle_inside == 0 or particles_overlap == 0:
            particle_inside   = 1
            particles_overlap = 1
            segments_P1 = []
            segments_P0 = []
            side = 0
            init_rosettes = copy.deepcopy(initial_rosettes)
            
            # Guess of the initial segment conformation:
            # 1 - each rod is placed inside the confining evironment
            # in a random position and with random orientation
            # 2 - possible clashes between generated rods are checked
            if fractional_radial_positions:
                if len(fractional_radial_positions) != len(chromosome_particle_numbers):
                    print "Please provide the desired fractional radial positions for all the chromosomes"
                    sys.exit()
                segments_P1 , segments_P0 = generate_rods_biased_conformation(rosettes_lengths, rosette_radius,
                                                                              confining_environment,
                                                                              fractional_radial_positions,
                                                                              max_number_of_temptative=100000)
            else:
                segments_P1 , segments_P0 = generate_rods_random_conformation(rosettes_lengths, rosette_radius,
                                                                              confining_environment,
                                                                              max_number_of_temptative=100000)

            # Roto-translation of the rosettes according to the segment position and orientation 
            final_rosettes = rosettes_rototranslation(init_rosettes, segments_P1, segments_P0)
            
            # Checking that the beads are all inside the confining environment and are not overlapping
            for rosette_pair in list(combinations(final_rosettes,2)):
                
                for x0,y0,z0 in zip(rosette_pair[0]['x'],rosette_pair[0]['y'],rosette_pair[0]['z']):
                    particle_inside = check_point_inside_the_confining_environment(x0, y0, z0,
                                                                                   particle_radius,
                                                                                   confining_environment)
                    
                    if particle_inside == 0: # 0 means that the particle is outside -> PROBLEM!!!
                        print "Particle",x0,y0,z0,"is out of the confining environment\n"
                        break

                    for x1,y1,z1 in zip(rosette_pair[1]['x'],rosette_pair[1]['y'],rosette_pair[1]['z']):
                        particles_overlap = check_particles_overlap(x0,y0,z0,x1,y1,z1,particle_radius)

                        if particles_overlap == 0: # 0 means that the particles are overlapping -> PROBLEM!!!
                            print "Particle",x0,y0,z0,"and",x1,y1,z1,"overlap\n"         
                            break
                    if particle_inside == 0 or particles_overlap == 0:
                        break
                if particle_inside == 0 or particles_overlap == 0:
                    break
            
        # Writing the final_rosettes conformation
        print "Succesfully generated conformation number %d\n" % (cnt+1)
        write_initial_conformation_file(final_rosettes,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        out_file=outfile)

##########
        
def generate_chromosome_rosettes_conformation_with_pbc ( chromosome_particle_numbers ,
                                                         fractional_radial_positions=None,
                                                         confining_environment=['cube',100.] ,
                                                         rosette_radius=12.0 , particle_radius=0.5 ,
                                                         seed_of_the_random_number_generator=1 ,
                                                         number_of_conformations=1,
                                                         outfile = "Initial_rosette_conformation_with_pbc.dat"):
    """
    Generates lammps initial conformation file by rosettes conformation
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param None fractional_radial_positions: list with fractional radial positions for all the chromosomes.
    :param ['cube',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
    :param 0.5 particle_radius: Radius of each particle.
    :param 1 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)    

    initial_rosettes , rosettes_lengths = generate_rosettes(chromosome_particle_numbers,
                                                            rosette_radius,
                                                            particle_radius)
    print rosettes_lengths

    
    # Constructing the rosettes conformations
    for cnt in xrange(number_of_conformations):

        particles_overlap = 0 # 0 means two particles are overlapping
        while particles_overlap == 0:
            particles_overlap = 1
            segments_P1 = []
            segments_P0 = []
            side = 0
            init_rosettes = copy.deepcopy(initial_rosettes)
            
            # Guess of the initial segment conformation:
            # 1 - each rod is placed in a random position and with random orientation
            # 2 - possible clashes between generated rods are checked taking into account pbc
            segments_P1 , segments_P0 = generate_rods_random_conformation_with_pbc (
                rosettes_lengths, 
                rosette_radius,
                confining_environment,
                max_number_of_temptative=100000)

            # Roto-translation of the rosettes according to the segment position and orientation 
            final_rosettes = rosettes_rototranslation(init_rosettes, segments_P1, segments_P0)
            
            # Checking that the beads once folded inside the simulation box (for pbc) are not overlapping
            folded_rosettes = copy.deepcopy(final_rosettes)
            for r in xrange(len(folded_rosettes)):
                particle = 0
                for x, y, z in zip(folded_rosettes[r]['x'],folded_rosettes[r]['y'],folded_rosettes[r]['z']):
                    #inside_1 = check_point_inside_the_confining_environment(x, y, z,
                    #                                                        particle_radius,
                    #                                                        confining_environment)
                    #if inside_1 == 0:
                    #    print inside_1, r, particle, x, y, z
                    
                    while x >  (confining_environment[1]*0.5):
                        x -=  confining_environment[1]
                    while x < -(confining_environment[1]*0.5):
                        x +=  confining_environment[1]
                            
                    while y >  (confining_environment[1]*0.5):
                        y -=  confining_environment[1]
                    while y < -(confining_environment[1]*0.5):
                        y +=  confining_environment[1]
                                    
                    while z >  (confining_environment[1]*0.5):
                        z -=  confining_environment[1]
                    while z < -(confining_environment[1]*0.5):
                        z +=  confining_environment[1]

                    #inside_2 = check_point_inside_the_confining_environment(x, y, z,
                    #                                                      particle_radius,
                    #                                                      confining_environment)
                    #if inside_2 == 1 and inside_1 == 0:
                    #    print inside_2, r, particle, x, y, z
                    folded_rosettes[r]['x'][particle] = x
                    folded_rosettes[r]['y'][particle] = y
                    folded_rosettes[r]['z'][particle] = z
                    particle += 1

            for rosette_pair in list(combinations(folded_rosettes,2)):
                
                for x0,y0,z0 in zip(rosette_pair[0]['x'],rosette_pair[0]['y'],rosette_pair[0]['z']):
                    for x1,y1,z1 in zip(rosette_pair[1]['x'],rosette_pair[1]['y'],rosette_pair[1]['z']):

                        particles_overlap = check_particles_overlap(x0,y0,z0,x1,y1,z1,particle_radius)

                        if particles_overlap == 0: # 0 means that the particles are overlapping -> PROBLEM!!!
                            print "Particle",x0,y0,z0,"and",x1,y1,z1,"overlap\n"         
                            break
                    if particles_overlap == 0:
                        break
                if particles_overlap == 0:
                    break
            
        # Writing the final_rosettes conformation
        print "Succesfully generated conformation number %d\n" % (cnt+1)
        write_initial_conformation_file(final_rosettes,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        out_file=outfile)

##########

def generate_rosettes(chromosome_particle_numbers, rosette_radius, particle_radius):
    # Genaration of the rosettes
    # XXXA. Rosa publicationXXX
    # List to contain the rosettes and the rosettes lengths
    rosettes = []
    rosettes_lengths  = []

    for number_of_particles in chromosome_particle_numbers:
        
        # Variable to build the chain
        phi = 0.0

        # Dictory of lists to contain the rosette
        rosette      = {}
        rosette['x'] = []
        rosette['y'] = []
        rosette['z'] = []        
        
        # Position of the first particle (x_0, 0.0, 0.0)
        rosette['x'].append(rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * cos(phi))
        rosette['y'].append(rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * sin(phi))
        rosette['z'].append(phi / (2.0 * pi))
        #print "First bead is in position: %f %f %f" % (rosette['x'][0], rosette['y'][0], rosette['z'][0])

        # Building the chain: The rosette is growing along the positive z-axes
        for particle in xrange(1,number_of_particles):

            distance = 0.0
            while distance < (particle_radius*2.0): 
                phi   = phi + 0.001
                x_tmp = rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * cos(phi)
                y_tmp = rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * sin(phi)
                z_tmp = phi / (2.0 * pi)     
                distance  = sqrt((x_tmp - rosette['x'][-1])*(x_tmp - rosette['x'][-1]) +
                                 (y_tmp - rosette['y'][-1])*(y_tmp - rosette['y'][-1]) +
                                 (z_tmp - rosette['z'][-1])*(z_tmp - rosette['z'][-1]))

            rosette['x'].append(x_tmp)
            rosette['y'].append(y_tmp)
            rosette['z'].append(z_tmp)
            if distance > ((particle_radius*2.0)*1.2):
                print "%f %d %d %d" % (distance, particle-1, particle)
            
        rosettes.append(rosette)
        rosettes_lengths.append(rosette['z'][-1]-rosette['z'][0])
        
    return rosettes , rosettes_lengths

##########

def generate_rods_biased_conformation(rosettes_lengths, rosette_radius,
                                      confining_environment,
                                      fractional_radial_positions,
                                      max_number_of_temptative=100000):
    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []

    if confining_environment[0] != 'sphere':
        print "ERROR: Biased chromosome positioning is currently implemented"
        print "only for spherical confinement. If you need other shapes, please"
        print "contact the developers"
    
    for length , target_radial_position in zip(rosettes_lengths,fractional_radial_positions):
        tentative            = 0
        clashes              = 0 # 0 means that there is an clash -> PROBLEM
        best_radial_position = 1.0
        best_radial_distance = 1.0
        best_segment_P0      = []
        best_segment_P1      = []
        
        # Positioning the rods
        while tentative < 100000 and best_radial_distance > 0.00005:                

            print "Length = %f" % length

            print "Trying to position terminus 0"
            segment_P0_tmp = []
            segment_P0_tmp = draw_point_inside_the_confining_environment(confining_environment,
                                                                         rosette_radius)
            print "Successfully positioned terminus 0: %f %f %f" % (segment_P0_tmp[0], segment_P0_tmp[1], segment_P0_tmp[2])
            
            print "Trying to position terminus 1"
            segment_P1_tmp = []                            
            segment_P1_tmp = draw_second_extreme_of_a_segment_inside_the_confining_environment(segment_P0_tmp[0],
                                                                                               segment_P0_tmp[1],
                                                                                               segment_P0_tmp[2],
                                                                                               length,
                                                                                               rosette_radius,
                                                                                               confining_environment)
            print "Successfully positioned terminus 1: %f %f %f" % (segment_P1_tmp[0], segment_P1_tmp[1], segment_P1_tmp[2])

            # Check clashes with the previously positioned rods
            clashes = 1
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes(segment_P1,
                                                 segment_P0,
                                                 segment_P1_tmp,
                                                 segment_P0_tmp,
                                                 rosette_radius)
                if clashes == 0:
                    break                

            if clashes == 1:
                # Check whether the midpoint of the segment is close to the target radial position
                segment_midpoint = []
                segment_midpoint.append((segment_P0_tmp[0] + segment_P1_tmp[0])*0.5)
                segment_midpoint.append((segment_P0_tmp[1] + segment_P1_tmp[1])*0.5)
                segment_midpoint.append((segment_P0_tmp[2] + segment_P1_tmp[2])*0.5)

                radial_position = sqrt( ( segment_midpoint[0] * segment_midpoint[0] +
                                          segment_midpoint[1] * segment_midpoint[1] +
                                          segment_midpoint[2] * segment_midpoint[2] ) /
                                        (confining_environment[1]*confining_environment[1]))

                radial_distance = fabs(radial_position-target_radial_position)

                print radial_position , target_radial_position , radial_distance , best_radial_distance , tentative
                
                # If the midpoint of the segment is closer to the target radial position than the
                # previous guesses. Store the points as the best guesses!
                if radial_distance < best_radial_distance:                    
                    best_radial_distance = radial_distance
                    best_radial_position = radial_position
                    best_tentative       = tentative+1 # The variable tentative starts from 0

                    best_segment_P0 = []
                    best_segment_P1 = []
                    for component_P0 , component_P1 in zip(segment_P0_tmp,segment_P1_tmp):
                        best_segment_P0.append(component_P0)
                        best_segment_P1.append(component_P1)
                    
                tentative = tentative + 1
                
        if best_segment_P0 == []:
            print "Valid placement not found for chromosome rosette after %d tentatives" % tentative
            sys.exit()

        print "Successfully positioned chromosome of length %lf at tentative %d of %d tentatives" % (length, best_tentative, tentative)        
        segments_P0.append(best_segment_P0)
        segments_P1.append(best_segment_P1)

    print "Successfully generated rod conformation!"
    return segments_P1 , segments_P0
    
##########

def generate_rods_random_conformation(rosettes_lengths, rosette_radius,
                                      confining_environment,
                                      max_number_of_temptative=100000):
    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []
    
    for length in rosettes_lengths:
        tentative = 0
        clashes   = 0
        # Random positioning of the rods
        while tentative < 100000 and clashes == 0:                

            tentative += 1
            clashes    = 1
            #print "Length = %f" % length

            print "Trying to position terminus 0"
            #pick uniformly within the confining environment using the rejection method 
            first_point = []
            first_point = draw_point_inside_the_confining_environment(confining_environment,
                                                                      rosette_radius)

            print "Successfully positioned terminus 0: %f %f %f" % (first_point[0], first_point[1], first_point[2])
            
            print "Trying to position terminus 1"
            #pick from P0 another point one the sphere of radius length inside the confining environment
            last_point = []
            last_point = draw_second_extreme_of_a_segment_inside_the_confining_environment(first_point[0],
                                                                                           first_point[1],
                                                                                           first_point[2],
                                                                                           length,
                                                                                           rosette_radius,
                                                                                           confining_environment)
            
            print "Successfully positioned terminus 1: %f %f %f" % (last_point[0], last_point[1], last_point[2])
                
            # Check clashes with the previously positioned rods
            clashes = 1 
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes(segment_P1,
                                                 segment_P0,
                                                 last_point,
                                                 first_point,
                                                 rosette_radius)
                if clashes == 0:
                    break                

            #print clashes
        print "Successfully positioned chromosome of length %lf at tentative %d\n" % (length, tentative)        
        segments_P1.append(last_point)
        segments_P0.append(first_point)            

    print "Successfully generated rod conformation!"
    return segments_P1 , segments_P0

##########

def generate_rods_random_conformation_with_pbc(rosettes_lengths, rosette_radius,
                                               confining_environment,
                                               max_number_of_temptative=100000):

    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []
    
    for length in rosettes_lengths:
        tentative = 0
        clashes   = 0
        # Random positioning of the rods
        while tentative < 100000 and clashes == 0:                

            tentative += 1
            clashes    = 1
            #print "Length = %f" % length

            print "Trying to position terminus 0"
            #pick uniformly within the confining environment using the rejection method 
            first_point = []
            first_point = draw_point_inside_the_confining_environment(confining_environment,
                                                                      rosette_radius)

            print "Successfully positioned terminus 0: %f %f %f" % (first_point[0], first_point[1], first_point[2])
            
            print "Trying to position terminus 1"
            #pick from P0 another point one the sphere of radius length inside the confining environment
            last_point = []
            last_point = draw_second_extreme_of_a_segment(first_point[0],
                                                          first_point[1],
                                                          first_point[2],
                                                          length,
                                                          rosette_radius)            
            
            print last_point
            # Check clashes with the previously positioned rods
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes_with_pbc(segment_P1,
                                                          segment_P0,
                                                          last_point,
                                                          first_point,
                                                          rosette_radius,
                                                          confining_environment)
                if clashes == 0:
                    break                

            #print clashes
        print "Successfully positioned chromosome of length %lf at tentative %d\n" % (length, tentative)        
        segments_P1.append(last_point)
        segments_P0.append(first_point)            

    print "Successfully generated rod conformation!"
    return segments_P1 , segments_P0

##########

def generate_random_walks(chromosome_particle_numbers,
                          particle_radius,
                          confining_environment,
                          center,
                          pbc=False):
    # Construction of the random walks initial conformation 
    random_walks = []
    
    for number_of_particles in chromosome_particle_numbers:
        #print "Trying to position random walk"
        random_walk      = {}
        random_walk['x'] = []
        random_walk['y'] = []
        random_walk['z'] = []        


        #print "Positioning first particle"            
        particle_overlap = 0
        while particle_overlap == 0:
            particle_overlap = 1
            first_particle = []
            first_particle = draw_point_inside_the_confining_environment(confining_environment,
                                                                         particle_radius)

            # Check if the particle is overlapping with any other particle in the system
            for rand_walk in random_walks:
                if pbc:
                    particle_overlap = check_particle_vs_all_overlap(first_particle[0],
                                                                     first_particle[1],
                                                                     first_particle[2],
                                                                     rand_walk,
                                                                     2.0*particle_radius)
                else:
                    particle_overlap = check_particle_vs_all_overlap(first_particle[0],
                                                                     first_particle[1],
                                                                     first_particle[2],
                                                                     rand_walk,
                                                                     2.0*particle_radius)
                    
                if particle_overlap == 0:
                    break

        random_walk['x'].append(first_particle[0])        
        random_walk['y'].append(first_particle[1])
        random_walk['z'].append(first_particle[2])

        for particle in xrange(1,number_of_particles):
            #print "Positioning particle %d" % (particle+1)
            particle_overlap = 0 # 0 means that there is an overlap -> PROBLEM
            overlapCounter = -1
            maxIter = 1000
            while particle_overlap == 0:
                overlapCounter += 1
                if overlapCounter > maxIter:
                    # raise error so log file is created to avoid k_seed
                    errorName = 'ERROR: Initial conformation non ending loop'
                    raise InitalConformationError(errorName)
                particle_overlap = 1
                new_particle = []
                if pbc:
                    new_particle = draw_second_extreme_of_a_segment(
                        random_walk['x'][-1],
                        random_walk['y'][-1],
                        random_walk['z'][-1],
                        2.0*particle_radius,
                        2.0*particle_radius)
                else:
                    new_particle = draw_second_extreme_of_a_segment_inside_the_confining_environment(
                        random_walk['x'][-1],
                        random_walk['y'][-1],
                        random_walk['z'][-1],
                        2.0*particle_radius,
                        2.0*particle_radius,
                        confining_environment)

                # Check if the particle is overlapping with any other particle in the system
                for rand_walk in random_walks:
                    particle_overlap = check_particle_vs_all_overlap(new_particle[0],
                                                                     new_particle[1],
                                                                     new_particle[2],
                                                                     rand_walk,
                                                                     2.0*particle_radius)
                    if particle_overlap == 0:
                        break
                if particle_overlap == 0:
                    continue
                
                # The current random walk is not yet in the list above
                particle_overlap = check_particle_vs_all_overlap(new_particle[0],
                                                                 new_particle[1],
                                                                 new_particle[2],
                                                                 random_walk,
                                                                 2.0*particle_radius)
                if particle_overlap == 0:
                    continue
                
            random_walk['x'].append(new_particle[0])        
            random_walk['y'].append(new_particle[1])
            random_walk['z'].append(new_particle[2])
                    
        #print "Successfully positioned random walk of %d particles" % number_of_particles
        random_walks.append(random_walk)

    #print "Successfully generated random walk conformation!"
    if center:
        for random_walk in random_walks:
            x_com, y_com, z_com = (0.0,0.0,0.0)
            cnt = 0
            for (x,y,z) in zip(random_walk['x'],random_walk['y'],random_walk['z']):
                x_com += x
                y_com += y
                z_com += z
                cnt += 1
            x_com, y_com, z_com = (x_com/cnt,y_com/cnt,z_com/cnt)
            print "#Old COM ",x_com,y_com,z_com

            for i in xrange(len(random_walk['x'])):
                random_walk['x'][i] -= x_com
                random_walk['y'][i] -= y_com
                random_walk['z'][i] -= z_com
            
            x_com, y_com, z_com = (0.0,0.0,0.0)
            cnt = 0
            for (x,y,z) in zip(random_walk['x'],random_walk['y'],random_walk['z']):
                x_com += x
                y_com += y
                z_com += z
                cnt += 1
            x_com, y_com, z_com = (x_com/cnt,y_com/cnt,z_com/cnt)
            print "#New COM ",x_com,y_com,z_com
            
    return random_walks

##########

def check_particle_vs_all_overlap(x,y,z,chromosome,overlap_radius):    
    particle_overlap = 1

    for x0, y0, z0 in zip(chromosome['x'],chromosome['y'],chromosome['z']):
        particle_overlap = check_particles_overlap(x0,y0,z0,x,y,z,overlap_radius)
        if particle_overlap == 0:
            return particle_overlap
        
    return particle_overlap
            
##########

def draw_second_extreme_of_a_segment_inside_the_confining_environment(x0, y0, z0, 
                                                                      segment_length, 
                                                                      object_radius, 
                                                                      confining_environment):
    inside = 0
    while inside == 0:
        particle = []
        temp_theta  = arccos(2.0*random()-1.0)
        temp_phi    = 2*pi*random()
        particle.append(x0 + segment_length * cos(temp_phi) * sin(temp_theta))
        particle.append(y0 + segment_length * sin(temp_phi) * sin(temp_theta))
        particle.append(z0 + segment_length * cos(temp_theta))
        # Check if the particle is inside the confining_environment
        inside = check_point_inside_the_confining_environment(particle[0],
                                                              particle[1],
                                                              particle[2],
                                                              object_radius,
                                                              confining_environment)

    return particle

##########

def draw_second_extreme_of_a_segment(x0, y0, z0, 
                                     segment_length, 
                                     object_radius):
    particle = []
    temp_theta  = arccos(2.0*random()-1.0)
    temp_phi    = 2*pi*random()
    particle.append(x0 + segment_length * cos(temp_phi) * sin(temp_theta))
    particle.append(y0 + segment_length * sin(temp_phi) * sin(temp_theta))
    particle.append(z0 + segment_length * cos(temp_theta))
    
    return particle

##########

def draw_point_inside_the_confining_environment(confining_environment, object_radius):
    #pick a point uniformly within the confining environment using the rejection method 

    if confining_environment[0] == 'cube':
        dimension_x = confining_environment[1] * 0.5
        dimension_y = confining_environment[1] * 0.5
        dimension_z = confining_environment[1] * 0.5        
        if len(confining_environment) > 2:
            print "# WARNING: Defined a cubical confining environment with reduntant paramenters."
            print "# Only 2 are needed the identifier and the side"
        
    if confining_environment[0] == 'sphere':
        dimension_x = confining_environment[1]
        dimension_y = confining_environment[1]
        dimension_z = confining_environment[1]
        if len(confining_environment) > 2:
            print "# WARNING: Defined a spherical confining environment with reduntant paramenters."
            print "# Only 2 are needed the identifier and the radius"
        
    if confining_environment[0] == 'ellipsoid':
        if len(confining_environment) < 4:
            print "# ERROR: Defined an ellipsoidal confining environment without the necessary paramenters."
            print "# 4 are needed the identifier, the x-semiaxes, the y-semiaxes, and the z-semiaxes"
            sys.exit()
        dimension_x = confining_environment[1]
        dimension_y = confining_environment[2]
        dimension_z = confining_environment[3]

    if confining_environment[0] == 'cylinder':
        if len(confining_environment) < 3:
            print "# WARNING: Defined a cylindrical confining environment without the necessary paramenters."
            print "# 3 are needed the identifier, the basal radius, and the height"
            sys.exit()
        dimension_x = confining_environment[1]
        dimension_y = confining_environment[1]
        dimension_z = confining_environment[2]
            
    inside = 0
    while inside == 0:
        particle = []
        particle.append((2.0*random()-1.0)*(dimension_x - object_radius))
        particle.append((2.0*random()-1.0)*(dimension_y - object_radius))
        particle.append((2.0*random()-1.0)*(dimension_z - object_radius))
        # Check if the particle is inside the confining_environment
        inside = check_point_inside_the_confining_environment(particle[0],
                                                              particle[1],
                                                              particle[2],
                                                              object_radius,
                                                              confining_environment)                        
    
    return particle
    
##########        

def check_point_inside_the_confining_environment(Px, Py, Pz,
                                                 object_radius,
                                                 confining_environment):
    # The shapes are all centered in the origin
    # - sphere    : radius r
    # - cube      : side
    # - cylinder  : basal radius , height
    # - ellipsoid : semi-axes a , b , c ;

    if confining_environment[0] == 'sphere':
        radius = confining_environment[1] - object_radius
        if ((Px*Px)/(radius*radius) + (Py*Py)/(radius*radius) + (Pz*Pz)/(radius*radius)) < 1.0 : return 1

    if confining_environment[0] == 'ellipsoid':
        a = confining_environment[1] - object_radius
        b = confining_environment[2] - object_radius
        c = confining_environment[3] - object_radius
        if ((Px*Px)/(a*a) + (Py*Py)/(b*b) + (Pz*Pz)/(c*c)) < 1.0 : return 1

    if confining_environment[0] == 'cube':
        hside = confining_environment[1] * 0.5 - object_radius
        if (((Px*Px)/(hside*hside)) < 1.0) and (((Py*Py)/(hside*hside)) < 1.0) and (((Pz*Pz)/(hside*hside)) < 1.0) : return 1

    if confining_environment[0] == 'cylinder':
        radius      = confining_environment[1]     - object_radius
        half_height = confining_environment[2]*0.5 - object_radius
        if (((Px*Px)/(radius*radius) + (Py*Py)/(radius*radius)) < 1.0) and (((Pz*Pz)/(half_height*half_height)) < 1.0): return 1
            
    return 0

##########

def check_segments_clashes(s1_P1, s1_P0, s2_P1, s2_P0, rosette_radius):

    # Check steric clashes without periodic boundary conditions
    if distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0) < 2.0*rosette_radius:
        # print "Clash between segments",s1_P1,s1_P0,"and",s2_P1_tmp,s2_P0_tmp,"at distance", distance
        return 0

    return 1

##########

def check_segments_clashes_with_pbc(s1_P1, s1_P0, s2_P1, s2_P0, 
                                    rosette_radius,
                                    confining_environment):

    # Check steric clashes with periodic boundary conditions
    if distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0) < 2.0*rosette_radius:
        # print "Clash between segments",s1_P1,s1_P0,"and",s2_P1_tmp,s2_P0_tmp,"at distance", distance
        return 0

    return 1

##########

def distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0):

    # Inspiration: http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm 
    # Copyright 2001, softSurfer (www.softsurfer.com)
    # This code may be freely used and modified for any purpose
    # providing that this copyright notice is included with it.
    # SoftSurfer makes no warranty for this code, and cannot be held
    # liable for any real or imagined damage resulting from its use.
    # Users of this code must verify correctness for their application.

    u  = []
    v  = []
    w  = []
    dP = []

    for c_s1_P1,c_s1_P0,c_s2_P1,c_s2_P0 in zip(s1_P1, s1_P0, s2_P1, s2_P0):        
        u.append(c_s1_P1 - c_s1_P0)
        v.append(c_s2_P1 - c_s2_P0)
        w.append(c_s1_P0 - c_s2_P0)
    
    a  = scalar_product(u, u)
    b  = scalar_product(u, v)
    c  = scalar_product(v, v)
    d  = scalar_product(u, w)
    e  = scalar_product(v, w)

    D  = a*c - b*b
    sD = tD = D
        
    if D < (1.0e-7):
        # Segments almost parallel 
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else:
        # Get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if (sN < 0.0):            
            # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD: # sc > 1 => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0: # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # Recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a        
    
    elif tN > tD: # tc > 1 => the t=1 edge is visible
        tN = tD
        # Recompute sc for this edge
        if (-d + b) < 0.0:
            sN = 0
        elif (-d + b) > a:
            sN = sD;
        else:
            sN = (-d + b)
            sD = a

    # Finally do the division to get sc and tc
    if abs(sN) < (1.0e-7):
        sc = 0.0
    else:
        sc = sN / sD
        
    if abs(tN) < (1.0e-7):
        tc = 0.0
    else:
        tc = tN / tD
     
    # Get the difference of the two closest points
    for i in xrange(len(w)):    
        dP.append(w[i] + ( sc * u[i] ) - ( tc * v[i] )) # = S1(sc) - S2(tc)
    
    return norm(dP)   # return the closest distance

##########

def rosettes_rototranslation(rosettes, segments_P1, segments_P0):

    for i in xrange(len(segments_P1)):
        vector = []
        theta  = []

        for component_P1 , component_P0 in zip(segments_P1[i], segments_P0[i]):
            vector.append(component_P1-component_P0)
            
        # Rotation Angles
        theta.append(atan2(vector[1],vector[2]))
      
        x_temp_2 =  vector[0]
        y_temp_2 =  cos(theta[0]) * vector[1] - sin(theta[0]) * vector[2]                
        z_temp_2 =  sin(theta[0]) * vector[1] + cos(theta[0]) * vector[2]        
        theta.append(atan2(x_temp_2,z_temp_2))
        
        x_temp_1 =  cos(theta[1]) * x_temp_2 - sin(theta[1]) * z_temp_2
        y_temp_1 =  y_temp_2
        z_temp_1 =  sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2
        
        if(z_temp_1 < 0.0):            
            z_temp_1 = -z_temp_1
            theta.append(pi)
        else:
            theta.append(0.0)
        #print x_temp_1 , y_temp_1 , z_temp_1 
        
        # Chromosome roto-translations
        for particle in xrange(len(rosettes[i]['x'])):

            x_temp_2 =   rosettes[i]['x'][particle]
            y_temp_2 =   cos(theta[2]) * rosettes[i]['y'][particle] + sin(theta[2]) * rosettes[i]['z'][particle]
            z_temp_2 = - sin(theta[2]) * rosettes[i]['y'][particle] + cos(theta[2]) * rosettes[i]['z'][particle]
            
            x_temp_1 =   cos(theta[1]) * x_temp_2 + sin(theta[1]) * z_temp_2
            y_temp_1 =   y_temp_2
            z_temp_1 = - sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2

            x =   x_temp_1;
            y =   cos(theta[0]) * y_temp_1 + sin(theta[0]) * z_temp_1;
            z = - sin(theta[0]) * y_temp_1 + cos(theta[0]) * z_temp_1;
            
            # Chromosome translations
            rosettes[i]['x'][particle] = segments_P0[i][0] + x;
            rosettes[i]['y'][particle] = segments_P0[i][1] + y;
            rosettes[i]['z'][particle] = segments_P0[i][2] + z;
    return rosettes

##########

def scalar_product(a, b):

    scalar = 0.0
    for c_a,c_b in zip(a,b):
        scalar = scalar + c_a*c_b 

    return scalar

##########

def norm(a):

    return sqrt(scalar_product(a, a))

##########

def write_initial_conformation_file(chromosomes,
                                    chromosome_particle_numbers,
                                    confining_environment,
                                    out_file="Initial_conformation.dat"):
    # Choosing the appropriate xlo, xhi...etc...depending on the confining environment
    xlim = []
    ylim = []
    zlim = []
    if confining_environment[0] == 'sphere':
        radius = confining_environment[1] + 1.0
        xlim.append(-radius)
        xlim.append(radius)
        ylim.append(-radius)
        ylim.append(radius)
        zlim.append(-radius)
        zlim.append(radius)
        
    if confining_environment[0] == 'ellipsoid':
        a = confining_environment[1] + 1.0
        b = confining_environment[2] + 1.0
        c = confining_environment[3] + 1.0
        xlim.append(-a)
        xlim.append(a)
        ylim.append(-b)
        ylim.append(b)
        zlim.append(-c)
        zlim.append(c)

    if confining_environment[0] == 'cube':
        hside = confining_environment[1] * 0.5 + 1.0
        xlim.append(-hside)
        xlim.append(hside)
        ylim.append(-hside)
        ylim.append(hside)
        zlim.append(-hside)
        zlim.append(hside)
        
    if confining_environment[0] == 'cylinder':
        radius      = confining_environment[1]   + 1.0
        hheight = confining_environment[2] * 0.5 + 1.0
        xlim.append(-radius)
        xlim.append(radius)
        ylim.append(-radius)
        ylim.append(radius)
        zlim.append(-hheight)
        zlim.append(hheight)
    
    fileout = open(out_file,'w')
    angle_types = 1
    bond_types  = 1
    atom_types  = 1
    n_chr=len(chromosomes)
    n_atoms=0
    for n in chromosome_particle_numbers:
        n_atoms+=n    
        
    fileout.write("LAMMPS input data file \n\n")
    fileout.write("%9d atoms\n" % (n_atoms))
    fileout.write("%9d bonds\n" % (n_atoms-n_chr))
    fileout.write("%9d angles\n\n" % (n_atoms-2*n_chr))
    fileout.write("%9s atom types\n" % atom_types)
    fileout.write("%9s bond types\n" % bond_types)
    fileout.write("%9s angle types\n\n" % angle_types)
    fileout.write("%6.3lf    %6.3lf     xlo xhi\n" % (xlim[0], xlim[1]))
    fileout.write("%6.3lf    %6.3lf     ylo yhi\n" % (ylim[0], ylim[1]))
    fileout.write("%6.3lf    %6.3lf     zlo zhi\n" % (zlim[0], zlim[1]))
  
    fileout.write("\n Atoms \n\n")
    particle_number = 1
    for chromosome in chromosomes:
        for x,y,z in zip(chromosome['x'],chromosome['y'],chromosome['z']):          
            fileout.write("%-8d %s %s %7.4lf %7.4lf %7.4lf\n" % (particle_number, "1", "1", x, y, z))
            particle_number += 1
            
    # for(i = 0; i < N_NUCL; i++)
    # {
    #    k++;
    #    fileout.write("%5d %s %s %7.4lf %7.4lf %7.4lf \n", k, "1", "1", P[i][0], P[i][1], P[i][2]);
    # }
  
    fileout.write("\n Bonds \n\n")
    bond_number          = 1
    first_particle_index = 1
    for chromosome in chromosomes:
        for i in xrange(len(chromosome['x'])-1):
            fileout.write("%-4d %s %4d %4d\n" % (bond_number, "1", first_particle_index, first_particle_index+1))
            bond_number          += 1
            first_particle_index += 1
        first_particle_index += 1 # I have to go to the end of the chromosome!
          
    fileout.write("\n Angles \n\n")
    angle_number         = 1
    first_particle_index = 1
    for chromosome in chromosomes:        
        for i in xrange(len(chromosome['x'])-2):
            fileout.write("%-4d %s %5d %5d %5d\n" % (angle_number, "1", first_particle_index, first_particle_index+1, first_particle_index+2))
            angle_number         += 1    
            first_particle_index += 1
        first_particle_index += 2 # I have to go to the end of the chromosome!

    fileout.close()

##########

def distance(x0,y0,z0,x1,y1,z1):
    return sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))

##########

def check_particles_overlap(x0,y0,z0,x1,y1,z1,overlap_radius):
    if distance(x0,y0,z0,x1,y1,z1) < overlap_radius:
        #print "Particle %f %f %f and particle %f %f %f are overlapping\n" % (x0,y0,z0,x1,y1,z1)
        return 0
    return 1

##########

def store_conformation_with_pbc(xc, result, confining_environment):
    # Reconstruct the different molecules and store them separatelly
    ix    , iy    , iz     = (0, 0, 0)
    ix_tmp, iy_tmp, iz_tmp = (0, 0, 0)
    x_tmp , y_tmp , z_tmp  = (0, 0, 0)

    molecule_number = 0 # We start to count from molecule number 0

    particles = []
    particles.append({})
    particles[molecule_number]['x'] = []
    particles[molecule_number]['y'] = []
    particles[molecule_number]['z'] = []

    particle_counts        = []
    particle_counts.append({}) # Initializing the particle counts for the first molecule
    
    max_bond_length        = (1.5*1.5) # This is the default polymer-based bond length

    for i in xrange(0,len(xc),3):
        particle = int(i/3)

        x = xc[i]   + ix * confining_environment[1] 
        y = xc[i+1] + iy * confining_environment[1] 
        z = xc[i+2] + iz * confining_environment[1] 
        
        # A - Check whether the molecule is broken because of pbc
        # or if we are changing molecule
        if particle > 0:             
        
            # Compute the bond_length
            bond_length  = (particles[molecule_number]['x'][-1] - x)* \
                           (particles[molecule_number]['x'][-1] - x)+ \
                           (particles[molecule_number]['y'][-1] - y)* \
                           (particles[molecule_number]['y'][-1] - y)+ \
                           (particles[molecule_number]['z'][-1] - z)* \
                           (particles[molecule_number]['z'][-1] - z)
            
            # Check whether the bond is too long. This could mean:
            # 1 - Same molecule disjoint by pbc
            # 2 - Different molecules
            if bond_length > max_bond_length:
                min_bond_length = bond_length                
                x_tmp , y_tmp , z_tmp = (x, y, z)

                # Check if we are in case 1: the same molecule continues 
                # in a nearby box
                indeces_sets = product([-1, 0, 1],
                                       [-1, 0, 1],
                                       [-1, 0, 1])
                
                for (l, m, n) in indeces_sets:
                    # Avoid to check again the same periodic copy
                    if (l, m, n) == (0, 0, 0):
                        continue

                    # Propose a new particle position
                    x = xc[i]   + (ix + l) * confining_environment[1] 
                    y = xc[i+1] + (iy + m) * confining_environment[1] 
                    z = xc[i+2] + (iz + n) * confining_environment[1] 
                    
                    # Check the new bond length
                    bond_length  = (particles[molecule_number]['x'][-1] - x)* \
                                   (particles[molecule_number]['x'][-1] - x)+ \
                                   (particles[molecule_number]['y'][-1] - y)* \
                                   (particles[molecule_number]['y'][-1] - y)+ \
                                   (particles[molecule_number]['z'][-1] - z)* \
                                   (particles[molecule_number]['z'][-1] - z)
                    
                    # Store the periodic copy with the minimum bond length
                    if bond_length < min_bond_length:
                        #print bond_length
                        x_tmp , y_tmp , z_tmp  = (x   , y   , z   )
                        ix_tmp, iy_tmp, iz_tmp = (ix+l, iy+m, iz+n)
                        min_bond_length = bond_length
                
                # If the minimum bond length is yet too large
                # we are dealing with case 2
                if min_bond_length > 10.:
                    # Start another molecule
                    molecule_number += 1

                    particles.append({})
                    particles[molecule_number]['x'] = []
                    particles[molecule_number]['y'] = []
                    particles[molecule_number]['z'] = []


                    particle_counts.append({}) # Initializing the particle counts for the new molecule

                # If the minimum bond length is sufficiently short
                # we are dealing with case 2
                ix, iy, iz = (ix_tmp, iy_tmp, iz_tmp)
                x , y , z  = (x_tmp , y_tmp , z_tmp)

        # To fullfill point B (see below), we have to count how many
        # particle we have of each molecule for each triplet
        # (ix, iy, iz)                
        try:
            particle_counts[molecule_number][(ix, iy, iz)] += 1.0
        except:
            particle_counts[molecule_number][(ix, iy, iz)] = 0.0
            particle_counts[molecule_number][(ix, iy, iz)] += 1.0
            
        particles[molecule_number]['x'].append(x)
        particles[molecule_number]['y'].append(y)
        particles[molecule_number]['z'].append(z)

    # B - Store in the final arrays each molecule in the periodic copy
    # with more particle in the primary cell (0, 0, 0)
    for molecule in xrange(molecule_number+1):
        max_number = 0
        # Get the periodic box with more particles
        for (l, m, n) in particle_counts[molecule]:
            if particle_counts[molecule][(l, m, n)] > max_number:
                ix, iy, iz = (l, m, n)
                max_number = particle_counts[molecule][(l, m, n)]

        # Translate the molecule to include the largest portion of particles
        # in the (0, 0, 0) image
        for (x, y, z) in zip(particles[molecule]['x'],particles[molecule]['y'],particles[molecule]['z']):
            x = x - ix * confining_environment[1]
            y = y - iy * confining_environment[1]
            z = z - iz * confining_environment[1]

            result['x'].append(x)
            result['y'].append(y)
            result['z'].append(z)


##### Loop extrusion dynamics functions #####
def read_target_loops_input(input_filename, chromosome_length, percentage):
    # Open input file
    fp_input = open(input_filename, "r")

    loops = []
    target_loops = []
    # Get each loop per line and fill the output list of loops
    for line in fp_input.readlines():

        if line.startswith('#'):
            continue

        splitted = line.strip().split()
        loop = []
        loop.append(int(splitted[1]))
        loop.append(int(splitted[2]))

        loops.append(loop)

    #ntarget_loops = int(len(loops)*percentage/100.)    
    ntarget_loops = int(len(loops))    
    shuffle(loops)
    target_loops = loops[0:ntarget_loops]

    return target_loops
        
##########

def draw_loop_extrusion_starting_points(target_loops, chromosome_length):
    initial_loops = []
    # Scroll the target loops and draw a point between each start and stop
    for target_loop in target_loops:

        random_particle =  randint(target_loop[0], target_loop[1])

        loop = []
        loop.append(random_particle)
        loop.append(random_particle+1)

        initial_loops.append(loop)

    return initial_loops

##########

def get_maximum_number_of_extruded_particles(target_loops, initial_loops):
    # The maximum is the maximum sequence distance between a target start/stop particle of a loop
    # and the initial random start/stop of a loop
    maximum = 0

    for target_loop,initial_loop in zip(target_loops,initial_loops):
        #print initial_loop,target_loop
        
        l = abs(target_loop[0]-initial_loop[0])+1
        if l > maximum:
            maximum = l

        l = abs(target_loop[1]-initial_loop[1])+1
        if l > maximum:
            maximum = l

    return maximum

##########

def compute_particles_distance(xc):
    
    particles = []
    distances = {}

    # Getting the coordiantes of the particles
    for i in xrange(0,len(xc),3):
        x = xc[i]  
        y = xc[i+1]
        z = xc[i+2]
        particles.append((x, y, z))

    # Checking whether the restraints are satisfied
    for pair in combinations(xrange(len(particles)), 2):
        dist = distance(particles[pair[0]][0],
                        particles[pair[0]][1],
                        particles[pair[0]][2],
                        particles[pair[1]][0],
                        particles[pair[1]][1],
                        particles[pair[1]][2])
        distances[pair] = dist

    return distances

##########

def compute_the_percentage_of_satysfied_restraints(input_file_name,
                                                   restraints,
                                                   output_file_name,
                                                   time_point,
                                                   timesteps_per_k_change):

    ### Change this function to use a posteriori the out.colvars.traj file similar to the obj funct calculation ###
    infile  = open(input_file_name , "r")
    outfile = open(output_file_name, "w")
    if os.path.getsize(output_file_name) == 0:
        outfile.write("#%s %s %s %s\n" % ("timestep","satisfied", "satisfiedharm", "satisfiedharmLowBound"))

    #restraints[pair] = [time_dependent_restraints[time_point+1][pair][0], # Restraint type -> Is the one at time point time_point+1
    #time_dependent_restraints[time_point][pair][1]*10.,                   # Initial spring constant 
    #time_dependent_restraints[time_point+1][pair][1]*10.,                 # Final spring constant 
    #time_dependent_restraints[time_point][pair][2],                       # Initial equilibrium distance 
    #time_dependent_restraints[time_point+1][pair][2],                     # Final equilibrium distance 
    #int(time_dependent_steering_pairs['timesteps_per_k_change']*0.5)]     # Number of timesteps for the gradual change

    # Write statistics on the restraints
    nharm = 0
    nharmLowBound = 0        
    ntot  = 0
    for pair in restraints:
        for i in xrange(len(restraints[pair][0])):
            if restraints[pair][0][i] == "Harmonic":
                nharm += 1
                ntot  += 1
            if restraints[pair][0][i] == "HarmonicLowerBound":
                nharmLowBound += 1
                ntot  += 1
    outfile.write("#NumOfRestraints = %s , Harmonic = %s , HarmonicLowerBound = %s\n" % (ntot, nharm, nharmLowBound))

    # Memorizing the restraint
    restraints_parameters = {}
    for pair in restraints:
        for i in xrange(len(restraints[pair][0])):
            #E_hlb_pot_p_106_189
            if restraints[pair][0][i] == "Harmonic":
                name  = "E_h_pot_%d_%d_%d" % (i, int(pair[0])+1, int(pair[1])+1)
            if restraints[pair][0][i] == "HarmonicLowerBound":
                name  ="E_hlb_pot_%d_%d_%d" % (i, int(pair[0])+1, int(pair[1])+1)
            restraints_parameters[name] = [restraints[pair][0][i],
                                           restraints[pair][1][i],
                                           restraints[pair][2][i],
                                           restraints[pair][3][i],
                                           restraints[pair][4][i],
                                           restraints[pair][5][i]]
    #print restraints_parameters
    
    # Checking whether the restraints are satisfied
    columns_to_consider = {}
    for line in infile.readlines():
        nsatisfied             = 0.
        nsatisfiedharm         = 0.
        nsatisfiedharmLowBound = 0.
        ntot                   = 0.
        ntotharm               = 0.
        ntotharmLowBound       = 0.

        line = line.strip().split()        
        
        # Checking which columns contain the pairwise distance
        if line[0][0] == "#":            
            for column in xrange(2,len(line)):
                # Get the columns with the distances
                if "_pot_" not in line[column]:
                    columns_to_consider[column-1] = line[column]
                    #print columns_to_consider
        else:
            for column in xrange(1,len(line)):
                if column in columns_to_consider:
                    if column >= len(line):
                        continue
                    dist = float(line[column])
                    
                    # Get which restraints are between the 2 particles
                    for name in ["E_h_pot_"+columns_to_consider[column], "E_hlb_pot_"+columns_to_consider[column]]:
                        if name not in restraints_parameters:
                            #print "Restraint %s not present" % name
                            continue
                        else:
                            pass
                            #print name, restraints_parameters[name] 
                    
                        restrainttype = restraints_parameters[name][0]
                        restraintd0   = float(restraints_parameters[name][3]) + float(line[0])/float(restraints_parameters[name][5])*(float(restraints_parameters[name][4]) - float(restraints_parameters[name][3]))
                        restraintk    = float(restraints_parameters[name][1]) + float(line[0])/float(restraints_parameters[name][5])*(float(restraints_parameters[name][2]) - float(restraints_parameters[name][1]))
                        sqrt_k = sqrt(restraintk)                    
                        limit1 = restraintd0 - 2./sqrt_k
                        limit2 = restraintd0 + 2./sqrt_k

                        if restrainttype == "Harmonic":
                            if dist >= limit1 and dist <= limit2:
                                nsatisfied     += 1.0
                                nsatisfiedharm += 1.0
                                #print "#ESTABLISHED",time_point,name,restraints_parameters[name],limit1,dist,limit2
                            else:
                                pass
                                #print "#NOESTABLISHED",time_point,name,restraints_parameters[name],limit1,dist,limit2
                            ntotharm += 1.0
                        if restrainttype == "HarmonicLowerBound":
                            if dist >= restraintd0:
                                nsatisfied             += 1.0
                                nsatisfiedharmLowBound += 1.0
                                #print "#ESTABLISHED",time_point,name,restraints_parameters[name],dist,restraintd0
                            else:
                                pass
                                #print "#NOESTABLISHED",time_point,name,restraints_parameters[name],dist,restraintd0
                            ntotharmLowBound += 1.0
                        ntot += 1.0
                        #print int(line[0])+(time_point)*timesteps_per_k_change, nsatisfied, ntot, nsatisfiedharm, ntotharm, nsatisfiedharmLowBound, ntotharmLowBound
            if ntotharm         == 0.:
                ntotharm         = 1.0
            if ntotharmLowBound == 0.:
                ntotharmLowBound = 1.0


            outfile.write("%d %lf %lf %lf\n" % (int(line[0])+(time_point)*timesteps_per_k_change, nsatisfied/ntot*100., nsatisfiedharm/ntotharm*100., nsatisfiedharmLowBound/ntotharmLowBound*100.))
    infile.close()
    outfile.close()

##########

def read_objective_function(fname):
    
    obj_func=[]
    fhandler = open(fname)
    line = fhandler.next()
    try:
        while True:
            if line.startswith('#'):
                line = fhandler.next()
                continue
            line = line.strip()
            if len(line) == 0:
                continue
            line_vals = line.split()
            obj_func.append(float(line_vals[1]))
            line = fhandler.next()
    except StopIteration:
        pass
    fhandler.close()        
            
    return obj_func      
##########

def compute_the_objective_function(input_file_name,
                                   output_file_name,
                                   time_point,
                                   timesteps_per_k_change):
    
    
    infile  = open(input_file_name , "r")
    outfile = open(output_file_name, "w")
    if os.path.getsize(output_file_name) == 0:
        outfile.write("#Timestep obj_funct\n")

    columns_to_consider = []

    # Checking which columns contain the energies to sum
    for line in infile.readlines():
        line = line.strip().split()        
        
        # Checking which columns contain the energies to sum
        if line[0][0] == "#":            
            for column in xrange(len(line)):
                if "_pot_" in line[column]:
                    columns_to_consider.append(column-1)
        else:
            obj_funct = 0.0
            for column in columns_to_consider:
                if column < len(line):
                    obj_funct += float(line[column])
            outfile.write("%d %s\n" % (int(line[0])+timesteps_per_k_change*(time_point), obj_funct))

    infile.close()
    outfile.close()
