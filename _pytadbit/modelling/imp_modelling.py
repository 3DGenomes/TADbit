"""
05 Jul 2013


"""

from math            import fabs
from cPickle         import load, dump
from sys             import stdout
from os.path         import exists
import multiprocessing as mu
from scipy           import polyfit

from pytadbit.modelling.IMP_CONFIG       import CONFIG, NROUNDS, STEPS, LSTEPS
from pytadbit.modelling.structuralmodels import StructuralModels
from pytadbit.modelling.impmodel         import IMPmodel
from pytadbit.modelling.restraints       import HiCBasedRestraints

#Local application/library specific imports
import IMP.core
import IMP.algebra
import IMP.display
from IMP.container import ListSingletonContainer
from IMP import Model
from IMP import FloatKey


IMP.set_check_level(IMP.NONE)
IMP.set_log_level(IMP.SILENT)


def generate_3d_models(zscores, resolution, nloci, start=1, n_models=5000,
                       n_keep=1000, close_bins=1, n_cpus=1, keep_all=False,
                       verbose=0, outfile=None, config=None,
                       values=None, experiment=None, coords=None, zeros=None,
                       first=None, container=None, use_HiC=True,
                       use_confining_environment=True, use_excluded_volume=True,
                       single_particle_restraints=None):
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
    :param None coords: a dictionary or a list of dictionaries like:
       ::

         {'crm'  : '19',
          'start': 14637,
          'end'  : 15689}

    :param 5000 n_models: number of models to generate
    :param 1000 n_keep: number of models used in the final analysis (usually
       the top 20% of the generated models). The models are ranked according to
       their objective function value (the lower the better)
    :param False keep_all: whether or not to keep the discarded models (if
       True, models will be stored under StructuralModels.bad_models)
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
    :param True  use_HiC
    :param True  use_confining_environment
    :param True  use_excluded_volume
    :param: None single_particle_restraints: a list containing restraints to single particles.
            Each restraint in the list is itself a list with the following information:
                [bin, [position_x, position_y, position_z], type, kforce, radius]
                bin: bin number of the particle to restraint
                [position_x, position_y, position_z](nm): center of the sphere of the restraint
                type: 'Harmonic', 'HarmonicLowerBound', 'HarmonicUpperBound'
                kforce: weigth of the restraint
                radius (nm): radius of the sphere

    :returns: a StructuralModels object

    """

    # Main config parameters
    global CONFIG
    
    # Setup CONFIG
    if isinstance(config, dict):
        CONFIG.update(config)
    elif config:
        raise Exception('ERROR: "config" must be a dictionary')
    
    CONFIG['resolution'] = resolution
    
    # scale factor
    global SCALE
    SCALE = float(resolution * CONFIG['scale'])
    
    # scale maxdist
    CONFIG['maxdist'] = CONFIG['maxdist'] / SCALE
    
    # Setup and scale CONFIG['container']
    try:
        CONFIG['container'] = {'shape' : container[0],
                               'radius': container[1] / SCALE,
                               'height': container[2] / SCALE,
                               'cforce': container[3]}
    except:
        CONFIG['container'] = {'shape' : None,
                               'radius': None,
                               'height': None,
                               'cforce': None}

    # print "Used",CONFIG,'\n'
    # print "Input",config,'\n'

    # Particles initial radius
    global RADIUS

    RADIUS = 0.5

    global LOCI
    # if z-scores are generated outside TADbit they may not start at zero
    if first == None:
        first = min([int(j) for i in zscores for j in zscores[i]] +
                    [int(i) for i in zscores])
    LOCI  = range(first, nloci + first)

    # random inital number
    global START
    START = start
    # verbose
    global VERBOSE
    VERBOSE = verbose
    #VERBOSE = 3

    HiCRestraints = HiCBasedRestraints(nloci,RADIUS,CONFIG,resolution,zscores,
                 chromosomes=coords, close_bins=close_bins,first=first)

    models, bad_models = multi_process_model_generation(
        n_cpus, n_models, n_keep, keep_all, HiCRestraints,
        use_HiC=use_HiC, use_confining_environment=use_confining_environment,
        use_excluded_volume=use_excluded_volume,
        single_particle_restraints=single_particle_restraints)

    try:
        xpr = experiment
        crm = xpr.crm
        description = {'identifier'        : xpr.identifier,
                       'chromosome'        : coords['crm'] if isinstance(coords,dict) else [c['crm'] for c in coords],
                       'start'             : xpr.resolution * (coords['start'] - 1) if isinstance(coords,dict) else [xpr.resolution * (c['start'] - 1) for c in coords],
                       'end'               : xpr.resolution * coords['end'] if isinstance(coords,dict) else [xpr.resolution *c['end'] for c in coords],
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
        for i, m in enumerate(models.values() + bad_models.values()):
            m['index'] = i
            m['description'] = description
    except AttributeError:  # case we are doing optimization
        description = None
        for i, m in enumerate(models.values() + bad_models.values()):
            m['index'] = i
    if outfile:
        if exists(outfile):
            old_models, old_bad_models = load(open(outfile))
        else:
            old_models, old_bad_models = {}, {}
        models.update(old_models)
        bad_models.update(old_bad_models)
        out = open(outfile, 'w')
        dump((models, bad_models), out)
        out.close()
    else:
        hicrestraints = HiCRestraints._get_restraints()
        hicrestraints = dict((r,(hicrestraints[r][0],hicrestraints[r][1]*SCALE, hicrestraints[r][2])) 
                             for r in hicrestraints)
        return StructuralModels(
            len(LOCI), models, bad_models, resolution, original_data=values,
            zscores=zscores, config=CONFIG, experiment=experiment, zeros=zeros,
            restraints=hicrestraints, description=description)

def multi_process_model_generation(n_cpus, n_models, n_keep, keep_all,HiCRestraints, use_HiC=True,
                                   use_confining_environment=True, use_excluded_volume=True,
                                   single_particle_restraints=None):
    """
    Parallelize the
    :func:`pytadbit.modelling.imp_model.StructuralModels.generate_IMPmodel`.

    :param n_cpus: number of CPUs to use
    :param n_models: number of models to generate
    """

    pool = mu.Pool(n_cpus, maxtasksperchild=1)
    jobs = {}
    for rand_init in xrange(START, n_models + START):
        jobs[rand_init] = pool.apply_async(generate_IMPmodel,
                                           args=(rand_init,HiCRestraints, use_HiC,
                                                 use_confining_environment, use_excluded_volume,
                                                 single_particle_restraints))

    pool.close()
    pool.join()

    results = []
    for rand_init in xrange(START, n_models + START):
        results.append((rand_init, jobs[rand_init].get()))

    models = {}
    bad_models = {}
    for i, (_, m) in enumerate(
        sorted(results, key=lambda x: x[1]['objfun'])[:n_keep]):
        models[i] = m
    if keep_all:
        for i, (_, m) in enumerate(
        sorted(results, key=lambda x: x[1]['objfun'])[n_keep:]):
            bad_models[i+n_keep] = m
    return models, bad_models



def generate_IMPmodel(rand_init, HiCRestraints,use_HiC=True, use_confining_environment=True,
                      use_excluded_volume=True, single_particle_restraints=None):
    """
    Generates one IMP model

    :param rand_init: random number kept as model key, for reproducibility.

    :returns: a model, that is a dictionary with the log of the objective
       function value optimization, and the coordinates of each particles.

    """
    verbose = VERBOSE
    # Set the IMP.random_number_generator.seed to get a different reproducible model
    IMP.random_number_generator.seed(rand_init)
    #print IMP.random_number_generator()

    log_energies = []
    model = {'radius'     : IMP.FloatKey("radius"),
             'model'      : Model(),
             'particles'  : None,
             'restraints' : None} # 2.6.1 compat
    model['particles'] = ListSingletonContainer(model['model'],
        IMP.core.create_xyzr_particles(model['model'], len(LOCI), RADIUS, 100000/SCALE))  # last number is box size

    try:
        model['restraints'] = IMP.RestraintSet(model['model']) # 2.6.1 compat
    except:
        pass

    # OPTIONAL:Set the name of each particle
    for i in range(0, len(LOCI)):
        p = model['particles'].get_particle(i)
        p.set_name(str(LOCI[i]))
        #print p.get_name()

    # Separated function for the confining environment restraint
    if use_confining_environment:
        # print "\nEnforcing the confining environment restraint"
        add_confining_environment(model)

    # Separated function for the bending rigidity restraints
    if CONFIG['kbending'] > 0.0:
        # print "\nEnforcing the bending rigidity restraint"
        theta0 = 0.0
        bending_kforce = CONFIG['kbending']
        add_bending_rigidity_restraint(model, theta0, bending_kforce)

    # Add restraints on single particles
    if single_particle_restraints:
        for ap in single_particle_restraints:
            ap[1] = [(c / SCALE) for c in ap[1]]
            ap[4] /= SCALE
        # This function is specific for IMP
        add_single_particle_restraints(model, single_particle_restraints)

    # Separated function fot the HiC-based restraints
    if use_HiC:
        # print "\nEnforcing the HiC-based Restraints"
        HiCbasedRestraints = HiCRestraints.get_hicbased_restraints()
        add_hicbased_restraints(model, HiCbasedRestraints)

    # Separated function for the excluded volume restraint
    if use_excluded_volume:
        # print "\nEnforcing the excluded_volume_restraints"
        excluded_volume_kforce = CONFIG['kforce']
        evr = add_excluded_volume_restraint(model, model['particles'], excluded_volume_kforce)

    if verbose == 1:
        try:
            print "Total number of restraints: %i" % (
                model['model'].get_number_of_restraints())
            if use_HiC:
                print len(HiCbasedRestraints)
        except:
            print "Total number of restraints: %i" % (
                model['restraints'].get_number_of_restraints()) # 2.6.1 compat
            if use_HiC:
                print len(HiCbasedRestraints)

    # Separated function for the Conjugate gradient optimization
    if verbose == 1:
        print ("\nPerforming the optimization of the Hi-C-based restraints "
               "using conjugate gradient optimisation\n")

    conjugate_gradient_optimization(model, log_energies)

    #try:
    #    log_energies.append(model['model'].evaluate(False))
    #except:
    #    log_energies.append(model['restraints'].evaluate(False)) # 2.6.1 compat
    if verbose >=1:
        if verbose >= 2 or not rand_init % 100:
            print 'Model %s IMP Objective Function: %s' % (
                rand_init, log_energies[-1])
    x, y, z, radius = (FloatKey("x"), FloatKey("y"),
                       FloatKey("z"), FloatKey("radius"))
    result = IMPmodel({'log_objfun' : log_energies,
                       'objfun'     : log_energies[-1],
                       'x'          : [],
                       'y'          : [],
                       'z'          : [],
                       'radius'     : None,
                       'cluster'    : 'Singleton',
                       'rand_init'  : str(rand_init)})
    for part in model['particles'].get_particles():
        result['x'].append(part.get_value(x) * SCALE)
        result['y'].append(part.get_value(y) * SCALE)
        result['z'].append(part.get_value(z) * SCALE)
        if verbose == 3:
            print (part.get_name(), part.get_value(x), part.get_value(y),
                   part.get_value(z), part.get_value(radius))
    # gets radius from last particle, assuming that all are the same
    # include in the loop when radius changes... should be a list then
    result['radius'] = part.get_value(radius)
    result['radius'] = SCALE / 2 
    #for log_energy in log_energies:
    #    print "%.30f" % log_energy
    #print rand_init, log_energies[-1]
    #stdout.flush()
    return result # rand_init, result

#Functions to add Centromeric, Connectivity, Hi-C-based, and Imaging-based restraints
def add_excluded_volume_restraint(model, particle_list, kforce): #, restraints):
    evr = IMP.core.ExcludedVolumeRestraint(particle_list, kforce)

    evr.set_name("ExcludedVolumeRestraint")
    #print "ExcludedVolume", particle_list.get_particles(), "%.30f" % CONFIG['kforce']
    try:
        model['model'].add_restraint(evr)
    except:
        model['restraints'].add_restraint(evr) # 2.6.1 compat
    #restraints.append(evr)
    return evr

def add_bending_rigidity_restraint(model, theta0, bending_kforce): #, restraints):

    harmonic = IMP.core.Harmonic(theta0, bending_kforce)
    for particle in xrange(0,len(LOCI)-2):
        p1  = model['particles'].get_particle(particle)
        p2  = model['particles'].get_particle(particle+1)
        p3  = model['particles'].get_particle(particle+2)

        try: # 2.6.1 compat
            brr = IMP.core.AngleRestraint(harmonic, p1, p2, p3)
        except TypeError:
            brr = IMP.core.AngleRestraint(model['model'], harmonic, p1, p2, p3)
        # print "Adding an Harmonic bending rigidity restraint between particles %s, %s, and %s using parameters theta0 %f and k %f" % (p1,p2,p3,theta0,bending_kforce)
        #restraints.append(brr)

        brr.set_name("BendingRigidityRestraint%s%s%s" % (p1, p2, p3))
        try:
            model['model'].add_restraint(brr)
        except:
            model['restraints'].add_restraint(brr) # 2.6.1 compat

def add_confining_environment(model): #, restraints):
    model['container'] = CONFIG['container']

    if model['container']['shape'] == 'cylinder':
        # define a segment of length equal to the cylinder height
        segment = IMP.algebra.Segment3D(
            IMP.algebra.Vector3D(0,0,0),
            IMP.algebra.Vector3D(model['container']['height'],0,0))
        bb = IMP.algebra.get_bounding_box(segment)
        # define a restraint to keep all the model['particles'] at a distance lower or equal to the cylinder radius
        confining_environment = IMP.container.SingletonsRestraint(
            IMP.core.BoundingBox3DSingletonScore(
                IMP.core.HarmonicUpperBound(model['container']['radius'],
                                            model['container']['cforce']), bb),
            model['particles'])

        confining_environment.set_name("ConfiningEnvironmentRestraint")
        try:
            model['model'].add_restraint(confining_environment)
        except AttributeError:  # 2.6.1 compat
            model['restraints'].add_restraint(confining_environment)
        confining_environment.evaluate(False)
        # print "Adding a confining environment restraint as a HarmonicUpperBound restraint"
        # print "of kforce %d: a cylinder of base radius %f and height %f" % (model['container']['cforce'], model['container']['radius'], model['container']['height'])
        #restraints.append(confining_environment)
    # else:
    #     print "ERROR the shape",model['container']['shape'],"is currently not defined!"

def add_single_particle_restraints(model, single_particle_restraints):

    for restraint in single_particle_restraints:

        if int(restraint[0]) >= len(LOCI):
            continue

        p1 = model['particles'].get_particle(int(restraint[0]))
        pos = IMP.algebra.Vector3D(restraint[1])

        kforce = float(restraint[3])
        dist   = float(restraint[4])

        if restraint[2]  == 'Harmonic':
            rb = IMP.core.Harmonic(dist, kforce)
        elif restraint[2] == 'HarmonicUpperBound':
            rb = IMP.core.HarmonicUpperBound(dist, kforce)
        elif restraint[2] == 'HarmonicLowerBound':
            rb = IMP.core.HarmonicLowerBound(dist, kforce)
        else:
            print "ERROR: RestraintType",restraint[2],"does not exist!"
            return

        ss = IMP.core.DistanceToSingletonScore(rb, pos)
        ar = IMP.core.SingletonRestraint(model['model'], ss, p1)
        ar.set_name("%sSingleDistanceRestraint%s" % (restraint[2], restraint[0]))
        try:
            model['model'].add_restraint(ar)
        except:
            model['restraints'].add_restraint(ar) # 2.6.1 compat

def add_hicbased_restraints(model, HiCbasedRestraints): #, restraints):
    # Add the restraints contained in HiCbasedRestraints
    #print HiCbasedRestraints
    for restraint in HiCbasedRestraints:
        p1 = model['particles'].get_particle(int(restraint[0]))
        p2 = model['particles'].get_particle(int(restraint[1]))

        kforce = float(restraint[3])
        dist   = float(restraint[4])

        if restraint[2] in ['Harmonic', 'NeighborHarmonic']:
            # print "Adding an HarmonicRestraint between particles %s and %s using parameters R0 %f and k %f" % (p1,p2,dist,kforce)
            add_harmonic_restraint(model, p1, p2, dist, kforce) #, restraints)

        elif restraint[2] in ['HarmonicUpperBound', 'NeighborHarmonicUpperBound']:
            # print "Adding an HarmonicLowerBoundRestraint between particles %s and %s using parameters R0 %f and k %f" % (p1,p2,dist,kforce)
            add_harmonic_upperbound_restraint(model, p1, p2, dist, kforce) #, restraints)

        elif restraint[2] == 'HarmonicLowerBound':
            # print "Adding an HarmonicUpperBoundRestraint between particles %s and %s using parameters R0 %f and k %f" % (p1,p2,dist,kforce)
            add_harmonic_lowerbound_restraint(model, p1, p2, dist, kforce) #, restraints)
        else:
            print "ERROR: RestraintType",restraint[2],"does not exist!"

def add_harmonic_restraint(model, p1, p2, dist, kforce, verbose=None): #, restraints, verbose=None):
    try:
        hr = IMP.core.DistanceRestraint(
            IMP.core.Harmonic(dist, kforce),
            p1,
            p2)
    except:
        hr = IMP.core.DistanceRestraint(
            model['model'],
            IMP.core.Harmonic(dist, kforce),
            p1,
            p2)

    #print p1, p2, "Harmonic", "%.30f" % kforce, "%.30f" % dist
    hr.set_name("HarmonicDistanceRestraint%s%s" % (p1, p2))
    try:
        model['model'].add_restraint(hr)
    except:
        model['restraints'].add_restraint(hr) # 2.6.1 compat
    #restraints.append(hr)

    if verbose == 3:
        try:
            print "Total number of restraints: %i" % (
                model['model'].get_number_of_restraints())
        except:
            print "Total number of restraints: %i" % (
                model['restraints'].get_number_of_restraints()) # 2.6.1 compat


def add_harmonic_upperbound_restraint(model, p1, p2, dist, kforce): #, restraints):
    try:
        hubr = IMP.core.DistanceRestraint(
            IMP.core.HarmonicUpperBound(dist, kforce),
            p1,
            p2) # older versions
    except:
        hubr = IMP.core.DistanceRestraint(
            model['model'],
            IMP.core.HarmonicUpperBound(dist, kforce),
            p1,
            p2)

    #print p1, p2, "HarmonicUpperBound", "%.30f" % kforce, "%.30f" % dist
    hubr.set_name("HarmonicUpperBoundDistanceRestraint%s%s" % (p1, p2))
    try:
        model['model'].add_restraint(hubr)
    except:
        model['restraints'].add_restraint(hubr) # 2.6.1 compat
    #restraints.append(hubr)



def add_harmonic_lowerbound_restraint(model, p1, p2, dist, kforce): #, restraints):
    try:
        hlbr = IMP.core.DistanceRestraint(
            IMP.core.HarmonicLowerBound(dist, kforce),
            p1,
            p2) # older versions
    except:
        hlbr = IMP.core.DistanceRestraint(
            model['model'],
            IMP.core.HarmonicLowerBound(dist, kforce),
            p1,
            p2)

    #print p1, p2, "HarmonicLowerBound", "%.30f" % kforce, "%.30f" % dist
    hlbr.set_name("HarmonicLowerBoundDistanceRestraint%s%s" % (p1, p2))
    try:
        model['model'].add_restraint(hlbr)
    except:
        model['restraints'].add_restraint(hlbr) # 2.6.1 compat
    #restraints.append(hlbr)

#Function to perform the scoring function optimization ConjugateGradient
def conjugate_gradient_optimization(model, log_energies):
    restraints = []
    # IMP COMMAND: Call the Conjugate Gradient optimizer
    try:
        restraints.append(model['restraints'])
        scoring_function = IMP.core.RestraintsScoringFunction(restraints)
    except:
        pass

    try:
        lo = IMP.core.ConjugateGradients()
        lo.set_model(model['model'])
    except: # since version 2.5, IMP goes this way
        lo = IMP.core.ConjugateGradients(model['model'])
    try:
        lo.set_scoring_function(scoring_function) # 2.6.1 compat
    except:
        pass

    #print LSTEPS, NROUNDS, STEPS
    o = IMP.core.MonteCarloWithLocalOptimization(lo, LSTEPS)
    try:
        o.set_scoring_function(scoring_function) # 2.6.1 compat
    except:
        pass

    o.set_return_best(True)
    fk = IMP.core.XYZ.get_xyz_keys()
    ptmp = model['particles'].get_particles()

    x, y, z, radius = (FloatKey("x"), FloatKey("y"),
                           FloatKey("z"), FloatKey("radius"))

    mov = IMP.core.NormalMover(ptmp, fk, 0.25 / SCALE)
    o.add_mover(mov)
    # o.add_optimizer_state(log)

    # Start optimization and save a VRML after 100 MC moves
    try:
        log_energies.append(model['model'].evaluate(False))
    except:
        log_energies.append(model['restraints'].evaluate(False)) # 2.6.1 compat

    #"""simulated_annealing: perform simulated annealing for at most nrounds
    # iterations. The optimization stops if the score does not change more than
    #    a value defined by endLoopValue and for stopCount iterations.
    #   @param endLoopCount = Counter that increments if the score of two models
    # did not change more than a value
    #   @param stopCount = Maximum values of iteration during which the score
    # did not change more than a specific value
    #   @paramendLoopValue = Threshold used to compute the value  that defines
    # if the endLoopCounter should be incremented or not"""
    # IMP.fivec.simulatedannealing.partial_rounds(m, o, nrounds, steps)
    endLoopCount = 0
    stopCount = 10
    endLoopValue = 0.00001

    # alpha is a parameter that takes into account the number of particles in
    # the model (len(LOCI)).
    # The multiplier (in this case is 1.0) is used to give a different weight
    # to the number of particles
    alpha = 1.0 * len(LOCI)

    # During the firsts hightemp iterations, do not stop the optimization
    hightemp = int(0.025 * NROUNDS)
    for i in range(0, hightemp):
        temperature = alpha * (1.1 * NROUNDS - i) / NROUNDS
        o.set_kt(temperature)
        log_energies.append(o.optimize(STEPS))
        #if i == 1:
        #    for p in ptmp:
        #        print p,p.get_value(x),p.get_value(y),p.get_value(z)
        if VERBOSE == 3:
            print i, log_energies[-1], o.get_kt()
    # After the firsts hightemp iterations, stop the optimization if the score
    # does not change by more than a value defined by endLoopValue and
    # for stopCount iterations
    lownrj = log_energies[-1]
    for i in range(hightemp, NROUNDS):
        temperature = alpha * (1.1 * NROUNDS - i) / NROUNDS
        o.set_kt(temperature)
        log_energies.append(o.optimize(STEPS))
        if VERBOSE == 3:
            print i, log_energies[-1], o.get_kt()
        # Calculate the score variation and check if the optimization
        # can be stopped or not
        if lownrj > 0:
            deltaE = fabs((log_energies[-1] - lownrj) / lownrj)
        else:
            deltaE = log_energies[-1]

        if (deltaE < endLoopValue and endLoopCount == stopCount):
            break
        elif (deltaE < endLoopValue and endLoopCount < stopCount):
            endLoopCount += 1
            lownrj = log_energies[-1]
        else:
            endLoopCount = 0
            lownrj = log_energies[-1]
    #"""simulated_annealing_full: preform simulated annealing for nrounds
    # iterations"""
    # # IMP.fivec.simulatedannealing.full_rounds(m, o, nrounds, steps)
    # alpha = 1.0 * len(LOCI)
    # for i in range(0,nrounds):
    #    temperature = alpha * (1.1 * nrounds - i) / nrounds
    #    o.set_kt(temperature)
    #    e = o.optimize(steps)
    #    print str(i) + " " + str(e) + " " + str(o.get_kt())