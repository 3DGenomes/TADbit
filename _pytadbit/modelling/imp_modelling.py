"""
05 Jul 2013


"""

from math            import fabs, pow as power
from cPickle         import load, dump
from sys             import stdout
from os.path         import exists
import multiprocessing as mu
from scipy           import polyfit

from pytadbit.modelling.IMP_CONFIG       import CONFIG, NROUNDS, STEPS, LSTEPS
from pytadbit.modelling.structuralmodels import StructuralModels
from pytadbit.modelling.impmodel         import IMPmodel

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
                       first=None, container=None):
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
          'dmel_01': {
              # Paramaters for the Hi-C dataset from:
              'reference' : 'victor corces dataset 2013',

              # Force applied to the restraints inferred to neighbor particles
              'kforce'    : 5,

              # Maximum experimental contact distance
              'maxdist'   : 600, # OPTIMIZATION: 500-1200

              # Maximum threshold used to decide which experimental values have to be
              # included in the computation of restraints. Z-score values greater than upfreq
              # and less than lowfreq will be included, while all the others will be rejected
              'upfreq'    : 0.3, # OPTIMIZATION: min/max Z-score

              # Minimum thresholds used to decide which experimental values have to be
              # included in the computation of restraints. Z-score values bigger than upfreq
              # and less that lowfreq will be include, whereas all the others will be rejected
              'lowfreq'   : -0.7 # OPTIMIZATION: min/max Z-score

              # Space occupied by a nucleotide (nm)
              'scale'     : 0.005

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

    :returns: a StructuralModels object

    """

    # Main config parameters
    global CONFIG
    CONFIG = config or CONFIG['dmel_01']
    CONFIG['kforce'] = CONFIG.get('kforce', 5)

    # setup container
    try:
        CONFIG['container'] = {'shape' : container[0],
                               'radius': container[1],
                               'height': container[2],
                               'cforce': container[3]}
    except:
        CONFIG['container'] = {'shape' : None,
                               'radius': None,
                               'height': None,
                               'cforce': None}
    # Particles initial radius
    global RADIUS
    RADIUS = float(resolution * CONFIG['scale']) / 2
    CONFIG['lowrdist'] = RADIUS * 2.


    if CONFIG['lowrdist'] > CONFIG['maxdist']:
        raise TADbitModelingOutOfBound(
            ('ERROR: we must prevent you from doing this for the safe of our' +
             'universe...\nIn this case, maxdist must be higher than %s\n' +
             '   -> resolution times scale -- %s*%s)') % (
                CONFIG['lowrdist'], resolution, CONFIG['scale']))

    # get SLOPE and regression for all particles of the z-score data
    global SLOPE, INTERCEPT
    zsc_vals = [zscores[i][j] for i in zscores for j in zscores[i]
                if abs(int(i) - int(j)) > 1] # condition is to avoid
                                             # taking into account selfies
                                             # and neighbors
    SLOPE, INTERCEPT   = polyfit([min(zsc_vals), max(zsc_vals)],
                                 [CONFIG['maxdist'], CONFIG['lowrdist']], 1)
    # get SLOPE and regression for neighbors of the z-score data
    global NSLOPE, NINTERCEPT
    xarray = [zscores[i][j] for i in zscores for j in zscores[i]
              if abs(int(i) - int(j)) <= (close_bins + 1)]
    yarray = [RADIUS * 2 for _ in xrange(len(xarray))]
    NSLOPE, NINTERCEPT = polyfit(xarray, yarray, 1)

    global LOCI
    # if z-scores are generated outside TADbit they may not start at zero
    if first == None:
        first = min([int(j) for i in zscores for j in zscores[i]] +
                    [int(i) for i in zscores])
    LOCI  = range(first, nloci + first)

    # Z-scores
    global PDIST
    PDIST = zscores
    # random inital number
    global START
    START = start
    # verbose
    global VERBOSE
    VERBOSE = verbose

    models, bad_models = multi_process_model_generation(
        n_cpus, n_models, n_keep, keep_all)

    try:
        xpr = experiment
        crm = xpr.crm
        description = {'identifier'        : xpr.identifier,
                       'chromosome'        : coords['crm'],
                       'start'             : xpr.resolution * coords['start'],
                       'end'               : xpr.resolution * coords['end'],
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
    except AttributeError: # case we are doing optimization
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
        return StructuralModels(
            len(LOCI), models, bad_models, resolution, original_data=values,
            zscores=zscores, config=CONFIG, experiment=experiment, zeros=zeros,
            restraints=_get_restraints(),
            description=description)


def _get_restraints():
    """
    Same function as addAllHarmonic but just to get restraints
    """
    restraint_names = {'None'               : None,
                       'Harmonic'           : 'a',
                       'NeighborHarmonic'   : 'n',
                       'HarmonicUpperBound' : 'u',
                       'NeighborHarmonicUpperBound' : 'u',
                       'HarmonicLowerBound' : 'l'}

    model = {'radius'     : IMP.FloatKey("radius"),
             'model'      : Model(),
             'restraints' : None, # 2.6.1 compat
             'particles'  : None}

    # set container
    try:
        model['restraints'] = IMP.RestraintSet(model['model']) # 2.6.1 compat
    except:
        pass

    model['particles'] = ListSingletonContainer(IMP.core.create_xyzr_particles(
        model['model'], len(LOCI), RADIUS, 100000))

    restraints = {}
    for i, j, RestraintType, kforce, dist in get_hicbased_restraints(model, CONFIG['kforce']):
        restraints[tuple(sorted((i, j)))] = restraint_names[RestraintType], dist, kforce
    return restraints


def multi_process_model_generation(n_cpus, n_models, n_keep, keep_all):
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
                                           args=(rand_init,))

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



def generate_IMPmodel(rand_init, use_HiC=True, use_confining_environment=True,
                      use_bending_rigidity=False, use_excluded_volume=True):
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

    #restraints   = []
    log_energies = []
    model = {'radius'     : IMP.FloatKey("radius"),
             'model'      : Model(),
             'particles'  : None,
             'restraints' : None} # 2.6.1 compat
    model['particles'] = ListSingletonContainer(
        IMP.core.create_xyzr_particles(model['model'], len(LOCI), RADIUS, 100000))  # last number is box size

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
    if use_bending_rigidity:
        # print "\nEnforcing the bending rigidity restraint"
        theta0 = 0.0
        bending_kforce=CONFIG['kforce']
        add_bending_rigidity_restraint(model, theta0, bending_kforce)

    # Separated function fot the HiC-based restraints
    if use_HiC:
        # print "\nEnforcing the HiC-based Restraints"
        HiCbasedRestraints = get_hicbased_restraints(model, CONFIG['kforce'])
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
            print len(HiCbasedRestraints)
        except:
            print "Total number of restraints: %i" % (
                model['restraints'].get_number_of_restraints()) # 2.6.1 compat
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
        result['x'].append(part.get_value(x))
        result['y'].append(part.get_value(y))
        result['z'].append(part.get_value(z))
        if verbose == 3:
            print (part.get_name(), part.get_value(x), part.get_value(y),
                   part.get_value(z), part.get_value(radius))
    # gets radius from last particle, assuming that all are the same
    # include in the loop when radius changes... should be a list then
    result['radius'] = part.get_value(radius)
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



def get_hicbased_restraints(model, nnkforce):

    # HiCbasedRestraints is a list of restraints returned by this function.
    # Each entry of the list is a list of 5 elements describing the details of the restraint:
    # 0 - particle_i
    # 1 - particle_j
    # 2 - type_of_restraint = Harmonic or HarmonicLowerBound or HarmonicUpperBound
    # 3 - the kforce of the restraint
    # 4 - the equilibrium (or maximum or minimum respectively) distance associated to the restraint

    HiCbasedRestraints = []

    for i in range(len(LOCI)):
        for j in range(i + 1, len(LOCI)):

            # Compute the sequence separation (in particles) depending on it the restraint changes
            seqdist = abs(j - i)

            # 1 - CASE OF TWO CONSECUTIVE LOCI (NEAREST NEIGHBOR PARTICLES)
            if seqdist == 1:
                RestraintType, dist = get_nearest_neighbors_restraint_distance(model, i, j)
                kforce = nnkforce

            # 2 - CASE OF 2 SECOND NEAREST NEIGHBORS SEQDIST = 2
            if seqdist == 2:
                RestraintType, dist = get_second_nearest_neighbors_restraint_distance(model, i, j)
                kforce = nnkforce

            # 3 - CASE OF TWO NON-CONSECUTIVE PARTICLES SEQDIST > 2
            if seqdist >  2:
                RestraintType, kforce, dist = get_long_range_restraints_kforce_and_distance(i, j)
                if RestraintType == "None":
                    #print "No HiC-based restraint between particles %d and %d" % (i,j)
                    continue
            
            HiCbasedRestraints.append([i, j, RestraintType, kforce, dist])

    return HiCbasedRestraints




#Functions to add restraints: HarmonicRestraints , HarmonicUpperBoundRestraints , HarmonicLowerBoundRestraints
#addNearestNeighborsRestraint , addSecondNearestNeighborsRestraint , addLongRangeRestraints
def get_nearest_neighbors_restraint_distance(model, i, j):
    x=str(i)
    y=str(j)

    if x in PDIST and y in PDIST[x] and PDIST[x][y] > CONFIG['upfreq']:
        # When p1 and p2 have a contact propensity larger that upfreq
        # their spatial proximity and a partial overlap between them is enforced
        # The equilibrium distance of the spring is inferred from the 3C based Z-score
        RestraintType = "NeighborHarmonic"
        dist = distance(PDIST[x][y],NSLOPE,NINTERCEPT)
        #print "Distance = ", dist
    else:
        # When p1 and p2 have a contact propensity lower than upfreq they are simply connected to each other
        p1 = model['particles'].get_particle(i)
        p2 = model['particles'].get_particle(j)
        RestraintType = "NeighborHarmonicUpperBound"
        dist = p1.get_value(model['radius']) + p2.get_value(model['radius'])
    return RestraintType , dist



def get_second_nearest_neighbors_restraint_distance(model, i, j):
    # IMP COMMAND: Consider the particles i, j and the particle between i and j
    p1      = model['particles'].get_particle(i)
    p2      = model['particles'].get_particle(j)
    pmiddle = model['particles'].get_particle(j-1)

    # The equilibrium distance is the sum of the radii of particles p1 and p2, and of the diameter of particle pmiddle
    RestraintType = "HarmonicUpperBound"
    dist = p1.get_value(model['radius']) + p2.get_value(model['radius']) + 2.0 * pmiddle.get_value(model['radius'])
    #print p1.get_value(model['radius']) , p2.get_value(model['radius']) , pmiddle.get_value(model['radius'])

    #print RestraintType , dist
    return RestraintType , dist



def get_long_range_restraints_kforce_and_distance(i, j):
    x = str(i)
    y = str(j)

    Zscore = float('nan')

    # For non consecutive particles the kforce is a function of the *C based Zscore
    # First we define the The kforce of the harmonic restraint. It is different for 3 scenarios...

    RestraintType = "None"
    kforce        = 0.0
    dist          = 0.0

    # 1 - If the Z-score between i and j is defined
    if x in PDIST and y in PDIST[x]:
        # Get the Zscore between particles p1 and p2
        Zscore = PDIST[x][y]
        kforce = k_force(Zscore)

    # 2 - If the Z-score is defined only for particle i (In the Hi-C matrix you could encounter zero values next to very high entries)
    elif x in PDIST:
        prevy = str(j - 1)
        posty = str(j + 1)
        # The Zscore is compute as the average of the Z-scores of p2 nearest neighbor particles with p1
        Zscore = (PDIST[x].get(prevy, PDIST[x].get(posty, float('nan'))) +
                  PDIST[x].get(posty, PDIST[x].get(prevy, float('nan')))) / 2
        kforce = 0.5 * k_force(Zscore)

    # 3 - If the Z-score is defined only for particle j
    else:
        prevx = str(i - 1)
        postx = str(i + 1)
        prevx = prevx if prevx in PDIST else postx
        postx = postx if postx in PDIST else prevx
        try:
            Zscore = (PDIST[prevx].get(y, PDIST[postx].get(y, float('nan'))) +
                      PDIST[postx].get(y, PDIST[prevx].get(y, float('nan')))) / 2
            # For non consecutive particles the kforce is a function of the *C based Zscore
        except KeyError:
            pass
        kforce = 0.5 * k_force(Zscore)


    # If the ZSCORE > UPFREQ the spatial proximity of particles p1 and p2 is favoured
    if Zscore > CONFIG['upfreq']:
        RestraintType = "Harmonic"
        dist = distance(Zscore, SLOPE, INTERCEPT)

    # If the ZSCORE < LOWFREQ the particles p1 and p2 are restrained to be far from each other.
    elif Zscore < CONFIG['lowfreq']:
        RestraintType = "HarmonicLowerBound"
        dist = distance(Zscore, SLOPE, INTERCEPT)

    return RestraintType, kforce, dist



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

    mov = IMP.core.NormalMover(ptmp, fk, 0.25)
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



#Function to translate the Zscore value into distances and kforce values
def distance(Zscore, slope, intercept):
    """
    Function mapping the Z-scores into distances for neighbor and non-neighbor fragments (slope, intercept) are different
    """
    return (slope * Zscore) + intercept

def k_force(Zscore):
    """
    Function to assign to each restraint a force proportional to the underlying
    experimental value.
    """
    return power(fabs(Zscore), 0.5)


class TADbitModelingOutOfBound(Exception):
    pass
