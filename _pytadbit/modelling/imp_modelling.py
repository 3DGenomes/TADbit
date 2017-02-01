"""
05 Jul 2013


"""
from pytadbit.modelling.IMP_CONFIG           import CONFIG, NROUNDS, STEPS, LSTEPS
from pytadbit.modelling.structuralmodels import StructuralModels
from pytadbit.modelling.impmodel         import IMPmodel
from scipy                         import polyfit
from math                          import fabs, pow as power
from cPickle                       import load, dump
from sys                           import stdout
from os.path                       import exists
import multiprocessing as mu

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
    #VERBOSE = 3

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
    model = {'rk'    : IMP.FloatKey("radius"),
             'model' : Model(),
             'rs'    : None, # 2.6.1 compat 
             'ps'    : None}
    model['ps'] = ListSingletonContainer(IMP.core.create_xyzr_particles(
        model['model'], len(LOCI), RADIUS, 100000))
    model['ps'].set_name("")

    # set container
    try:
        model['rs'] = IMP.RestraintSet(model['model']) # 2.6.1 compat 
    except:
        pass
    model['container'] = CONFIG['container']
    if model['container']['shape'] == 'cylinder':
         # define a segment of a given size
        segment = IMP.algebra.Segment3D(
            IMP.algebra.Vector3D(0,0,0),
            IMP.algebra.Vector3D(model['container']['height'],0,0))
        bb = IMP.algebra.get_bounding_box(segment)
        rb = IMP.container.SingletonsRestraint(
            IMP.core.BoundingBox3DSingletonScore(
                IMP.core.HarmonicUpperBound(model['container']['radius'],
                                            model['container']['cforce']), bb),
            model['ps'])
	try:
	       model['model'].add_restraint(rb)
	except:
	       model['rs'].add_restraint(rb) # 2.6.1 compat
        rb.evaluate(False)
    # elif model['container']['shape']:
    #     raise noti
    
    for i in range(0, len(LOCI)):
        p = model['ps'].get_particle(i)
        p.set_name(str(LOCI[i]))
        #p.set_value(model['rk'], RADIUS)
    restraints = {}
    for i in range(len(LOCI)):
        p1 = model['ps'].get_particle(i)
        x = p1.get_name()
        if model['container']['shape'] == 'sphere':
            ub  = IMP.core.HarmonicUpperBound(
                model['container']['properties'][0], CONFIG['kforce'] * 10)
            ss  = IMP.core.DistanceToSingletonScore(
                ub, model['container']['center'])
            rss = IMP.core.SingletonRestraint(ss, p1)
            try:
	            model['model'].add_restraint(rss) 
            except:
	            model['rs'].add_restraint(rss) # 2.6.1 compat
            rss.evaluate(False)
        for j in range(i+1, len(LOCI)):
            p2 = model['ps'].get_particle(j)
            y = p2.get_name()
            typ, dist, frc = addHarmonicPair(model, p1, p2, x, y, j, dry=True)
            if VERBOSE >= 1:
                stdout.write('%s\t%s\t%s\t%s\t%s\n' % (typ, x, y, dist, frc))
            if typ[-1] == 'a':
                typ = 'H'
            elif typ[-1] == 'l':
                typ = 'L'
            elif typ[-1] == 'u':
                typ = 'U'
            elif typ[-1] == 'n':
                typ = 'C'
            else:
                continue
            restraints[tuple(sorted((x, y)))] = typ[-1], dist, frc
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


def generate_IMPmodel(rand_init):
    """
    Generates one IMP model
    
    :param rand_init: random number kept as model key, for reproducibility.

    :returns: a model, that is a dictionary with the log of the objective
       function value optimization, and the coordinates of each particles.

    """
    verbose = VERBOSE
    IMP.random_number_generator.seed(rand_init)

    log_energies = []
    model = {'rk'    : IMP.FloatKey("radius"),
             'model' : Model(),
             'rs'    : None, # 2.6.1 compat
             'ps'    : None,
             'pps'   : None}
    model['ps'] = ListSingletonContainer(IMP.core.create_xyzr_particles(
        model['model'], len(LOCI), RADIUS, 100000))
    model['ps'].set_name("")

    # initialize each particles
    for i in range(0, len(LOCI)):
        p = model['ps'].get_particle(i)
        p.set_name(str(LOCI[i]))
        # computed following the relationship with the 30nm vs 40nm fiber
        #p.set_value(model['rk'], RADIUS)

    # Restraints between pairs of LOCI proportional to the PDIST
    try:
        model['pps']  = IMP.kernel.ParticlePairsTemp()
    except:
        model['pps']  = IMP.ParticlePairsTemp() # 2.6.1 compat

    # CALL BIG FUNCTION
    if rand_init == START and verbose == 0.5:
        verbose = 1
        stdout.write("# Harmonic\tpart1\tpart2\tdist\tkforce\n")

    # set container
    try:
        model['rs'] = IMP.RestraintSet(model['model']) # 2.6.1 compat
    except:
        pass
    restraints = [] # 2.6.1 compat
    model['container'] = CONFIG['container']
    if model['container']['shape'] == 'cylinder':
         # define a segment of a given size
        segment = IMP.algebra.Segment3D(
            IMP.algebra.Vector3D(0,0,0),
            IMP.algebra.Vector3D(model['container']['height'],0,0))
        bb = IMP.algebra.get_bounding_box(segment)
        rb = IMP.container.SingletonsRestraint(
            IMP.core.BoundingBox3DSingletonScore(
                IMP.core.HarmonicUpperBound(model['container']['radius'],
                                            model['container']['cforce']), bb),
            model['ps'])
        try:
	        model['model'].add_restraint(rb)
        except:
	        model['rs'].add_restraint(rb) # 2.6.1 compat
        rb.evaluate(False)
        
    # elif model['container']['shape']:
    #     raise noti

    addAllHarmonics(model)

    # Setup an excluded volume restraint between a bunch of particles
    # with radius
    r = IMP.core.ExcludedVolumeRestraint(model['ps'], CONFIG['kforce'])
    try:
        model['model'].add_restraint(r)
    except:
        model['rs'].add_restraint(r) # 2.6.1 compat
        restraints.append(model['rs'])
        scoring_function = IMP.core.RestraintsScoringFunction(restraints)

    if verbose == 3:
	   try:
	       "Total number of restraints: %i" % (
		      model['model'].get_number_of_restraints())
	   except:
	       "Total number of restraints: %i" % (
		      model['rs'].get_number_of_restraints()) # 2.6.1 compat

    # Set up optimizer
    try:
        lo = IMP.core.ConjugateGradients()
        lo.set_model(model['model'])
    except: # since version 2.5, IMP goes this way
        lo = IMP.core.ConjugateGradients(model['model'])
    try:
        lo.set_scoring_function(scoring_function) # 2.6.1 compat
    except:
        pass
    o = IMP.core.MonteCarloWithLocalOptimization(lo, LSTEPS)
    try:
        o.set_scoring_function(scoring_function) # 2.6.1 compat
    except:
        pass
    o.set_return_best(True)
    fk = IMP.core.XYZ.get_xyz_keys()
    ptmp = model['ps'].get_particles()
    mov = IMP.core.NormalMover(ptmp, fk, 0.25)
    o.add_mover(mov)
    # o.add_optimizer_state(log)

    # Optimizer's parameters
    if verbose == 3:
         "nrounds: %i, steps: %i, lsteps: %i" % (NROUNDS, STEPS, LSTEPS)

    # Start optimization and save an VRML after 100 MC moves
    try:
	     log_energies.append(model['model'].evaluate(False))
    except:
	     log_energies.append(model['rs'].evaluate(False)) # 2.6.1 compat
    if verbose == 3:
         "Start", log_energies[-1]

    #"""simulated_annealing: preform simulated annealing for at most nrounds
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
        if verbose == 3:
             i, log_energies[-1], o.get_kt()
    # After the firsts hightemp iterations, stop the optimization if the score
    # does not change by more than a value defined by endLoopValue and
    # for stopCount iterations
    lownrj = log_energies[-1]
    for i in range(hightemp, NROUNDS):
        temperature = alpha * (1.1 * NROUNDS - i) / NROUNDS
        o.set_kt(temperature)
        log_energies.append(o.optimize(STEPS))
        if verbose == 3:
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

    try:
        log_energies.append(model['model'].evaluate(False))
    except:
        log_energies.append(model['rs'].evaluate(False)) # 2.6.1 compat
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
    for part in model['ps'].get_particles():
        result['x'].append(part.get_value(x))
        result['y'].append(part.get_value(y))
        result['z'].append(part.get_value(z))
        if verbose == 3:
            print (part.get_name(), part.get_value(x), part.get_value(y),
                   part.get_value(z), part.get_value(radius))
    # gets radius from last particle, assuming that all are the same
    result['radius'] = part.get_value(radius)
    return result # rand_init, result


def addAllHarmonics(model):
    """
    Add harmonics to all pair of particles.
    """
    for i in range(len(LOCI)):
        p1 = model['ps'].get_particle(i)
        x = p1.get_name()
        for j in range(i+1, len(LOCI)):
            p2 = model['ps'].get_particle(j)
            y = p2.get_name()
            addHarmonicPair(model, p1, p2, x, y, j)


def addHarmonicPair(model, p1, p2, x, y, j, dry=False):
    """
    add harmonic to a given pair of particles

    :param model: a model dictionary that contains IMP model, singleton
       containers...
    :param p1: first particle
    :param p2: second particle
    :param x: first particle name
    :param y: second particle name
    :param j: id of second particle
    :param num_loci1: index of the first particle
    :param num_loci2: index of the second particle
    """
    num_loci1, num_loci2 = int(x), int(y)
    seqdist = num_loci2 - num_loci1
    restraint = ('no', 0, 0)
    freq = float('nan')
    # SHORT RANGE DISTANCE BETWEEN TWO CONSECUTIVE LOCI
    if seqdist == 1:
        kforce = CONFIG['kforce']
        if x in PDIST and y in PDIST[x] and PDIST[x][y] > CONFIG['upfreq']:
            dist = distConseq12(PDIST[p1.get_name()][p2.get_name()])
            if not dry:
                addHarmonicNeighborsRestraints(model, p1, p2, dist, kforce)
            else:
                return ("addHn", dist, kforce)
        else:
            dist = (p1.get_value(model['rk']) + p2.get_value(model['rk']))
            # dist = (p1.get_value(rk) + p2.get_value(rk))
            if not dry:
                addHarmonicUpperBoundRestraints(model, p1, p2, dist, kforce)
            else:
                return ("addHu", dist, kforce)
    # SHORT RANGE DISTANCE BETWEEN TWO SEQDIST = 2
    elif seqdist == 2:
        p3 = model['ps'].get_particle(j-1)
        kforce = CONFIG['kforce']
        dist = (p1.get_value(model['rk']) + p2.get_value(model['rk'])
                + 2.0 * p3.get_value(model['rk']))
        # dist = (p1.get_value(rk) + p2.get_value(rk))
        if not dry:
            addHarmonicUpperBoundRestraints(model, p1, p2, dist, kforce)
        else:
            return ("addHu", dist, kforce)
    # LONG RANGE DISTANCE DISTANCE BETWEEN TWO NON-CONSECUTIVE LOCI
    elif x in PDIST and y in PDIST[x]:
        freq = PDIST[x][y]
        kforce = kForce(freq)
    # X IN PDIST BUT Y NOT IN PDIST[X]
    elif x in PDIST:
        prevy = str(num_loci2 - 1)
        posty = str(num_loci2 + 1)
        # mean dist to prev and next part are used with half weight
        freq = (PDIST[x].get(prevy, PDIST[x].get(posty, float('nan'))) +
                PDIST[x].get(posty, PDIST[x].get(prevy, float('nan')))) / 2

        kforce = 0.5 * kForce(freq)
    # X NOT IN PDIST
    else:
        prevx = str(num_loci1 - 1)
        postx = str(num_loci1 + 1)
        prevx = prevx if prevx in PDIST else postx
        postx = postx if postx in PDIST else prevx
        try:
            freq = (PDIST[prevx].get(y, PDIST[postx].get(y, float('nan'))) +
                    PDIST[postx].get(y, PDIST[prevx].get(y, float('nan')))) / 2
        except KeyError:
            pass
        kforce = 0.5 * kForce(freq)

    # FREQUENCY > UPFREQ
    if freq > CONFIG['upfreq']:
        if not dry:
            addHarmonicRestraints(model, p1, p2, distance(freq), kforce)
        else:
            return ("addHa", distance(freq), kforce)
    # FREQUENCY > LOW THIS HAS TO BE THE THRESHOLD FOR
    # "PHYSICAL INTERACTIONS"
    elif freq < CONFIG['lowfreq']:
        if not dry:
            addHarmonicLowerBoundRestraints(model, p1, p2, distance(freq), kforce)
        else:
            return ("addHl", distance(freq), kforce)
    if dry:
        return restraint

def distConseq12(freq):
    """
    Function mapping the Z-scores into distances for neighbor fragments
    """
    return (NSLOPE * freq) + NINTERCEPT

def distance(freq):
    """
    Function mapping the Z-scores into distances for non-neighbor fragments
    """
    return (SLOPE * freq) + INTERCEPT

def addHarmonicNeighborsRestraints(model, p1, p2, dist, kforce):
    p = IMP.ParticlePair(p1, p2)
    model['pps'].append(p)
    try:
        dr = IMP.core.DistanceRestraint(
            model['model'], IMP.core.Harmonic(dist, kforce),p1, p2)
    except TypeError:
        dr = IMP.core.DistanceRestraint(
            IMP.core.Harmonic(dist, kforce),p1, p2) # older versions
    try:
        model['model'].add_restraint(dr)
    except:
        model['rs'].add_restraint(dr) # 2.6.1 compat

def addHarmonicUpperBoundRestraints(model, p1, p2, dist, kforce):
    p = IMP.ParticlePair(p1, p2)
    model['pps'].append(p)
    try:
        dr = IMP.core.DistanceRestraint(
            model['model'], IMP.core.HarmonicUpperBound(dist, kforce), p1, p2)
    except TypeError:
        dr = IMP.core.DistanceRestraint(
            IMP.core.HarmonicUpperBound(dist, kforce), p1, p2) # older versions
    try:
        model['model'].add_restraint(dr)
    except:
        model['rs'].add_restraint(dr) # 2.6.1 compat

def addHarmonicRestraints(model, p1, p2, dist, kforce):
    p = IMP.ParticlePair(p1, p2)
    model['pps'].append(p)
    try:
        dr = IMP.core.DistanceRestraint(
            model['model'], IMP.core.Harmonic(dist, kforce), p1, p2)
    except TypeError:
        dr = IMP.core.DistanceRestraint(
            IMP.core.Harmonic(dist, kforce), p1, p2) # older versions
    try:
        model['model'].add_restraint(dr)
    except:
        model['rs'].add_restraint(dr) # 2.6.1 compat

def addHarmonicLowerBoundRestraints(model, p1, p2, dist, kforce):
    p = IMP.ParticlePair(p1, p2)
    model['pps'].append(p)
    try:
        dr = IMP.core.DistanceRestraint(
            model['model'], IMP.core.HarmonicLowerBound(dist, kforce), p1, p2)
    except TypeError:
        dr = IMP.core.DistanceRestraint(
            IMP.core.HarmonicLowerBound(dist, kforce), p1, p2) # older versions
    try:
        model['model'].add_restraint(dr)
    except:
        model['rs'].add_restraint(dr) # 2.6.1 compat


def kForce(freq):
    """
    Function to assign to each restraint a force proportional to the underlying
    experimental value.
    """
    return power(fabs(freq), 0.5 )


class TADbitModelingOutOfBound(Exception):
    pass
