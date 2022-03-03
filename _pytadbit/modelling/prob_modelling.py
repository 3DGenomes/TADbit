"""
23 Jul 2020


"""
from __future__ import print_function

from future import standard_library
from scipy.sparse.base import issparse
standard_library.install_aliases()

from os import path
from pickle import dump
from itertools import product
from scipy.spatial.distance import pdist,squareform
from scipy.ndimage import gaussian_filter
from scipy.sparse.linalg import eigsh
import multiprocessing as mu
import numpy as np

sparse_avail=False
try:
    from scipy.sparse import issparse, csr_matrix, triu
    sparse_avail=True
except ImportError:
    
    pass

import pytadbit.modelling.globals as globals
from pytadbit import HiC_data
from pytadbit.modelling.IMP_CONFIG       import CONFIG
from pytadbit.modelling.structuralmodels import StructuralModels
from pytadbit.modelling.imp_modelling    import generate_IMPmodel,\
    add_single_particle_restraints
from pytadbit.modelling.impmodel         import IMPmodel
from pytadbit.modelling.restraints       import ProbabilityBasedRestraintsList
from pytadbit.modelling.restraints       import ProbabilityBasedRestraints



def generate_3d_models(hic_data, beg, end, dist, tf_model, 
                       binsAround, resolution, short_restr_dist=1.5e6,
                       tf_model_far=None, resolution_far=None, far_restr_dist=1.5e6,
                       start=1, n_models=5000, n_keep=1000,
                       n_cpus=1, keep_all=False,
                       verbose=0, outfile=None, config=None,
                       experiment=None, coords=None, zeros=None,
                       container=None, use_HiC=True,
                       use_confining_environment=False,
                       single_particle_restraints=None,
                       save_restraints_folder=None,
                       include_restraints_in_models=True):

    """
    This function generates three-dimensional models from a probabilistic distribution
    of distances derived from a machine learning (ML) model trained with microscopy
    data and feed from Hi-C data.
    The final analysis will be performed on the n_keep top models.

    :param hic_data: the dictionary of the of Hi-C pairwise interactions
    :param dist: scipy.stats distribution for which the ML model has been trained
    :param tf_model: tensorflow model which input will be hic submatrices
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
    :param False verbose: if set to True, information about the distance and force
       between particles will be printed. If verbose is 0.5 than
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
    
    globals.init()
    np.random.seed(seed=start)
    # Setup CONFIG
    if isinstance(config, dict):
        CONFIG.update(config)
    elif config:
        raise Exception('ERROR: "config" must be a dictionary')

    CONFIG['resolution'] = resolution

    globals.SCALE = float(resolution * CONFIG['scale'])

    # Setup and scale CONFIG['container']
    try:
        CONFIG['container'] = {'shape' : container[0],
                               'radius': container[1] / globals.SCALE,
                               'height': container[2] / globals.SCALE,
                               'cforce': container[3]}
    except:
        CONFIG['container'] = {'shape' : None,
                               'radius': None,
                               'height': None,
                               'cforce': None}

    nloci = int((end - beg + 1))
    globals.LOCI  = list(range(nloci))

    globals.START = start
    globals.VERBOSE = verbose

    mat_width = mat_height = binsAround*2+1
    if sparse_avail and issparse(hic_data):
        if not isinstance(hic_data, csr_matrix):
            raise ValueError('Matrix given must be of CSR format.')
        fullmat = hic_data.copy()
        fullmat.setdiag(np.NaN)
        max_mat = np.log1p(fullmat.max())
        min_mat = np.log1p(fullmat.min())
        fullmat=fullmat.tolil()
    else:
        fullmat = np.asarray(hic_data)
        np.fill_diagonal(fullmat,np.nan)
        max_mat = np.log1p(np.nanmax(fullmat))
        min_mat = np.log1p(np.nanmin(fullmat))
    region_zeros = []
    if zeros is not None:
        zeros_models = tuple([zeros[i] for i in range(beg,end+1)])
        region_zeros = [i-beg for i in range(beg,end+1) if not zeros[i]]
    else:
        zeros_models = tuple([True for i in range(beg,end+1)])

    pairs = []
    submatrices = []
    pair_dists = []
    for i in range(beg, end+1):
        if i-beg in region_zeros:
            continue
        for j in range(i+1, end+1):
            if j-beg in region_zeros:
                continue
            if abs(i-j)*resolution>short_restr_dist:
                continue
            submat = fullmat[max(0,i-binsAround):i+binsAround+1,
                             max(0,j-binsAround):j+binsAround+1]
            if sparse_avail and issparse(submat):
                submat = submat.toarray()
            submat = np.array(submat)
            if i-binsAround < 0:
                pad_v = mat_width-submat.shape[0]
                submat = np.pad(submat,((pad_v, 0), (0, 0)), mode='reflect')
            if j-binsAround < 0:
                pad_v = mat_height-submat.shape[1]
                submat = np.pad(submat,((0, 0), (pad_v, 0)), mode='reflect')
            if i+binsAround+1 >= fullmat.shape[0]:
                pad_v = mat_width-submat.shape[0]
                submat = np.pad(submat,((0, pad_v), (0, 0)), mode='reflect')
            if j+binsAround+1 >= fullmat.shape[1]:
                pad_v = mat_height-submat.shape[1]
                submat = np.pad(submat,((0, 0), (0, pad_v)), mode='reflect')
            hic_freq = (np.log1p(fullmat[i,j])-min_mat)/(max_mat-min_mat)
            submat = submat.reshape((mat_width,mat_height,1))
            submat = np.log1p(submat)
            submat = (submat-np.nanmin(submat))/(np.nanmax(submat)-np.nanmin(submat))
            submat[np.isnan(submat)]=0
            
            pairi= int(i-beg)
            pairj= int(j-beg)
            pairs.append((pairi,pairj))
            submatrices.append(np.asarray(submat))
            pair_dists.append([abs(i-j),hic_freq])

    submatrices = np.asarray(submatrices)
    pair_dists = np.asarray(pair_dists)
    
    (pair_predictions_K,
     pair_predictions_loc,
     pair_predictions_scale) = tf_model.predict([pair_dists, submatrices])
    
    pair_predictions=[]
    valid_pairs=[]
    invalid_pairs=[]
    for i in range(len(pair_predictions_K)):
        pred_list_i = [pair_predictions_K[i],
                       pair_predictions_loc[i],
                       pair_predictions_scale[i]]
        if all(pi >= 0.2 for pi in pred_list_i[1:]):
            valid_pairs.append(pairs[i])
            # correct too high K values
            pred_list_i[0] = 1.5 if pred_list_i[0] > 1.5 else pred_list_i[0]
            pred_list_i[0] = 0.5 if pred_list_i[0] < 0.5 else pred_list_i[0]
            pair_predictions.append(pred_list_i)
        else:
            invalid_pairs.append(pairs[i])

    nbr_short_pairs = len(valid_pairs)
    bins_group = distance_mats  = None
    if tf_model_far:
        if globals.VERBOSE > 0:
            print('Generating restraints at %d resolution\n'%(resolution_far))
        bins_group, distance_mats, far_positions = generate_parent_locations(hic_data, region_zeros,
                                                              beg, end,
                                                              dist, resolution, tf_model_far,
                                                              resolution_far,
                                                              0,n_cpus=n_cpus,
                                                              n_models=n_models)
    if globals.VERBOSE > 0:
        print('Generating restraints at %d resolution\n'%(resolution))
    
    
    
    ProbRestraintsList = ProbabilityBasedRestraintsList(nloci, globals.RADIUS,
                                                        globals.SCALE, CONFIG['kforce'],
                                                        n_models, dist, valid_pairs,
                                                        invalid_pairs, pair_predictions,
                                                        resolution, region_zeros=region_zeros,
                                                        short_restr_dist=short_restr_dist,
                                                        far_restr_dist=far_restr_dist,
                                                        chromosomes=coords, n_cpus=n_cpus,
                                                        random_seed=start, bins_group=bins_group,
                                                        distance_mats=distance_mats,
                                                        save_restraints_folder=save_restraints_folder)
    
    if save_restraints_folder:
        paramsfile = path.join(save_restraints_folder,'_tmp_common_params.pickle')
        tmp_params = open(paramsfile, 'wb')
        dump(nloci, tmp_params)
        dump(CONFIG, tmp_params)
        tmp_params.close()
        return
    
    globals.LOCI  = list(range(nloci))
    models, bad_models = multi_process_model_generation(
        n_cpus, n_models, n_keep, keep_all, ProbRestraintsList,
        use_HiC=use_HiC, use_confining_environment=use_confining_environment,
        use_excluded_volume=False,
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
        for i, m in enumerate(list(models.values()) + list(bad_models.values())):
            m['index'] = i
            m['description'] = description
    except AttributeError:  # case we are doing optimization
        description = None
        for i, m in enumerate(list(models.values()) + list(bad_models.values())):
            m['index'] = i
    if outfile:
        if exists(outfile):
            old_models, old_bad_models = load(open(outfile, 'rb'))
        else:
            old_models, old_bad_models = {}, {}
        models.update(old_models)
        bad_models.update(old_bad_models)
        out = open(outfile, 'wb')
        dump((models, bad_models), out)
        out.close()
    else:
        probrestraints = None
        if include_restraints_in_models:
            probrestraints = ProbRestraintsList._get_restraints(globals.SCALE)
        if sparse_avail and issparse(fullmat):
            fullmat = fullmat.tocsr()
            fullmat.data = np.log10(fullmat.data)
            original_data = fullmat
        else:
            original_data = np.where(fullmat != 0, np.log10(fullmat), 0)
        original_data = original_data[beg:end+1,beg:end+1]
        return StructuralModels(
            nloci, models, bad_models, resolution, original_data=original_data,
            config=CONFIG, experiment=experiment, zeros=zeros_models,
            restraints=probrestraints, description=description)

def generate_parent_locations(hic_data, region_zeros, beg, end, dist, resolution,
                              tf_model_far, resolution_far, far_distance,
                              include_bins=[], end_distance=100e6, n_models=5000,
                              coords=None, n_cpus=1, use_confining_environment=False):
    
    reso_ratio = resolution_far / resolution
    if sparse_avail and issparse(hic_data):
        if not isinstance(hic_data, csr_matrix):
            raise ValueError('Matrix given must be of CSR format.')
        fullmat = hic_data[beg:end+1,beg:end+1].copy()
        fullmat=fullmat.todense()
    else:
        fullmat = np.asarray(hic_data)[beg:end+1,beg:end+1]
        
    fullmat_far = rebin(hic_data, resolution, resolution_far)
    np.fill_diagonal(fullmat_far,np.nan)    
    fullmat_far = np.log1p(fullmat_far)
    fullmat_far = (fullmat_far-np.nanmin(fullmat_far))/(np.nanmax(fullmat_far)-np.nanmin(fullmat_far))
    fullmat_far[np.isnan(fullmat_far)]=0

    nloci = int((end - beg + 1))
    beg_far = round((beg * resolution)/resolution_far)
    end_far = round((end * resolution)/resolution_far)

    include_far_bins = [round((i * resolution)/resolution_far) for i in include_bins]
    pairs = []
    bins_group = []
    pair_freqs = []
    for i in range(beg_far, end_far+1):
        pairi= (max(0, int(i*reso_ratio)-beg),
                min(end-beg+1, int((i+1)*reso_ratio)-beg))
        bins_group.append(pairi)
        if all(p in region_zeros for p in range(pairi[0],pairi[1]+1)):
            continue
        for j in range(i+1, end_far+1):
            pairj= (max(0, int(j*reso_ratio)-beg),
                    min(end-beg, int((j+1)*reso_ratio)-beg))
            if all(p in region_zeros for p in range(pairj[0],pairj[1]+1)):
                continue
            if not (i in include_far_bins or j in include_far_bins)  \
                and (abs(i-j)*resolution_far < far_distance \
                or abs(i-j)*resolution_far > end_distance):
                continue
            hic_freq = fullmat_far[i,j]
            if hic_freq==0:
                continue
            pairi_far= int(i-beg_far)
            pairj_far= int(j-beg_far)
            pairs.append((pairi_far,pairj_far))
            pair_freqs.append([abs(i-j),hic_freq])
            
    pair_freqs = np.asarray(pair_freqs)
    
    if len(pair_freqs) > 0:
        (pair_predictions_far_K,
         pair_predictions_far_loc,
         pair_predictions_far_scale) = tf_model_far.predict(pair_freqs)
        pair_predictions_far=[]
        valid_pairs=[]
        for i in range(len(pair_predictions_far_K)):
            pred_list_i = [float(pair_predictions_far_K[i]),
                           float(pair_predictions_far_loc[i]),
                           float(pair_predictions_far_scale[i])]
            if all(pi >= 0.5 for pi in pred_list_i[1:]):    
                valid_pairs.append(pairs[i])
                pred_list_i[0] = 1.5 if pred_list_i[0] > 1.5 else pred_list_i[0]
                pred_list_i[0] = 0.5 if pred_list_i[0] < 0.5 else pred_list_i[0]
                pair_predictions_far.append(pred_list_i)
    else:
        return None
    
    nloci = int((end_far - beg_far + 1))
    globals.LOCI  = list(range(nloci))
    ProbRestraintsList = ProbabilityBasedRestraintsList(nloci, globals.RADIUS*reso_ratio,
                                                    globals.SCALE, CONFIG['kforce'],
                                                    n_models*2, dist, valid_pairs,
                                                    [], pair_predictions_far,
                                                    resolution_far, short_restr_dist=0,
                                                    far_restr_dist=0,
                                                    chromosomes=coords, n_cpus=n_cpus,
                                                    random_seed=globals.START,
                                                    perc_restraints=0.6)
    if not CONFIG['container']['shape'] and len(ProbRestraintsList.bounding_boxes) > 0:
        CONFIG['container'] = {'shape' : 'cylinder',
                               'radius': (max(ProbRestraintsList.bounding_boxes)/2) / globals.SCALE,
                               'height': 0,
                               'cforce': CONFIG['kforce']*1000}
        
        use_confining_environment=False
    
    globals.LOCI  = list(range(nloci))
    models, bad_models = multi_process_model_generation(
        n_cpus, n_models*2, n_models, False, ProbRestraintsList,
        use_HiC=True, use_confining_environment=use_confining_environment,
        use_excluded_volume=False)

    positions = []
    for m in list(models.values()):
        positions.append(list(zip(m['x'],m['y'],m['z'])))
    
    xyzs = np.array(positions)
    
    distance_mats=np.array(list(map(squareform,list(map(pdist,xyzs)))))
    
    return bins_group, distance_mats, positions   

def multi_process_model_generation(n_cpus, n_models, n_keep, keep_all, ProbRestraintsList,
                                   bounding_restraints=None, use_HiC=True, use_confining_environment=True,
                                   use_excluded_volume=True, single_particle_restraints=None):
    """
    Parallelize the
    :func:`pytadbit.modelling.imp_model.StructuralModels.generate_IMPmodel`.

    :param n_cpus: number of CPUs to use
    :param n_models: number of models to generate
    """
    
    def update_progress(res):
        if globals.VERBOSE > 0:
            print('Generated model %s'%(res['rand_init']))
    pool = mu.Pool(n_cpus)
    jobs = {}
    for rand_init in range(globals.START, n_models + globals.START):
        ProbRestraints =  ProbabilityBasedRestraints(
                ProbRestraintsList.get_probability_restraints(rand_init-globals.START))
        jobs[rand_init] = pool.apply_async(generate_IMPmodel,
                                           args=(rand_init, ProbRestraints,
                                                 use_HiC,use_confining_environment, use_excluded_volume,
                                                 single_particle_restraints),
                                           callback=update_progress)

    pool.close()
    pool.join()

    results = []
    for rand_init in range(globals.START, n_models + globals.START):
        
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

def rebin(a, resolution, resolution_far):

    M, N = a.shape
    ratio = resolution_far / resolution
    m, n = int(M//ratio)+1, int(N//ratio)+1
    b = np.zeros(shape=(m,n))
    if sparse_avail and issparse(a):
        a_s = triu(a, k=1).tolil()
        for i in range(m):
            subi_start = int(max(0,i*ratio-(ratio/2)))
            subi_end = int(min(M,i*ratio+(ratio/2)))
            for j in range(i+1,n):
                subj_start = int(max(0,j*ratio-(ratio/2)))
                subj_end = int(min(N,j*ratio+(ratio/2)))
                b[i,j] = np.nanmean(a_s[subi_start:subi_end,
                                        subj_start:subj_end].todense())
                b[j,i] = b[i,j]
    else:
        a_t = np.triu(a, k=1)
        for i in range(m):
            subi_start = int(max(0,i*ratio-(ratio/2)))
            subi_end = int(min(M,i*ratio+(ratio/2)))
            for j in range(i+1,n):
                subj_start = int(max(0,j*ratio-(ratio/2)))
                subj_end = int(min(N,j*ratio+(ratio/2)))
                b[i,j] = np.nanmean(a_t[subi_start:subi_end,
                                        subj_start:subj_end])
                b[j,i] = b[i,j]

    return b
