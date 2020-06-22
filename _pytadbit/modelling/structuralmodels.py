"""
19 Jul 2013
"""
from cPickle                          import load, dump, HIGHEST_PROTOCOL
from subprocess                       import Popen, PIPE
from math                             import acos, degrees, pi, sqrt
from warnings                         import warn
from string                           import uppercase as uc, lowercase as lc
from random                           import random, randint
from os.path                          import exists, isdir
from os                               import system 
from itertools                        import combinations
from uuid                             import uuid5, UUID
from hashlib                          import md5

from numpy                            import exp as np_exp
from numpy                            import median as np_median
from numpy                            import mean as np_mean
from numpy                            import std as np_std, log2
from numpy                            import array, cross, dot, ma, isnan
from numpy                            import histogram, linspace
from numpy                            import nanmin, nanmax
from numpy                            import seterr
from numpy.linalg                     import norm

from scipy.optimize                   import curve_fit
from scipy.stats                      import spearmanr, pearsonr, chisquare
from scipy.stats                      import linregress
from scipy.stats                      import normaltest, norm as sc_norm
from scipy.cluster.hierarchy          import linkage, fcluster

from pytadbit                         import get_dependencies_version
from pytadbit.utils.three_dim_stats   import calc_consistency, mass_center
from pytadbit.utils.three_dim_stats   import dihedral, calc_eqv_rmsd
from pytadbit.utils.three_dim_stats   import get_center_of_mass, distance
from pytadbit.utils.tadmaths          import calinski_harabasz, nozero_log_list
from pytadbit.utils.tadmaths          import mean_none
from pytadbit.utils.extraviews        import plot_3d_model, setup_plot
from pytadbit.utils.extraviews        import chimera_view, tadbit_savefig
from pytadbit.utils.extraviews        import augmented_dendrogram, plot_hist_box
from pytadbit.utils.extraviews        import tad_coloring
from pytadbit.utils.extraviews        import tad_border_coloring
from pytadbit.utils.extraviews        import color_residues
from pytadbit.modelling.impmodel      import IMPmodel
from pytadbit.modelling.lammpsmodel   import LAMMPSmodel
from pytadbit.centroid                import centroid_wrapper
from pytadbit.aligner3d               import aligner3d_wrapper
from pytadbit.squared_distance_matrix import squared_distance_matrix_calculation_wrapper

try:
    from matplotlib import pyplot as plt
    from matplotlib.cm import jet, bwr
except ImportError:
    warn('matplotlib not found\n')


def R2_vs_L(L, P):
    """
    Calculates the persistence length (Lp) of given section of the model.
    Persistence length is calculated according to [Bystricky2004]_ :

    .. math::

    <R^2> = 2 \\times Lp^2 \\times (\\frac{Lc}{Lp} - 1 + e^{\\frac{-Lc}{Lp}})

    with the contour length as :math:`Lc = \\frac{d}{c}` where :math:`d` is
    the genomic dstance in bp and :math:`c` the linear mass density of the
    chromatin (in bp/nm).

    :returns: persistence length, or 2 times the Kuhn length
    """
    return 2.0 * P * ( L - P * ( 1.0 - np_exp( - L / P ) ) )

def load_structuralmodels(path_f):
    """
    Loads :class:`pytadbit.modelling.structuralmodels.StructuralModels` from a file
    (generated with
    :class:`pytadbit.modelling.structuralmodels.StructuralModels.save_models`).

    :param path: to the pickled StructuralModels object.

    :returns: a :class:`pytadbit.modelling.imp_model.StructuralModels`.
    """
    svd = load(open(path_f))
    try:
        return StructuralModels(
            nloci=svd['nloci'], models=svd['models'], bad_models=svd['bad_models'],
            resolution=svd['resolution'], original_data=svd['original_data'],
            clusters=svd['clusters'], config=svd['config'], zscores=svd['zscore'],
            zeros=svd['zeros'], restraints=svd.get('restraints', None),
            description=svd.get('description', None), stages=svd.get('stages', None),
            models_per_step=svd.get('models_per_step', 0))
    except KeyError:  # old version
        return StructuralModels(
            nloci=svd['nloci'], models=svd['models'], bad_models=svd['bad_models'],
            resolution=svd['resolution'], original_data=svd['original_data'],
            clusters=svd['clusters'], config=svd['config'], zscores=svd['zscore'],
            restraints=svd.get('restraints', None), stages=svd.get('stages', None),
            models_per_step=svd.get('models_per_step', 0))


class StructuralModels(object):
    """
    This class contains three-dimensional models generated from a single Hi-C
    data. They can be reached either by their index (integer representing their
    rank according to objective function value), or by their IMP random intial
    number (as string).

    :param nloci: number of particles in the selected region
    :param models: a dictionary containing the generated
       :class:`pytadbit.modelling.impmodel.IMPmodel` to be used as 'best models'
    :param bad_models: a dictionary of :class:`pytadbit.modelling.impmodel.IMPmodel`,
       these model will not be used, just stored in case the set of
       'best models' needs to be extended later-on (
       :func:`pytadbit.modelling.structuralmodels.StructuralModels.define_best_models`
       ).
    :param resolution: number of nucleotides per Hi-C bin. This will be the
       number of nucleotides in each model's particle.
    :param None original_data: a list of list (equivalent to a square matrix) of
       the normalized Hi-C data, used to build this list of models.
    :param None clusters: a dictionary of type
       :class:`pytadbit.modelling.structuralmodels.ClusterOfModels`
    :param None config: a dictionary containing the parameter to be used for the
       generation of three dimensional models.

    """

    def __init__(self, nloci, models, bad_models, resolution,
                 original_data=None, zscores=None, clusters=None,
                 config=None, experiment=None, zeros=None, restraints=None,
                 description=None, stages=None, models_per_step=0):

        self.__models       = models
        self._bad_models    = bad_models
        self.nloci          = nloci
        self.clusters       = clusters or ClusterOfModels()
        self.resolution     = float(resolution)
        self._original_data = original_data  # only used for correlation
        self._zscores       = zscores        # only used for plotting
        self._zeros         = zeros or {}    # filtered out columns
        self._config        = config or {}
        self.experiment     = experiment
        self._restraints    = restraints
        self.description    = description
        self.stages         = stages or {}
        self.models_per_step = models_per_step

    def __getitem__(self, nam):
        if isinstance(nam, str):
            for m in self.__models:
                if self.__models[m]['rand_init'] == nam:
                    return self.__models[m]
            else:
                raise KeyError('Model %s not found\n' % (nam))
        try:
            return self.__models[nam]
        except TypeError:
            for i, key in self.__models:
                if nam == i:
                    return self.__models[key]
            raise KeyError('Model %s not found\n' % (i))

    def __iter__(self):
        for m in self.__models:
            yield self.__models[m]

    def __len__(self):
        return len(self.__models)

    def __repr__(self):
        model_type = 'IMP'
        if isinstance(self.__models[0],LAMMPSmodel):
            model_type = 'LAMMPS'
            
        return ('StructuralModels with %s models of %s particles\n' +
                '   (objective function range: %s - %s)\n' +
                '   (corresponding to the best models out of %s models).\n' +
                '  %s modeling used this parameters:\n' +
                '%s\n' +
                '  Models where clustered into %s clusters') % (
            len(self.__models),
            self.nloci,
            int(self.__models[0]['objfun']),
            int(self.__models[len(self.__models) - 1]['objfun']),
            len(self.__models) + len(self._bad_models),
            model_type,
            '\n'.join(['   - %-12s: %s' % (k, v)
                       for k, v in self._config.iteritems()]),
            len(self.clusters))

    def _extend_models(self, models, stages=None):
        """
        add new models to structural models
        """
        if isinstance(models, StructuralModels):
            stages = models.stages
            models = models._StructuralModels__models
        if not isinstance(models, dict):
            warn('ERROR: models has to be a StructuralModels object '
                     'or a dictionary')
            return
        nbest = len(self.__models)
        nall  = len(self.__models) + len(self._bad_models)
        self.define_best_models(nall)
        ids = set(self.__models[m]['rand_init'] for m in self.__models)
        for m in models.keys():
            if models[m]['rand_init'] in ids:
                warn('WARNING: found model with same random seed number, '
                     'CHANGING rand_init')
                models[m]['rand_init'] += '-' + str(randint(1,10000))
                #del(models[m]) # why?
        new_models = {}
        if len(self.stages) > 1:
            if stages and len(set(stages.keys()) & set(self.stages.keys())) == len(self.stages):
                new_stages = {}
                offset= 0
                for stg in self.stages:
                    new_stages[stg] = []
                    stage_models = [self.__models[m] for m in set(self.stages[stg])] + [models[m] for m in set(stages[stg])] 
                    for i, m in enumerate(stage_models):
                        new_models[i+offset] = m
                        new_stages[stg].append(i+offset)
                    offset += len(stage_models)
            else:
                warn('WARNING: we need the same number of stages '
                     'to extend the structural models')
                return
            self.stages = new_stages
            self.__models = new_models
        else:
            for i, m in enumerate(sorted(models.values() + self.__models.values(),
                                         key=lambda x: x['objfun'])):
                new_models[i] = m
                new_models[i]['index'] = i
            self.__models = new_models
            # keep the same number of best models if best models were not all
            if len(self._bad_models) == 0:
                nbest = len(new_models) 
            self.define_best_models(nbest)

    def align_models(self, models=None, cluster=None, in_place=False,
                     reference_model=None, **kwargs):
        """
        Three-dimensional aligner for structural models.

        :param None models: if None (default) the average model will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the average model only for the models in the
           cluster number 'cluster'
        :param False in_place: if True, the result of the alignment will REPLACE
           the coordinates of each model. Default is too yield new coordinates of
           each model
        :param None reference_model: align given model to reference model
        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        ref_model = models[0] if reference_model is None else reference_model
        firstx, firsty, firstz = (self[ref_model]['x'],
                                  self[ref_model]['y'],
                                  self[ref_model]['z'])
        aligned = []
        for sec in models[1 if reference_model is None else 0:]:
            coords = aligner3d_wrapper(firstx, firsty, firstz,
                                       self[sec]['x'],
                                       self[sec]['y'],
                                       self[sec]['z'],
                                       self._zeros,
                                       self.nloci)
            if in_place:
                self[sec]['x'], self[sec]['y'], self[sec]['z'] = coords
            else:
                aligned.append(coords)

        if in_place:
            mass_center(self[ref_model]['x'], self[ref_model]['y'],
                        self[ref_model]['z'], self._zeros)
        else:
            x, y, z = (self[ref_model]['x'][:], self[ref_model]['y'][:],
                       self[ref_model]['z'][:])
            mass_center(x, y, z, self._zeros)
            aligned.insert(ref_model, (x, y, z))
            return aligned

    def fetch_model_by_rand_init(self, rand_init, all_models=False):
        """
        Models are stored according to their objective function value (first
        best), but in order to reproduce a model, we need its initial random
        number. This method helps to fetch the model corresponding to a given
        initial random number stored under
        StructuralModels.models[N]['rand_init'].

        :param rand_init: the wanted rand_init number.
        :param False all_models: whether or not to use 'bad' models

        :returns: index of 3d model
        """
        for i, model in enumerate(self):
            if model['rand_init'] == rand_init:
                return i
        if all_models:
            for m in self._bad_models:
                if self._bad_models[m]['rand_init'] == rand_init:
                    return m
        raise IndexError(('Model with initial random number: %s, not found\n' +
                          '') % (rand_init))

    def centroid_model(self, models=None, cluster=None, verbose=False):
        """
        Estimates and returns the centroid model of a given group of models.

        :param None models: if None (default) the centroid model will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the centroid model only for the models in the
           cluster number 'cluster'
        :param False verbose: prints the distance of each model to average model
           (in stderr)

        :returns: the centroid model of a given group of models (the most model
           representative)
        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        # remove particles with zeros from calculation
        x = []
        y = []
        z = []
        for model in xrange(len(models)):
            x.append([self[model]['x'][i] for i in xrange(self.nloci)
                      if self._zeros[i]])
            y.append([self[model]['y'][i] for i in xrange(self.nloci)
                      if self._zeros[i]])
            z.append([self[model]['z'][i] for i in xrange(self.nloci)
                      if self._zeros[i]])
        zeros = tuple([True for _ in xrange(len(x[0]))])
        idx = centroid_wrapper(x, y, z, zeros, len(x[0]), len(models),
                               int(verbose), 0)
        return models[idx]

    def average_model(self, models=None, cluster=None, verbose=False):
        """
        Builds and returns an average model representing a given group of models

        :param None models: if None (default) the average model will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the average model only for the models in the
           cluster number 'cluster'
        :param False verbose: prints the distance of each model to average model
           (in stderr)

        :returns: the average model of a given group of models (a new and
           ARTIFICIAL model)

        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        # remove particles with zeros from calculation
        x = []
        y = []
        z = []
        for model in xrange(len(models)):
            x.append([self[model]['x'][i] for i in xrange(self.nloci)])
            y.append([self[model]['y'][i] for i in xrange(self.nloci)])
            z.append([self[model]['z'][i] for i in xrange(self.nloci)])
        idx = centroid_wrapper(x, y, z, self._zeros, len(x[0]), len(models),
                               int(verbose), 1)
        avgmodel = IMPmodel((('x', idx[0]), ('y', idx[1]), ('z', idx[2]),
                             ('rand_init', 'avg'), ('objfun', None),
                             ('radius', float(self.resolution *
                                              self._config['scale']) / 2)))
        return avgmodel

    def cluster_models(self, fact=0.75, dcutoff=None, method='mcl',
                       mcl_bin='mcl', tmp_file=None, verbose=True, n_cpus=1,
                       mclargs=None, external=False, what='score'):
        """
        This function performs a clustering analysis of the generated models
        based on structural comparison. The result will be stored in
        StructuralModels.clusters

        Clustering is done according to a score of pairwise comparison
        calculated as:

       .. math::

         score_i = eqvs_i \\times \\frac{dRMSD_i / max(dRMSD)}
                                         {RMSD_i / max(RMSD)}

       where :math:`eqvs_i` is the number of equivalent position for the ith
       pairwise model comparison.


        :param 0.75 fact: factor to define the percentage of equivalent
           positions to be considered in the clustering
        :param None dcutoff: distance threshold (nm) to determine if two
           particles are in contact, default is 1.5 times resolution times scale
        :param 'mcl' method: clustering method to use, which can be either
           'mcl' or 'ward'. MCL method is recommended. WARD method uses a scipy
           implementation of this hierarchical clustering, and selects the best
           number of clusters using the
           :func:`pytadbit.utils.tadmaths.calinski_harabasz` function.
        :param 'mcl' mcl_bin: path to the mcl executable file, in case of the
           'mcl is not in the PATH' warning message
        :param None tmp_file: path to a temporary file created during
           the clustering computation. Default will be created in /tmp/ folder
        :param True verbose: same as print StructuralModels.clusters
        :param 1 n_cpus: number of cpus to use in MCL clustering
        :param mclargs: list with any other command line argument to be passed
           to mcl (i.e,: mclargs=['-pi', '10', '-I', '2.0'])
        :param False external: if True returns the cluster found instead of
           storing it as StructuralModels.clusters
        :param 'score' what: Statistic used for clustering. Can be one of
           'score', 'rmsd', 'drmsd' or 'eqv'.

        """
        tmp_file = '/tmp/tadbit_tmp_%s.txt' % (
            ''.join([(uc + lc)[int(random() * 52)] for _ in xrange(4)]))
        if not dcutoff:
            dcutoff = int(1.5 * self.resolution * self._config['scale'])
            #dcutoff = int(1.5) # * self.resolution * self._config['scale'])
        scores = calc_eqv_rmsd(self.__models, self.nloci, self._zeros, dcutoff,
                               what=what, normed=True)
        from distutils.spawn import find_executable
        if not find_executable(mcl_bin):
            print('\nWARNING: MCL not found in path using WARD clustering\n')
            method = 'ward'
        # Initialize cluster definition of models:
        for model in self:
            model['cluster'] = 'Singleton'
        new_singles = 0
        if method == 'ward':

            matrix = [[0.0 for _ in xrange(len(self))]
                      for _ in xrange(len(self))]
            for (i, j), score in scores.iteritems():
                matrix[i][j] = score if score > fact * self.nloci else 0.0
            clust = linkage(matrix, method='ward')
            # score each possible cut in hierarchical clustering
            solutions = {}
            for k in clust[:, 2]:
                clusters = ClusterOfModels()
                [clusters.setdefault(j, []).append(i) for i, j in
                 enumerate(fcluster(clust, k, criterion='distance'))]
                solutions[k] = {'out': clusters}
                solutions[k]['score'] = calinski_harabasz(scores, clusters)
            # take best cluster according to calinski_harabasz score
            clusters = [solutions[s] for s in sorted(
                solutions, key=lambda x: solutions[x]['score'])
                if solutions[s]['score'] > 0][-1]['out']
            # sort clusters, the more populated, the first.
            clusters = dict([(i + 1, j) for i, j in
                             enumerate(sorted(clusters.values(),
                                              key=len, reverse=True))])
            if external:
                return clusters
            self.clusters = ClusterOfModels()
            for cluster in clusters:
                self.clusters[cluster] = []
                for model in clusters[cluster]:
                    self[model]['cluster'] = cluster
                    self.clusters[cluster].append(str(self[model]['rand_init']))
                self.clusters[cluster].sort(
                    key=lambda x: self[str(x)]['objfun'])
        else:
            out_f = open(tmp_file, 'w')
            uniqs = list(set([tuple(sorted((m1, m2))) for m1, m2 in scores]))
            cut = fact * (self.nloci - self._zeros.count(False))
            for md1, md2 in uniqs:
                score = scores[(md1, md2)]
                if score >= cut:
                    out_f.write('model_%s\tmodel_%s\t%s\n' % (md1, md2, score))
            out_f.close()
            Popen('%s %s --abc -te %s -V all -o %s.mcl %s' % (
                mcl_bin, tmp_file, n_cpus, tmp_file, ' '.join(
                    mclargs or [])), stdout=PIPE, stderr=PIPE,
                  shell=True).communicate()
            clusters = ClusterOfModels()
            if not exists(tmp_file + '.mcl'):
                raise Exception('Problem with clustering, try increasing ' +
                                '"dcutoff", now: %s\n' % (dcutoff))
            for cluster, line in enumerate(open(tmp_file + '.mcl')):
                models = line.split()
                if len(models) == 1:
                    new_singles += 1
                else:
                    clusters[cluster + 1] = []
                    for model in models:
                        model = int(model.split('_')[1])
                        if not external:
                            self[model]['cluster'] = cluster + 1
                        clusters[cluster + 1].append(
                            str(self[model]['rand_init']))
                    clusters[cluster + 1].sort(
                        key=lambda x: self[str(x)]['objfun'])
            if external:
                return clusters
            self.clusters = clusters
        if verbose:
            singletons = len([1 for m in self
                              if m['cluster'] == 'Singleton'])
            print ('Number of singletons excluded from clustering: %s (total' +
                   ' singletons: %s)') % (singletons - new_singles, singletons)
            print self.clusters

    def _build_distance_matrix(self, n_best_clusters):
        """
        """
        clusters = sorted(self.clusters.keys())[:n_best_clusters]
        matrix = [[0.0 for _ in clusters] for _ in clusters]
        clust_count = dict([(c, len([m for m in self if m['cluster'] == c]))
                            for c in clusters])
        objfun = dict([(c, [m for m in self if m['cluster'] == c][0]['objfun'])
                       for c in clusters])
        for i, cl1 in enumerate(clusters):
            md1 = md2 = None
            # find model with lowest energy for each cluster
            for md1 in self:
                if md1['cluster'] == cl1:
                    # the first one found is the best :)
                    break
            for j, cl2 in enumerate(clusters[i + 1:]):
                # find model with lowest energy for each cluster
                for md2 in self:
                    if md2['cluster'] == cl2:
                        # the first one found is the best :)
                        break
                matrix[i][j + i + 1] = calc_eqv_rmsd({0: md1, 1: md2},
                                                     self.nloci,
                                                     self._zeros, one=True)
        return clust_count, objfun, matrix

    def cluster_analysis_dendrogram(self, n_best_clusters=None, color=False,
                                    axe=None, savefig=None, **kwargs):
        """
        Representation of the clustering results. The length of the leaves if
        proportional to the final objective function value of each model. The
        branch widths are proportional to the number of models in a given
        cluster (or group of clusters, if it is an internal branch).

        :param None n_best_clusters: number of clusters to represent (by
           default all clusters will be shown)
        :param False color: color the dendrogram based on the significance of
           the clustering (basically it depends of the internal branch lengths)
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param 10.0 width_factor: multiplicator for the width of the line
           representing the number of models in a given cluster.
        :param 8 fontsize: size of the smallest font represented in the plot
        :param (8,8) figsize: a tuple of width and height, to set the size of
           the plot.
        """

        if not self.clusters:
            warn('WARNING: no clusters found, clustering with default' +
                 ' parameters\n')
            self.cluster_models()
        if not n_best_clusters:
            n_best_clusters = len(self.clusters)
        if n_best_clusters <= 1:
            warn("Need at least 2 clusters to display...")
            return None
        clust_count, objfun, matrix = self._build_distance_matrix(n_best_clusters)
        z = linkage(matrix)
        minnrj = min(objfun.values()) - 1
        maxnrj = max(objfun.values()) - 1
        val = (maxnrj - minnrj)
        maxz = max([i[2] for i in z])
        for i in z:
            i[2] = i[2] / maxz * val

        dads = {}
        i = max(clust_count)
        for a, b, _, _ in z:
            i += 1
            clust_count[i] = clust_count[a + 1] + clust_count[b + 1]
            dads[a + 1] = i
            dads[b + 1] = i

        d = augmented_dendrogram(clust_count, dads, objfun, color,
                                 axe, savefig, z, **kwargs)
        return d

    def get_contact_matrix(self, models=None, cluster=None, 
                           stage=None, cutoff=None,
                           distance=False, show_bad_columns=True):
        """
        Returns a matrix with the number of interactions observed below a given
        cutoff distance.

        :param None models: if None (default) the contact matrix will be
           computed using all the models. A list of numbers corresponding to a
           given set of models can be passed
        :param None cluster: compute the contact matrix only for the models in
           the cluster number 'cluster'
        :param None stage: compute the contact matrix only for the models in
            stage number 'stage'
        :param None cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not, default is 2 times resolution, times scale.
           Cutoff can also be a list of values, in wich case the returned object
           will be a dictionnary of matrices (keys being square cutoffs)
        :param False distance: returns the distance matrix of all_angles against
           all_angles particles instead of a contact_map matrix using the cutoff
        :param True show_bad_columns: show bad columns in contact map

        :returns: matrix frequency of interaction
        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        elif stage > -1 and stage in self.stages:
            models = [m for m in self.stages[stage]]
        else:
            models = [m for m in self.__models]
        if not cutoff:
            cutoff = [int(2 * self.resolution * self._config['scale'])]
            #cutoff = [int(2)] # * self.resolution * self._config['scale'])]
        cutoff_list = True
        if not isinstance(cutoff, list):
            cutoff = [cutoff]
            cutoff_list = False
        cutoff.sort(reverse=True)
        cutoff = [c**2 for c in cutoff]
        matrix = dict([(c, [[0. for _ in xrange(self.nloci)]
                            for _ in xrange(self.nloci)]) for c in cutoff])
        # remove (or not) interactions from bad columns
        if show_bad_columns:
            wloci = [i for i in xrange(self.nloci) if self._zeros[i]]
        else:
            wloci = [i for i in xrange(self.nloci)]
        models = [self[mdl] for mdl in models]

        frac = 1.0 / len(models)
        #print "#Frac",frac

        all_matrix = []
        for model in models:
            #print model
            squared_distance_matrix = squared_distance_matrix_calculation_wrapper(
                model['x'], model['y'], model['z'], self.nloci)

            #print model, len(x), len(y), len(z)
            for c in cutoff:
                #print "#Cutoff",c
                for i, j in combinations(wloci, 2):
                    if squared_distance_matrix[i][j] <= c:
                        matrix[c][i][j] += frac  # * 100
                        matrix[c][j][i] += frac  # * 100

        if cutoff_list:
            return matrix
        return matrix.values()[0]

    def define_best_models(self, nbest):
        """
        Defines the number of top models (based on the objective function) to
        keep. If keep_all is set to True in
        :func:`pytadbit.modelling.imp_model.generate_3d_models` or in
        :func:`pytadbit.experiment.Experiment.model_region`, then the full set
        of models (n_models parameter) will be used, otherwise only the n_keep
        models will be available.

        :param nbest: number of top models to keep (usually 20% of the
           generated models).
        """
        tmp_models = self.__models
        tmp_models.update(self._bad_models)
        self.__models = dict([(i, tmp_models[i]) for i in xrange(nbest)])
        self._bad_models = dict([(i, tmp_models[i]) for i in
                                 xrange(nbest, len(tmp_models))])

    def deconvolve(self, fact=0.75, dcutoff=None, method='mcl',
                   mcl_bin='mcl', tmp_file=None, verbose=True, n_cpus=1,
                   mclargs=None, what='score', n_best_clusters=10,
                   savefig=None, represent_models=False, figsize=(11, 11),
                   clusters=None, **kwargs):
        """
        This function performs a deconvolution analysis of a given froup of models.
        It first clusters models based on structural comparison (dRMSD), and
        then, performs a differential contact map between each possible pair
        of cluster.

        .. note::

          Clusters defined here are different from the one defined when using
          :func:`pytadbit.modelling.structuralmodels.StructuralModels.cluster_models`.
          They are also not stored into StructuralModels.clusters

        :param 0.75 fact: factor to define the percentage of equivalent
           positions to be considered in the clustering
        :param None dcutoff: distance threshold (nm) to determine if two
           particles are in contact, default is 1.5 times resolution times scale
        :param 'mcl' method: clustering method to use, which can be either
           'mcl' or 'ward'. MCL method is recommended. WARD method uses a scipy
           implementation of this hierarchical clustering, and selects the best
           number of clusters using the
           :func:`pytadbit.utils.tadmaths.calinski_harabasz` function.
        :param 'mcl' mcl_bin: path to the mcl executable file, in case of the
           'mcl is not in the PATH' warning message
        :param None tmp_file: path to a temporary file created during
           the clustering computation. Default will be created in /tmp/ folder
        :param True verbose: same as print StructuralModels.clusters
        :param 1 n_cpus: number of cpus to use in MCL clustering
        :param mclargs: list with any other command line argument to be passed
           to mcl (i.e,: mclargs=['-pi', '10', '-I', '2.0'])
        :param 10 n_best_clusters: number of clusters to represent
        :param None clusters: provide clusters as a dictionary with keys=cluster
           number, or name, and values list of model numbers.
        :param False represent_models: To generate an interactive visualization
           of a representative model for each cluster. Representative model
           depends on the value passed to this option, it can be either
           'centroid' or 'best' (this last standing for the model with lowest
           IMP objective function value).
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param (11,11) figsize: dimension of the plot
        """
        #fact /= self.nloci
        if not dcutoff:
            dcutoff = int(1.5 * self.resolution * self._config['scale'])
        if not clusters:
            clusters = self.cluster_models(fact=fact, dcutoff=dcutoff,
                                           mcl_bin=mcl_bin,
                                           method=method, tmp_file=tmp_file,
                                           n_cpus=n_cpus, mclargs=mclargs,
                                           external=True, what=what)
        if len(clusters) <= 1:
            raise Exception('ERROR: did not found clusters to be compared ' +
                            '(try different clustering parameters).\n')
        if verbose:
            print clusters
        n_best_clusters = min(len(clusters), n_best_clusters)
        add = 1 if represent_models else 0
        fig, axes = plt.subplots(n_best_clusters - 1 + add,
                                 n_best_clusters - 1 + add,
                                 sharex=True, sharey=True, figsize=figsize)
        # pre-calculate contact-matrices
        cmatrices = [self.get_contact_matrix(
            [str(m) for m in clusters[i + 1]], cutoff=dcutoff)
            for i in xrange(n_best_clusters)]
        for i in xrange(n_best_clusters - 1 + add):
            for j in xrange(n_best_clusters - 1 + add):
                try:
                    axes[i, j].set(adjustable='box-forced', aspect=1)
                    axes[i, j].set_visible(False)
                except TypeError:
                    axes.set(adjustable='box-forced', aspect=1)
                    axes.set_visible(False)
        # doing the plot
        for i in xrange(n_best_clusters - 1):
            for j in xrange(1, n_best_clusters):
                if j < i + 1:
                    continue
                try:
                    axes[i + add, j - 1].set_visible(True)
                    axes[i + add, j - 1].set(adjustable='box-forced', aspect=1)
                except TypeError:
                    axes.set_visible(True)
                    axes.set(adjustable='box-forced', aspect=1)
                matrix3 = [[cmatrices[i][k][l] - cmatrices[j][k][l]
                            for l in xrange(self.nloci)]
                           for k in xrange(self.nloci)]
                try:
                    ims = axes[i + add, j - 1].imshow(
                        matrix3, origin='lower', cmap=bwr,
                        interpolation="nearest", vmin=-1, vmax=1,
                        extent=(0.5, len(matrix3) + 0.5,
                                0.5, len(matrix3) + 0.5))
                    axes[i + add, j - 1].grid()
                except TypeError:
                    ims = axes.imshow(
                        matrix3, origin='lower', cmap=bwr,
                        interpolation="nearest", vmin=-1, vmax=1,
                        extent=(0.5, len(matrix3) + 0.5,
                                0.5, len(matrix3) + 0.5))
                    axes.grid()
                if not i and not represent_models:
                    try:
                        axes[i + add, j - 1].set_title('Cluster #%s' % (j + 1),
                                                       color='blue')
                    except TypeError:
                        axes.set_title('Cluster #%s' % (j + 1),
                                       color='blue')
                if j != i + 1:
                    try:
                        axes[i + add, j - 1].yaxis.set_ticks_position('none')
                    except TypeError:
                        axes.yaxis.set_ticks_position('none')
                else:
                    try:
                        plt.setp(axes[i + add, j - 1].get_yticklabels(), visible=True)
                        axes[i + add, j - 1].yaxis.set_ticks_position('left')
                        for item in axes[i + add, j - 1].get_yticklabels():
                            item.set_fontsize(9)
                    except TypeError:
                        plt.setp(axes.get_yticklabels(), visible=True)
                        axes.yaxis.set_ticks_position('left')
                        for item in axes.get_yticklabels():
                            item.set_fontsize(9)
                if i != j - 1:
                    try:
                        axes[i + add, j - 1].xaxis.set_ticks_position('none')
                    except TypeError:
                        axes.xaxis.set_ticks_position('none')
                else:
                    try:
                        plt.setp(axes[i + add, j - 1].get_xticklabels(), visible=True)
                        axes[i + add, j - 1].xaxis.set_ticks_position('bottom')
                        for item in axes[i + add, j - 1].get_xticklabels():
                            item.set_fontsize(9)
                    except TypeError:
                        plt.setp(axes.get_xticklabels(), visible=True)
                        axes.xaxis.set_ticks_position('bottom')
                        for item in axes.get_xticklabels():
                            item.set_fontsize(9)
                if j == n_best_clusters - 1 and not represent_models:
                    try:
                        axes[i + add, j - 1].yaxis.set_label_position('right')
                        axes[i + add, j - 1].set_ylabel('Cluster #%s' % (i + 1),
                                                        rotation=-90,
                                                        fontsize='large',
                                                        color='red',
                                                        va='bottom')
                    except TypeError:
                        axes.yaxis.set_label_position('right')
                        axes.set_ylabel('Cluster #%s' % (i + 1),
                                        rotation=-90, fontsize='large',
                                        color='red', va='bottom')
                try:
                    axes[i + add, j - 1].set_xlim((0.5, len(matrix3) + 0.5))
                    axes[i + add, j - 1].set_ylim((0.5, len(matrix3) + 0.5))
                except TypeError:
                    axes.set_xlim((0.5, len(matrix3) + 0.5))
                    axes.set_ylim((0.5, len(matrix3) + 0.5))
        # new axe for the color bar
        cell = fig.add_axes([0.125, 0.1, 0.01, 0.25])
        cbar = fig.colorbar(ims, cax=cell, cmap=jet)
        cbar.set_ticks([float(k) / 100
                        for k in xrange(-100, 150, 50)])
        cbar.set_ticklabels(['%3s%% ' % (p)
                             for p in xrange(-100, 150, 50)])
        for item in cell.get_yticklabels():
            item.set_fontsize(10)
        cbar.ax.set_ylabel('\n\nPercentage of models with a\n' +
                           'given pair of particles closer\n' +
                           'than the %s nm cutoff\n' % dcutoff +
                           '(red clusters minus blue)\n\n' +
                           (''.join([' ' * 10 +
                                     'Cluster #%-2s: %3s models\n' % (
                                         i + 1, len(clusters[i + 1]))
                                     for i in xrange(n_best_clusters)])),
                           rotation=0, ha='left', va='center')
        plt.suptitle(('Deconvolution analysis for the %s top clusters ' +
                      '(cutoff=%s nm)') % (
                          n_best_clusters, dcutoff), size='x-large')
        if represent_models:
            for i in range(1, n_best_clusters):
                ax = fig.add_subplot(n_best_clusters, n_best_clusters,
                                     i, projection='3d')
                ax.set_title('Cluster #%s' % (i + 1), color='blue')
                if represent_models == 'centroid':
                    mdl = self[self.centroid_model(
                        models=[m for m in clusters[i + 1]])]['index']
                else:
                    if represent_models != 'best':
                        warn("WARNING: represent_model value should be one of" +
                             "'centroid' or 'best' not %s\n"  % (
                                 represent_models) + "Showing best model.")
                self.view_models(models=[self[m]['index']
                                         for m in clusters[i + 1]],
                                 tool='plot', axe=ax, **kwargs)
                for item in [ax]:
                    item.patch.set_visible(False)
            for i in range(n_best_clusters - 1):
                ax = fig.add_subplot(n_best_clusters, n_best_clusters,
                                     n_best_clusters * (i + 2), projection='3d')
                self.view_models(models=[self[m]['index']
                                         for m in clusters[i + 1]],
                                 tool='plot', axe=ax, **kwargs)
                ax.yaxis.set_label_position('top')
                ax.set_title('Cluster #%s' % (i + 1), rotation=-90,
                             fontsize='large', color='red', position=(1, .5),
                             va='center', ha='left')
                for item in [ax]:
                    item.patch.set_visible(False)
            axes[0, n_best_clusters - 1].set_visible(False)
        if savefig:
            tadbit_savefig(savefig)
        else:
            plt.show()

    def contact_map(self, models=None, cluster=None, dynamics=False, stage=None, 
                    cutoff=None, axe=None, savefig=None, savedata=None,
                    cmap='viridis'):
        """
        Plots a contact map representing the frequency of interaction (defined
        by a distance cutoff) between two particles.

        :param None models: if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param None dynamics: compute the contact map for all the stages
        :param None stage: compute the contact map only for the models in
            stage number 'stage'
        :param None cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not, default is 2 times resolution, times scale.
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param None savedata: path to a file where to save the contact map data
           generated, in three columns format (particle1, particle2, percentage
           of models where these two particles are in contact)
        :param viridis cmap: The Colormap instance

        """
        if dynamics: 
            if not (savefig or savedata):
                raise Exception('ERROR: dynamics should only be called ' +
                                'with savefig or savedata option.\n')
                return
            if (savefig and not isdir(savefig)) or (savedata and not isdir(savedata)):
                raise Exception('ERROR: savefig or savedata should ' +
                                'be a folder with dynamics option.\n')
                return
        if not cutoff:
            cutoff = int(2 * self.resolution * self._config['scale'])
        matrices = []
        if dynamics:
            for stg in self.stages:
                matrices.append(self.get_contact_matrix(stage=stg, cutoff=cutoff))
        else:
            matrices.append(self.get_contact_matrix(models, cluster, stage=stage, cutoff=cutoff))
        show = False
        if savedata:
            for nbr, matrix in enumerate(matrices):
                if dynamics:
                    out = open(savedata+'/stage_'+str(nbr)+'.txt', 'w')
                else:
                    out = open(savedata, 'w')
                out.write('#Particle1\tParticle2\tModels_percentage\n')
                for i in xrange(len(matrix)):
                    for j in xrange(i + 1, len(matrix)):
                        out.write('%s\t%s\t%s\n' % (i, j, matrix[i][j]))
                out.close()
        if not savefig and not show and not axe:
            return  # stop here, we do not want to display anything
        cbar = None
        for nbr, matrix in enumerate(matrices):
            if not axe:
                fig = plt.figure(figsize=(8, 6))
                axe = fig.add_subplot(111)
                show = True
            else:
                fig = axe.get_figure()
            cmap = plt.get_cmap(cmap)
            cmap.set_bad('darkgrey', 1)
            ims = axe.imshow(matrix, origin='lower', interpolation="nearest",
                             vmin=0, vmax=1, cmap=cmap,
                             extent=(0.5, self.nloci + 0.5, 0.5, self.nloci + 0.5))
            axe.set_ylabel('Particle')
            axe.set_xlabel('Particle')
            if not cbar:
                cbar = axe.figure.colorbar(ims)
                oldlabels = cbar.ax.get_yticklabels()
                newlabels = map(lambda x: str(int(100 * float(x.get_text())))+'%', oldlabels)
                cbar.ax.set_yticklabels(newlabels)
                cbar.ax.set_ylabel('Percentage of models with particles at <' +
                                   '%s nm' % (cutoff))
            if dynamics:
                axe.set_title('Contact map stage %s' % str(nbr))
            else:
                axe.set_title('Contact map')
            if savefig:
                if dynamics:
                    tadbit_savefig(savefig+'/contact_map_stage_'+str(nbr)+'.png')
                else:
                    tadbit_savefig(savefig)
            elif show and not dynamics:
                plt.show()
        if dynamics:
            try:
                system('ffmpeg -r 50 -i '+savefig+'/contact_map_stage_%d.png -vb 20M -y '+savefig+'/contact_map_all_stages.mpeg')
            except:
                pass

    def accessibility(self, radius, models=None, cluster=None, nump=100,
                      superradius=200, savefig=None, savedata=None, axe=None,
                      plot=True, error=True, steps=(1, )):
        """
        Calculates a mesh surface around the model (distance equal to input
        **radius**) and checks if each point of this mesh could be replaced by
        an object (i.e. a protein) of a given **radius**

        Outer part of the model can be excluded from the estimation of
        accessible surface, as the occupancy outside the model is unkown (see
        superradius option).

        :param radius: radius of the object we want to fit in the model.
        :param 100 nump: number of points to draw around a given particle. This
           number also sets the number of points drawn around edges, as each
           point occupies a given surface (see maths below). *Note that this
           number is considerably lowered by the occupancy of edges, depending
           of the angle formed by the edges surrounding a given particle, only
           10% to 50% of the ``nump`` will be drawn in fact.*
        :param None savefig: path where to save chimera image
        :param (1, ) steps: how many particles to group for the
           estimation. By default 1 curve is drawn
        :param 200 superradius: radius of an object used to exclude outer
           surface of the model. Superradius must be higher than radius.

        This function will first define a mesh around the chromatin,
        representing all possible position of the center of the object we want
        to fit. This mesh will be at a distance of *radius* from the chromatin
        strand. All dots in the mesh represents an equal area (*a*), the whole
        surface of the chromatin strand being: :math:`A=n \\times a` (*n* being
        the total number of dots in the mesh).

        The mesh consists of spheres around particles of the model, and
        cylinders around edges joining particles (no overlap is allowed between
        sphere and cylinders or cylinder and cylinder when they are
        consecutive).

        If we want that all dots of the mesh representing the surface of the
        chromatin, corresponds to an equal area (:math:`a`)

        .. math::

          a = \\frac{4\pi r^2}{s} = \\frac{2\pi r N_{(d)}}{c}

        with:

        * :math:`r` radius of the object to fit (as the input parameter **radius**)
        * :math:`s` number of points in sphere
        * :math:`c` number of points in circle (as the input parameter **nump**)
        * :math:`N_{(d)}` number of circles in an edge of length :math:`d`

        According to this, when the distance between two particles is equal
        to :math:`2r` (:math:`N=2r`), we would have :math:`s=c`.

        As :

        .. math::

          2\pi r = \sqrt{4\pi r^2} \\times \sqrt{\pi}

        It is fair to state the number of dots represented along a circle as:

        .. math::

          c = \sqrt{s} \\times \sqrt{\pi}

        Thus the number of circles in an edge of length :math:`d` must be:

        .. math::

          N_{(d)}=\\frac{s}{\sqrt{s}\sqrt{\pi}}\\times\\frac{d}{2r}

        :returns: a list of *1-* the number of dots in the mesh that could be
           occupied by an object of the given radius *2-* the total number of
           dots in the mesh *3-* the estimated area of the mesh (in square
           micrometers) *4-* the area of the mesh of a virtually straight strand
           of chromatin defined as
           :math:`contour\\times 2\pi r + 4\pi r^2` (also in
           micrometers) *5-* a list of number of (accessibles, inaccessible) for
           each particle (percentage burried can be infered afterwards by
           accessible/(accessible+inaccessible) )
        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        acc = []
        for model in models:
            acc_vs_inacc = self[model].accessible_surface(
                radius, nump=nump, superradius=superradius,
                include_edges=False)[-1]
            acc.append([(float(j) / (j + k))  if (j + k) else 0.0
                        for _, j, k in acc_vs_inacc])
        accper, errorn, errorp = self._windowize(zip(*acc), steps, average=True)
        if savedata:
            out = open(savedata, 'w')
            out.write('# Particle\t%s\n' % ('\t'.join([
                str(c) + '\t' + '2*stddev(%d)' % c for c in steps])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (part + 1, '\t'.join(
                    ['nan\tnan' if part >= len(accper[c]) else
                     (str(round(accper[c][part], 3)) + '\t' +
                      str(round(errorp[c][part] - accper[c][part], 3)))
                     if accper[c][part] else 'nan\tnan'
                     for c in steps])))
            out.close()
        if not plot:
            return  # stop here, we do not want to display anything
        # plot
        ylabel = 'Accessibility to an object with a radius of %s nm ' % (radius)
        xlabel = 'Particle number'
        title = 'Accesibility per particle'
        self._generic_per_particle_plot(steps, accper, error, errorp,
                                        errorn, savefig, axe, xlabel=xlabel,
                                        ylabel=ylabel, title=title, ylim=(0, 1))
        #if savefig:
        #    tadbit_savefig(savefig)
        #elif not axe:
        #    plt.show()
        #plt.close('all')


    def _get_density(self, models, interval, use_mass_center):
        dists = [[None] * len(models)] * interval
        for p in range(interval, self.nloci - interval):
            part1, part2, part3 = p - interval, p, p + interval
            if use_mass_center:
                subdists = []
                for m in models:
                    try:
                        coord1 = get_center_of_mass(
                            self[m]['x'][part1:part2],
                            self[m]['y'][part1:part2],
                            self[m]['z'][part1:part2],
                            self._zeros)
                        coord2 = get_center_of_mass(
                            self[m]['x'][part2:part3],
                            self[m]['y'][part2:part3],
                            self[m]['z'][part2:part3],
                            self._zeros)
                        subdists.append(distance(coord1, coord2))
                    except ZeroDivisionError:  # part1==part2 or part2==part3
                        subdists.append(float('nan'))
                dists.append([float(interval * self.resolution) / d for d in subdists])
            else:
                dist1 = self.median_3d_dist(part1 + 1, part2 + 1, models,
                                            plot=False, median=False)
                dist2 = self.median_3d_dist(part2 + 1, part3 + 1, models,
                                            plot=False, median=False)
                dist = [(d1 + d2) for d1, d2 in zip(dist1, dist2)]
                dists.append([float(interval * self.resolution * 2) / d
                              for d in dist])
        return dists

    def density_plot(self, models=None, cluster=None, steps=(1, 2, 3, 4, 5),
                     interval=1, use_mass_center=False, error=False, axe=None,
                     savefig=None, savedata=None, plot=True):
        """
        Plots the number of nucleotides per nm of chromatin vs the modeled
        region bins.

        :param None models: if None (default) the density plot will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the density plot only for the models in the
           cluster number 'cluster'
        :param (1, 2, 3, 4, 5) steps: how many particles to group for the
           estimation. By default 5 curves are drawn
        :param False error: represent the error of the estimates
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param 1 interval: distance are measure with this given interval
           between two bins.
        :param False use_mass_center: if interval is higher than one, calculates the
           distance between the center of mass of the particles *n* to
           *n+interval* and the center of mass of the particles *n+interval* and
           *n+2interval*
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param None savedata: path to a file where to save the density data
           generated (1 column per step + 1 for particle number).
        :param True plot: e.g. if False, only saves data. No plotting done

        """
        if isinstance(steps, int):
            steps = (steps, )

        models = self._get_models(models, cluster)
        dists = self._get_density(models, interval, use_mass_center)
        distsk, errorn, errorp = self._windowize(dists, steps, interval=interval,
                                                 average=False)

        # write consistencies to file
        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle\t%s\n' % ('\t'.join([
                str(c) + '\t' + '2*stddev(%d)' % c for c in steps])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (part + 1, '\t'.join(
                    ['nan\tnan' if part >= len(distsk[c]) else
                     (str(round(distsk[c][part], 3)) + '\t' +
                      str(round(errorp[c][part] - round(distsk[c][part], 3))))
                     if distsk[c][part] else 'nan\tnan'
                     for c in steps])))
            out.close()
        if plot:
            xlabel = 'Particle number'
            ylabel = 'Density (bp / nm)'
            title  = 'Chromatin density'
            # self._generic_per_particle_plot(steps, distsk, error, errorp, errorn,
            #                                 xlabel=xlabel, ylabel=ylabel, title=title)
            self._generic_per_particle_plot(steps, distsk, error, errorp,
                                            errorn, savefig, axe, xlabel=xlabel,
                                            ylabel=ylabel, title=title)

    def _get_interactions(self, models, cutoff):
        interactions = [[] for _ in xrange(self.nloci)]
        if not cutoff:
            cutoff = int(2 * self.resolution * self._config['scale'])
        cutoff2 = cutoff**2
        for i in xrange(self.nloci):
            for m in models:
                val = 0
                mdl = self[m]
                for j in xrange(self.nloci):
                    if i == j:
                        continue
                    val += ((mdl['x'][i] - mdl['x'][j])**2 +
                            (mdl['y'][i] - mdl['y'][j])**2 +
                            (mdl['z'][i] - mdl['z'][j])**2) < cutoff2
                    # val += self.__square_3d_dist(i + 1, j + 1, models=[m])[0] < cutoff2
                interactions[i].append(val)
        return interactions

    def interactions(self, models=None, cluster=None, cutoff=None,
                     steps=(1, 2, 3, 4, 5), axe=None, error=False,
                     savefig=None, savedata=None, average=True, plot=True):
        """
        Plots, for each particle, the number of interactions (particles closer
        than the given cut-off). The value given is the average for all models.

        :param None models: if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param None cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not, default is 2 times resolution, times scale.
        :param (1, 2, 3, 4, 5) steps: how many particles to group for the
           estimation. By default 5 curves are drawn
        :param False error: represent the error of the estimates
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param None savedata: path to a file where to save the contact map data
           generated, in three columns format (particle1, particle2, percentage
           of models where these two particles are in contact)
        :param True average: calculate average interactions along models,
           otherwise, the median.
        :param True plot: e.g. only saves data. No plotting done

        """
        if isinstance(steps, int):
            steps = (steps, )

        models = self._get_models(models, cluster)

        interactions = self._get_interactions(models, cutoff)

        distsk, errorn, errorp = self._windowize(interactions, steps,
                                                 average=average)
        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle\t%s\n' % (
                '\t'.join(['%s_interactions(%s)\t2*stddev(%s)' % (
                    'Average' if average else 'Median', k, k) for k in steps])))
            for i in xrange(self.nloci):
                out.write('%s\t%s\n' % (i + 1, '\t'.join(
                    ['%s\t%s' % (str('nan' if (len(distsk[k]) <= i or
                                               distsk[k][i] is None)
                                     else round(distsk[k][i], 2)), (
                        str('nan' if (len(distsk[k]) <= i or
                                      distsk[k][i] is None)
                            else round(errorp[k][i] - distsk[k][i], 2))))
                     for k in steps])))
            out.close()
        if plot:
            ylabel = 'Number of particles closer than %s nm' % (cutoff)
            xlabel = 'Particle number'
            title = 'Interactions per particle'
            self._generic_per_particle_plot(steps, distsk, error, errorp,
                                            errorn, savefig, axe, xlabel=xlabel,
                                            ylabel=ylabel, title=title)

    def model_consistency(self, cutoffs=None, models=None,
                          cluster=None, axe=None, savefig=None, savedata=None,
                          plot=True):
        """
        Plots the particle consistency, over a given set of models, vs the
        modeled region bins. The consistency is a measure of the variability
        (or stability) of the modeled region (the higher the consistency value,
        the higher stability).

        :param None cutoffs: list of distance cutoffs (nm) used to compute the
           consistency. Two particle are considered consistent if their distance
           is less than the given cutoff, default is a tuple of 0.5, 1, 1.5 and
           2 times resolution, times scale.
        :param None models:  if None (default) the consistency will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the consistency only for the models in the
           cluster number 'cluster'
        :param '/tmp/tmp_cons' tmp_path: location of the input files for
           TM-score program
        :param '' tmsc: path to the TMscore_consistency script (assumed to be
           installed by default)
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param None savedata: path to a file where to save the consistency data
           generated (1 column per cutoff + 1 for particle number).
        :param True plot: e.g. only saves data. No plotting done

        """
        models = self._get_models(models, cluster)
        models = [self.__models[m] for m in models]

        if not cutoffs:
            cutoffs = (int(0.5 * self.resolution * self._config['scale']),
                       int(1.0 * self.resolution * self._config['scale']),
                       int(1.5 * self.resolution * self._config['scale']),
                       int(2.0 * self.resolution * self._config['scale']))
        consistencies = {}
        for cut in cutoffs:
            consistencies[cut] = calc_consistency(models, self.nloci,
                                                  self._zeros, cut)
        # write consistencies to file
        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle\t%s\n' % ('\t'.join([str(c) for c in cutoffs])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (str(part + 1), '\t'.join(
                    [str(round(consistencies[c][part], 3)) for c in cutoffs])))
            out.close()
        if not plot:
            return
        # plot
        show = False if axe else True
        axe = setup_plot(axe)
        plots = []
        self._plot_polymer(axe)
        scat1 = plt.Line2D((0, 1), (0, 0), color=(0.15, 0.15, 0.15), marker='o',
                           linestyle='')
        scat2 = plt.Line2D((0, 1), (0, 0), color=(0.7 , 0.7 , 0.7 ), marker='o',
                           linestyle='')
        for i, cut in enumerate(cutoffs[::-1]):
            plots += axe.plot(range(1, self.nloci + 1),
                              consistencies[cut], color='darkred',
                              alpha=1 - i / float(len(cutoffs)))
        try:
            axe.legend(plots + [scat1, scat2],
                       ['%s nm' % (k) for k in cutoffs[::-1]] + ['particles with restraints',
                                                                 'particles without restraints'],
                       numpoints=1, fontsize='small', loc='center left',
                       bbox_to_anchor=(1, 0.5))
        except TypeError:
            axe.legend(plots + [scat1, scat2],
                       ['%s nm' % (k) for k in cutoffs[::-1]] + ['particles with restraints',
                                                                 'particles without restraints'],
                       numpoints=1, loc='center left',
                       bbox_to_anchor=(1, 0.5))
        axe.set_xlim((1, self.nloci))
        axe.set_ylim((0, 100))
        axe.set_xlabel('Particle')
        axe.set_ylabel('Consistency (%)')
        if cluster:
            axe.set_title('Cluster %s' % (cluster))
        elif len(models) == len(self):
            axe.set_title('All clusters')
        else:
            axe.set_title('Selected models')
        plt.subplots_adjust(left=0.1, right=0.77)
        if savefig:
            tadbit_savefig(savefig)
        elif show:
            plt.show()
        plt.close('all')

    def walking_dihedral(self, models=None, cluster=None, steps=(1, 3),
                         span=(-2, 1, 0, 1, 3), error=False,
                         plot=True, savefig=None, axe=None, savedata=None):
        """
        Plots the dihedral angle between successive planes. A plane is formed by
        3 successive loci.

        :param None models: if None (default) the dihedral angle will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the dihedral angle only for the models in the
           cluster number 'cluster'
        :param (1, 3) steps: how many particles to group for the estimation.
           By default 2 curves are drawn
        :param True signed: whether to compute the sign of the angle according
           to a normal plane, or not.
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param (-3, -1, 0, 1, 2) span: by default computes angle between a plane
           comprising the current residue (*0*), the residue after (*1*) and the
           residue 3 positions before (*-3*), and the plane with the current
           residue, the residue directly after (*1*) and the residue 3 positions
           after (*3*). e.g. In the scheme bellow it's plane ACD vs plane CDF.
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param False error: represent the error of the estimates
        :param True plot: e.g. if False, only saves data. No plotting done
        :param None savedata: path to a file where to save the angle data
           generated (1 column per step + 1 for particle number).


        ::

                                C..........D
                             ...            ...
                          ...                 ...
                       ...                       ...
            A..........B                            .E
           ..                                        .
          .                                          .
                                                     .                 .
                                                     .                .
                                                     F...............G


        """
        models = self._get_models(models, cluster)
        rads = {}
        if span[0] > 0:
            raise ValueError('ERROR: first element of span should be negative')
        if span[-1] < 0:
            raise ValueError('ERROR: last element of span should be negative')

        rads = [[None] * len(models)] * (-span[0])
        for res in xrange(-span[0], self.nloci - span[-1]):
            subrad = self.dihedral_angle(res + span[0], res + span[1],
                                         res + span[2],
                                         res + span[3], res + span[4], models)
            # rads.append([None] * abs(span[0]) + subrad + [None] * span[3])
            rads.append(subrad)
        rads += [[None] * len(models)] * (span[-1])
        radsk, errorn, errorp = self._windowize(rads, steps, interval=0,
                                                average=False, minerr=-360)
        if plot:
            xlabel = 'Particle number'
            ylabel = 'Dihedral angle in degrees'
            title  = 'Dihedral angle between consecutive loci'
            self._generic_per_particle_plot(steps, radsk, error, errorp,
                                            errorn, savefig, axe, xlabel=xlabel,
                                            ylabel=ylabel, title=title)

        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle\t%s\n' % ('\t'.join([
                str(c) + '\t' + '2*stddev(%d)' % c for c in steps])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (part + 1, '\t'.join(
                    ['nan\tnan' if part >= len(radsk[c]) else
                     (str(round(radsk[c][part], 3)) + '\t' +
                      str(round(errorp[c][part] - round(radsk[c][part], 3))))
                     if radsk[c][part] else 'nan\tnan'
                     for c in steps])))
            out.close()

        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()
        plt.close('all')
        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()
        plt.close('all')

    def walking_angle(self, models=None, cluster=None, steps=(1, 3), signed=True,
                      savefig=None, savedata=None, axe=None, plot=True,
                      error=False):
        """
        Plots the angle between successive loci in a given model or set of
        models. In order to limit the noise of the measure angle is calculated
        between 3 loci, between each are two other loci. E.g. in the scheme
        bellow, angle are calculated between loci A, D and G.

        :param None models: if None (default) all models will be used for
           computation. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the angle only for the models in the
           cluster number 'cluster'
        :param (1, 3) steps: how many particles to group for the estimation.
           By default 2 curves are drawn
        :param True signed: whether to compute the sign of the angle according
           to a normal plane, or not.
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param False error: represent the error of the estimates
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param True plot: e.g. if False, only saves data. No plotting done
        :param None savedata: path to a file where to save the angle data
           generated (1 column per step + 1 for particle number).


        ::

                                C..........D
                             ...            ...
                          ...                 ...
                       ...                       ...
            A..........B                            .E
           ..                                        .
          .                                          .
                                                     .                 .
                                                     .                .
                                                     F...............G


        """
        if not isinstance(steps, tuple):
            steps = (steps,)
        models = self._get_models(models, cluster)
        rads = []
        ones = array([1., 1., 1.])
        sign = 1
        for model in models:
            mox = self[model]['x']
            moy = self[model]['y']
            moz = self[model]['z']
            subrad = [None] * 3
            for res in xrange(1, self.nloci - 5):
                subrad.append(self.angle_between_3_particles(
                    res, res + 3, res + 6, models=[model], all_angles=False))
                if signed:
                    res1 = array([mox[res - 1], moy[res - 1], moz[res - 1]])  # faster than get_particle coordinate
                    res2 = array([mox[res + 2], moy[res + 2], moz[res + 2]])
                    res3 = array([mox[res + 5], moy[res + 5], moz[res + 5]])
                    vec1 = res1 - res2 / norm(res1 - res2)
                    vec2 = res1 - res3 / norm(res1 - res3)
                    sign = dot(ones, cross(vec1, vec2))
                    sign = -1 if sign < 0 else 1
                subrad[-1] *= sign
            subrad += [None] * 3
            rads.append(subrad)

        radsk, errorn, errorp = self._windowize(zip(*rads), steps, interval=0,
                                                average=False, minerr=-360)
        if plot:
            xlabel = 'Particle number'
            ylabel = 'Angle in degrees'
            title  = 'Angle between consecutive loci'
            self._generic_per_particle_plot(steps, radsk, error, errorp,
                                            errorn, savefig, axe, xlabel=xlabel,
                                            ylabel=ylabel, title=title)
        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle\t%s\n' % ('\t'.join([
                str(c) + '\t' + '2*stddev(%d)' % c for c in steps])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (part + 1, '\t'.join(
                    ['nan\tnan' if part >= len(radsk[c]) else
                     (str(round(radsk[c][part], 3)) + '\t' +
                      str(round(errorp[c][part] - round(radsk[c][part], 3))))
                     if radsk[c][part] else 'nan\tnan'
                     for c in steps])))
            out.close()

        #if savefig:
        #    tadbit_savefig(savefig)
        #elif not axe:
        #    plt.show()
        #plt.close('all')

    def zscore_plot(self, axe=None, savefig=None, do_normaltest=False,
                    stage=0, cmap='Reds'):
        """
        Generate 3 plots. Two heatmaps of the Z-scores used for modeling, one
        of which is binary showing in red Z-scores higher than upper cut-off;
        and in blue Z-scores lower than lower cut-off. Last plot is an histogram
        of the distribution of Z-scores, showing selected regions. Histogram
        also shows the fit to normal distribution.
        :param None axe: a matplotlib.axes.Axes object to define the plot
            appearance
        :param None savefig: path to a file where to save the image generated;
            if None, the image will be shown using matplotlib GUI (the extension
            of the file name will determine the desired format).
        :param False do_normaltest: to display the result of a test of normality
            (D'Agostino and Pearson's test, that combines skew and kurtosis).
        :param viridis cmap: The Colormap instance
        """

        if stage > -1 and stage in self.stages:
            stage_zscore = self._zscores[stage]
        elif len(self.stages) == 0:
            stage_zscore = self._zscores
        else:
            raise ValueError('ERROR: please specify a correct stage')
        zsc_mtrx = reduce(lambda x, y: x + y, [[k] + stage_zscore[k].keys()
                                                for k in stage_zscore.keys()])
        max_bin = max([int(i) for i in zsc_mtrx])
        zsc_mtrx = [[float('nan') for _ in xrange(max_bin)]
                    for _ in xrange(max_bin)]
        for i in xrange(max_bin):
            for j in xrange(max_bin):
                try:
                    zsc_mtrx[i][j] = stage_zscore[str(i)][str(j)]
                except KeyError:
                    try:
                        zsc_mtrx[i][j] = stage_zscore[str(j)][str(i)]
                    except KeyError:
                        zsc_mtrx[i][j] = float('Nan')

        # suppress warnings due to nans in following lines
        prevErrorSet = seterr()
        seterr(invalid='ignore')

        masked_array = ma.array (zsc_mtrx, mask=isnan(zsc_mtrx))
        masked_array_top = ma.array (
            masked_array, mask=masked_array < self._config['upfreq'])
        masked_array_bot = ma.array (
            masked_array, mask=self._config['lowfreq'] < masked_array)
        
        # get back to previous error warnings
        seterr(invalid=prevErrorSet['invalid'])

        # color coding
        if cmap == 'viridis':
            lowf = plt.get_cmap('Purples')
            upf = plt.get_cmap('summer')
        else:
            lowf = plt.get_cmap('Blues')
            upf = plt.get_cmap('Reds')
        cmap =  plt.get_cmap(cmap)
        cmap.set_bad('white', 0.)
        lowf.set_bad('white', 0.)
        upf.set_bad('white', 0.)
        
        if not axe:
            fig = plt.figure(figsize=(25, 5.5))
        else:
            fig = axe.get_figure()
        ax = fig.add_subplot(131)
        ims = ax.imshow(zsc_mtrx, origin='lower',
                        interpolation="nearest", cmap=cmap)
        ax.set_ylabel('Particles')
        ax.set_xlabel('Particles')
        ax.set_title('Z-scores of the normalized Hi-C count')
        cbar = ax.figure.colorbar(ims, cmap=cmap)
        cbar.ax.set_ylabel('Z-score value')

        ax = plt.axes([.38, 0.11, .28, .61])
        zdata = sorted(reduce(lambda x, y: x + y,
                                [stage_zscore[v].values()
                                for v in stage_zscore.keys()]))
        try:
            _, _, patches = ax.hist(zdata, bins=25, linewidth=1,
                                    facecolor='none', edgecolor='k', density=True)
        except AttributeError:
            _, _, patches = ax.hist(zdata, bins=25, linewidth=1,
                                    facecolor='none', edgecolor='k')
        k2, pv = normaltest(zdata)
        normfit = sc_norm.pdf(zdata, np_mean(zdata), np_std(zdata))
        normplot = ax.plot(zdata, normfit, ':o', color='grey', ms=3, alpha=.4)
        try:
            ax.hist(
                reduce(lambda x, y: x + y, [stage_zscore[v].values()
                                            for v in stage_zscore.keys()]),
                bins=25, linewidth=2, facecolor='none', edgecolor='k',
                histtype='stepfilled', density=True)
        except AttributeError:
            ax.hist(
                reduce(lambda x, y: x + y, [stage_zscore[v].values()
                                            for v in stage_zscore.keys()]),
                bins=25, linewidth=2, facecolor='none', edgecolor='k',
                histtype='stepfilled')
        height1 = height2 = 0
        minv = nanmin(masked_array)
        maxv = nanmax(masked_array)
        for thispatch in patches:
            beg = thispatch.get_x()
            end = thispatch.get_x() + thispatch.get_width()
            ax.fill_betweenx([0] + [thispatch.get_height()] * 100,
                                max(beg, self._config['upfreq'])
                                if end > self._config['upfreq' ] else beg,
                                min(end, self._config['lowfreq' ])
                                if beg < self._config['lowfreq'] else end,
                                color=(
                                    lowf(beg / minv) if beg < self._config['lowfreq'] else
                                    upf (end / maxv) if end > self._config['upfreq']
                                    else 'w'))
            if end > self._config['lowfreq' ] and beg < self._config['lowfreq']:
                height1 = thispatch.get_height()
            elif beg < self._config['upfreq'] and end > self._config['upfreq' ]:
                height2 = thispatch.get_height()
        labels = []
        p1 = plt.Rectangle((0, 0), 1, 1, fc=lowf(0.7), color='k')
        labels.append('< %.2f (force particles apart)' % (
            self._config['lowfreq']))
        p2 = plt.Rectangle((0, 0), 1, 1, fc="w", color='k')
        labels.append('Not used (no constraints)')
        p3 = plt.Rectangle((0, 0), 1, 1, fc=upf(0.7), color='k')
        labels.append('> %.2f (force particles together)' % (
            self._config['upfreq']))
        try:
            ax.legend([p1, p2, p3, normplot[0]], labels + [
                "Fitted normal distribution" +
                (("\n D'Agostino Pearson's normality test $K^2$=%.2f pv=%.3f"
                    % (k2, pv)) if do_normaltest else '')],
                fontsize='small', frameon=False,
                bbox_to_anchor=(1.013, 1.3 + (.03 if do_normaltest else 0)))
        except TypeError:
            ax.legend([p1, p2, p3, normplot[0]], labels + [
                "Fitted normal distribution" +
                (("\n D'Agostino Pearson's normality test $K^2$=%.2f pv=%.3f"
                    % (k2, pv)) if do_normaltest else '')],
                frameon=False,
                bbox_to_anchor=(1.013, 1.3 + (.03 if do_normaltest else 0)))
        ax.set_xlabel('Z-scores')
        ax.set_ylabel('Proportion of particles')
        ax.vlines(self._config['lowfreq'], 0, height1, color='k',
                    linestyle='-', lw=2, alpha=1)
        ax.vlines(self._config['upfreq'] , 0, height2, color='k',
                    linestyle='-', lw=2, alpha=1)

        ax = fig.add_subplot(133)
        _ = ax.imshow(masked_array_top, origin='lower',
                        interpolation="nearest", cmap=upf)
        _ = ax.imshow(masked_array_bot, origin='lower',
                        interpolation="nearest", cmap=lowf)
        ax.set_ylabel('Particles')
        ax.set_xlabel('Particles')
        ax = plt.axes([.42, 0.11, .48, .79])
        ax.set_title('Binary representation of Z-scores used in the ' +
                        'computation of restraints')
        ax.set_title('Binary representation of Z-scores used in the ' +
                        'computation of restraints')
        ax.axison = False
        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()
        plt.close('all')


    def correlate_with_real_data(self, models=None, cluster=None,
                                 stage=None, index=0,
                                 dynamics=False, cutoff=None,
                                 off_diag=1, plot=False, axe=None, savefig=None,
                                 corr='spearman', midplot='hexbin',
                                 log_corr=True, contact_matrix=None,
                                 cmap='viridis', show_bad_columns=True):
        """
        Plots the result of a correlation between a given group of models and
        original Hi-C data.

        :param None models: if None (default) the correlation will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the correlation only for the models in the
           cluster number 'cluster'
        :param None dynamics: compute the correlation for all the stages
        :param None cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not, default is 2 times resolution, times scale.
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param False plot: to display the plot
        :param True log_corr: log plot for correlation
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None contact_matrix: input a contact matrix instead of computing
           it from the models
        :param 'viridis' cmap: The Colormap instance
        :param True show_bad_columns: Wether to hide or not bad columns in the 
            contact map

        :returns: correlation coefficient rho, between the two
           matrices. A rho value greater than 0.7 indicates a very good
           correlation
        """
        if dynamics:
            if not savefig:
                raise Exception('ERROR: dynamics should only be called ' +
                                'with savefig option.\n')
                return
            if not isdir(savefig):
                raise Exception('ERROR: savefig should ' +
                                'be a folder with dynamics option.\n')
                return
        elif stage is not None and stage not in self.stages:
            raise Exception('ERROR: stage ' +
                            'not found in stages.\n')
            return
        if not cutoff:
            cutoff = int(2 * self.resolution * self._config['scale'])
        if contact_matrix:
            all_original_data = [0]
            all_model_matrix = [contact_matrix]
        else:
            if dynamics:
                all_model_matrix = []
                all_original_data = []
                for st in range(0,int((len(self.stages)-1)/self.models_per_step)+1):
                    all_original_data.append(st)
                    all_model_matrix.append(self.get_contact_matrix(stage=int(st*self.models_per_step), cutoff=cutoff, show_bad_columns=show_bad_columns))
            elif stage is not None:
                all_original_data = [index]
                all_model_matrix = [self.get_contact_matrix(stage=stage,cutoff=cutoff)]
            else:
                all_original_data = [index]
                all_model_matrix = [self.get_contact_matrix(models=models, cluster=cluster,
                                                   cutoff=cutoff, show_bad_columns=show_bad_columns)]
        correl = {}
        for model_matrix, od in zip(all_model_matrix,all_original_data):
            oridata = []
            moddata = []
            if len(model_matrix) == 0:
                correl[od] = 'Nan'
                continue
            if dynamics:
                original_data = self._original_data[od]
            elif stage is not None or len(self.stages) > 0:
                original_data = self._original_data[od]
            else:
                original_data = self._original_data
            for i in xrange(len(original_data)):
                for j in xrange(i + off_diag, len(original_data)):
                    if not original_data[i][j] > 0:
                        continue
                    oridata.append(original_data[i][j])
                    moddata.append(model_matrix[i][j])
            if corr == 'spearman':
                correl[od] = spearmanr(moddata, oridata)
            elif corr == 'pearson':
                correl[od] = pearsonr(moddata, oridata)
            elif corr == 'logpearson':
                correl[od] = pearsonr(nozero_log_list(moddata), nozero_log_list(oridata))
            elif corr == 'chi2':
                tmpcorr = chisquare(array(moddata), array(oridata))
                tmpcorr = 1. / tmpcorr[0], tmpcorr[1]
                correl[od] = tmpcorr
            else:
                raise NotImplementedError('ERROR: %s not implemented, must be one ' +
                                          'of spearman, pearson or frobenius\n')
        if not plot and not savefig:
            if len(correl) < 2:
                return correl[next(iter(correl))]
            return correl
        cbar = None
        for model_matrix, od in zip(all_model_matrix,all_original_data):
            oridata = []
            moddata = []
            if correl[od] == 'Nan':
                continue
            if dynamics:
                original_data = self._original_data[od]
                stage_label = ' for stage %d'%int(od*self.models_per_step)
                hic_label = ' %d'%od
            elif stage is not None or len(self.stages) > 0:
                original_data = self._original_data[od]
                stage_label = ' for stage %d'%stage
                hic_label = ' %d'%od
            else:
                original_data = self._original_data
                stage_label = ''
                hic_label = ''
            for i in xrange(len(original_data)):
                for j in xrange(i + off_diag, len(original_data)):
                    if not original_data[i][j] > 0:
                        continue
                    oridata.append(original_data[i][j])
                    moddata.append(model_matrix[i][j])
            if not axe:
                fig = plt.figure(figsize=(20, 4.5))
            else:
                fig = axe.get_figure()
            fig.suptitle('Correlation between normalized-real%s and modeled '%stage_label +
                         'contact maps%s (correlation=%.4f)' % (hic_label, correl[od][0]),
                         size='x-large')
            ax = fig.add_subplot(131)
            # imshow of the modeled data
            cmap = plt.get_cmap(cmap)
            cmap.set_bad('darkgrey', 1)
            ims = ax.imshow(model_matrix, origin='lower', interpolation="nearest",
                             vmin=0, vmax=1, cmap=cmap,
                             extent=(0.5, self.nloci + 0.5, 0.5, self.nloci + 0.5))
            ax.set_ylabel('Particle')
            ax.set_xlabel('Particle')
            if not cbar:
                cbar = ax.figure.colorbar(ims)
                cbar.ax.set_yticklabels(['%3s%%' % (p) for p in range(0, 110, 10)])
                cbar.ax.set_ylabel('Percentage of models with particles at <' +
                                   '%s nm' % (cutoff))
            ax.set_title('Contact map')        # correlation
    
            ax = fig.add_subplot(132)
            try:
                if log_corr:
                    minmoddata = float(min([m for m in moddata if m]))
                    minoridata = float(min([m for m in oridata if m]))
                    moddata, oridata = (log2([(m if m else minmoddata / 2) * 100 for m in moddata]),
                                        log2([m if m else minoridata / 2 for m in oridata]))
            except:
                warn('WARNING: unable to log for correlation with real data...')
            slope, intercept, r_value, p_value, _ = linregress(moddata, oridata)
            # slope, intercept, r_value, p_value, std_err = linregress(moddata, oridata)
            if midplot == 'classic':
                lnr = ax.plot(moddata, intercept + slope * array (moddata), color='k',
                              ls='--', alpha=.7)
                ax.legend(lnr, ['p-value: %.3f, R: %.3f' % (p_value, r_value)])
                ax.plot(moddata, oridata, 'ro', alpha=0.5)
                ax.set_xlabel('Modelled data')
                ax.set_ylabel('Real data')
            elif midplot == 'hexbin':
                hb = ax.hexbin(moddata, oridata, mincnt=1,
                               gridsize=50, cmap=plt.cm.Spectral_r)
                lnr = ax.plot(moddata, intercept + slope * array (moddata), color='k',
                              ls='--', alpha=.7)
                ax.set_xlabel(
                    '%sroportion of models with a particle pair closer than cutoff' % (
                        'Log p' if log_corr else 'P'))
                ax.set_ylabel('%sormalized Hi-C count for a particle pair' % (
                    'Log n' if log_corr else 'N'))
                cbaxes = fig.add_axes([0.41, 0.42, 0.005, 0.45])
                cbar = plt.colorbar(hb, cax=cbaxes)  # orientation='horizontal')
                cbar.set_label('Number of particle pairs')
            elif midplot == 'triple':
                maxval = max(oridata)
                minval = min(oridata)
                ax.set_visible(False)
                axleft = fig.add_axes([0.42, 0.18, 0.1, 0.65])
                axleft.spines['right'].set_color('none')
                axleft.spines['bottom'].set_color('none')
                axleft.spines['left'].set_smart_bounds(True)
                axleft.spines['top'].set_smart_bounds(True)
                axleft.xaxis.set_ticks_position('top')
                axleft.yaxis.set_ticks_position('left')
                axleft.set_ylabel('Normalized Hi-C count for a particle pair')
                axleft.patch.set_visible(False)
                axbott = fig.add_axes([0.44, 0.13, 0.17, 0.5])
                axbott.spines['left'].set_color('none')
                axbott.spines['top'].set_color('none')
                axbott.spines['left'].set_smart_bounds(True)
                axbott.spines['bottom'].set_smart_bounds(True)
                axbott.xaxis.set_ticks_position('bottom')
                axbott.yaxis.set_ticks_position('right')
                axbott.patch.set_visible(False)
                axbott.set_xlabel('Proportion of models with a particle pair ' +
                                  ' interacting')
                axmidl = fig.add_axes([0.44, 0.18, 0.17, 0.65])
                axbott.hist(moddata, bins=20, alpha=.2)
                x, _  = histogram([i if str(i) != '-inf' else 0. for i in oridata],
                                  bins=20)
                axleft.barh(linspace(minval, maxval, 20), x,
                            height=(maxval - minval) / 20, alpha=.2)
                axleft.set_ylim((minval -
                                 (maxval - minval) / 20, maxval +
                                 (maxval - minval) / 20))
                axmidl.plot(moddata, oridata, 'k.', alpha=.3)
                axmidl.plot(moddata, intercept + slope * array (moddata), color='k',
                            ls='--', alpha=.7)
                axmidl.set_ylim(axleft.get_ylim())
                axmidl.set_xlim(axbott.get_xlim())
                axmidl.axis('off')
                # axmidl.patch.set_visible(False)
            ax.set_title('Real versus modelled data')
            ax = fig.add_subplot(133)
            cmap = plt.get_cmap(cmap)
            cmap.set_bad('darkgrey', 1)
            ims = ax.imshow(log2(original_data), origin='lower',
                            interpolation="nearest", cmap=cmap,
                            extent=(0.5, self.nloci + 0.5, 0.5, self.nloci + 0.5))
            ax.set_ylabel('Genomic bin')
            ax.set_xlabel('Genomic bin')
            ax.set_title('Normalized Hi-C count')
            cbar = ax.figure.colorbar(ims)
            cbar.ax.set_ylabel('Log2 (normalized Hi-C data)')
            if savefig:
                if dynamics:
                    tadbit_savefig(savefig+'/correlation_plot_stage_'+str(od)+'.png')
                else:
                    tadbit_savefig(savefig)
                
            elif not axe:
                plt.show()
            plt.close('all')
        return correl


    def view_centroid(self, **kwargs):
        """
        shortcut for
        view_models(tool='plot', show='highlighted', highlight='centroid')

        :param kwargs: any parameters to be passed to view_models (i.e.:
           view_centroid(azimuth=30, elevation=10, show_axe=True, label=True))
        """
        self.view_models(tool='plot', show='highlighted', highlight='centroid',
                         **kwargs)

    def view_models(self, models=None, cluster=None, stage=None, dynamics=None,
                    tool='chimera', show='all', highlight='centroid',
                    savefig=None, cmd=None, color='index', align=True, **kwargs):
        """
        Visualize a selected model in the three dimensions (either with Chimera
        or through matplotlib).

        :param None models:  if None (default) the visualization will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the visualization only for the models in the
           cluster number 'cluster'
        :param None stage: compute the visualization only for the models in
            stage number 'stage'
        :param None dynamics: compute the visualization for all the stages of the
            replica number 'dynamics'
        :param 'chimera' tool: path to the external tool used to visualize the
           model. Can also be 'plot', to use matplotlib.
        :param None savefig: path to a file where to save the image OR movie
           generated (depending on the extension; accepted formats are png, mov
           and webm). if set to None, the image or movie will be shown using
           the default GUI.
        :param 'index' color: can be:

             * a string as:
                 * '**index**' to color particles according to their position in the
                   model (:func:`pytadbit.utils.extraviews.color_residues`)
                 * '**tad**' to color particles according to the TAD they belong to
                   (:func:`pytadbit.utils.extraviews.tad_coloring`)
                 * '**border**' to color particles marking borders. Color according to
                   their score (:func:`pytadbit.utils.extraviews.tad_border_coloring`)
                   coloring function like.
             * a function, that takes as argument a model and any other parameter
               passed through the kwargs.
             * a list of (r, g, b) tuples (as long as the number of particles).
               Each r, g, b between 0 and 1.
        :param 'centroid' highlight: higlights a given model, or group of models.
           Can be either 'all', 'centroid' or 'best' ('best' being the model
           with the lowest IMP objective function value
        :param 'all' show: models to be displayed. Can be either 'all', 'grid'
           or 'highlighted'.

        :param None cmd: list of commands to be passed to the viewer. The chimera list is:

           ::

             focus
             set bg_color white
             windowsize 800 600
             bonddisplay never #0
             represent wire
             shape tube #0 radius 5 bandLength 100 segmentSubdivisions 1 followBonds on
             clip yon -500
             ~label
             set subdivision 1
             set depth_cue
             set dc_color black
             set dc_start 0.5
             set dc_end 1
             scale 0.8

           Followed by the movie command to record movies:

           ::

             movie record supersample 1
             turn y 3 120
             wait 120
             movie stop
             movie encode output SAVEFIG

           Or the copy command for images:

           ::

             copy file SAVEFIG png

           Passing as the following list as 'cmd' parameter:
           ::

             cmd = ['focus', 'set bg_color white', 'windowsize 800 600', 'bonddisplay never #0', 'shape tube #0 radius 10 bandLength 200 segmentSubdivisions 100 followBonds on', 'clip yon -500', '~label', 'set subdivision 1', 'set depth_cue', 'set dc_color black', 'set dc_start 0.5', 'set dc_end 1', 'scale 0.8']

           will return the default image (other commands can be passed to
           modified the final image/movie).
        :param True align: show aligned models
        :param 15 radius: radius for the chimera particles
        :param kwargs: see :func:`pytadbit.utils.extraviews.plot_3d_model` or
           :func:`pytadbit.utils.extraviews.chimera_view` for other arguments
           to pass to this function. See also coloring function


        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        elif stage > -1 and stage in self.stages:
            models = [m for m in self.stages[stage]]
        elif dynamics > -1:
            models = self.stages[0] + [self.stages[s+1][dynamics] for s in xrange(len(self.stages)-1)]
        else:
            models = [m for m in self.__models]
        models = [m['rand_init'] if 'IMPmodel' in str(type(m))
                  else m for m in models]
        if color in ['tad', 'border'] and 'tads' not in kwargs:
            start = (float(self[models[0]]['description']['start']) /
                     self[models[0]]['description']['resolution'] - 1)
            end   = (float(self[models[0]]['description']['end'  ]) /
                     self[models[0]]['description']['resolution'])
            kwargs.update((('tads', self.experiment.tads),
                           ('mstart', start ), ('mend', end)))
        centroid_model = models[0]
        if 'centroid' in [show, highlight] and len(models) > 1:
            centroid_model = self.centroid_model(models)
        if highlight == 'centroid':
            mdl = centroid_model
        elif highlight == 'best':
            mdl = self[sorted(models, key=lambda x:
                              self[x]['objfun'])[0]]['index']
        else:
            if highlight != 'all':
                warn("WARNING: represent_model value should be one of" +
                     "'centroid', 'best' or 'all' not %s\n"  % (
                         highlight) + "Highlighting no models.")
            mdl = 'all'
        # View with Matplotlib
        if tool == 'plot':
            pltshow = 'axe' not in kwargs
            model_coords = []
            if len(models) > 1 and align:
                for model in self.align_models(models, **kwargs):
                    model_coords.append(model)
            elif kwargs.get('reference_model', None) is not None and align:
                model_coords.append(self.align_models(models, **kwargs)[0])
            else:
                for model in models:
                    model_coords.append((
                        self[model]['x'], self[model]['y'], self[model]['z']))
            if show in ['all', 'highlighted']:
                if 'axe' not in kwargs:
                    fig = plt.figure(figsize=kwargs.get('figsize', (8, 8)))
                    kwargs['axe'] = fig.add_subplot(1, 1, 1, projection='3d')
                for i in models:
                    if show == 'all' or i == mdl or mdl == 'all':
                        plot_3d_model(
                            *model_coords[models.index(i)], color=color,
                            thin=False if highlight == 'all' else (i != mdl),
                            **kwargs)
                if pltshow:
                    try:
                        kwargs['axe'].set_title('Model %s highlighted as %s' % (
                            self[mdl]['rand_init'], highlight))
                    except KeyError:
                        kwargs['axe'].set_title('All models highlighted')
            else:
                sqrmdl = sqrt(len(models))
                cols = int(round(sqrmdl +
                                 (0.0 if int(sqrmdl) == sqrmdl else .5)))
                rows = int(sqrmdl + .5)
                if pltshow:
                    fig = plt.figure()
                for i in range(cols):
                    for j in range(rows):
                        if i * rows + j >= len(models):
                            break
                        this = self[models[i * rows + j]]['index']
                        if pltshow:
                            kwargs['axe'] = fig.add_subplot(
                                rows, cols, i * rows + j + 1, projection='3d')
                        plot_3d_model(
                            *model_coords[i * rows + j], color=color,
                            thin=False if highlight == 'all' else (this != mdl),
                            **kwargs)
                        if pltshow:
                            kwargs['axe'].set_title(
                                'Model %s' % self[this]['rand_init'])
            if savefig:
                tadbit_savefig(savefig)
            elif pltshow:
                plt.show()
            return
        # View with Chimera
        cmm_files = []
        radius = 10
        for model_num in models:
            if show in ['all', 'grid'] or model_num == mdl or mdl == 'all':
                self.write_cmm('/tmp/', model_num=model_num, color=color, **kwargs)
                cmm_files.append('/tmp/model.%s.cmm' % (self[model_num]['rand_init']))
                radius = self[model_num]['radius']
        chimera_view(cmm_files,
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd,
                     highlight=(0 if (show == 'highlighted' and mdl != 'all')
                                else models.index(mdl) if mdl != 'all' else mdl),
                     align=align, grid=show == 'grid', radius=radius)

    def angle_between_3_particles(self, parta, partb, partc,
                                  models=None, cluster=None,
                                  radian=False, all_angles=False):
        """
        Calculates the angle between 3 particles.


        Given three particles A, B and C, the angle g (angle ACB, shown below):

        ::


                              A
                             /|
                            /i|
                          c/  |
                          /   |
                         /    |
                        B )g  |b
                         \    |
                          \   |
                          a\  |
                            \h|
                             \|
                              C

        is given by the theorem of Al-Kashi:

        .. math::

          b^2 = a^2 + c^2 - 2ac\cos(g)

        :param parta: A particle number
        :param partb: A particle number
        :param partc: A particle number
        :param None models:  if None (default) the angle will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the angle only for the models in the
           cluster number 'cluster'
        :param False radian: if True, return value in radians (in degrees
           otherwise)

        :returns: an angle, either in degrees or radians. If all_angles is true
           returns a list of the angle g, h, i (see picture above)

        """
        # WARNING: here particle numbers are +1, they will be reduced
        # inside median_3d_dist
        a2 = np_median(self.__square_3d_dist(partb, partc, models=models,
                                             cluster=cluster))
        c2 = np_median(self.__square_3d_dist(parta, partb, models=models,
                                             cluster=cluster))
        b2 = np_median(self.__square_3d_dist(parta, partc, models=models,
                                             cluster=cluster))
        a = a2**0.5
        c = c2**0.5

        try:
            g = acos((a2 - b2 + c2) / (2 * a * c))
        except ValueError:
            g = 0.

        if not all_angles:
            return g if radian else degrees(g)

        b = b2**0.5
        try:
            h = acos((a2 + b2 - c2) / (2 * a * b))
        except ValueError:
            h = 0.

        i = pi - g - h

        return (g, h, i)

    def particle_coordinates(self, part, models=None, cluster=None):
        """
        Returns the mean coordinate of a given particle in a group of models.

        :param part: the index number of a particle
        :param None models:  if None (default) the angle will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the angle only for the models in the
           cluster number 'cluster'
        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]

        part -= 1
        xis = 0.
        yis = 0.
        zis = 0.
        for mod in models:
            xis += self[mod]['x'][part]
            yis += self[mod]['y'][part]
            zis += self[mod]['z'][part]
        xis /= len(models)
        yis /= len(models)
        zis /= len(models)

        return [xis, yis, zis]

    def dihedral_angle(self, pa, pb, pc, pd, pe, models):
        """
        Calculates the dihedral angle between 2 planes formed by 5 particles
           (one common to both planes).

        :param None models:  if None (default) the angle will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the angle only for the models in the
           cluster number 'cluster'
        """
        pa -= 1
        pb -= 1
        pc -= 1
        pd -= 1
        pe -= 1
        dangles = []
        for m in models:
            parta = array([self[m]['x'][pa], self[m]['y'][pa], self[m]['z'][pa]])
            partb = array([self[m]['x'][pb], self[m]['y'][pb], self[m]['z'][pb]])
            partc = array([self[m]['x'][pc], self[m]['y'][pc], self[m]['z'][pc]])
            partd = array([self[m]['x'][pd], self[m]['y'][pd], self[m]['z'][pd]])
            parte = array([self[m]['x'][pe], self[m]['y'][pe], self[m]['z'][pe]])
            dangles.append(dihedral(parta, partb, partc, partd, parte))
        return dangles

    def median_3d_dist(self, part1, part2, models=None, cluster=None,
                       plot=True, median=True, axe=None, savefig=None):
        """
        Computes the median distance between two particles over a set of models

        :param part1: number corresponding to the first particle
        :param part2: number corresponding to the second particle
        :param None models:  if None (default) the distance will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the distance only for the models in the
           cluster number 'cluster'
        :param True plot: if True, display a histogram and a box-plot of the
           distribution of the calculated distances. If False, return either
           the full list of the calculated distances or their median value
        :param True median: return either the full list of the calculated
           distances (False) or their median value (True), when 'plot' is set
           to False

        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :returns: if 'plot' is False, return either the full list of the
           calculated distances or their median value distances, either the
           list of distances.
        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        models = [self[mdl] for mdl in models]
        dists = [mdl.distance(part1, part2) for mdl in models]
        if not plot:
            if median:
                return np_median(dists)
            else:
                return dists
        plot_hist_box(dists, part1, part2, axe, savefig)

    def __square_3d_dist(self, part1, part2, models=None, cluster=None):
        """
        same as median_3d_dist, but return the square of the distance instead
        """
        part1 -= 1
        part2 -= 1
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        models = [self[mdl] for mdl in models]
        return [(mdl['x'][part1] - mdl['x'][part2])**2 +
                (mdl['y'][part1] - mdl['y'][part2])**2 +
                (mdl['z'][part1] - mdl['z'][part2])**2
                for mdl in models]

    def __fast_square_3d_dist(self, part1, part2, models):
        """
        same as median_3d_dist, but return the square of the distance instead
        we know what we are doing.
        """
        return [(mdl['x'][part1] - mdl['x'][part2])**2 +
                (mdl['y'][part1] - mdl['y'][part2])**2 +
                (mdl['z'][part1] - mdl['z'][part2])**2
                for mdl in models]

    def objective_function_model(self, model, log=False, smooth=True, axe=None,
                                 savefig=None):
        """
        This function plots the objective function value per each Monte-Carlo
        step

        :param model: the number of the model to plot
        :param False log: log plot
        :param True smooth: curve smoothing
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        """
        self[model].objective_function(log=log, smooth=smooth, axe=axe,
                                       savefig=savefig)

    def write_cmm(self, directory, model_num=None, models=None, cluster=None,
                  color='index', rndname=True, **kwargs):
        """
        Save a model in the cmm format, read by Chimera
        (http://www.cgl.ucsf.edu/chimera).

        **Note:** If none of model_num, models or cluster parameter are set,
        ALL the models will be written.

        :param directory: location where the file will be written (note: the
           name of the file will be model_1.cmm if model number is 1)
        :param None model_num: the number of the model to save
        :param None models: a list of numbers corresponding to a given set of
           models to save
        :param None cluster: save the models in the cluster number 'cluster'
        :param True rndname: If True, file names will be formatted as:
           model.RND.cmm, where RND is the random number feed used by IMP to
           generate the corresponding model. If False, the format will be:
           model_NUM_RND.cmm where NUM is the rank of the model in terms of
           objective function value
        :param 'index' color: can be:

             * a string as:
                 * '**index**' to color particles according to their position in the
                   model (:func:`pytadbit.utils.extraviews.color_residues`)
                 * '**tad**' to color particles according to the TAD they belong to
                   (:func:`pytadbit.utils.extraviews.tad_coloring`)
                 * '**border**' to color particles marking borders. Color according to
                   their score (:func:`pytadbit.utils.extraviews.tad_border_coloring`)
                   coloring function like.
             * a function, that takes as argument a model and any other parameter
               passed through the kwargs.
             * a list of (r, g, b) tuples (as long as the number of particles).
               Each r, g, b between 0 and 1.
        :param kwargs: any extra argument will be passed to the coloring
           function
        """
        if model_num > -1:
            models = [model_num]
        elif models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        for model_num in models:
            try:
                model = self[model_num]
            except KeyError:
                model = self._bad_models[model_num]
            model.write_cmm(directory, color=color, rndname=rndname,
                            model_num=model_num, **kwargs)

    def write_json(self, filename, color='index', models=None, cluster=None,
                   title=None, **kwargs):
        """
        Save a model in the json format, read by TADkit.

        **Note:** If none of model_num, models or cluster parameter are set,
        ALL the models will be written.

        :param directory: location where the file will be written (note: the
           name of the file will be model_1.cmm if model number is 1)
        :param None model_num: the number of the model to save
        :param True rndname: If True, file names will be formatted as:
           model.RND.cmm, where RND is the random number feed used by IMP to
           generate the corresponding model. If False, the format will be:
           model_NUM_RND.cmm where NUM is the rank of the model in terms of
           objective function value
        :param None filename: overide the default file name writing
        :param 'index' color: can be:

             * a string as:
                 * '**index**' to color particles according to their position in the
                   model (:func:`pytadbit.utils.extraviews.color_residues`)
                 * '**tad**' to color particles according to the TAD they belong to
                   (:func:`pytadbit.utils.extraviews.tad_coloring`)
                 * '**border**' to color particles marking borders. Color according to
                   their score (:func:`pytadbit.utils.extraviews.tad_border_coloring`)
                   coloring function like.
             * a function, that takes as argument a model and any other parameter
               passed through the kwargs.
             * a list of (r, g, b) tuples (as long as the number of particles).
               Each r, g, b between 0 and 1.
        :param kwargs: any extra argument will be passed to the coloring
           function
        """
        if isinstance(color, str):
            if color == 'index':
                color = color_residues(self, **kwargs)
            elif color == 'tad':
                if 'tads' not in kwargs:
                    raise Exception('ERROR: missing TADs\n   ' +
                                    'pass an Experiment.tads disctionary\n')
                color = tad_coloring(self, **kwargs)
            elif color == 'border':
                if 'tads' not in kwargs:
                    raise Exception('ERROR: missing TADs\n   ' +
                                    'pass an Experiment.tads disctionary\n')
                color = tad_border_coloring(self, **kwargs)
            else:
                raise NotImplementedError(('%s type of coloring is not yet ' +
                                           'implemeted\n') % color)
        elif hasattr(color, '__call__'):  # it's a function
            color = color(self, **kwargs)
        elif not isinstance(color, list):
            raise TypeError('one of function, list or string is required\n')
        form = '''
{
        "metadata" : {
                "version"  : 1.0,
                "type"     : "dataset",
                "generator": "TADbit"
                },
        "object": {\n%(descr)s
                   "uuid": "%(sha)s",
                   "title": "%(title)s",
                   "bp_per_nm": %(scale)s,
                   "datatype": "xyz",
                   "components": 3,
                   "source": "local",
                   "dependencies": %(dep)s
                  },
        "models":
                 [\n%(xyz)s
                 ],
        "clusters":%(cluster)s,
        "centroids":%(centroid)s,
        "restraints": %(restr)s,
        "hic_data": { "data": {'''
        form_end = '''}, "n": %(len_hic_data)i , "tads": [%(tad_def)s]}
}
'''
        fil = {}
        fil['title']   = title or "Sample TADbit data"
        fil['scale']   = str(self._config['scale'])
        versions = get_dependencies_version(dico=True)
        fil['dep'] = str(dict([(k.strip(), v.strip()) for k, v in versions.items()
                               if k in ['  TADbit', 'IMP', 'MCL']])).replace("'", '"')
        ukw = 'UNKNOWN'

        def tocamel(key):
            key = ''.join([(w[0].upper() + w[1:]) if i else w
                           for i, w in enumerate(key.split('_'))])
            key = ''.join([(w[0].upper() + w[1:]) if i else w
                           for i, w in enumerate(key.split(' '))])
            return key
        try:
            try:
                my_descr = dict(self.description)
            except TypeError:
                my_descr = {}
            if not my_descr.get('start', 0):
                my_descr['start'] = 0
            if not my_descr.get('end', 0):
                my_descr['end'  ] = self.nloci
            my_descr['chrom'] = my_descr['chromosome'] if 'chromosome' in my_descr and isinstance(my_descr['chromosome'], list) else ["%s" % (my_descr.get('chromosome', 'Chromosome'))]
            if 'chromosome' in my_descr:
                del my_descr['chromosome']
            if 'chrom_start' not in my_descr:
                warn("WARNING: chrom_start variable wasn't set, setting it to" +
                     " the position in the experiment matrix (%s)" % (
                         str([m for m in my_descr['start']])))
                my_descr['chrom_start'] = [m for m in my_descr['start']]
            if 'chrom_end' not in my_descr:
                warn("WARNING: chrom_end variable wasn't set, setting it to" +
                     " the position in the experiment matrix (%s)" % (
                         str([m for m in my_descr['end']])))
                my_descr['chrom_end'] = [m for m in my_descr['end']]
            if not my_descr['species']:
                warn("WARNING: species wasn't set, The resulting JSON will not work properly in TADkit.")
            # coordinates inside an array in case different models
            # from different places in the genome
            if not isinstance(my_descr['chrom_start'],list):
                my_descr['chrom_start'] = [my_descr['chrom_start']]
                my_descr['chrom_end'  ] = [my_descr['chrom_end'  ]]

            fil['descr']   = ',\n'.join([
                (' ' * 19) + '"%s" : %s' % (tocamel(k),
                                            ('"%s"' % (v))
                                            if not ((isinstance(v, int) and not isinstance(v, bool)) or
                                                    isinstance(v, list) or
                                                    isinstance(v, float))
                                            else str(v).replace("'", '"'))
                for k, v in my_descr.items()])
            if fil['descr']:
                fil['descr'] += ','
        except AttributeError:
            fil['descr']   = '"description": "Just some models"'
        
        if self.__models:
            aligned_coords = self.align_models(models=models, cluster=cluster)
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        fil['xyz'] = []
        for m_idx in xrange(len(models)):
            m = models[m_idx]
            model = {'rand_init':self[models[m_idx]]['rand_init'],'x':aligned_coords[m_idx][0],'y':aligned_coords[m_idx][1],'z':aligned_coords[m_idx][2]}
            fil['xyz'].append((' ' * 18) + '{"ref": %s,"data": [' % (
                model['rand_init']) + ','.join(
                    ['%.0f,%.0f,%.0f' % (model['x'][i],
                                         model['y'][i],
                                         model['z'][i])
                     for i in xrange(len(model['x']))]) + ']}')
        fil['xyz'] = ',\n'.join(fil['xyz'])
        # creates a UUID for this particular set of coordinates AND for TADbit version
        fil['sha'] = str(uuid5(UUID(md5(versions['  TADbit']).hexdigest()),
                               fil['xyz']))
        try:
            fil['restr']  = '[' + ','.join(['[%s,%s,"%s",%f]' % (
                k[0], k[1], self._restraints[k][0], self._restraints[k][2])
                for k in self._restraints]) + ']'
        except:
            fil['restr'] = '[]'
        fil['cluster'] = '[' + ','.join(['[' + ','.join(self.clusters[c]) + ']'
                                         for c in self.clusters]) + ']'
        fil['centroid'] = '[' + ','.join(
            [self[self.centroid_model(cluster=c)]['rand_init']
             for c in self.clusters]) + ']'
        fil['len_hic_data'] = len(self._original_data)
        try:
            fil['tad_def'] = ','.join(['['+','.join([str(i),str(self.experiment.tads[tad]['start']*self.resolution),
                                    str(self.experiment.tads[tad]['end']*self.resolution),
                                    str(self.experiment.tads[tad]['score'])])+']'
                                       for i,tad in enumerate(self.experiment.tads)
                                        if self.experiment.tads[tad]['start']*self.resolution >= my_descr['chrom_start'][0]
                                            and self.experiment.tads[tad]['end']*self.resolution <= my_descr['chrom_end'][0]])
        except:
            fil['tad_def'] = ''
        out_f = open(filename, 'w')
        out_f.write(form % fil)
        first = True
        for i, nrow in enumerate(self._original_data):
            for j, ncol in enumerate(nrow):
                if not isnan(ncol) and int(ncol) != 0:
                    if not first:
                        out_f.write(',')
                    first = False
                    if isinstance( ncol, ( int, long ) ):
                        out_f.write('"'+str((i*len(nrow))+j)+'":'+str(ncol))
                    else:
                        out_f.write('"'+str((i*len(nrow))+j)+'":'+"{:2.6f}".format(ncol))
        out_f.write(form_end % fil)
        out_f.close()

    def write_xyz(self, directory, model_num=None, models=None, cluster=None,
                  get_path=False, rndname=True):
        """
        Writes a xyz file containing the 3D coordinates of each particle in the
        model.

        .. note::

          If none of model_num, models or cluster parameter are set,
          ALL the models will be written.

        :param directory: location where the file will be written (note: the
           file name will be model.1.xyz, if the model number is 1)
        :param None model_num: the number of the model to save
        :param None models: a list of numbers corresponding to a given set of
           models to be written
        :param None cluster: save the models in the cluster number 'cluster'
        :param True rndname: If True, file names will be formatted as:
           model.RND.xyz, where RND is the random number feed used by IMP to
           generate the corresponding model. If False, the format will be:
           model_NUM_RND.xyz where NUM is the rank of the model in terms of
           objective function value
        :param False get_path: whether to return, or not, the full path where
           the file has been written
        """
        if model_num > -1:
            models = [model_num]
        elif models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        for model_num in models:
            try:
                model = self[model_num]
            except KeyError:
                model = self._bad_models[model_num]
            path_f = model.write_xyz(directory, model_num=model_num,
                                     get_path=get_path, rndname=rndname)
        if get_path:
            return path_f

    def write_xyz_babel(self, directory, model_num=None, models=None, cluster=None,
                        get_path=False, rndname=True):
        """
        Writes a xyz file containing the 3D coordinates of each particle in the
        model using a file format compatible with babel
        (http://openbabel.org/wiki/XYZ_%28format%29).
        .. note::
          If none of model_num, models or cluster parameter are set,
          ALL the models will be written.
        :param directory: location where the file will be written (note: the
           file name will be model.1.xyz, if the model number is 1)
        :param None model_num: the number of the model to save
        :param None models: a list of numbers corresponding to a given set of
           models to be written
        :param None cluster: save the models in the cluster number 'cluster'
        :param True rndname: If True, file names will be formatted as:
           model.RND.xyz, where RND is the random number feed used by IMP to
           generate the corresponding model. If False, the format will be:
           model_NUM_RND.xyz where NUM is the rank of the model in terms of
           objective function value
        :param False get_path: whether to return, or not, the full path where
           the file has been written
        """
        if model_num > -1:
            models = [model_num]
        elif models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = [m for m in self.__models]
        for model_num in models:
            try:
                model = self[model_num]
            except KeyError:
                model = self._bad_models[model_num]
            path_f = model.write_xyz_babel(directory, model_num=model_num,
                                           get_path=get_path, rndname=rndname)
        if get_path:
            return path_f
    
    def get_persistence_length(self, begin=0, end=None, axe=None, savefig=None, savedata=None,
                               plot=True):
        """
        Calculates the persistence length (Lp) of given section of the model.
        Persistence length is calculated according to [Bystricky2004]_ :

        .. math::

          <R^2> = 2 \\times Lp^2 \\times (\\frac{Lc}{Lp} - 1 + e^{\\frac{-Lc}{Lp}})

        with the contour length as :math:`Lc = \\frac{d}{c}` where :math:`d` is
        the genomic dstance in bp and :math:`c` the linear mass density of the
        chromatin (in bp/nm).

        :param 0  begin: starting particle of the region to consider
        :param None end: ending particle of the region to consider


        :returns: (float) persistence length
        """
        if not end:
            end = self.nloci

        wloci = [i for i in xrange(self.nloci)]

        # Maximum genomic distance
        max_gen_dist=end-begin
        # Quantities for all models
        R2_all  = [0]*max_gen_dist
        R4_all  = [0]*max_gen_dist
        cnt_all = [0]*max_gen_dist

        for model in xrange(len(self.__models)):

            # Quantities within each model
            R2  = [0]*max_gen_dist
            R4  = [0]*max_gen_dist
            cnt = [0]*max_gen_dist

            x = [] ; y = [] ; z = []
            for locus in xrange(begin,end):
                x.append(self[model]['x'][locus])
                y.append(self[model]['y'][locus])
                z.append(self[model]['z'][locus])

            #Compute the contact matrix
            squared_distance_matrix = squared_distance_matrix_calculation_wrapper(x, y, z, max_gen_dist)

            # Compute the average R2 per single model
            for i, j in combinations(wloci, 2):
                R2[abs(i-j)]  += squared_distance_matrix[i][j]
                R4[abs(i-j)]  += (squared_distance_matrix[i][j]*squared_distance_matrix[i][j])
                cnt[abs(i-j)] += 1

            for i in xrange(max_gen_dist):
                if cnt[i] != 0:
                    R2[i] = R2[i] / cnt[i]
                    R4[i] = R4[i] / cnt[i]
                else:
                    R2[i] = 0.0
                    R4[i] = 0.0

            for i in xrange(max_gen_dist):
                R2_all[i]  += R2[i]
                R4_all[i]  += R4[i]
                cnt_all[i] += 1

        avgs     = []
        std_devs = []
        for i in xrange(max_gen_dist):
            if cnt_all[i] != 0:
                avg     = R2_all[i]/cnt_all[i]
                avg2    = R4_all[i]/cnt_all[i]
                std_dev = sqrt(avg2-avg*avg)
                avgs.append(avg)
                std_devs.append(std_dev)
            else:
                avgs.append(0)
                std_devs.append(0)

        x = linspace(0.0 , float(max_gen_dist), num=max_gen_dist, endpoint=False)
        persistence_length, pcov = curve_fit(R2_vs_L, x, avgs)

        # write a 3 column file with genomic_distance | R2 | std_dev_R2
        if savedata:
            out = open(savedata, 'w')
            out.write('#gen_dist\tR2\tstd_dev_R2\n')
            out.write('#persistence_length = %f\n' % (persistence_length))
            for gen_dist in xrange(max_gen_dist):
                out.write('%f\t%f\t%f\n' % (gen_dist,avgs[gen_dist],std_devs[gen_dist]))
            out.close()

        if not plot:
            return persistence_length[0]

        # If plot we do the plot of R2 vs L
        show = False if axe else True
        axe = setup_plot(axe)
        plots = []
        plots = axe.plot(range(0, self.nloci),
                         avgs, color='black')

        axe.set_xlim((0, max_gen_dist))
        axe.set_xlabel('Genomic distance (particle)')
        axe.set_ylim((0, max(avgs)))
        axe.set_ylabel('<$R^{2}$> $(nm^2)$')

        if savefig:
            tadbit_savefig(savefig)
        elif show:
            plt.show()
        plt.close('all')

        return persistence_length[0]

    def save_models(self, outfile):
        """
        Saves all the models in pickle format (python object written to disk).

        :param path_f: path where to save the pickle file
        """

        out = open(outfile, 'w')
        dump(self._reduce_models(), out, protocol=HIGHEST_PROTOCOL)
        out.close()

    def _reduce_models(self, minimal=False):
        """
        reduce strural models objects to a dictionary to be saved

        :param False minimal: do not save info about log_objfun decay nor
           zscores

        :returns: this dictionary
        """
        to_save = {}

        if minimal:
            for m in self.__models:
                self.__models[m]['log_objfun'] = None
            to_save['models']    = self.__models
        else:
            to_save['models']    = self.__models
        to_save['bad_models']    = self._bad_models
        to_save['description']   = self.description
        to_save['nloci']         = self.nloci
        to_save['clusters']      = self.clusters
        to_save['resolution']    = self.resolution
        to_save['original_data'] = self._original_data
        to_save['config']        = self._config
        to_save['zscore']        = {} if minimal else self._zscores
        to_save['restraints']    = {} if minimal else self._restraints
        to_save['zeros']         = self._zeros
        to_save['stages']         = self.stages
        to_save['models_per_step']= self.models_per_step

        return to_save

    def _get_models(self, models, cluster, stage=None):
        """
        Internal function to transform cluster name, model name, or model list
        into proper list of models processable by StructuralModels functions
        """
        if models:
            models = [m if isinstance(m, int) else self[m]['index']
                      if isinstance(m, str) else m['index'] for m in models]
        elif cluster > -1 and len(self.clusters) > 0:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        elif stage > -1 and stage in self.stages:
            models = [m for m in self.stages[stage]]
        else:
            models = [m for m in self.__models]
        return models

    def _windowize(self, dists, steps, average=True, interval=0, minerr=0.):
        lmodels = len(dists[0])
        distsk = {1: dists}
        for k in steps[1:] if steps[0] == 1 else steps:
            distsk[k] = [None for _ in range(k / 2 + interval)]
            for i in range(interval, self.nloci - k - interval + 1):
                distsk[k].append(reduce(lambda x, y: x + y,
                                        [dists[i + j] for j in range(k)]))
                if k == 1:
                    continue
                # calculate the mean for steps larger than 1
                distsk[k][-1] = [mean_none([distsk[k][-1][i + lmodels * j]
                                            for j in xrange(k)])
                                 for i in xrange(lmodels)]
        new_distsk = {}
        errorp     = {}
        errorn     = {}
        for k, dists in distsk.iteritems():
            new_distsk[k] = []
            errorp[k] = []
            errorn[k] = []
            for part in dists:
                if not part:
                    new_distsk[k].append(None)
                    errorp[k].append(None)
                    errorn[k].append(None)
                    continue
                try:
                    mean_part = (np_mean([p for p in part if p is not None]) if
                                 average else
                                 np_median([p for p in part if p is not None]))
                except IndexError:  # bug in new version of numpy?
                    mean_part = float('nan')
                new_distsk[k].append(mean_part)
                try:
                    errorn[k].append(mean_part - 2 * np_std(
                        [p for p in part if p is not None]))
                    if errorn[k][-1] < minerr:
                        errorn[k][-1] = minerr
                    errorp[k].append(mean_part + 2 * np_std(
                        [p for p in part if p is not None]))
                    if errorp[k][-1] < minerr:
                        errorp[k][-1] = minerr
                except TypeError:
                    errorn[k].append(None)
                    errorp[k].append(None)
        return new_distsk, errorn, errorp

    def _plot_polymer(self, axe):
        where = axe.get_ylim()[0]
        axe.scatter(range(1, self.nloci + 1), [where] * self.nloci,
                    s=40,
                    color=[(0.15, 0.15, 0.15) if i else (0.7, 0.7, 0.7)
                           for i in self._zeros], clip_on=False,
                    zorder=100, edgecolor='k')
        axe.set_ylim((where, axe.get_ylim()[1]))

    def _generic_per_particle_plot(self, steps, distsk, error, errorp, errorn,
                                   savefig, axe, xlabel='', ylabel='', title='',
                                   colors=None, ylim=None):
        if not colors:
            colors = ['grey', 'darkgreen', 'darkblue', 'purple', 'darkorange',
                      'darkred'][-len(steps):]
        if len(steps) > 6:
            raise Exception('Sorry not enough colors to do this :)')
        ax = setup_plot(axe)
        plots = []
        for k in steps:
            plots += ax.plot(range(1, len(distsk[k]) + 1),
                             distsk[k], color=colors[steps.index(k)],
                             lw=steps.index(k) + 1, alpha=0.5)
        if error:
            for k in steps:
                plots += ax.plot(range(1, len(errorp[k]) + 1),
                                 errorp[k], color=colors[steps.index(k)], ls='--')
                ax.plot(range(1, len(errorp[k]) + 1), errorn[k],
                        color=colors[steps.index(k)], ls='--')
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        self._plot_polymer(ax)
        scat1 = plt.Line2D((0, 1), (0, 0), color=(0.15, 0.15, 0.15), marker='o',
                           linestyle='')
        scat2 = plt.Line2D((0, 1), (0, 0), color=(0.7,  0.7,  0.7 ), marker='o',
                           linestyle='')
        try:
            ax.legend(plots + [scat1, scat2],
                      ['Average for %s particle%s' % (k, 's' if k else '')
                       for k in steps] + (
                ['+/- 2 standard deviations'
                 for k in steps] if error else []) +
                ['particles with restraints', 'particles without restraints'],
                fontsize='small',
                numpoints=1, bbox_to_anchor=(1, 0.5), loc='center left')
        except TypeError:
            ax.legend(plots + [scat1, scat2],
                      ['Average for %s particle%s' % (k, 's' if k else '')
                       for k in steps] + (
                ['+/- 2 standard deviations'
                 for k in steps] if error else [] + ['particles with restraints',
                                                     'particles without restraints']),
                numpoints=1, bbox_to_anchor=(1, 0.5), loc='center left')
        ax.set_xlim((1, self.nloci))
        if ylim:
            a, b = ax.get_ylim()
            ax.set_ylim((max(0, a), min(1, b)))
        ax.set_title(title)
        plt.subplots_adjust(left=0.1, right=0.77)
        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()
        plt.close('all')


class ClusterOfModels(dict):
    def __str__(self):
        out1 = '   Cluster #%s has %s models [top model: %s]\n'
        out = 'Total number of clusters: %s\n%s' % (
            len(self),
            ''.join([out1 % (k, len(self[k]), self[k][0]) for k in self]))
        return out
