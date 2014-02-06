"""
19 Jul 2013


"""
from pytadbit.utils.three_dim_stats import calc_consistency, mass_center
from pytadbit.utils.three_dim_stats import dihedral, calc_eqv_rmsd
from pytadbit.utils.tadmaths        import calinski_harabasz
from pytadbit.utils.extraviews      import plot_3d_model
from pytadbit.utils.extraviews      import chimera_view, tadbit_savefig
from pytadbit.utils.extraviews      import augmented_dendrogram, plot_hist_box
from pytadbit.imp.impmodel          import IMPmodel
from pytadbit.centroid              import centroid_wrapper
from pytadbit.aligner3d             import aligner3d_wrapper
from cPickle                        import load, dump
from subprocess                     import Popen, PIPE
from math                           import acos, degrees, pi, sqrt
from numpy                          import median as np_median
from numpy                          import std as np_std, log2
from numpy                          import array, cross, dot, ma, isnan
from numpy.linalg                   import norm
from scipy.cluster.hierarchy        import linkage, fcluster
from scipy.stats                    import spearmanr
from warnings                       import warn
from string                         import uppercase as uc, lowercase as lc
from random                         import random
from os.path                        import exists

try:
    from matplotlib import pyplot as plt
    from matplotlib.cm import jet, bwr
except ImportError:
    warn('matplotlib not found\n')


def load_structuralmodels(path_f):
    """
    Loads :class:`pytadbit.imp.structuralmodels.StructuralModels` from a file
    (generated with
    :class:`pytadbit.imp.structuralmodels.StructuralModels.save_models`).
    
    :param path: to the pickled StructuralModels object.

    :returns: a :class:`pytadbit.imp.imp_model.StructuralModels`.
    """
    svd = load(open(path_f))
    return StructuralModels(
        nloci=svd['nloci'], models=svd['models'], bad_models=svd['bad_models'],
        resolution=svd['resolution'], original_data=svd['original_data'],
        clusters=svd['clusters'], config=svd['config'], zscores=svd['zscore'])


class StructuralModels(object):
    """
    This class contains three-dimensional models generated from a single Hi-C
    data. They can be reached either by their index (integer representing their
    rank according to objective function value), or by their IMP random intial
    number (as string).

    :param nloci: number of particles in the selected region
    :param models: a dictionary containing the generated
       :class:`pytadbit.imp.impmodel.IMPmodel` to be used as 'best models'
    :param bad_models: a dictionary of :class:`pytadbit.imp.impmodel.IMPmodel`,
       these model will not be used, just stored in case the set of
       'best models' needs to be extended later-on (
       :func:`pytadbit.imp.structuralmodels.StructuralModels.define_best_models`
       ).
    :param resolution: number of nucleotides per Hi-C bin. This will be the
       number of nucleotides in each model's particle.
    :param None original_data: a list of list (equivalent to a square matrix) of
       the normalized Hi-C data, used to build this list of models.
    :param None clusters: a dictionary of type
       :class:`pytadbit.imp.structuralmodels.ClusterOfModels`
    :param None config: a dictionary containing the parameter to be used for the
       generation of three dimensional models.

    """

    def __init__(self, nloci, models, bad_models, resolution,
                 original_data=None, zscores=None, clusters=None,
                 config=None, experiment=None):

        self.__models       = models
        self._bad_models    = bad_models
        self.nloci          = nloci
        self.clusters       = clusters or ClusterOfModels()
        self.resolution     = float(resolution)
        self._original_data = original_data # only used for correlation
        self._zscores       = zscores       # only used for plotting
        self._config        = config
        self.experiment     = experiment


    def __getitem__(self, nam):
        if type(nam) == str:
            nam = int(nam)
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
        return ('StructuralModels with %s models (objective function range: %s - %s)\n' +
                '   (corresponding to the best models out of %s models).\n' +
                '  IMP modeling used this parameters:\n' +
                '%s\n' +
                '  Models where clustered into %s clusters') % (
            len(self.__models),
            int(self.__models[0]['objfun']),
            int(self.__models[len(self.__models) - 1]['objfun']),
            len(self.__models) + len(self._bad_models),
            '\n'.join(['   - %-12s: %s' % (k, v)
                       for k, v in self._config.iteritems()]),
            len(self.clusters))


    def align_models(self, models=None, cluster=None, in_place=False):
        """
        Three-dimentional aligner for structural models.
        
        :param None models: if None (default) the average model will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the average model only for the models in the
           cluster number 'cluster'
        :param False in_place: if True, the result of the alignment will REPLACE
           the coordinates of each model. Default is too yield new coordinates of
           each model
        """
        if models:
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        firstx, firsty, firstz = (self[models[0]]['x'], self[models[0]]['y'],
                                  self[models[0]]['z'])
        for sec in models[1:]:
            coords = aligner3d_wrapper(firstx, firsty, firstz,
                                       self[sec]['x'],
                                       self[sec]['y'],
                                       self[sec]['z'],
                                       self.nloci)
            if in_place:
                models[sec]['x'], models[sec]['y'], models[sec]['z'] = coords
            else:
                yield coords

        if in_place:
            mass_center(self[models[0]]['x'], self[models[0]]['y'],
                        self[models[0]]['z'])
        else:
            x, y, z = (self[models[0]]['x'][:], self[models[0]]['y'][:],
                       self[models[0]]['z'][:])
            mass_center(x, y, z)
            yield x, y, z


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
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        idx = centroid_wrapper([self[x]['x'] for x in models],
                               [self[x]['y'] for x in models],
                               [self[x]['z'] for x in models],
                               self.nloci, len(models), int(verbose), 0)
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
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        idx = centroid_wrapper([self[x]['x'] for x in models],
                               [self[x]['y'] for x in models],
                               [self[x]['z'] for x in models],
                               self.nloci, len(models), int(verbose), 1)
        avgmodel = IMPmodel((('x', idx[0]), ('y', idx[1]), ('z', idx[2]),
                             ('rand_init', 'avg'), ('objfun', None),
                             ('radius',
                              self.resolution * self._config['scale'])))
        return avgmodel


    def cluster_models(self, fact=0.75, dcutoff=200, method='mcl',
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
        :param 200 dcutoff: distance threshold (nm) to determine if two
           particles are in contact
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
        scores = calc_eqv_rmsd(self.__models, self.nloci, dcutoff, what=what,
                               normed=True)
        from distutils.spawn import find_executable
        if not find_executable(mcl_bin):
            print('\nWARNING: MCL not found in path using WARD clustering\n')
            method = 'ward'
        # Initialize cluster definition of models:
        for model in self:
            model['cluster'] = 'Singleton'
        if method == 'ward':

            matrix = [[0.0 for _ in xrange(len(self))]
                      for _ in xrange(len(self))]
            for (i, j), score in scores.iteritems():
                matrix[i][j] = score if score > fact * self.nloci else 0.0
            clust = linkage(matrix, method='ward')
            # score each possible cut in hierarchical clustering
            solutions = {}
            for k in clust[:,2]:
                clusters = ClusterOfModels()
                _ = [clusters.setdefault(j, []).append(i) for i, j in
                     enumerate(fcluster(clust, k, criterion='distance'))]
                solutions[k] = {'out': clusters}
                solutions[k]['score'] = calinski_harabasz(scores, clusters)
            # take best cluster according to calinski_harabasz score
            clusters = [solutions[s] for s in sorted(
                solutions, key=lambda x: solutions[x]['score'])
                         if solutions[s]['score']>0][-1]['out']
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
                    self.clusters[cluster].append(self[model]['rand_init'])
                self.clusters[cluster].sort(
                    key=lambda x: self[str(x)]['objfun'])
        else:
            out_f = open(tmp_file, 'w')
            uniqs = list(set([tuple(sorted((m1, m2))) for m1, m2 in scores]))
            cut = fact * self.nloci
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
                raise Exception(
                    'Problem with clustering, try increasing "dcutoff"\n')
            new_singles = 0
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
                        clusters[cluster + 1].append(self[model]['rand_init'])
                    clusters[cluster + 1].sort(
                        key=lambda x: self[str(x)]['objfun'])
            if external:
                return clusters
            self.clusters = clusters
        if verbose:
            singletons = len([1 for model in self
                              if model['cluster'] == 'Singleton'])
            print ('Number of singletons excluded from clustering: %s (total' +
                   ' singletons: %s)') % (singletons - new_singles, singletons)
            print self.clusters


    def _build_distance_matrix(self, n_best_clusters):
        """
        """
        clusters = sorted(self.clusters.keys())[:n_best_clusters]
        matrix = [[0.0 for _ in clusters] for _ in clusters]
        clust_count = dict([(c, len([m for m in self if m['cluster']==c]))
                            for c in clusters])
        objfun = dict([(c, [m for m in self if m['cluster']==c][0]['objfun'])
                       for c in clusters])
        for i, cl1 in enumerate(clusters):
            md1 = md2 = None
            # find model with lowest energy for each cluster
            for md1 in self:
                if md1['cluster'] == cl1:
                    # the first on found is the best :)
                    break
            for j, cl2 in enumerate(clusters[i+1:]):
                # find model with lowest energy for each cluster
                for md2 in self:
                    if md2['cluster'] == cl2:
                        # the first on found is the best :)
                        break
                matrix[i][j+i+1] = calc_eqv_rmsd({0: md1, 1: md2}, self.nloci,
                                                 one=True)
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
        val = (maxnrj-minnrj)
        maxz = max([i[2] for i in z])
        for i in z:
            i[2] = i[2]/maxz * val

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


    def density_plot(self, models=None, cluster=None, steps=(1, 2, 3, 4, 5),
                     error=False, axe=None, savefig=None, savedata=None,
                     plot=True):
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
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param None savedata: path to a file where to save the density data
           generated (1 column per step + 1 for particle number).
        :param True plot: e.g. only saves data. No plotting done

        """
        if type(steps) == int:
            steps = (steps, )
        if len(steps) > 6:
            raise Exception('Sorry not enough colors to do this.\n')
        colors = ['grey', 'darkgreen', 'darkblue', 'purple', 'darkorange',
                  'darkred'][-len(steps):]
        dists = []
        for part1, part2 in zip(range(self.nloci - 1), range(1, self.nloci)):
            dists.append(self.median_3d_dist(part1 + 1, part2 + 1, models,
                                             cluster, plot=False, median=False))
        lmodels = len(dists[0])
        distsk = {1: dists}
        for k in (steps[1:] if steps[0]==1 else steps):
            distsk[k] = [None for _ in range(k/2)]
            for i in range(self.nloci - k):
                distsk[k].append(reduce(lambda x, y: x + y,
                                        [dists[i+j] for j in range(k)]))
                if k == 1:
                    continue
                # calculate the mean for steps larger than 1
                distsk[k][-1] = [float(sum([distsk[k][-1][i+lmodels*j]
                                            for j in xrange(k)])) / k
                                 for i in xrange(lmodels)]
        new_distsk = {}
        errorp    = {}
        errorn    = {}
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
                part = [self.resolution / p for p in part]
                new_distsk[k].append(np_median(part))
                try:
                    errorn[k].append(new_distsk[k][-1] - 2 * np_std(part))
                    errorn[k][-1] = errorn[k][-1] if errorn[k][-1] > 0 else 0.0
                    errorp[k].append(new_distsk[k][-1] + 2 * np_std(part))
                    errorp[k][-1] = errorp[k][-1] if errorp[k][-1] > 0 else 0.0
                except TypeError:
                    errorn[-1].append(None)
                    errorp[-1].append(None)
        distsk = new_distsk
        # write consistencies to file
        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle\t%s\n' % ('\t'.join([str(c) + '\t' +
            '2*stddev(%d)' % c for c in steps])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (part + 1, '\t'.join(
                    ['None\tNone' if part >= len(distsk[c]) else
                    (str(round(distsk[c][part], 3)) + '\t' +
                     str(round(errorp[c][part], 3)))
                     if distsk[c][part] else 'None\tNone'
                     for c in steps])))
            out.close()
        if not plot:
            return
        # plot
        if axe:
            ax = axe
            fig = ax.get_figure()
        else:
            fig = plt.figure(figsize=(11, 5))
            ax = fig.add_subplot(111)
            ax.patch.set_facecolor('lightgrey')
            ax.patch.set_alpha(0.4)
            ax.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
            ax.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
            ax.set_axisbelow(True)
            ax.minorticks_on() # always on, not only for log
            # remove tick marks
            ax.tick_params(axis='both', direction='out', top=False, right=False,
                           left=False, bottom=False)
            ax.tick_params(axis='both', direction='out', top=False, right=False,
                           left=False, bottom=False, which='minor')
        plots = []
        for k in steps:
            plots += ax.plot(range(1, len(distsk[k]) + 1), distsk[k],
                             color=colors[steps.index(k)],
                             lw=steps.index(k) + 1, alpha=0.5)
        if error:
            for k in steps:
                plots += ax.plot(range(1, len(errorp[k]) + 1), errorp[k],
                                 color=colors[steps.index(k)], ls='--')
                ax.plot(range(1, len(errorp[k]) + 1), errorn[k],
                        color=colors[steps.index(k)], ls='--')
        ax.set_ylabel('Density (bp / nm)')
        ax.set_xlabel('Particle number')
        try:
            ax.legend(plots, ['Average for %s particle%s' % (k, 's' if k else '')
                              for k in steps] + (
                          ['+/- 2 standard deviations'
                           for k in steps] if error else []), fontsize='small',
                      bbox_to_anchor=(1, 0.5), loc='center left')
        except TypeError:
            ax.legend(plots, ['Average for %s particle%s' % (k, 's' if k else '')
                              for k in steps] + (
                          ['+/- 2 standard deviations'
                           for k in steps] if error else []), 
                      bbox_to_anchor=(1, 0.5), loc='center left')
        ax.set_xlim((1, self.nloci))
        ax.set_title('Chromatin density')
        plt.subplots_adjust(left=0.1, right=0.78)
        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()


    def get_contact_matrix(self, models=None, cluster=None, cutoff=150):
        """
        Returns a matrix with the number of interactions observed below a given
        cutoff distance.

        :param None models: if None (default) the contact matrix will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact matrix only for the models in the
           cluster number 'cluster'
        :param 150 cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not

        :returns: matrix frequency of interaction
        """
        if models:
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        matrix = [[float('nan') for _ in xrange(self.nloci)]
                  for _ in xrange(self.nloci)]
        cutoff = cutoff**2
        for i in xrange(self.nloci):
            for j in xrange(i + 1, self.nloci):
                val = len([k for k in self.__square_3d_dist(
                    i + 1, j + 1, models=models)
                           if k < cutoff])
                matrix[i][j] = matrix[j][i] = float(val) / len(models)# * 100
        return matrix


    def define_best_models(self, nbest):
        """
        Defines the number of top models (based on the objective function) to
        keep. If keep_all is set to True in
        :func:`pytadbit.imp.imp_model.generate_3d_models` or in
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


    def deconvolve(self, fact=0.75, dcutoff=200, method='mcl',
                   mcl_bin='mcl', tmp_file=None, verbose=True, n_cpus=1,
                   mclargs=None, what='dRMSD', n_best_clusters=10,
                   savefig=None, represent_models=False, figsize=(11, 11),
                   **kwargs):
        """
        This function performs a clustering analysis of the generated models
        based on structural comparison (dRMSD).
        Then, performs a differential contact map between each possible pair
        of cluster.

        .. note::

          Clusters defined here are different from the one defined when using
          :func:`pytadbit.imp.structuralmodels.StructuralModels.cluster_models`.
          They are also not stored into StructuralModels.clusters

        :param 0.75 fact: factor to define the percentage of equivalent
           positions to be considered in the clustering
        :param 200 dcutoff: distance threshold (nm) to determine if two
           particles are in contact
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
        fact /= self.nloci
        clusters = self.cluster_models(fact=fact, dcutoff=dcutoff, mcl_bin=mcl_bin,
                                       method=method, tmp_file=tmp_file,
                                       n_cpus=n_cpus, mclargs=mclargs,
                                       external=True, what=what)
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
                axes[i,j].set(adjustable='box-forced', aspect=1)
                axes[i,j].set_visible(False)
        # doing the plot
        for i in xrange(n_best_clusters - 1):
            for j in xrange(1, n_best_clusters):
                if j < i+1:
                    continue
                axes[i+add,j-1].set_visible(True)
                axes[i+add,j-1].set(adjustable='box-forced', aspect=1)
                matrix3 = [[cmatrices[i][k][l] - cmatrices[j][k][l]
                            for l in xrange(self.nloci)]
                           for k in xrange(self.nloci)]
                ims = axes[i+add,j-1].imshow(matrix3, origin='lower', cmap=bwr,
                                         interpolation="nearest", vmin=-1, vmax=1,
                                         extent=(0.5, len(matrix3) + 0.5,
                                                 0.5, len(matrix3) + 0.5))
                axes[i+add, j-1].grid()
                if not i and not represent_models:
                    axes[i+add,j-1].set_title('Cluster #%s' % (j + 1),
                                              color='blue')
                if j != i+1:
                    axes[i+add,j-1].yaxis.set_ticks_position('none')
                else:
                    plt.setp(axes[i+add,j-1].get_yticklabels(), visible=True)
                    axes[i+add,j-1].yaxis.set_ticks_position('left')
                    for item in axes[i+add,j-1].get_yticklabels():
                        item.set_fontsize(9)
                if i != j-1:
                    axes[i+add,j-1].xaxis.set_ticks_position('none')
                else:
                    plt.setp(axes[i+add,j-1].get_xticklabels(), visible=True)
                    axes[i+add,j-1].xaxis.set_ticks_position('bottom')
                    for item in axes[i+add,j-1].get_xticklabels():
                        item.set_fontsize(9)
                if j == n_best_clusters - 1 and not represent_models:
                    axes[i+add,j-1].yaxis.set_label_position('right')
                    axes[i+add,j-1].set_ylabel('Cluster #%s' % (i + 1),
                                               rotation=-90, fontsize='large',
                                               color='red', va='bottom')
                axes[i+add,j-1].set_xlim((0.5, len(matrix3) + 0.5))
                axes[i+add,j-1].set_ylim((0.5, len(matrix3) + 0.5))
        # new axe for the color bar
        cell = fig.add_axes([0.125, 0.1, 0.01, 0.25])
        cbar = fig.colorbar(ims, cax=cell, cmap=jet)
        cbar.set_ticks([float(k)/100
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
                               i+1, len(clusters[i+1]))
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
                if represent_models=='centroid':
                    mdl = str(self.centroid_model(
                        models=[m for m in clusters[i+1]]))
                else:
                    if represent_models != 'best':
                        warn("WARNING: represent_model value should be one of" +
                             "'centroid' or 'best' not %s\n"  % (
                                 represent_models) + "Showing best model.")
                    mdl = str(clusters[i+1][0])
                self[mdl].view_model(tool='plot', axe=ax, **kwargs)
                for item in [ax]:
                    item.patch.set_visible(False)                
            for i in range(n_best_clusters - 1):
                ax = fig.add_subplot(n_best_clusters, n_best_clusters,
                                     n_best_clusters * (i + 2), projection='3d')
                self[str(clusters[i+1][0])].view_model(tool='plot', axe=ax,
                                                     **kwargs)
                ax.yaxis.set_label_position('top')
                ax.set_title('Cluster #%s' % (i + 1), rotation=-90,
                             fontsize='large', color='red', position=(1,.5),
                             va='center', ha='left')
                for item in [ax]:
                    item.patch.set_visible(False)                
            axes[0, n_best_clusters-1].set_visible(False)
        if savefig:
            tadbit_savefig(savefig)
        else:
            plt.show()


    def contact_map(self, models=None, cluster=None, cutoff=150, axe=None,
                    savefig=None, savedata=None):
        """
        Plots a contact map representing the frequency of interaction (defined
        by a distance cutoff) between two particles.

        :param None models: if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param 150 cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param None savedata: path to a file where to save the contact map data
           generated, in three columns format (particle1, particle2, percentage
           of models where these two particles are in contact)

        """
        matrix = self.get_contact_matrix(models, cluster, cutoff=cutoff)
        show = False
        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle1\tParticle2\tModels_percentage\n')
            for i in xrange(len(matrix)):
                for j in xrange(i+1, len(matrix)):
                    out.write('%s\t%s\t%s\n' % (i, j, matrix[i][j]))
            out.close()
            return # stop here, we do not want to display anything
        if not axe:
            fig = plt.figure(figsize=(8, 6))
            axe = fig.add_subplot(111)
            show=True
        else:
            fig = axe.get_figure()
        ims = axe.imshow(matrix, origin='lower', interpolation="nearest",
                         vmin=0, vmax=1,
                         extent=(0.5, self.nloci + 0.5, 0.5, self.nloci + 0.5))
        axe.set_ylabel('Particle')
        axe.set_xlabel('Particle')
        cbar = axe.figure.colorbar(ims)
        cbar.ax.set_yticklabels(['%3s%%' % (p) for p in range(0, 110, 10)])
        cbar.ax.set_ylabel('Percentage of models with particles closer than ' +
                           '%s nm' % (cutoff))
        axe.set_title('Contact map')
        if savefig:
            tadbit_savefig(savefig)
        elif show:
            plt.show()


    def zscore_plot(self, axe=None, savefig=None):
        """
        Generate 3 plots. Two heatmaps of the Z-scores used for modeling, one
        of which is binary showing in red Z-scores higher than upper cut-off;
        and in blue Z-scores lower than lower cut-off. Last plot is an histogram
        of the distribution of Z-scores, showing selected regions.

        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).

        """

        zsc_mtrx = reduce(lambda x, y: x + y, [[k] + self._zscores[k].keys()
                                               for k in self._zscores.keys()])
        max_bin = max([int(i) for i in zsc_mtrx])
        zsc_mtrx = [[float('nan') for _ in xrange(max_bin)] for _ in xrange(max_bin)]
        for i in xrange(max_bin):
            for j in xrange(max_bin):
                try:
                    zsc_mtrx[i][j] = self._zscores[str(i)][str(j)]
                except KeyError:
                    continue
        for i in xrange(max_bin):
            for j in xrange(max_bin):
                if zsc_mtrx[i][j]:
                    zsc_mtrx[j][i] = zsc_mtrx[i][j]
        for i in xrange(max_bin):
            for j in xrange(max_bin):
                if self._config['lowfreq'] < zsc_mtrx[i][j] < self._config['upfreq']:
                    zsc_mtrx[j][i] = float('nan')
                    zsc_mtrx[i][j] = float('nan')
        masked_array = ma.array (zsc_mtrx, mask=isnan(zsc_mtrx))
        cmap = jet
        cmap.set_bad('w', 1.)
        if not axe:
            fig = plt.figure(figsize=(25, 5.5))
        else:
            fig = axe.get_figure()
        ax = fig.add_subplot(131)
        ims = ax.imshow(masked_array, origin='lower',
                        interpolation="nearest", cmap=cmap)
        ax.set_ylabel('Particles')
        ax.set_xlabel('Particles')
        ax.set_title('Z-scores of the observed Hi-C count')
        cbar = ax.figure.colorbar(ims, cmap=cmap)
        cbar.ax.set_ylabel('Z-score value')
        #
        ax = fig.add_subplot(132)
        _, _, patches = ax.hist(
            reduce(lambda x, y: x+y, [self._zscores[v].values()
                                      for v in self._zscores.keys()]),
            bins=50)
        for thispatch in patches:
            color = ('grey' if self._config['lowfreq'] < thispatch.get_x() + thispatch.get_width() < self._config['upfreq'] 
                     else 'green')
            thispatch.set_facecolor(color)
            thispatch.set_alpha(0.7)
        ax.set_title('Histogram of Z-scores')
        ax.vlines(self._config['lowfreq'], 1, ax.get_ylim()[1], color='red',
                  linestyle='--')
        ax.vlines(self._config['upfreq'] , 1, ax.get_ylim()[1], color='red',
                  linestyle='--')
        #
        ax = fig.add_subplot(133)
        masked_array = ma.array (zsc_mtrx, mask=isnan(zsc_mtrx))
        for i in masked_array:
            for j in xrange(len(i)):
                try:
                    i[j].mask
                    continue
                except AttributeError:
                    pass
                if i[j] > self._config['upfreq']:
                    i[j] = 1
                elif i[j] < self._config['lowfreq']:
                    i[j] = -1
        ims = ax.imshow(masked_array, origin='lower', vmin=-1.2, vmax=1.2,
                        interpolation="nearest", cmap=cmap)
        ax.set_ylabel('Particles')
        ax.set_xlabel('Particles')
        ax.set_title('Binary representation of Z-scores\nred: > ' +
                     '%.2f; blue: < %.2f' % (self._config['upfreq'],
                                             self._config['lowfreq']))
        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()


    def correlate_with_real_data(self, models=None, cluster=None, cutoff=200,
                                 plot=False, axe=None, savefig=None):
        """
        Plots the result of a correlation between a given group of models and
        original Hi-C data.
        
        :param None models: if None (default) the correlation will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the correlation only for the models in the
           cluster number 'cluster'
        :param 200 cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param False plot: to display the plot
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance

        :returns: correlation coefficient rho, between the two
           matrices. A rho value greater than 0.7 indicates a very good
           correlation

        """
        model_matrix = self.get_contact_matrix(models=models, cluster=cluster,
                                               cutoff=cutoff)
        oridata = []
        moddata = []
        for i in xrange(len(self._original_data)):
            for j in xrange(i + 1, len(self._original_data)):
                if not self._original_data[i][j] > 0:
                    continue
                oridata.append(self._original_data[i][j])
                moddata.append(model_matrix[i][j])
        # corr = spearmanr(model_matrix, self._original_data, axis=None)
        corr = spearmanr(moddata, oridata)
        # corr = spearmanr(reduce(lambda x, y, : x + y, model_matrix),
        #                  reduce(lambda x, y, : x + y, self._original_data))
        # corr = corrcoef(moddata, oridata)[1]
        if not plot and not savefig:
            return corr
        if not axe:
            fig = plt.figure(figsize=(15, 5.5))
        else:
            fig = axe.get_figure()
        fig.suptitle('Correlation between normalized-real and modeled '
                     + 'contact maps (correlation=%.4f)' % (corr[0]),
                     size='x-large')
        ax = fig.add_subplot(121)
        self.contact_map(models, cluster, cutoff, axe=ax)
        ax = fig.add_subplot(122)
        ims = ax.imshow(log2(self._original_data), origin='lower',
                        interpolation="nearest",
                        extent=(0.5, self.nloci + 0.5, 0.5, self.nloci + 0.5))
        ax.set_ylabel('Particles')
        ax.set_xlabel('Particles')
        ax.set_title('Normalized Hi-C count')
        cbar = ax.figure.colorbar(ims)
        cbar.ax.set_ylabel('Log2 (normalized Hi-C data)')

        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()
        return corr


    def model_consistency(self, cutoffs=(50, 100, 150, 200), models=None,
                          cluster=None, axe=None, savefig=None, savedata=None,
                          plot=True):
        """
        Plots the particle consistency, over a given set of models, vs the
        modeled region bins. The consistency is a measure of the variability
        (or stability) of the modeled region (the higher the consistency value,
        the higher stability).

        :param (50,100,150,200) cutoffs: list of distance cutoffs (nm) used to
           compute the consistency. Two particle are considered consistent if
           their distance is less than the given cutoff
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
        if models:
            models = dict([(i, self[m]) for i, m in enumerate(models)])
        elif cluster > -1:
            models = dict([(i, self[str(m)]) for i, m in
                           enumerate(self.clusters[cluster])])
        else:
            models = self.__models
        consistencies = {}
        for cut in cutoffs:
            consistencies[cut] = calc_consistency(models, self.nloci, cut)
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
        show = False
        if not axe:
            fig = plt.figure(figsize=(11, 5))
            axe = fig.add_subplot(111)
            show=True
            axe.patch.set_facecolor('lightgrey')
            axe.patch.set_alpha(0.4)
            axe.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
            axe.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
            axe.set_axisbelow(True)
            axe.minorticks_on() # always on, not only for log
            # remove tick marks
            axe.tick_params(axis='both', direction='out', top=False,
                            right=False, left=False, bottom=False)
            axe.tick_params(axis='both', direction='out', top=False,
                            right=False, left=False, bottom=False,
                            which='minor')
        else:
            fig = axe.get_figure()
        # plot!
        # colors = ['grey', 'darkgreen', 'darkblue', 'purple', 'darkorange', 'darkred'][-len(cutoffs):]
        plots = []
        for i, cut in enumerate(cutoffs[::-1]):
            plots += axe.plot(range(1, self.nloci + 1),
                              consistencies[cut], color='darkred',
                              alpha= 1 - i / float(len(cutoffs)))
        try:
            axe.legend(plots, ['%s nm' % (k) for k in cutoffs[::-1]],
                       fontsize='small', loc='center left',
                       bbox_to_anchor=(1, 0.5))
        except TypeError:
            axe.legend(plots, ['%s nm' % (k) for k in cutoffs[::-1]],
                       loc='center left',
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
        if savefig:
            tadbit_savefig(savefig)
        elif show:
            plt.show()


    def view_centroid(self, **kwargs):
        """
        shortcut for
        models.view_models(tool='plot', show='stressed', stress='centroid')

        :param kwargs: any parameters to be passed to view_models (i.e.:
           models.view_centroid(azimuth=30, elevation=10, show_axe=True, label=True))
        """
        self.view_models(tool='plot', show='stressed', stress='centroid',
                         **kwargs)
        

    def view_models(self, models=None, cluster=None, tool='chimera',
                    show='all', stress='centroid', savefig=None,
                    cmd=None, color='index', align=True, **kwargs):
        """
        Visualize a selected model in the three dimensions (either with Chimera
        or through matplotlib).

        :param None models:  if None (default) the visualization will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the visualization only for the models in the
           cluster number 'cluster'
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
        :param 'centroid' stress: higlights a given model, or group of models.
           Can be either 'all', 'centroid' or 'best' ('best' being the model
           with the lowest IMP objective function value
        :param 'all' show: models to be displayed. Can be either 'all', 'grid'
           or 'stressed'.

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
        :param kwargs: see :func:`pytadbit.utils.extraviews.plot_3d_model` or
           :func:`pytadbit.utils.extraviews.chimera_view` for other arguments
           to pass to this function. See also coloring function

           
        """
        if models:
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        models = [m['rand_init'] if 'IMPmodel' in str(type(m))
                  else m for m in models]
        if color in ['tad', 'border'] and not 'tads' in kwargs:
            start = (float(self[models[0]]['description']['start']) /
                     self[models[0]]['description']['resolution'] - 1)
            end   = (float(self[models[0]]['description']['end'  ]) /
                     self[models[0]]['description']['resolution'])
            kwargs.update((('tads', self.experiment.tads),
                           ('mstart', start ), ('mend', end)))
        centroid_model = 0
        if 'centroid' in [show, stress] and len(models) > 1:
            centroid_model = self.centroid_model(models)
        if stress=='centroid':
            mdl = centroid_model
        elif stress=='best':
            mdl = self[sorted(models, key=lambda x:
                              self[x]['objfun'])[0]]['index']
        else:
            if stress != 'all':
                warn("WARNING: represent_model value should be one of" +
                     "'centroid', 'best' or 'all' not %s\n"  % (
                         stress) + "Stressing no models.")
            mdl = 'all'
        ## View with Matplotlib
        if tool == 'plot':
            model_coords = []
            if len(models) > 1 and align:
                for model in self.align_models(models):
                    model_coords.append(model)
            else:
                for model in models:
                    model_coords.append((
                        self[model]['x'], self[model]['y'],self[model]['z']))
            if show in ['all', 'stressed']:
                fig = plt.figure()
                axe = fig.add_subplot(1,1,1, projection='3d')
                for i in models:
                    if show=='all' or i==mdl or mdl=='all':
                        plot_3d_model(*model_coords[models.index(i)],
                                      axe=axe, color=color,
                                      thin=False if stress=='all' else (i!=mdl),
                                      **kwargs)
                try:
                    axe.set_title('Model %s stressed as %s' % (
                        self[mdl]['rand_init'], stress))
                except ValueError:
                    axe.set_title('All models stressed')
            else:
                sqrmdl = sqrt(len(models))
                cols = int(round(sqrmdl + (0.0 if int(sqrmdl)==sqrmdl else .5)))
                rows = int(sqrmdl+.5)
                fig = plt.figure()
                for i in range(cols):
                    for j in range(rows):
                        if i * rows + j >= len(models):
                            break
                        this = self[models[i * rows + j]]['index']
                        axe = fig.add_subplot(rows, cols, i * rows + j+1,
                                              projection='3d')
                        plot_3d_model(
                            *model_coords[i * rows + j], axe=axe, color=color,
                            thin=False if stress=='all' else (this!=mdl),
                            **kwargs)
                        axe.set_title(
                            'Model %s' % self[this]['rand_init'])
            if savefig:
                tadbit_savefig(savefig)
            else:
                plt.show()
            return
        ## View with Chimera
        cmm_files = []
        for model_num in models:
            if show in ['all', 'grid'] or model_num==mdl or mdl=='all':
                self.write_cmm('/tmp/', model_num=model_num, color=color, **kwargs)
                cmm_files.append('/tmp/model.%s.cmm' % (self[model_num]['rand_init']))
        chimera_view(cmm_files,
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd,
                     stress=(0 if (show=='stressed' and mdl!='all')
                             else models.index(mdl) if mdl!='all' else mdl),
                     align=align, grid=show=='grid')


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
        a = self.median_3d_dist(partb, partc, models=models,
                                cluster=cluster, plot=False)
        c = self.median_3d_dist(parta, partb, models=models,
                                cluster=cluster, plot=False)
        b = self.median_3d_dist(parta, partc, models=models,
                                cluster=cluster, plot=False)

        try:
            g = acos((a**2 - b**2 + c**2) / (2 * a * c))
        except ValueError:
            g = 0.

        if not all_angles:
            return g if radian else degrees(g)

        try:
            h = acos((a**2 + b**2 - c**2) / (2 * a * b))
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
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models

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


    def dihedral_angle(self, parta, partb, partc, partd, models=None,
                       cluster=None):
        """
        Calculates the dihedral angle between 2 planes formed by 4 particles.
        
        :param None models:  if None (default) the angle will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the angle only for the models in the
           cluster number 'cluster'
        """
        parta = array(self.particle_coordinates(parta, models, cluster))
        partb = array(self.particle_coordinates(partb, models, cluster))
        partc = array(self.particle_coordinates(partc, models, cluster))
        partd = array(self.particle_coordinates(partd, models, cluster))
        return dihedral(parta, partb, partc, partd)


    def walking_dihedral(self, models=None, cluster=None, steps=(1,3),
                         plot=True, savefig=None, axe=None):
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
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
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
        # plot
        if axe:
            ax = axe
            fig = ax.get_figure()
        else:
            fig = plt.figure(figsize=(11, 5))
            ax = fig.add_subplot(111)
            ax.patch.set_facecolor('lightgrey')
            ax.patch.set_alpha(0.4)
            ax.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
            ax.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
            ax.set_axisbelow(True)
            ax.minorticks_on() # always on, not only for log
            # remove tick marks
            ax.tick_params(axis='both', direction='out', top=False, right=False,
                           left=False, bottom=False)
            ax.tick_params(axis='both', direction='out', top=False, right=False,
                           left=False, bottom=False, which='minor')
        colors = ['grey', 'darkgreen', 'darkblue', 'purple', 'darkorange',
                  'darkred'][-len(steps):]
        #
        rads = {}
        rads[1] = []
        for res in xrange(self.nloci - 6):
            rads[1].append(self.dihedral_angle(res + 1, res + 4,
                                               res + 5, res + 7,
                                               models=models, cluster=cluster))
        lmodels = len(rads[1])
        for k in (steps[1:] if steps[0]==1 else steps):
            rads[k] = [None for _ in range(k/2)]
            for i in range(1, self.nloci - k - 5):
                rads[k].append(reduce(lambda x, y: x + y,
                                      [rads[1][i+j] for j in range(k)]) / k)
                if k == 1:
                    continue
        plots = []
        for k in steps:
            plots += ax.plot(range(1, len(rads[k]) + 1), rads[k],
                             color=colors[steps.index(k)],
                             lw=steps.index(k) + 1, alpha=0.5)

        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()


    def walking_angle(self, models=None, cluster=None, steps=(1,3), signed=True,
                      savefig=None, savedata=None, axe=None):
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
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
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
        # plot
        if axe:
            ax = axe
            fig = ax.get_figure()
        else:
            fig = plt.figure(figsize=(11, 5))
            ax = fig.add_subplot(111)
            ax.patch.set_facecolor('lightgrey')
            ax.patch.set_alpha(0.4)
            ax.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
            ax.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
            ax.set_axisbelow(True)
            ax.minorticks_on() # always on, not only for log
            # remove tick marks
            ax.tick_params(axis='both', direction='out', top=False, right=False,
                           left=False, bottom=False)
            ax.tick_params(axis='both', direction='out', top=False, right=False,
                           left=False, bottom=False, which='minor')
        colors = ['grey', 'darkgreen', 'darkblue', 'purple', 'darkorange',
                  'darkred'][-len(steps):]
        #
        if not type(steps) == tuple:
            steps = (steps,)
        rads = {}
        rads[1] = []
        sign = 1
        for res in xrange(self.nloci - 6):
            rads[1].append(self.angle_between_3_particles(res + 1, res + 4,
                                                          res + 7,
                                                          models=models,
                                                          cluster=cluster))
            if signed:
                res1 = self.particle_coordinates(res+1)
                res2 = self.particle_coordinates(res+4)
                res3 = self.particle_coordinates(res+7)
                vec1 = array(res1) - array(res2) / norm(array(res1) - array(res2))
                vec2 = array(res1) - array(res3) / norm(array(res1) - array(res3))
                sign = dot(array([1.,1.,1.]), cross(vec1, vec2))
                sign = -1 if sign < 0 else 1
            rads[1][-1] *= sign
        for k in (steps[1:] if steps[0]==1 else steps):
            rads[k] = [None for _ in range(k/2)]
            for i in range(1, self.nloci - k - 5):
                rads[k].append(reduce(lambda x, y: x + y,
                                      [rads[1][i+j] for j in range(k)]) / k)
                if k == 1:
                    continue
        plots = []
        for k in steps:
            plots += ax.plot(range(1, len(rads[k]) + 1), rads[k],
                             color=colors[steps.index(k)],
                             lw=steps.index(k) + 1, alpha=0.5)
        if savedata:
            out = open(savedata, 'w')
            out.write('#Particle\t' +
                      '\t'.join(['angle(step:%s)' % (s) for s in steps]) + '\n')
            for p in xrange(len(rads[1])):
                out.write(str(p+1))
                for s in steps:
                    try:
                        out.write('\t%s' % rads[s][p])
                    except IndexError:
                        out.write('\tNone')
                out.write('\n')
            out.close()

        ax.set_ylabel('Angle in degrees')
        ax.set_xlabel('Particle number')
        try:
            ax.legend(plots, ['Average for %s angle%s' % (k, 's' if k else '')
                              for k in steps], fontsize='small',
                      bbox_to_anchor=(1, 0.5), loc='center left')
        except TypeError:
            ax.legend(plots, ['Average for %s angle%s' % (k, 's' if k else '')
                              for k in steps],
                      bbox_to_anchor=(1, 0.5), loc='center left')
        ax.set_xlim((1, self.nloci))
        ax.set_title('Angle between consecutive loci')
        plt.subplots_adjust(left=0.1, right=0.8)

        if savefig:
            tadbit_savefig(savefig)
        elif not axe:
            plt.show()


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
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
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
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        models = [self[mdl] for mdl in models]
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
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        for model_num in models:
            try:
                model = self[model_num]
            except KeyError:
                model = self._bad_models[model_num]
            model.write_cmm(directory, color=color, rndname=rndname,
                            model_num=model_num, **kwargs)


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
            models = [m if type(m) is int else self[m]['index'] for m in models]
        elif cluster > -1:
            models = [self[str(m)]['index'] for m in self.clusters[cluster]]
        else:
            models = self.__models
        for model_num in models:
            try:
                model = self[model_num]
            except KeyError:
                model = self._bad_models[model_num]
            path_f = model.write_xyz(directory, model_num=model_num,
                                     get_path=get_path, rndname=rndname)
        if get_path:
            return path_f


    def save_models(self, outfile):
        """
        Saves all the models in pickle format (python object written to disk).

        :param path_f: path where to save the pickle file
        """
        to_save = {}

        to_save['models']        = self.__models
        to_save['bad_models']    = self._bad_models
        to_save['nloci']         = self.nloci
        to_save['clusters']      = self.clusters
        to_save['resolution']    = self.resolution
        to_save['original_data'] = self._original_data
        to_save['config']        = self._config
        to_save['zscore']        = self._zscores

        out = open(outfile, 'w')
        dump(to_save, out)
        out.close()


class ClusterOfModels(dict):
    def __str__(self):
        out1 = '   Cluster #%s has %s models [top model: %s]\n'
        out = 'Total number of clusters: %s\n%s' % (
            len(self),
            ''.join([out1 % (k, len(self[k]), self[k][0]) for k in self]))
        return out
