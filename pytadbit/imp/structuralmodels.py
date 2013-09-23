"""
19 Jul 2013


"""
from pytadbit.utils.tadmaths   import calinski_harabasz, calc_eqv_rmsd
from pytadbit.utils.tadmaths   import calc_consistency, dihedral
from pytadbit.utils.extraviews import color_residues, chimera_view
from pytadbit.utils.extraviews import augmented_dendrogram, plot_hist_box
from cPickle                   import load, dump
from subprocess                import Popen, PIPE
from math                      import sqrt, acos, degrees, pi
from numpy                     import median as np_median
from numpy                     import std as np_std, log2
from numpy                     import array, cross, dot
from numpy.linalg              import norm
from scipy.cluster.hierarchy   import linkage, fcluster
from scipy.stats               import spearmanr
from warnings                  import warn
from string                    import uppercase as uc, lowercase as lc
from random                    import random

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')


def load_structuralmodels(path_f):
    """
    
    :param path: to the pickled StructuralModels object.

    :returns: a :class:`pytadbit.imp.imp_model.StructuralModels`.
    """
    svd = load(open(path_f))
    return StructuralModels(
        nloci=svd['nloci'], models=svd['models'], bad_models=svd['bad_models'],
        resolution=svd['resolution'], original_data=svd['original_data'],
        clusters=svd['clusters'], config=svd['config'])


class StructuralModels(object):
    """
    This function generates three-dimensional models starting from Hi-C data. 
    The final analysis will be performed on the n_keep top models.

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
                 original_data=None, clusters=None, config=None):
        
        self.__models       = models
        self._bad_models    = bad_models
        self.nloci          = nloci
        self.clusters       = clusters or ClusterOfModels()
        self.resolution     = float(resolution)
        self._original_data = original_data # only used for correlation
        self._config        = config
        

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


    def cluster_models(self, fact=0.75, dcutoff=200, var='score', method='mcl',
                       mcl_bin='mcl', tmp_file=None, verbose=True):
        """
        This function performs a clustering analysis of the generated models 
        based on structural comparison. The result will be stored in 
        StructuralModels.clusters
        
        :param 0.75 fact: factor to define the percentage of equivalent 
           positions to be considered in the clustering
        :param 200 dcutoff: distance threshold (nm) to determine if two 
           particles are in contact
        :param 'score' var: value to return, which can be either (i) 'drmsd' 
           (symmetry independent: mirrors will show no differences), or (ii) 
           'score', defined as:
           
           ::
           
                                    dRMSD[i] / max(dRMSD)
              score[i] = eqvs[i] * -----------------------
                                     RMSD[i] / max(RMSD)
           
           where eqvs[i] is the number of equivalent position for the ith 
           pairwise model comparison
        :param 'mcl' method: clustering method to use, which can be either 
           'mcl' or 'ward'. The last one uses scipy implementation, and is 
           NOT RECOMMENDED.
        :param 'mcl' mcl_bin: path to the mcl executable file, in case of the 
           'mcl is not in the PATH' warning message
        :param None tmp_file: path to a temporary file created during
           the clustering computation. Default will be created in /tmp/ folder
        :param True verbose: same as print StructuralModels.clusters
        
        """
        tmp_file = '/tmp/tadbit_tmp_%s.txt' % (
            ''.join([(uc + lc)[int(random() * 52)] for _ in xrange(4)]))
        scores = calc_eqv_rmsd(self.__models, self.nloci, dcutoff, var)
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
            solutions = {}
            # x = []
            # y = []
            for k in clust[:,2]:
                clusters = ClusterOfModels()
                _ = [clusters.setdefault(j, []).append(i) for i, j in
                     enumerate(fcluster(clust, k, criterion='distance'))]
                solutions[k] = {'out': clusters}
                solutions[k]['score'] = calinski_harabasz(scores, clusters)
                # print len(a), k, solutions[k]['score']
                # x.append(solutions[k]['score'])
                # y.append(len(a))
            # plt.plot(y,x)
            # plt.show()
            clusters = [solutions[s] for s in sorted(
                solutions, key=lambda x: solutions[x]['score'])
                         if solutions[s]['score']>0][-1]['out']
            for cluster in clusters:
                for model in clusters[cluster]:
                    self[model]['cluster'] = cluster
            self.clusters = clusters
        else:
            out_f = open(tmp_file, 'w')
            uniqs = list(set([tuple(sorted((m1, m2))) for m1, m2 in scores]))
            for md1, md2 in uniqs:
                score = scores[(md1, md2)]
                if score >= fact * self.nloci:
                    out_f.write('model_%s\tmodel_%s\t%s\n' % (md1, md2, score))
            out_f.close()
            Popen('%s %s --abc -V all -o %s.mcl' % (
                mcl_bin, tmp_file, tmp_file), stdout=PIPE, stderr=PIPE,
                  shell=True).communicate()
            self.clusters = ClusterOfModels()
            for cluster, line in enumerate(open(tmp_file + '.mcl')):
                self.clusters[cluster] = []
                for model in line.split():
                    model = int(model.split('_')[1])
                    self[model]['cluster'] = cluster
                    self.clusters[cluster].append(self[model]['rand_init'])
                self.clusters[cluster].sort(
                    key=lambda x: self[str(x)]['objfun'])
            # sort clusters according to their lowest energy
            # for clt in clusters:
            #     clusters[clt].sort()
            # for i, clt in enumerate(sorted(
            #     clusters, key=lambda x: self[clusters[x][0]]['energy'])):
            #     self.clusters[i] = clusters[clt]
            #     for model in self.clusters[i]:
            #         self.__models[model]['cluster'] = i
        if verbose:
            print 'Number of Singletons excluded from clustering: %s' % (
                len([1 for model in self if model['cluster'] == 'Singleton']))
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
                                                 var='drmsd', one=True)
        return clust_count, objfun, matrix


    def cluster_analysis_dendrogram(self, n_best_clusters=None, color=False,
                                    axe=None, savefig=None):
        """
        Representation of the clustering results. The length of the leaves if 
        proportional to the final objective function value of each model. The 
        branch widths are proportional to the number of models in a given 
        cluster (or group of clusters, if it is an internal branch).
           
        :param None n_best_clusters: number of clusters to represent (by 
           default all clusters will be shown)
        :param False color: color the dendrogram based on the significance of
           the clustering (basically it depends of the internal branch lengths)
        """

        if not self.clusters:
            self.cluster_models()
        if not n_best_clusters:
            n_best_clusters = len(self.clusters)
        if n_best_clusters <= 1:
            raise Exception("Need at least 2 clusters to display...")
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
            clust_count[i] = clust_count[a] + clust_count[b]
            dads[a] = i
            dads[b] = i

        d = augmented_dendrogram(clust_count, dads, objfun, color, 
                                 axe, savefig, z)
        return d


    def density_plot(self, models=None, cluster=None, steps=(1, 2, 3, 4, 5),
                     error=False, axe=None, savefig=None, savedata=None):
        """
        Plots the number of nucleotides per nm of chromatin vs the modeled
        region bins.

        :param None models: if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param (1, 2, 3, 4, 5) steps: how many particles to group for the
           estimation. By default 5 curves are drawn
        :param False error: represent the error of the estimates
        :param None savedata: path to a file where to save the density data
           generated (1 column per step + 1 for particle number).
        
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
                     str(errorp[c][part]))
                     if distsk[c][part] else 'None\tNone'
                     for c in steps])))
            out.close()
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
        ax.legend(plots, ['Average for %s particle%s' % (k, 's' if k else '')
                          for k in steps] + (
                      ['+/- 2 standard deviations'
                       for k in steps] if error else []), fontsize='small',
                  bbox_to_anchor=(1, 0.5), loc='center left')
        ax.set_xlim((1, self.nloci))
        ax.set_title('Chromatin density')
        plt.subplots_adjust(left=0.1, right=0.78)
        if savefig:
            fig.savefig(savefig)
        elif not axe:
            plt.show()


    def get_contact_matrix(self, models=None, cluster=None, cutoff=150):
        """
        Returns a matrix with the number of interactions observed below a given 
        cutoff distance.

        :param None models: if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param 150 cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not

        :returns: matrix frequency of interaction
        """
        if models:
            models = models
        elif cluster > -1:
            models = [str(m) for m in self.clusters[cluster]]
        else:
            models = self.__models
        matrix = [[float('nan') for _ in xrange(self.nloci)]
                  for _ in xrange(self.nloci)]
        for i in xrange(self.nloci):
            for j in xrange(i + 1, self.nloci):
                val = len([k for k in self.median_3d_dist(
                    i + 1, j + 1, plot=False, median=False, models=models)
                           if k < cutoff])
                matrix[i][j] = matrix[j][i] = float(val) / len(models) * 100
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
           if None, the image will be shown using matplotlib GUI
        :param None savedata: path to a file where to save the contact map data
           generated, in three columns format (particle1, particle2, percentage
           of models where these two particles are in contact)
           
        """
        matrix = self.get_contact_matrix(models, cluster, cutoff)
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
                         extent=(0.5, self.nloci + 0.5, 0.5, self.nloci + 0.5))
        axe.set_ylabel('Particle')
        axe.set_xlabel('Particle')
        cbar = axe.figure.colorbar(ims)
        cbar.ax.set_yticklabels(['%3s%%' % (p) for p in range(0, 110, 10)])
        cbar.ax.set_ylabel('Percentage of models with particles closer than ' +
                           '%s nm' % (cutoff))
        axe.set_title('Contact map')
        if savefig:
            fig.savefig(savefig)
        elif show:
            plt.show()


    def correlate_with_real_data(self, models=None, cluster=None, cutoff=200,
                                 plot=False, axe=None, savefig=None):
        """
        :param None models: if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param 150 cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not
        :param None savefig: path to a file where to save the image generated; 
           if None, the image will be shown using matplotlib GUI
        :param False plot: to display the plot
        :param None axe: a matplotlib.axes.Axes object to define the plot 
           appearance

        :returns: Spearman correlation rho and p-value, between the two
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
        if not plot:
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
        ax.set_title('Z-scores of the observed Hi-C count')
        cbar = ax.figure.colorbar(ims)
        cbar.ax.set_ylabel('Log2 (normalized Hi-C data)')
        
        if savefig:
            fig.savefig(savefig)
        elif not axe:
            plt.show()


    def model_consistency(self, cutoffs=(50, 100, 150, 200), models=None,
                          cluster=None, axe=None, savefig=None, savedata=None):
        """
        Plots the particle consistency, over a given set of models, vs the 
        modeled region bins. The consistency is a measure of the variability
        (or stability) of the modeled region (the higher the consistency value,
        the higher stability).  

        :param (50,100,150,200) cutoffs: list of distance cutoffs (nm) used to
           compute the consistency. Two particle are considered consistent if 
           their distance is less than the given cutoff
        :param None models:  if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param '/tmp/tmp_cons' tmp_path: location of the input files for 
           TM-score program
        :param '' tmsc: path to the TMscore_consistency script (assumed to be 
           installed by default)
        :param None savedata: path to a file where to save the consistency data
           generated (1 column per cutoff + 1 for particle number).
        
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
        axe.legend(plots, ['%s nm' % (k) for k in cutoffs[::-1]],
                   fontsize='small', loc='center left', bbox_to_anchor=(1, 0.5))
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
            fig.savefig(savefig)
        elif show:
            plt.show()


    def view_model(self, model_num, tool='chimera', savefig=None, cmd=None):
        """
        Visualize a selected model in the three dimensions.

        :param model_num: model to visualize
        :param 'chimera' tool: path to the external tool used to visualize the 
           model
        :param None savefig: path to a file where to save the image OR movie
           generated (depending on the extension; accepted formats are png, mov
           and webm). if set to None, the image or movie will be shown using 
           the default GUI.
        :param None cmd: list of commands to be passed to the viewer. The chimera list is:

           ::

             focus
             set bg_color white
             windowsize 800 600
             bonddisplay never #0
             shape tube #0 radius 10 bandLength 200 segmentSubdivisions 100 followBonds on
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
        
        """
        self.write_cmm('/tmp/', model_num=model_num)
        chimera_view('/tmp/model.%s.cmm' % (self[model_num]['rand_init']),
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd)
    

    def measure_angle_3_particles(self, parta, partb, partc,
                                  models=None, cluster=None,
                                  radian=False, all_angles=False):
        """
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
        :param None models:  if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
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
        """
        if models:
            models=models
        elif cluster > -1:
            models = [str(m) for m in self.clusters[cluster]]
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
        """
        parta = array(self.particle_coordinates(parta, models, cluster))
        partb = array(self.particle_coordinates(partb, models, cluster))
        partc = array(self.particle_coordinates(partc, models, cluster))
        partd = array(self.particle_coordinates(partd, models, cluster))
        return dihedral(parta, partb, partc, partd)


    def walking_dihedral(self, models=None, cluster=None, steps=(1,3),
                         plot=True, savefig=None, axe=None):
        """
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
            fig.savefig(savefig)
        elif not axe:
            plt.show()



    def walking_angle(self, models=None, cluster=None, steps=(1,3), signed=True,
                      plot=True, savefig=None, savedata=None, axe=None):
        """
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
            rads[1].append(self.measure_angle_3_particles(res + 1, res + 4,
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
        ax.legend(plots, ['Average for %s angle%s' % (k, 's' if k else '')
                          for k in steps], fontsize='small',
                  bbox_to_anchor=(1, 0.5), loc='center left')
        ax.set_xlim((1, self.nloci))
        ax.set_title('Angle between consecutive loci')
        plt.subplots_adjust(left=0.1, right=0.8)

        if savefig:
            fig.savefig(savefig)
        elif not axe:
            plt.show()


    def median_3d_dist(self, part1, part2, models=None, cluster=None,
                       plot=True, median=True, axe=None, savefig=None):
        """
        Computes the median distance between two particles over a set of models
        
        :param part1: number corresponding to the first particle
        :param part2: number corresponding to the second particle
        :param None models:  if None (default) the contact map will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the contact map only for the models in the
           cluster number 'cluster'
        :param True plot: if True, display a histogram and a box-plot of the 
           distribution of the calculated distances. If False, return either
           the full list of the calculated distances or their median value
        :param True median: return either the full list of the calculated 
           distances (False) or their median value (True), when 'plot' is set 
           to False
        
        :returns: if 'plot' is False, return either the full list of the 
           calculated distances or their median value distances, either the 
           list of distances.
        """
        part1 -= 1
        part2 -= 1
        dists = []
        if models:
            models=models
        elif cluster > -1:
            models = [str(m) for m in self.clusters[cluster]]
        else:
            models = self.__models
        for mdl in models:
            dists.append(
                sqrt(
                    (self[mdl]['x'][part1] - self[mdl]['x'][part2])**2 + 
                    (self[mdl]['y'][part1] - self[mdl]['y'][part2])**2 +
                    (self[mdl]['z'][part1] - self[mdl]['z'][part2])**2)
                )
        if not plot:
            if median:
                return np_median(dists)
            else:
                return dists
        plot_hist_box(dists, part1 + 1, part2 + 1, axe, savefig)


    def objective_function_model(self, model, log=False, smooth=True, axe=None,
                                 savefig=None):
        """
        This function plots the objective function value per each Monte-Carlo 
        step

        :param model: the number of the model to plot
        :param False log: log plot
        :param True smooth: curve smoothing
        """
        self[model].objective_function(log=log, smooth=smooth, axe=axe,
                                       savefig=savefig)
        

    def write_cmm(self, directory, model_num=None, models=None, cluster=None,
                  color=color_residues, rndname=True):
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
        :param color_residues color: either a coloring function like
           :func:`pytadbit.imp.imp_model.color_residues` or a list of (r, g, b)
           tuples (as long as the number of particles) 
        """
        if model_num > -1:
            models = [model_num]
        elif models:
            models = models
        elif cluster > -1:
            models = [str(m) for m in self.clusters[cluster]]
        else:
            models = self.__models        
        for model_num in models:
            try:
                model = self[model_num]
            except KeyError:
                model = self._bad_models[model_num]
            if type(color) != list:
                color = color(self.nloci)
            out = '<marker_set name=\"%s\">\n' % (model['rand_init'])
            form = ('<marker id=\"%s\" x=\"%s\" y=\"%s\" z=\"%s\"' +
                    ' r=\"%s\" g=\"%s\" b=\"%s\" ' +
                    'radius=\"' +
                    str(self.resolution * self._config['scale']) +
                    '\" note=\"%s\"/>\n')
            for n in xrange(self.nloci):
                out += form % (n + 1,
                               model['x'][n], model['y'][n], model['z'][n],
                               color[n][0], color[n][1], color[n][2], n + 1)
            form = ('<link id1=\"%s\" id2=\"%s\" r=\"1\" ' +
                    'g=\"1\" b=\"1\" radius=\"' +
                    str(1) +
                    '\"/>\n')
            for n in xrange(1, self.nloci):
                out += form % (n, n + 1)
            out += '</marker_set>\n'

            if rndname:
                out_f = open('%s/model.%s.cmm' % (directory,
                                                  model['rand_init']), 'w')
            else:
                out_f = open('%s/model_%s_rnd%s.cmm' % (
                    directory, model_num, model['rand_init']), 'w')
            out_f.write(out)
            out_f.close()


    def write_xyz(self, directory, model_num=None, models=None, cluster=None,
                  get_path=False, rndname=True):
        """
        Writes a xyz file containing the 3D coordinates of each particle in the
        model.

        **Note:** If none of model_num, models or cluster parameter are set,
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
            models = models
        elif cluster > -1:
            models = [str(m) for m in self.clusters[cluster]]
        else:
            models = self.__models        
        for model_num in models:
            try:
                model = self[model_num]
            except KeyError:
                model = self._bad_models[model_num]
            if rndname:
                path_f = '%s/model.%s.xyz' % (directory, model['rand_init'])
            else:
                path_f = '%s/model_%s_rnd%s.xyz' % (directory, model_num,
                                                    model['rand_init'])
            out = ''
            form = "%12s%12s%12.3f%12.3f%12.3f\n"
            for n in xrange(self.nloci):
                out += form % ('p' + str(n + 1), n + 1, round(model['x'][n], 3),
                               round(model['y'][n], 3), round(model['z'][n], 3))
            out_f = open(path_f, 'w')
            out_f.write(out)
            out_f.close()
            if get_path:
                return path_f


    def save_models(self, path_f):
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
        
        out = open(path_f, 'w')
        dump(to_save, out)
        out.close()


class ClusterOfModels(dict):
    def __str__(self):
        out1 = '   Cluster #%s has %s models [top model: %s]\n'
        out = 'Total number of clusters: %s\n%s' % (
            len(self), 
            ''.join([out1 % (k, len(self[k]), self[k][0]) for k in self]))
        return out
