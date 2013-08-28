"""
19 Jul 2013


"""
from pytadbit.utils.tadmaths   import calinski_harabasz, calc_eqv_rmsd
from pytadbit.utils.tadmaths   import calc_consistency
from pytadbit.utils.extraviews import color_residues, chimera_view
from pytadbit.utils.extraviews import augmented_dendrogram, plot_hist_box
from cPickle                   import load, dump
from subprocess                import Popen, PIPE
from math                      import sqrt, acos, degrees
from numpy                     import median as np_median
from numpy                     import std as np_std, log2
from scipy.cluster.hierarchy   import linkage, fcluster
from scipy.stats               import spearmanr
from warnings                  import warn


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
    To generate three-dimensional models from Hi-C data (z-scores of the
    interactions). A given number of models is kept, and can be used to draw
    some statistics or search for some specific interactions.

    :param nloci: length of the chromatin fragment modelled (in number of
       particles)
    :param models: a dictionary contatining 
    :param resolution: of the Hi-C experiment, this will be the number of
       nucleotides in each particle of the models

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
                '  IMP modelling used this parameters:\n' +
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
        Models are stored according to their objctive function value (first
        best), but in order to reproduce a model, we need its initial random
        number. This method helps to fetch the model corresponding to a given
        intitial random number stored under
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
                       mcl_bin='mcl', tmp_file='/tmp/tmp.xyz', verbose=True):
        """
        Runs a clustering analysis over models. Clusters found will be stored
           in StructuralModels.clusters
        
        :param 0.75 fact: Factor for equivalent positions
        :param 200 dcutoff: distance in nanometer from which it is considered
           that two particles are separated.
        :param 'score' var: value to return, can be either (i) 'drmsd' (symmetry
           independent: mirrors will show no differences) (ii) 'score' that is:
           
           ::
           
                                    dRMSD[i] / max(dRMSD)
              score[i] = eqvs[i] * -----------------------
                                     RMSD[i] / max(RMSD)
           
           where eqvs[i] is the number of equivalent position for the ith
           pairwise model comparison.
        :param 'mcl' method: clustering method to use, can be either 'mcl' or
           'ward'. Last one uses scipy implementation, and is NOT RECOMENDED.
        :param 'mcl' mcl_bin: path to MCL executable file, in case 'mcl is not
           in the PATH'
        :param '/tmp/tmp.xyz' tmp_file: path to a temporary file
        :param True verbose: same as print StructuralModels.clusters
        
        """
        scores = calc_eqv_rmsd(self.__models, self.nloci, dcutoff, var)
        from distutils.spawn import find_executable
        if not find_executable(mcl_bin):
            print('\nWARNING: MCL not found in path using WARD clustering\n')
            method = 'ward'
            
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
            for (md1, md2), score in scores.iteritems():
                out_f.write('model_%s\tmodel_%s\t%s\n' % (
                    md1, md2, score if score > fact * self.nloci else 0.0))
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
                matrix[i][j+i+1] = calc_eqv_rmsd({1: md1, 2: md2}, self.nloci,
                                                 var='drmsd', one=True)
        return clust_count, objfun, matrix


    def cluster_analysis_dendrogram(self, n_best_clusters=None, color=False,
                                    axe=None, savefig=None):
        """
        Representation of clusters of models. The length of the leaves if
           proportional to the final objective function value of each model.
           Branch widths are proportional to the number of models in a given
           cluster (or group of clusters, ifit is an internal branch)
           
        :param None n_best_clusters: number of clusters to represent, by default
           all clusters will be shown
        :param False color: display colors for the significance of the
           clustering (basically dependent of internal branch lengths)
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
                     error=False, axe=None, savefig=None, outfile=None):
        """
        Represents the number of nucletotide base pairs can be found in 1 nm of
           chromatin, this along strand modelled.

        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param (1, 2, 3, 4, 5) steps: how many particles to group for the
           estimation. By default 5 curves are drawn.
        :param False error: represent the error of the estimates.
        :param None outfile: path to a file where to save the density data
           generated (1 column per step + 1 for particle number).
        
        """
        if type(steps) == int:
            steps = (steps, )
        if len(steps) > 6:
            raise Exception('Sorry not enough colors to do this.\n')
        colors = ['grey', 'darkgreen', 'darkblue', 'purple', 'darkorange', 'darkred'][-len(steps):]
        dists = []
        for part1, part2 in zip(range(self.nloci - 1), range(1, self.nloci)):
            dists.append(self.median_3d_dist(part1, part2, models, cluster,
                                              plot=False, median=False))
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
        if outfile:
            out = open(outfile, 'w')
            out.write('#Particle\t%s\n' % ('\t'.join([str(c) for c in steps])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (part + 1, '\t'.join(
                    [('-' if not part in distsk[c] else str(distsk[c][part])
                      if distsk[c][part] else 'None')
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
            plots += ax.plot(distsk[k], color=colors[steps.index(k)],
                             lw=steps.index(k) + 1, alpha=0.5)
        if error:
            for k in steps:
                plots += ax.plot(errorp[k], color=colors[steps.index(k)], ls='--')
                ax.plot(errorn[k], color=colors[steps.index(k)], ls='--')
        ax.set_ylabel('Density (bp / nm)')
        ax.set_xlabel('Particle number')
        ax.legend(plots, ['Average for %s particle%s' % (k, 's' if k else '')
                          for k in steps] + (
                      ['+/- 2 standard deviations for %s' % (k)
                       for k in steps] if error else []), fontsize='small',
                  bbox_to_anchor=(1, 0.5), loc='center left')
        ax.set_xlim((0, self.nloci))
        ax.set_title('Chromatin density')
        if savefig:
            fig.savefig(savefig)
        elif not axe:
            plt.show()


    def get_contact_matrix(self, models=None, cluster=None, cutoff=150):
        """
        Draws a heatmap representing the proportion of times two particles are
           closer than a given cutoff.

        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param 150 cutoff: distance in nanometer from which it is considered
           that two particles are separated.

        :returns: matrix of contact counts
        """
        if models:
            models = models
        elif cluster > -1:
            models = [str(m) for m in self.clusters[cluster]]
        else:
            models = self.__models
        matrix = [[float('nan') for _ in xrange(self.nloci)] for _ in xrange(self.nloci)]
        for i in xrange(self.nloci):
            for j in xrange(i + 1, self.nloci):
                val = len([k for k in self.median_3d_dist(
                    i, j, plot=False, median=False, models=models)
                           if k < cutoff])
                matrix[i][j] = matrix[j][i] = float(val) / len(models) * 100
        return matrix
        

    def define_best_models(self, nbest):
        """
        Define the number of best models to keep. If keep_all was True in
        :func:`pytadbit.imp.imp_model.generate_3d_models` or in
        :func:`pytadbit.experiment.Experiment.model_region` the full set of models
        (n_models parameter) will be available, otherwise only the n_keep models.

        :param nbest: number of models to consider as best models.
           Usually 20% of the models generated are kept.
        """
        tmp_models = self.__models
        tmp_models.update(self._bad_models)
        self.__models = dict([(i, tmp_models[i]) for i in xrange(nbest)])
        self._bad_models = dict([(i, tmp_models[i]) for i in
                                 xrange(nbest, len(tmp_models))])


    def contact_map(self, models=None, cluster=None, cutoff=150, axe=None,
                    savefig=None, outfile=None):
        """
        Draws a heatmap representing the proportion of times two particles are
           closer than a given cutoff.

        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param 150 cutoff: distance in nanometer from which it is considered
           that two particles are separated.
        :param None axe: a matplotlib.axes.Axes object, using it allows to
           redefine figure size, background etc...
        :param None savefig: path to a file where to save the image generated
           if None, the image will be shown using matplotlib GUI.
        :param None outfile: path to a file where to save the contact map data
           generated (1 column particle1, 1 for particle2, one with the
           percentage of models where this particles are close in space
           
        """
        matrix = self.get_contact_matrix(models, cluster, cutoff)
        show = False
        if outfile:
            out = open(outfile, 'w')
            out.write('#Particle1\tParticle2\tModels_percentage\n')
            for i in xrange(len(matrix)):
                for j in xrange(i+1, len(matrix)):
                    out.write('%s\t%s\t%s\n' % (i, j, matrix[i][j]))
            out.close()
        if not axe:
            fig = plt.figure(figsize=(8, 6))
            axe = fig.add_subplot(111)
            show=True
        else:
            fig = axe.get_figure()
        ims = axe.imshow(matrix, origin='lower', interpolation="nearest")
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
        :param hic_matrix: a matrix representing the normalized count of Hi-C
           interactions, used to generated these models.
        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param 150 cutoff: distance in nanometer from which it is considered
           that two particles are separated.
        :returns: spearman correlation rho and p-value, between the two
           matrices. Good correlation may have a Rho value upper than 0.7
        :param None savefig: path to a file where to save the image generated
           if None, the image will be shown using matplotlib GUI.
        
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
        fig.suptitle('Correlation bettween normalized-real and modelled '
                     + 'contact maps (correlation=%.4f)' % (corr[0]),
                     size='x-large')
        ax = fig.add_subplot(121)
        self.contact_map(models, cluster, cutoff, axe=ax)
        ax = fig.add_subplot(122)
        ims = ax.imshow(log2(self._original_data), origin='lower',
                        interpolation="nearest")
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
                          cluster=None, axe=None, savefig=None, outfile=None):
        """
        Plots the consistency of a given set of models, along the chromatin
           fragment modelled. This plot can also be viewed as how well defined,
           or how stable, is a given portion of the chromatin model.

        :param (50,100,150,200) cutoffs: list of cutoff values (in nanometer)
           to plot. These distances are used to know when to consider two
           particles as diferent.
        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param '/tmp/tmp_cons' tmp_path: where to write input files for TM-score
           program
        :param '' tmsc: path to TMscore_consistency, by default it assumes that
           it is installed
        :param None outfile: path to a file where to save the consistency data
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
        if outfile:
            out = open(outfile, 'w')
            out.write('#Particle\t%s' % ('\t'.join([c for c in cutoffs])))
            for part in xrange(self.nloci):
                out.write('%s\t%s\n' % (part + 1, '\t'.join(
                    [consistencies[c][part] for c in cutoffs])))
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
            plots += axe.plot(consistencies[cut], color='darkred',
                              alpha= 1 - i / float(len(cutoffs)))
        axe.legend(plots, ['%s nm' % (k) for k in cutoffs[::-1]],
                   fontsize='small', loc='center left', bbox_to_anchor=(1, 0.5))
        axe.set_xlim((0, self.nloci))
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
        Visualize a given model in three dimensions

        :param model_num: number of the model to visualize
        :param 'chimera' tool: path to external tool to visualize the model
        :param None savefig: path to a file where to save the image OR movie
           generated (depending on the extension, accepted formats are png,
           mov and webm) if None, the image or movie will be shown using
           default GUI.
        :param None cmd: a list of commands to be passed to the viewer. The
           default is (using chimera tool):

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

           Followed by the movie record for movies:

           ::
          
             movie record supersample 1
             turn y 3 120
             wait 120
             movie stop
             movie encode output SAVEFIG

           Or the copy for images:

           ::

             copy file SAVEFIG png

           Thus If you pass as 'cmd' parameter this list:
           ::

             cmd = ['focus', 'set bg_color white', 'windowsize 800 600', 'bonddisplay never #0', 'shape tube #0 radius 10 bandLength 200 segmentSubdivisions 100 followBonds on', 'clip yon -500', '~label', 'set subdivision 1', 'set depth_cue', 'set dc_color black', 'set dc_start 0.5', 'set dc_end 1', 'scale 0.8']
             
           you would obtain the same result as with default (do not forget to
           add commands to save figure/movie if you which).
        
        """
        self.write_cmm('/tmp/', model_num=model_num)
        chimera_view('/tmp/model.%s.cmm' % (self[model_num]['rand_init']),
                     savefig=savefig, chimera_bin=tool, chimera_cmd=cmd)
    

    def measure_angle_3_particles(self, parta, partc, partb,
                                  models=None, cluster=None,
                                  radian=False):
        """
        Given three particles (A, C and B) the angle g shown below:
        
        ::


                              A
                             /|
                            / |
                          b/  |
                          /   |
                         /    |
                        C )g  |c
                         \    |
                          \   |
                          a\  |
                            \ |
                             \|
                              B

        is given by the theorem of Al-Kashi:
        
        .. math::

          c^2 = a^2 + b^2 - 2 \cos(g)

        :param part1: A particle number
        :param part2: A particle number
        :param part3: A particle number
        :param None models: If None (default) will do calculate the distance
            along all models. A list of numbers corresponding to a given set of
            models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param False radian: return value in radian, in degree if false

        :returns: an angle,  either in degrees or radians
        """

        a = self.median_3d_dist(partc, partb, models=models,
                                cluster=cluster, plot=False)
        b = self.median_3d_dist(parta, partc, models=models,
                                cluster=cluster, plot=False)
        c = self.median_3d_dist(parta, partb, models=models,
                                cluster=cluster, plot=False)

        g = acos((a**2 + b**2 - c**2) / (2 * a * b))

        return g if radian else degrees(g)


    def median_3d_dist(self, part1, part2, models=None, cluster=None,
                       plot=True, median=True, axe=None, savefig=None):
        """
        Computes the distance between two particles. This is done by calculating
           the median value corresponding to the set of given models.
        
        :param part1: number corresponding to the first particle
        :param part2: number corresponding to the second particle
        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param True plot: if True, will display a histogram and a boxplot
           representing the distribution of distances calculated. Else it will
           return either the median value of these distances or the full list
        :param True median: if plot is False this option set the return to,
           either the median value of these distances or the full list
        
        :returns: if plot is false, returns either the median values of the
           distances, either the list of distances.
        """
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
        plot_hist_box(dists, part1, part2, axe, savefig)


    def objective_function_model(self, model, log=False, smooth=True, axe=None,
                                 savefig=None):
        """
        Plots the fall in objective function through the Monte-Carlo search.

        :param model: the number of the model to plot
        :param False log: to plot in log scale
        :param True smooth: to smooth the curve
        """
        self[model].objective_function(log=log, smooth=smooth, axe=axe,
                                       savefig=savefig)
        

    def write_cmm(self, directory, model_num=None, models=None, cluster=None,
                  color=color_residues, rndname=True):
        """
        Writes cmm file read by Chimera (http://www.cgl.ucsf.edu/chimera).

        **Note:** If none of model_num, models or cluster parameter are set,
        ALL models will be writen.
        
        :param directory: where to write the file (note: the name of the file
           will be model_1.cmm if model number is 1)
        :param None model_num: the number of the model to write.
        :param None models: A list of numbers corresponding to a given set of
           models to be written.
        :param None cluster: A number can be passed in order to write models
           corresponding to the cluster number 'cluster'
        :param True rndname: If True, file names will be formated as:
           model.RND.cmm, where RND is the random initial value used by IMP to
           generate this model. If False, the format will be:
           model_NUM_RND.cmm where NUM is the rank of the model in terms of
           objective function value.
        :param color_residues color: either a coloring function like
           :func:`pytadbit.imp.imp_model.color_residues` or a list of (r, g, b)
           tuples (as long as the number of particles). 
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
                model = self.__models[model_num]
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
        Writes xyz file containing the 3D coordinates of each particles.

        **Note:** If none of model_num, models or cluster parameter are set,
        ALL models will be writen.
        
        :param directory: where to write the file (note: the name of the file
           will be model.1.xyz if model number is 1)
        :param None model_num: the number of the model to write.
        :param None models: A list of numbers corresponding to a given set of
           models to be written.
        :param None cluster: A number can be passed in order to write models
           corresponding to the cluster number 'cluster'
        :param True rndname: If True, file names will be formated as:
           model.RND.xyz, where RND is the random initial value used by IMP to
           generate this model. If False, the format will be:
           model_NUM_RND.xyz where NUM is the rank of the model in terms of
           objective function value.
        :param False get_path: whether to return, or not the full path where
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
                out += form % ('p' + str(n + 1), n + 1, model['x'][n],
                               model['y'][n], model['z'][n])
            out_f = open(path_f, 'w')
            out_f.write(out)
            out_f.close()
            if get_path:
                return path_f


    def save_models(self, path_f):
        """
        Saves all models in pickle format (python object written to disk).
        
        :param path_f: path to file where to pickle
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
    def __repr__(self):
        out1 = '   Cluster #%s has %s models [top model: %s]\n'
        out = 'Total number of clusters: %s\n%s' % (
            len(self), 
            ''.join([out1 % (k, len(self[k]), self[k][0]) for k in self]))
        return out
