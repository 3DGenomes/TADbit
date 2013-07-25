"""
19 Jul 2013


"""
from pytadbit.utils import color_residues, calc_eqv_rmsd
from pytadbit.utils import augmented_dendrogram, plot_hist_box
from cPickle import load, dump
from subprocess import Popen, PIPE
from math import sqrt
from numpy import median as np_median
from numpy import std as np_std


def load_threedeemodels(path):
    """
    
    :param path: to the pickled ThreeDeeModels object.

    :returns: a :class:`pytadbit.imp.imp_model.ThreeDeeModels`.
    """
    return load(open(path))


class ThreeDeeModels(object):
    """
    To generate three-dimensional models from Hi-C data (z-scores of the
    interactions). A given number of models is kept, and can be used to draw
    some statistics or search for some specific interactions.

    :param nloci: length of the chromatine fragment modelled (in number of
       particles)
    :param models: a dictionary contatining 
    :param resolution: of the Hi-C experiment, this will be the number of
       nucleotides in each particle of the models

    """

    def __init__(self, nloci, models, bad_models, resolution):
        
        self.__models = models
        self._bad_models = bad_models
        self.nloci = nloci
        self.clusters = {}
        self.resolution = float(resolution)


    def __getitem__(self, nam):
        try:
            return self.__models[nam]
        except TypeError:
            for i, key in self.__models:
                if nam == i:
                    return self.__models[key]
            raise KeyError('Model {} not found\n'.format(i))


    def __repr__(self):
        return ('ThreeDeeModels with {} models (energy range: {}-{})\n' +
                '   (corresponding to the best models out of {} models).\n' +
                '  Models where clustered into {} clusters').format(
            len(self.__models),
            self.__models[0]['energy'],
            self.__models[len(self.__models) - 1]['energy'],
            len(self.__models) + len(self._bad_models), len(self.clusters))


    def fetch_model_by_rand_init(self, rand_init, all_models=False):
        """
        Models are stored according to their energy value (first best), but in
        order to reproduce a model, we need its initial random number. This
        method helps to fetch the model corresponding to a given intitial random
        number stored under ThreeDeeModels.models[N]['rand_init'].

        :param rand_init: the wanted rand_init number.
        :param False all_models: whether or not to use 'bad' models

        :returns: index of 3d model
        """
        for m in self.__models:
            if self.__models[m]['rand_init'] == rand_init:
                return m
        if all_models:
            for m in self._bad_models:
                if self._bad_models[m]['rand_init'] == rand_init:
                    return m
        raise IndexError('Model {} not found\n'.format(rand_init))


    def cluster_models(self, fact=0.75, dcutoff=200,
                       mcl_bin='mcl', tmp_file='/tmp/tmp.xyz'):
        """
        Runs a clustering analysis avoer models. Clusters found will be stored
           in ThreeDeeModels.clusters
        
        :param 0.75 fact: Factor for equivalent positions
        :param 200 dcutoff: Distance Cut-off for TMScore.
        :param 'mcl' mcl_bin: path to MCL executable file
        :param '/tmp/tmp.xyz' tmp_file: path to a temporary file
        
        """
        scores = calc_eqv_rmsd(self.__models, self.nloci, dcutoff, fact,
                               tmp_file=tmp_file)
        # this may disappear if we use ward clustering
        out_f = open(tmp_file, 'w')
        for score in scores:
            out_f.write('{}\t{}\t{}\n'.format(*score))
        out_f.close()
        Popen('{0} {1} --abc -V all -o {1}.mcl'.format(
            mcl_bin, tmp_file, stdout=PIPE, stderr=PIPE),
              shell=True).communicate()
        clusters= {}
        for cluster, line in enumerate(open(tmp_file + '.mcl')):
            clusters[cluster] = []
            for model in line.split():
                model = int(model.split('_')[1])
                clusters[cluster].append(model)
        # sort clusters according to their lowest energy
        for clt in clusters:
            clusters[clt].sort()
        for i, clt in enumerate(sorted(
            clusters, key=lambda x: self[clusters[x][0]]['energy'])):
            self.clusters[i] = clusters[clt]
            for model in self.clusters[i]:
                self.__models[model]['cluster'] = i


    def _build_distance_matrix(self, n_best_clusters):
        """
        """
        clusters = sorted(self.clusters.keys())[:n_best_clusters]
        matrix = [[0.0 for _ in clusters] for _ in clusters]
        clust_count = dict([(c, len([m for m in self.__models if self.__models[m]['cluster'] == c])) for c in clusters])
        energy = dict([(c, self.__models[[m for m in self.__models if self.__models[m]['cluster'] == c][0]]['log_energies'][-1]) for c in clusters])
        for i, cl1 in enumerate(clusters):
            # find model with lowest energy for each cluster
            for md1 in self.__models:
                if self.__models[md1]['cluster'] == cl1:
                    # the first on found is the best :)
                    break
            for j, cl2 in enumerate(clusters[i+1:]):
                # find model with lowest energy for each cluster
                for md2 in self.__models:
                    if self.__models[md2]['cluster'] == cl2:
                        # the first on found is the best :)
                        break
                matrix[i][j+i+1] = float(calc_eqv_rmsd({1: self.__models[md1],
                                                        2: self.__models[md2]},
                                                       self.nloci, var='drmsd')[0])
        return clust_count, energy, matrix


    def cluster_analysis_dendrogram(self, n_best_clusters=None, color=False):
        """
        Representation of clusters of models. The length of the leaves if
           proportional to the final energy of each model. Branch widths are
           proportional to the number of models in a given cluster (or group of
           clusters, ifit is an internal branch)
           
        :param None n_best_clusters: number of clusters to represent, by default
           all clusters will be shown
        :param False color: display colors for the significance of the
           clustering (basically dependent of internal branch lengths)
        """

        if not self.clusters:
            self.cluster_models()
        from scipy.cluster.hierarchy import linkage
        if not n_best_clusters:
            n_best_clusters = len(self.clusters)
        clust_count, energy, matrix = self._build_distance_matrix(n_best_clusters)
        z = linkage(matrix)
        minnrj = min(energy.values()) - 1
        maxnrj = max(energy.values()) - 1
        val = (maxnrj-minnrj)
        maxz = max([i[2] for i in z])
        for i in z:
            i[2] = i[2]/maxz * val

        dad = {}
        i = max(clust_count)
        for a, b, _, _ in z:
            i += 1
            clust_count[i] = clust_count[a] + clust_count[b]
            dad[a] = i
            dad[b] = i

        my_count = {}
        for i in range(max(clust_count)+1):
            my_count[i] = clust_count[i]
        d = augmented_dendrogram(my_count, dad, energy, color, z)
        return d


    def density_plot(self, models=None, cluster=None, steps=(1, 2, 3, 4, 5),
                     error=False):
        """
        Represents the number of nucletotide base pairs can be found in 1 nm of
           chromatine, this along strand modelled.

        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param (1, 2, 3, 4, 5) steps: how many particles to group for the
           estimation. By default 5 curves are drawn.
        :param False error: represent the error of the estimates.
        """
        if len(steps) > 7:
            raise Exception('Sorry not enough colors to do this.\n')
        colors = ['grey', 'darkgreen', 'darkblue', 'purple', 'darkorange', 'darkred', 'red'][-len(steps):]
        dists = []
        for part1, part2 in zip(range(self.nloci - 1), range(1, self.nloci)):
            dists.append(self.average_3d_dist(part1, part2, models, cluster,
                                              plot=False, median=False))
        distsk = {1: dists}
        for k in (steps[1:] if steps[0]==1 else steps):
            distsk[k] = [None for _ in range(k/2)]
            for i in range(self.nloci - k):
                distsk[k].append(reduce(lambda x, y: x + y,
                                        [dists[i+j] for j in range(k)]))
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
        from matplotlib import pyplot as plt
        fig = plt.figure()
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
            plots += plt.plot(distsk[k], color=colors[steps.index(k)],
                              lw=steps.index(k) + 1, alpha=0.5)
        if error:
            for k in steps:
                plots += plt.plot(errorp[k], color=colors[steps.index(k)], ls='--')
                plt.plot(errorn[k], color=colors[steps.index(k)], ls='--')
        plt.ylabel('Density (bp / nm)')
        plt.xlabel('Particle number')
        plt.legend(plots, ['Average for {} particle{}'.format(k,
                                                             's' if k else '')
                           for k in steps] + (
                       ['+/- 2 standard deviations for {}'.format(k)
                           for k in steps] if error else []), fontsize='small')
        plt.xlim((0, self.nloci))
        plt.title('Chromatine density')
        plt.show()


    def contact_map_consistency(self, models=None, cluster=None, cutoff=150):
        """
        Draws a heatmap representing the proportion of times two particles are
           closer than a given cutoff.

        :param None models: If None (default) will do calculate the distance
           along all models. A list of numbers corresponding to a given set of
           models can also be passed.
        :param None cluster: A number can be passed in order to calculate the
           distance between particles in all models corresponding to the cluster
           number 'cluster'
        :param 150 cutoff: distance in nanometers from which it is considered
           that two particles are separated.
           
        """
        if models:
            models=models
        elif cluster > -1:
            models = self.clusters[cluster]
        else:
            models = self.__models
        matrix = [[100.0 for _ in xrange(self.nloci)] for _ in xrange(self.nloci)]
        for i in xrange(self.nloci):
            for j in xrange(i + 1, self.nloci):
                val = len([k for k in self.average_3d_dist(i, j, plot=False,
                                                           median=False)
                           if k < cutoff])
                matrix[i][j] = matrix[j][i] = float(val) / len(models) * 100
        from matplotlib import pyplot as plt
        ims = plt.imshow(matrix, origin='lower', interpolation="nearest")
        plt.ylabel('Particles')
        plt.xlabel('Particles')
        cbar = plt.colorbar(ims)
        cbar.ax.set_yticklabels(['{:>3}%'.format(p) for p in range(0, 110, 10)])
        cbar.ax.set_ylabel('Percentage of models with particles closer than {} nm'.format(cutoff))
        plt.title('Contact map consistency')
        plt.show()


    def average_3d_dist(self, part1, part2, models=None, cluster=None,
                        plot=True, median=True):
        """
        Computes the distance between two particles. This is done by averaging
           between a set of given models.
        
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
            models = self.clusters[cluster]
        else:
            models = self.__models
        for model in models:
            dists.append(sqrt((self.__models[model]['x'][part1] -
                               self.__models[model]['x'][part2])**2 + 
                              (self.__models[model]['y'][part1] -
                               self.__models[model]['y'][part2])**2 +
                              (self.__models[model]['z'][part1] -
                               self.__models[model]['z'][part1])**2))
        if not plot:
            if median:
                return np_median(dists)
            else:
                return dists
        plot_hist_box(dists, part1, part2)


    def objective_function_model(self, model, log=False, smooth=True):
        """
        Plots the fall in energy through the Monte-Carlo search.

        :param model: the number of the model to plot
        :param False log: to plot in log scale
        :param True smooth: to smooth the curve
        """
        self.__models[model].objective_function(log=log, smooth=smooth)
        

    def write_cmm(self, model_num, directory, color=color_residues):
        """
        Writes cmm file read by Chimera (http://www.cgl.ucsf.edu/chimera).
        
        :param model: the number of the model to write
        :param directory: where to write the file (note: the name of the file
           will be model_1.cmm if model number is 1)
        :param color_residues color: either a coloring function like
           :func:`pytadbit.imp.imp_model.color_residues` or a list of (r, g, b)
           tuples (as long as the number of particles). 
        """
        try:
            model = self.__models[model_num]
        except KeyError:
            model = self._bad_models[model_num]
        if type(color) != list:
            color = color(self.nloci)
        out = '<marker_set name=\"marker set $marker\">\n'
        form = ('<marker id=\"{0}\" x=\"{1}\" y=\"{2}\" z=\"{3}\" r=\"{4}\" ' +
                'g=\"{5}\" b=\"{6}\" radius=\"0.5\" note=\"{0}\"/>\n')
        for n in xrange(self.nloci):
            out += form.format(n + 1,
                               model['x'][n], model['y'][n], model['z'][n],
                               color[n][0], color[n][1], color[n][2])
        form = ('<link id1=\"{}\" id2=\"{}\" r=\"1\" ' +
                'g=\"1\" b=\"1\" radius=\"0.1\"/>\n')
        for n in xrange(1, self.nloci):
            out += form.format(n, n + 1)
        out += '</marker_set>\n'

        out_f = open('{}/model_{}_rnd{}.cmm'.format(
            directory, model_num, model['rand_init']), 'w')
        out_f.write(out)
        out_f.close()


    def write_xyz(self, model_num, directory):
        """
        Writes xyz file containing the 3D coordinates of each particles.
        
        :param model: the number of the model to write
        :param directory: where to write the file (note: the name of the file
           will be model_1.cmm if model number is 1)
        """
        try:
            model = self.__models[model_num]
        except KeyError:
            model = self._bad_models[model_num]
        out = ''
        form = "{:>12}{:>12}{:>12.3f}{:>12.3f}{:>12.3f}\n"
        for n in xrange(self.nloci):
            out += form.format('p' + str(n + 1), n + 1, model['x'][n],
                               model['y'][n], model['z'][n])
        out_f = open('{}/model_{}_rnd{}.cmm'.format(
            directory, model_num, model['rand_init']), 'w')
        out_f.write(out)
        out_f.close()


    def pickle_models(self, path):
        """
        Saves all models in pickle format (python object written to disk).
        
        :param path: path to file where to pickle
        """
        out = open(path, 'w')
        dump(self, out)
        out.close()

