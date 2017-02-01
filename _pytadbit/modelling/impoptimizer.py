"""
28 Aug 2013


"""
from pytadbit.modelling.imp_modelling    import generate_3d_models
from pytadbit.utils.extraviews     import plot_2d_optimization_result
from pytadbit.utils.extraviews     import plot_3d_optimization_result
from pytadbit.modelling.structuralmodels import StructuralModels
from cPickle                       import dump, load
from sys                           import stderr
import numpy           as np
import multiprocessing as mu


class IMPoptimizer(object):
    """
    This class optimizes a set of paramaters (scale, maxdist, lowfreq and
    upfreq) in order to maximize the correlation between the models generated 
    by IMP and the input data.

    :param experiment: an instance of the class pytadbit.experiment.Experiment
    :param start: first bin to model (bin number)
    :param end: last bin to model (bin number)
    :param 5000 n_models: number of models to generate
    :param 1000 n_keep: number of models used in the final analysis (usually 
       the top 20% of the generated models). The models are ranked according to
       their objective function value (the lower the better)
    :param 1 close_bins: number of particles away (i.e. the bin number 
       difference) a particle pair must be in order to be considered as 
       neighbors (e.g. 1 means consecutive particles)
    :param None cutoff: distance cutoff (nm) to define whether two particles
       are in contact or not, default is 2 times resolution, times scale.
    :param None container: restrains particle to be within a given object. Can 
       only be a 'cylinder', which is, in fact a cylinder of a given height to 
       which are added hemispherical ends. This cylinder is defined by a radius, 
       its height (with a height of 0 the cylinder becomes a sphere) and the 
       force applied to the restraint. E.g. for modeling E. coli genome (2 
       micrometers length and 0.5 micrometer of width), these values could be 
       used: ['cylinder', 250, 1500, 50], and for a typical mammalian nuclei
       (6 micrometers diameter): ['cylinder', 3000, 0, 50]
    """
    def __init__(self, experiment, start, end, n_models=500,
                 n_keep=100, close_bins=1, container=None):

        self.resolution = experiment.resolution
        (self.zscores,
         self.values, zeros) = experiment._sub_experiment_zscore(start, end)
        self.zeros = tuple([i not in zeros for i in xrange(end - start + 1)])
        self.nloci       = end - start + 1
        if not self.nloci == len(self.zeros):
            raise Exception('ERROR: in optimization, bad number of particles\n')
        self.n_models    = n_models
        self.n_keep      = n_keep
        self.close_bins  = close_bins

        self.scale_range   = []
        self.maxdist_range = []
        self.lowfreq_range = []
        self.upfreq_range  = []
        self.dcutoff_range = []
        self.container     = container
        self.results = {}


    def run_grid_search(self, 
                        upfreq_range=(0, 1, 0.1),
                        lowfreq_range=(-1, 0, 0.1),
                        maxdist_range=(400, 1500, 100),
                        scale_range=0.01,
                        dcutoff_range=2,
                        corr='spearman', off_diag=1,
                        savedata=None, n_cpus=1, verbose=True):
        """
        This function calculates the correlation between the models generated 
        by IMP and the input data for the four main IMP parameters (scale, 
        maxdist, lowfreq and upfreq) in the given ranges of values.
        
        :param n_cpus: number of CPUs to use
        :param (-1,0,0.1) lowfreq_range: range of lowfreq values to be 
           optimized. The last value of the input tuple is the incremental 
           step for the lowfreq values
        :param (0,1,0.1) upfreq_range: range of upfreq values to be optimized.
           The last value of the input tuple is the incremental step for the
           upfreq values
        :param (400,1400,100) maxdist_range: upper and lower bounds
           used to search for the optimal maximum experimental distance. The 
           last value of the input tuple is the incremental step for maxdist 
           values
        :param 0.01 scale_range: upper and lower bounds used to search for
           the optimal scale parameter (nm per nucleotide). The last value of
           the input tuple is the incremental step for scale parameter values
        :param 2 dcutoff_range: upper and lower bounds used to search for
           the optimal distance cutoff parameter (distance, in number of beads,
           from which to consider 2 beads as being close). The last value of the
           input tuple is the incremental step for scale parameter values
        :param None savedata: concatenate all generated models into a dictionary
           and save it into a file named by this argument
        :param True verbose: print the results to the standard output
        """
        if verbose:
            stderr.write('Optimizing %s particles\n' % self.nloci)
        if isinstance(maxdist_range, tuple):
            maxdist_step = maxdist_range[2]
            maxdist_arange = range(maxdist_range[0],
                                        maxdist_range[1] + maxdist_step,
                                        maxdist_step)
        else:
            if isinstance(maxdist_range, (float, int)):
                maxdist_range = [maxdist_range]
            maxdist_arange = maxdist_range
        #
        if isinstance(lowfreq_range, tuple):
            lowfreq_step = lowfreq_range[2]
            lowfreq_arange = np.arange(lowfreq_range[0],
                                            lowfreq_range[1] + lowfreq_step / 2,
                                            lowfreq_step)
        else:
            if isinstance(lowfreq_range, (float, int)):
                lowfreq_range = [lowfreq_range]
            lowfreq_arange = lowfreq_range
        #
        if isinstance(upfreq_range, tuple):
            upfreq_step = upfreq_range[2]
            upfreq_arange = np.arange(upfreq_range[0],
                                           upfreq_range[1] + upfreq_step / 2,
                                           upfreq_step)
        else:
            if isinstance(upfreq_range, (float, int)):
                upfreq_range = [upfreq_range]
            upfreq_arange = upfreq_range
        #
        if isinstance(scale_range, tuple):
            scale_step = scale_range[2]
            scale_arange = np.arange(scale_range[0],
                                          scale_range[1] + scale_step / 2,
                                          scale_step)
        else:
            if isinstance(scale_range, (float, int)):
                scale_range = [scale_range]
            scale_arange = scale_range
        #
        if isinstance(dcutoff_range, tuple):
            dcutoff_step = dcutoff_range[2]
            dcutoff_arange = np.arange(dcutoff_range[0],
                                          dcutoff_range[1] + dcutoff_step / 2,
                                          dcutoff_step)
        else:
            if isinstance(dcutoff_range, (float, int)):
                dcutoff_range = [dcutoff_range]
            dcutoff_arange = dcutoff_range

        # round everything
        if not self.maxdist_range:
            self.maxdist_range = [my_round(i) for i in maxdist_arange]
        else:
            self.maxdist_range = sorted([my_round(i) for i in maxdist_arange
                                         if not my_round(i) in self.maxdist_range] +
                                        self.maxdist_range)
        if not self.upfreq_range:
            self.upfreq_range  = [my_round(i) for i in upfreq_arange ]
        else:
            self.upfreq_range = sorted([my_round(i) for i in upfreq_arange
                                        if not my_round(i) in self.upfreq_range] +
                                       self.upfreq_range)
        if not self.lowfreq_range:
            self.lowfreq_range = [my_round(i) for i in lowfreq_arange]
        else:
            self.lowfreq_range = sorted([my_round(i) for i in lowfreq_arange
                                         if not my_round(i) in self.lowfreq_range] +
                                        self.lowfreq_range)
        if not self.scale_range:
            self.scale_range   = [my_round(i) for i in scale_arange  ]
        else:
            self.scale_range = sorted([my_round(i) for i in scale_arange
                                       if not my_round(i) in self.scale_range] +
                                      self.scale_range)
        if not self.dcutoff_range:
            self.dcutoff_range = [my_round(i) for i in dcutoff_arange]
        else:
            self.dcutoff_range = sorted([my_round(i) for i in dcutoff_arange
                                         if not my_round(i) in self.dcutoff_range] +
                                        self.dcutoff_range)
        # grid search
        models = {}
        count = 0
        if verbose:
            stderr.write('# %3s %6s %7s %7s %6s %7s %7s\n' % (
                                    "num", "upfrq", "lowfrq", "maxdist",
                                    "scale", "cutoff", "corr"))
        for scale in [my_round(i) for i in scale_arange]:
            for maxdist in [my_round(i) for i in maxdist_arange]:
                for upfreq in [my_round(i) for i in upfreq_arange]:
                    for lowfreq in [my_round(i) for i in lowfreq_arange]:
                        # check if this optimization has been already done
                        if (scale, maxdist, upfreq, lowfreq) in [
                            tuple(k[:4]) for k in self.results]:
                            k = [k for k in self.results
                                 if (scale, maxdist, upfreq,
                                     lowfreq) == tuple(k[:4])][0]
                            result = self.results[(scale, maxdist, upfreq,
                                                   lowfreq, k[-1])]
                            if verbose:
                                verb = '%5s %6s %7s %7s %6s %7s  ' % (
                                    'xx', upfreq, lowfreq, maxdist,
                                    scale, k[-1])
                                if verbose == 2:
                                    stderr.write(verb + str(round(result, 4))
                                                 + '\n')
                                else:
                                    print verb + str(round(result, 4))
                            continue
                        tmp = {'kforce'   : 5,
                               'lowrdist' : 100,
                               'maxdist'  : int(maxdist),
                               'upfreq'   : float(upfreq),
                               'lowfreq'  : float(lowfreq),
                               'scale'    : float(scale)}
                        try:
                            count += 1
                            tdm = generate_3d_models(
                                self.zscores, self.resolution,
                                self.nloci, n_models=self.n_models,
                                n_keep=self.n_keep, config=tmp,
                                n_cpus=n_cpus, first=0,
                                values=self.values, container=self.container,
                                close_bins=self.close_bins, zeros=self.zeros)
                            result = 0
                            cutoff = my_round(dcutoff_arange[0])

                            matrices = tdm.get_contact_matrix(
                                cutoff=[int(i * self.resolution * float(scale))
                                        for i in dcutoff_arange])
                            for m in matrices:
                                cut = int(m**0.5)
                                sub_result = tdm.correlate_with_real_data(
                                    cutoff=cut, corr=corr,
                                    off_diag=off_diag, contact_matrix=matrices[m])[0]
                                if result < sub_result:
                                    result = sub_result
                                    cutoff = my_round(float(cut) / self.resolution /
                                                      float(scale))
                        except Exception, e:
                            print '  SKIPPING: %s' % e
                            result = 0
                            cutoff = my_round(dcutoff_arange[0])
                        if verbose:
                            verb = '%5s %6s %7s %7s %6s %7s  ' % (
                                count, upfreq, lowfreq, maxdist,
                                scale, cutoff)
                            if verbose == 2:
                                stderr.write(verb + str(round(result, 4))
                                             + '\n')
                            else:
                                print verb + str(round(result, 4))
                        # store
                        self.results[(scale, maxdist,
                                      upfreq, lowfreq, cutoff)] = result
                        if savedata and result:
                            models[(scale, maxdist, upfreq, lowfreq, cutoff)
                                   ] = tdm._reduce_models(minimal=True)
        if savedata:
            out = open(savedata, 'w')
            dump(models, out)
            out.close()
        self.scale_range.sort(  key=float)
        self.maxdist_range.sort(key=float)
        self.lowfreq_range.sort(key=float)
        self.upfreq_range.sort( key=float)
        self.dcutoff_range.sort(key=float)


    def load_grid_search(self, filenames, corr='spearman', off_diag=1,
                         verbose=True, n_cpus=1):
        """
        Loads one file or a list of files containing pre-calculated Structural
        Models (keep_models parameter used). And correlate each set of models
        with real data. Usefull to run different correlation on the same data
        avoiding to re-calculate each time the models.

        :param filenames: either a path to a file or a list of paths.
        :param spearman corr: correlation coefficient to use
        'param 1 off_diag:
        :param True verbose: print the results to the standard output

        """
        if isinstance(filenames, str):
            filenames = [filenames]
        models = {}
        for filename in filenames:
            inf = open(filename)
            models.update(load(inf))
            inf.close()
        count = 0
        pool = mu.Pool(n_cpus, maxtasksperchild=1)
        jobs = {}
        for scale, maxdist, upfreq, lowfreq, dcutoff in models:
            svd = models[(scale, maxdist, upfreq, lowfreq, dcutoff)]
            jobs[(scale, maxdist, upfreq, lowfreq, dcutoff)] = pool.apply_async(
                _mu_correlate, args=(svd, corr, off_diag,
                                     scale, maxdist, upfreq, lowfreq, dcutoff,
                                     verbose, count))
            count += 1
        pool.close()
        pool.join()
        for scale, maxdist, upfreq, lowfreq, dcutoff in models:
            self.results[(scale, maxdist, upfreq, lowfreq, dcutoff)] = \
                                 jobs[(scale, maxdist, upfreq, lowfreq, dcutoff)].get()
            if not scale in self.scale_range:
                self.scale_range.append(scale)
            if not maxdist in self.maxdist_range:
                self.maxdist_range.append(maxdist)
            if not lowfreq in self.lowfreq_range:
                self.lowfreq_range.append(lowfreq)
            if not upfreq in self.upfreq_range:
                self.upfreq_range.append(upfreq)
            if not dcutoff in self.dcutoff_range:
                self.dcutoff_range.append(dcutoff)
        self.scale_range.sort(  key=float)
        self.maxdist_range.sort(key=float)
        self.lowfreq_range.sort(key=float)
        self.upfreq_range.sort( key=float)
        self.dcutoff_range.sort(key=float)


    def get_best_parameters_dict(self, reference=None, with_corr=False):
        """
        :param None reference: a description of the dataset optimized
        :param False with_corr: if True, returns also the correlation value

        :returns: a dict that can be used for modelling, see config parameter in
           :func:`pytadbit.experiment.Experiment.model_region`
           
        """
        if not self.results:
            stderr.write('WARNING: no optimization done yet\n')
            return
        best = ((None, None, None, None), 0.0)
        for (sca, mxd, ufq, lfq, cut), val in self.results.iteritems():
            if val > best[-1]:
                best = ((sca, mxd, ufq, lfq, cut), val)
        if with_corr:
            return (dict((('scale'  , float(best[0][0])),
                          ('maxdist', float(best[0][1])),
                          ('upfreq' , float(best[0][2])),
                          ('lowfreq', float(best[0][3])),
                          ('dcutoff', float(best[0][4])),
                          ('reference', reference or ''), ('kforce', 5))),
                    best[-1])
        else:
            return dict((('scale'  , float(best[0][0])),
                         ('maxdist', float(best[0][1])),
                         ('upfreq' , float(best[0][2])),
                         ('lowfreq', float(best[0][3])),
                         ('dcutoff', float(best[0][4])),
                         ('reference', reference or ''), ('kforce', 5)))
    

    def plot_2d(self, axes=('scale', 'maxdist', 'upfreq', 'lowfreq'),
                show_best=0, skip=None, savefig=None,clim=None):
        """
        A grid of heatmaps representing the result of the optimization.

        :param 'scale','maxdist','upfreq','lowfreq' axes: list of axes to be
           represented in the plot. The order will define which parameter will
           be placed on the x, y, z or w axe.
        :param 0 show_best: number of best correlation values to highlight in 
           the plot
        :param None skip: if passed (as a dictionary), fix a given axe,
           e.g.: {'scale': 0.001, 'maxdist': 500}
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).

        """
        results = self._result_to_array()
        plot_2d_optimization_result((('scale', 'maxdist', 'upfreq', 'lowfreq'),
                                     ([float(i) for i in self.scale_range],
                                      [float(i) for i in self.maxdist_range],
                                      [float(i) for i in self.upfreq_range],
                                      [float(i) for i in self.lowfreq_range]),
                                     results), axes=axes, show_best=show_best,
                                    skip=skip, savefig=savefig,clim=clim)


    def plot_3d(self, axes=('scale', 'maxdist', 'upfreq', 'lowfreq')):
        """
        A grid of heatmaps representing the result of the optimization.

        :param 'scale','maxdist','upfreq','lowfreq' axes: tuple of axes to be
           represented in the plot. The order will define which parameter will 
           be placed on the x, y, z or w axe.

        """
        results = self._result_to_array()
        plot_3d_optimization_result((('scale', 'maxdist', 'upfreq', 'lowfreq'),
                                     ([float(i) for i in self.scale_range],
                                      [float(i) for i in self.maxdist_range],
                                      [float(i) for i in self.upfreq_range],
                                      [float(i) for i in self.lowfreq_range]),
                                     results), axes=axes)


    def _result_to_array(self):
        results = np.empty((len(self.scale_range), len(self.maxdist_range),
                            len(self.upfreq_range), len(self.lowfreq_range)))
        for w, scale in enumerate(self.scale_range):
            for x, maxdist in enumerate(self.maxdist_range):
                for y, upfreq in enumerate(self.upfreq_range):
                    for z, lowfreq in enumerate(self.lowfreq_range):
                        try:
                            cut = [c for c in self.dcutoff_range
                                   if (scale, maxdist, upfreq, lowfreq, c)
                                   in self.results][0]
                        except IndexError:
                            results[w, x, y, z] = float('nan')
                            continue
                        try:
                            results[w, x, y, z] = self.results[
                                (scale, maxdist, upfreq, lowfreq, cut)]
                        except KeyError:
                            results[w, x, y, z] = float('nan')
        return results


    def write_result(self, f_name):
        """
        This function writes a log file of all the values tested for each 
        parameter, and the resulting correlation value.

        This file can be used to load or merge data a posteriori using 
        the function pytadbit.modelling.impoptimizer.IMPoptimizer.load_from_file
        
        :param f_name: file name with the absolute path
        """
        out = open(f_name, 'w')
        out.write(('## n_models: %s n_keep: %s ' +
                   'close_bins: %s\n') % (self.n_models, 
                                          self.n_keep, self.close_bins))
        out.write('# scale\tmax_dist\tup_freq\tlow_freq\tdcutoff\tcorrelation\n')
        for scale in self.scale_range:
            for maxdist in self.maxdist_range:
                for upfreq in self.upfreq_range:
                    for lowfreq in self.lowfreq_range:
                        try:
                            cut = sorted(
                                [c for c in self.dcutoff_range
                                 if (scale, maxdist, upfreq, lowfreq, c)
                                 in self.results],
                                key=lambda x: self.results[
                                    (scale, maxdist, upfreq, lowfreq, x)])[0]
                        except IndexError:
                            print 'Missing dcutoff', (scale, maxdist, upfreq, lowfreq)
                            continue
                        try:
                            result = self.results[(scale, maxdist,
                                                   upfreq, lowfreq, cut)]
                            out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
                                scale, maxdist, upfreq, lowfreq, cut, result))
                        except KeyError:
                            print 'KeyError', (scale, maxdist, upfreq, lowfreq, cut, result)
                            continue
        out.close()


    def load_from_file(self, f_name):
        """
        Loads the optimized parameters from a file generated with the function:
        pytadbit.modelling.impoptimizer.IMPoptimizer.write_result.
        This function does not overwrite the parameters that were already 
        loaded or calculated.

        :param f_name: file name with the absolute path
        """
        for line in open(f_name):
            # Check same parameters
            if line.startswith('##'):
                n_models, _, n_keep, _, close_bins = line.split()[2:]
                if ([int(n_models), int(n_keep), int(close_bins)]
                    != 
                    [self.n_models, self.n_keep, self.close_bins]):
                    raise Exception('Parameters does in %s not match: %s\n%s' %(
                        f_name,
                        [int(n_models), int(n_keep), int(close_bins)],
                        [self.n_models, self.n_keep, self.close_bins]))
            if line.startswith('#'):
                continue
            scale, maxdist, upfreq, lowfreq, dcutoff, result = line.split()
            scale, maxdist, upfreq, lowfreq, dcutoff = (
                float(scale), int(maxdist), float(upfreq), float(lowfreq),
                float(dcutoff))
            scale   = my_round(scale, val=5)
            maxdist = my_round(maxdist)
            upfreq  = my_round(upfreq)
            lowfreq = my_round(lowfreq)
            dcutoff = my_round(dcutoff)
            self.results[(scale, maxdist, upfreq, lowfreq, dcutoff)] = float(result)
            if not scale in self.scale_range:
                self.scale_range.append(scale)
            if not maxdist in self.maxdist_range:
                self.maxdist_range.append(maxdist)
            if not upfreq in self.upfreq_range:
                self.upfreq_range.append(upfreq)
            if not lowfreq in self.lowfreq_range:
                self.lowfreq_range.append(lowfreq)
            if not dcutoff in self.dcutoff_range:
                self.dcutoff_range.append(dcutoff)
        self.scale_range.sort(  key=float)
        self.maxdist_range.sort(key=float)
        self.lowfreq_range.sort(key=float)
        self.upfreq_range.sort( key=float)
        self.dcutoff_range.sort(key=float)


def my_round(num, val=4):
    num = round(float(num), val)
    return str(int(num) if num == int(num) else num)


def _mu_correlate(svd, corr, off_diag, scale, maxdist, upfreq, lowfreq,
                  dcutoff, verbose, count):
    tdm = StructuralModels(
        nloci=svd['nloci'], models=svd['models'],
        bad_models=svd['bad_models'],
        resolution=svd['resolution'],
        original_data=svd['original_data'],
        clusters=svd['clusters'], config=svd['config'],
        zscores=svd['zscore'])
    try:
        result = tdm.correlate_with_real_data(
            cutoff=dcutoff, corr=corr,
            off_diag=off_diag)[0]
        if verbose:
            verb = '%5s  %s %s %s %s %s' % (
                count, upfreq, lowfreq, maxdist, scale, dcutoff)
            if verbose == 2:
                stderr.write(verb + str(result) + '\n')
            else:
                print verb + str(result)
    except Exception, e:
        print 'ERROR %s' % e
    return result

