"""
20 Feb 2013


"""

from pytadbit.parsers.hic_parser         import read_matrix
from pytadbit.utils.extraviews           import nicer
from pytadbit.utils.tadmaths             import zscore
from pytadbit.utils.hic_filtering        import hic_filtering_for_modelling
from pytadbit.parsers.tad_parser         import parse_tads
from warnings                            import warn
from math                                import sqrt
from pytadbit.imp.CONFIG                 import CONFIG

try:
    from pytadbit.imp.imp_modelling          import generate_3d_models
    from pytadbit.imp.modelling_optimization import grid_search
except ImportError:
    warn('IMP not found, check PYTHONPATH\n')

class Experiment(object):
    """
    Hi-C experiment.

    :param name: name of the experiment
    :param resolution: the resolution of the experiment (size of a bin in
       bases)
    :param None hic_data: whether a file or a list of lists corresponding to
       the Hi-C data
    :param None tad_def: a file or a dict with precomputed TADs for this
       experiment
    :param None parser: a parser function that returns a tuple of lists 
       representing the data matrix and the length of a row/column. With
       the file example.tsv:

       ::
       
         chrT_001	chrT_002	chrT_003	chrT_004
         chrT_001	629	164	88	105
         chrT_002	164	612	175	110
         chrT_003	88	175	437	100
         chrT_004	105	110	100	278

       the output of parser('example.tsv') would be be:
       ``[([629, 164, 88, 105, 164, 612, 175, 110, 88, 175, 437, 100, 105,
       110, 100, 278]), 4]``
    :param None conditions: :py:func:`list` of experimental conditions, e.g. 
       the cell type, the enzyme... (i.e.: ['HindIII', 'cortex', 'treatment']).
       This parameter may be used to compare the effect of this conditions on
       the TADs
    :param True filter_columns: filter the columns with unexpectedly high 
       content of low values

    TODO: doc conditions
    TODO: normalization
    """


    def __init__(self, name, resolution, hic_data=None, tad_def=None,
                 parser=None, no_warn=False, weights=None,
                 conditions=None, filter_columns=True):
        self.name            = name
        self.resolution      = resolution
        self.crm             = None
        self._ori_resolution = resolution
        self.hic_data        = None
        self._ori_hic        = None
        self.conditions      = sorted(conditions) if conditions else []
        self.size            = None
        self.tads            = {}
        self.wght            = None
        self._normalization  = None
        self._zeros          = None
        self._zscores        = {}
        if hic_data:
            self.load_experiment(hic_data, parser,
                              filter_columns=filter_columns)
        if tad_def:
            self.load_tad_def(tad_def, weights=weights)
        elif not hic_data and not no_warn:
            warn('WARNING: this is an empty shell, no data here.\n')


    def __repr__(self):
        return 'Experiment %s (resolution: %s, TADs: %s, Hi-C rows: %s, normalized: %s)' % (
            self.name, nicer(self.resolution), len(self.tads) or None,
            self.size, self._normalization if self._normalization else 'None')


    def __add__(self, other):
        """
        sum Hi-C data of experiments into a new one.
        """
        reso1, reso2 = self.resolution, other.resolution
        if self.resolution == other.resolution:
            resolution = self.resolution
        else:
            resolution = max(reso1, reso2)
            self.set_resolution(resolution)
            other.set_resolution(resolution)
            
        xpr = Experiment(name='%s+%s' % (self.name, other.name),
                         resolution=resolution,
                         hic_data=tuple([i + j for i, j in zip(
                             self.hic_data[0], other.hic_data[0])]))
        self.set_resolution(reso1)
        other.set_resolution(reso2)
        xpr.crm = self.crm
        return xpr


    def set_resolution(self, resolution, keep_original=True):
        """
        Set a new value for the resolution. Copy the original data into
        Experiment._ori_hic and replace the Experiment.hic_data
        with the data corresponding to new data 
        (:func:`pytadbit.Chromosome.compare_condition`).

        :param resolution: an integer representing the resolution. This number
           must be a multiple of the original resolution, and higher than it
        :param True keep_original: either to keep or not the original data

        """
        if resolution < self._ori_resolution:
            raise Exception('New resolution might be higher than original.')
        if resolution % self._ori_resolution:
            raise Exception('New resolution might be a multiple original.\n' +
                            '  otherwise it is too complicated for me :P')
        if resolution == self.resolution:
            return
        # if we want to go back to original resolution
        if resolution == self._ori_resolution:
            self.hic_data   = self._ori_hic
            self.size       = self.resolution / self._ori_resolution * self.size
            self.resolution = self._ori_resolution
            return
        # if current resolution is the original one
        if self.resolution == self._ori_resolution:
            self._ori_hic = self.hic_data[:]
        self.resolution = resolution
        fact = self.resolution / self._ori_resolution
        # super for!
        size = int(sqrt(len(self._ori_hic[0])))
        self.hic_data = [[]]
        self.size     = size / fact
        rest = size % fact
        if rest:
            self.size += 1
        for i in xrange(0, size, fact):
            for j in xrange(0, size, fact):
                val = 0
                for k in xrange(fact):
                    if i + k >= size:
                        break
                    for l in  xrange(fact):
                        if j + l >= size:
                            break
                        val += self._ori_hic[0][(i + k) * size + j + l]
                self.hic_data[0].append(val)
        # hic_data needs always to be stored as tuple
        self.hic_data[0] = tuple(self.hic_data[0])
        if not keep_original:
            del(self._ori_hic)


    def load_experiment(self, hic_data, parser=None, resolution=None,
                        filter_columns=True):
        """
        Add a Hi-C experiment to the Chromosome object.
        
        :param None hic_data: whether a file or a list of lists corresponding to
           the Hi-C data
        :param name: name of the experiment
        :param False force: overwrite the experiments loaded under the same 
           name
        :param None parser: a parser function that returns a tuple of lists
           representing the data matrix and the length of a row/column. 
           With the file example.tsv:

           ::
           
             chrT_001	chrT_002	chrT_003	chrT_004
             chrT_001	629	164	88	105
             chrT_002	86	612	175	110
             chrT_003	159	216	437	105
             chrT_004	100	111	146	278
           
           the output of parser('example.tsv') would be:
           ``[([629, 86, 159, 100, 164, 612, 216, 111, 88, 175, 437, 146, 105,
           110, 105, 278]), 4]``
        :param None resolution: resolution of the experiment in the file; it
           will be adjusted to the resolution of the experiment. By default the
           file is expected to contain a Hi-C experiment with the same resolution
           as the :class:`pytadbit.Experiment` created, and no change is made
        :param True filter_columns: filter the columns with unexpectedly high content
           of low values
        
        """
        nums, size = read_matrix(hic_data, parser=parser)
        self.hic_data = nums
        self.size     = size
        resolution = resolution or self.resolution
        self.set_resolution(resolution, keep_original=False)
        # self._zeros   = [int(pos) for pos, raw in enumerate(
        #     xrange(0, self.size**2, self.size))
        #                  if sum(self.hic_data[0][raw:raw + self.size]) <= 100]
        if filter_columns:
            self._zeros = hic_filtering_for_modelling(self.get_hic_matrix())
        

    def load_tad_def(self, tad_def, weights=None):
        """
         Add the Topologically Associated Domains definition detection to Slice
        
        :param None tad_def: a file or a dict with precomputed TADs for this
           experiment
        :param None name: name of the experiment, if None f_name will be used
        :param None weights: Store information about the weights, corresponding
           to the normalization of the Hi-C data (see tadbit function
           documentation)
        
        """
        tads, wght = parse_tads(tad_def)
        self.tads = tads
        self.wght  = weights or wght
        

    def normalize_hic(self, method='sqrt'):
        """
        Normalize the Hi-C data. This normalization step does the same of
        the :func:`pytadbit.tadbit.tadbit` function (default parameters),

        It fills the Experiment.wght variable with the Hi-C values divided by
        the calculated weight.

        The weight of a given cell in column i and row j corresponds to the
        square root of the product of the sum of column i by the sum of row
        j.

        :param sqrt method: either 'sqrt' or 'visibility'. Depending on this
           parameter, the weight of the Hi-C count in row I, column J of the
           Hi-C matrix will be, under 'sqrt':
           ::
              
                                   _________________________________________
                         \        / N                    N                  |
                          \      / ___                  ___             
            weight(I,J) =  \    /  \                    \           
                            \  /   /__ (matrix(i, J)) * /__  (matrix(I, j))
                             \/    i=0                  j=0

           and under 'visibility': 
           ::
   
                             N                    N                 
                            ___                  ___                
                            \                    \                  
                            /__ (matrix(i, J)) * /__  (matrix(I, j))
                            i=0                  j=0                
            weight(I,J) =  -----------------------------------------         
                                     N     N                                 
                                    ___   ___                                
                                    \     \                                  
                                    /__   /__ (matrix(i, j))
                                    j=0   i=0                                
   
   
           
           with N being the number or rows/columns of the Hi-C matrix in both
           cases.

           Note that the default behavior (also used in
           :func:`pytadbit.tadbit.tadbit`)
           corresponds to method='sqrt'.
        """
        if not self.hic_data:
            raise Exception('ERROR: No Hi-C data loaded\n')
        if not method in ['sqrt', 'visibility']:
            raise LookupError('Only "sqrt" and "visibility" methods are implemented')
        if self.wght:
            warn('WARNING: removing previous weights\n')
        # removes columns where there is no data in the diagonal
        forbidden = [i for i in xrange(self.size)
                     if not self.hic_data[0][i*self.size+i]]
        rowsums = []
        for i in xrange(self.size):
            i *= self.size
            rowsums.append(0)
            if i in forbidden:
                continue
            for j in xrange(self.size):
                if not j in forbidden:
                    rowsums[-1] += self.hic_data[0][i + j]
        self.wght = [[0. for _ in xrange(self.size * self.size)]]
        if method == 'visibility':
            total = sum(rowsums)
            func = lambda x, y: float(rowsums[x] * rowsums[y]) / total
        elif method == 'sqrt':
            func = lambda x, y: sqrt(rowsums[x] * rowsums[y])
        for i in xrange(self.size):
            for j in xrange(self.size):
                if i in forbidden or j in forbidden:
                    self.wght[0][i * self.size + j] = 0.0
                else:
                    try:
                        self.wght[0][i * self.size + j] = (
                            self.hic_data[0][i * self.size + j] / func(i, j))
                    except ZeroDivisionError:
                        self.wght[0][i * self.size + j] = 0.0
        self._normalization = method


    def get_hic_zscores(self, normalized=True, zscored=True, remove_zeros=True):
        """
        Normalize the Hi-C raw data. The result will be stored into
        the private Experiment._zscore list.

        :param True normalized: whether to normalize the result using the
           weights (see :func:`normalize_hic`)
        :param True zscored: calculate the z-score of the data
        :param True remove_zeros: remove null interactions
        
        """
        values = []
        zeros = {}
        self._zscores = {}
        if normalized:
            for i in xrange(self.size):
                # zeros are rows or columns having a zero in the diagonal
                if i in self._zeros:
                    continue
                for j in xrange(i + 1, self.size):
                    if j in self._zeros:
                        continue
                    if (not self.hic_data[0][i * self.size + j] 
                        or not self.hic_data[0][i * self.size + j])\
                        and remove_zeros:
                        zeros[(i, j)] = None
                        continue
                    values.append(self.wght[0][i * self.size + j])
        else:
            for i in xrange(self.size):
                if i in self._zeros:
                    continue
                for j in xrange(i + 1, self.size):
                    if j in self._zeros:
                        continue
                    values.append(self.hic_data[0][i * self.size + j])
        # compute Z-score
        if zscored:
            zscore(values, self.size)
        iterval = values.__iter__()
        for i in xrange(self.size):
            if i in self._zeros:
                continue
            for j in xrange(i + 1, self.size):
                if j in self._zeros:
                    continue
                if (i, j) in zeros and remove_zeros:
                    continue
                zsc = iterval.next()
                self._zscores.setdefault(str(i), {})
                self._zscores[str(i)][str(j)] = zsc
                # self._zscores.setdefault(j, {})
                # self._zscores[j][i] = zsc


    def model_region(self, start, end, n_models=5000, n_keep=1000, n_cpus=1,
                     verbose=0, keep_all=False, close_bins=1, outfile=None,
                     config=CONFIG['dmel_01']):
        """

        :param start:  first bin to model (bin number)
        :param end: last bin to model (bin number)
        :param 5000 n_models: number of modes to generate
        :param 1000 n_keep: number of models used in the final analysis 
           (usually the top 20% of the generated models). The models are ranked
           according to their objective function value (the lower the better)
        :param False keep_all: whether or not to keep the discarded models (if
           True, models will be stored under tructuralModels.bad_models)
        :param 1 close_bins: number of particles away (i.e. the bin number
           difference) a particle pair must be in order to be considered as
           neighbors (e.g. 1 means consecutive particles)
        :param n_cpus: number of CPUs to use
        :param 0 verbose: the information printed can be: nothing (0), the
           objective function value the selected models (1), the objective
           function value of all the models (2), all the modeling 
           information (3)
        :param CONFIG['dmel_01'] a dictionary containing the standard
           parameters used to generate the models. The dictionary should
           contain the keys kforce, maxdist, upfreq and lowfreq.
           Examples can be seen by doing:
           
           ::
           
             from pytadbit.imp.CONFIG import CONFIG

           where CONFIG is a dictionarry of dictionnaries to be passed to this
           function:
           
           :::
           
             CONFIG = {
              'dmel_01': {
                  # use these paramaters with the Hi-C data from:
                  'reference' : 'victor corces dataset 2013',
             
                  # Force applied to the restraints inferred to neighbor particles
                  'kforce'    : 5,
             
                  # Maximum experimental contact distance
                  'maxdist'   : 600, # OPTIMIZATION: 500-1200
             
                  # Minimum and maximum thresholds used to decide which experimental values have to be
                  # included in the computation of restraints. Z-score values bigger than upfreq
                  # and less that lowfreq will be include, whereas all the others will be rejected
                  'upfreq'    : 0.3, # OPTIMIZATION: min/max Z-score
             
                  'lowfreq'   : -0.7 # OPTIMIZATION: min/max Z-score
             
                  }
              }

        """
        if self._normalization != 'visibility':
            warn('WARNING: normalizing according to visibility method')
            self.normalize_hic(method='visibility')
        zscores, values = self._sub_experiment_zscore(start, end)
        return generate_3d_models(zscores, self.resolution, values=values,
                                  n_models=n_models, outfile=outfile,
                                  n_keep=n_keep, n_cpus=n_cpus, verbose=verbose,
                                  keep_all=keep_all, close_bins=close_bins,
                                  config=config)


    def optimal_imp_parameters(self, start, end, n_models=500, n_keep=100,
                               n_cpus=1, upfreq_range=(0, 1, 0.1), close_bins=1,
                               lowfreq_range=(-1, 0, 0.1),
                               scale_range=[0.005][:],
                               maxdist_range=(400, 1400), cutoff=300,
                               outfile=None, verbose=True):
        """
        Find the optimal set of parameters to be used for the 3D modeling in
        IMP.

        :param start: first bin to model (bin number)
        :param end: last bin to model (bin number)
        :param 500 n_models: number of modes to generate
        :param 100 n_keep: number of models used in the final analysis (usually
           the top 20% of the generated models). The models are ranked
           according to their objective function value (the lower the better)
        :param 1 close_bins: number of particles away (i.e. the bin number
           difference) a particle pair must be in order to be considered as
           neighbors (e.g. 1 means consecutive particles)
        :param n_cpus: number of CPUs to use
        :param False verbose: if set to True, information about the distance,
           force and Z-score between particles will be printed
        :param (-1,0,0.1) lowfreq_range:  range of lowfreq values to be 
           optimized. The last value of the input tuple is the incremental step
           for the lowfreq values
        :param (0,1,0.1,0.1) upfreq_range: range of upfreq values to be
           optimized. The last value of the input tuple is the incremental step
           for the upfreq values
        :param (400,1400,100) maxdist_range: upper and lower bounds used to
           search for the optimal maximum experimental distance. The last value
           of the input tuple is the incremental step for maxdist values 
        :param [0.01] scale_range: upper and lower bounds used to search for
           the optimal scale parameter (nm per nucleotide). The last value of
           the input tuple is the incremental step for scale parameter values
        :param True verbose: print the results to the standard output

        .. note::
        
          Each of the *_range* parameters accept tuples in the form
           *(start, end, step)*, or a list with the list of values to test.

           E.g.:
             * scale_range=[0.001, 0.005, 0.006] will test these three values.
             * scale_range=(0.001, 0.005, 0.001) will test the values 0.001,
               0.002, 0.003, 0.004 and 0.005

        :returns: a tuple containing:

             - a 3D numpy array with the values of the correlations calculated
             - the range of scale used
             - the range of maxdist used
             - the range of upfreq used
             - the range of lowfreq used

        """
        zscores, values = self._sub_experiment_zscore(start, end)
        (matrix, scale_arange, max_dist_arange,
         upfreq_arange, lowfreq_arange) = grid_search(
            upfreq_range=upfreq_range, lowfreq_range=lowfreq_range,
            scale_range=scale_range, zscores=zscores,
            resolution=self.resolution, values=values,
            maxdist_range=maxdist_range, n_cpus=n_cpus, n_models=n_models,
            n_keep=n_keep, cutoff=cutoff,
            close_bins=close_bins, verbose=verbose)
        if outfile:
            out = open(outfile, 'w')
            out.write('# scale\tmax_dist\tup_freq\tlow_freq\tcorrelation\n')
            for h, hh in enumerate(scale_arange):
                for i, ii in enumerate(max_dist_arange):
                    for j, jj in enumerate(upfreq_arange):
                        for k, kk in enumerate(lowfreq_arange):
                            out.write('%s\t%s\t%s\t%s\t%s\n' % (
                                hh, ii, jj, kk, matrix[h, i, j, k]))
            out.close()
        return (('scale', 'maxdist', 'upfreq', 'lowfreq'),
                (scale_arange, max_dist_arange, upfreq_arange, lowfreq_arange),
                matrix)

    
    def _sub_experiment_zscore(self, start, end):
        """
        Get the z-score of a sub-region of an  experiment.

        TODO: find a nicer way to do this...

        :param start: first bin to model (bin number)
        :param end: first bin to model (bin number)

        :returns: z-score and raw values of the experiment
        """
        if self._normalization != 'visibility':
            warn('WARNING: normalizing according to visibility method')
            self.normalize_hic(method='visibility')
        from pytadbit import Chromosome
        matrix = self.get_hic_matrix()
        end += 1
        new_matrix = [[] for _ in range(end-start)]
        for i in xrange(start, end):
            for j in xrange(start, end):
                new_matrix[i - start].append(matrix[i][j])
                
        tmp = Chromosome('tmp')
        tmp.add_experiment('exp1', hic_data=[new_matrix],
                           resolution=self.resolution, filter_columns=False)
        exp = tmp.experiments[0]
        # We want the weights and zeros calculated in the full chromosome
        siz = self.size
        exp.wght = [[self.wght[0][i + siz * j] for i in xrange(start, end)
                     for j in xrange(start, end)]]
        exp._zeros = dict([(z - start, None) for z in self._zeros
                           if start <= z <= end])
        if len(exp._zeros) == (end + 1 - start):
            raise Exception('ERROR: no interaction found in selected regions')
        # ... but the z-scores in this particular region
        exp.get_hic_zscores(remove_zeros=True)
        values = [[float('nan') for _ in xrange(exp.size)]
                  for _ in xrange(exp.size)]
        for i in xrange(exp.size):
            # zeros are rows or columns having a zero in the diagonal
            if i in exp._zeros:
                continue
            for j in xrange(i + 1, exp.size):
                if j in exp._zeros:
                    continue
                if (not exp.hic_data[0][i * exp.size + j] 
                    or not exp.hic_data[0][i * exp.size + j]):
                    continue
                values[i][j] = exp.wght[0][i * exp.size + j]
                values[j][i] = exp.wght[0][i * exp.size + j]
        return exp._zscores, values
        

    def write_interaction_pairs(self, fname, normalized=True, zscored=True,
                                diagonal=False, cutoff=None, header=False,
                                true_position=False, uniq=True,
                                remove_zeros=False):
        """
        Creates a tab separated file with all the pairwise interactions.
        
        :param fname: file name where to write the  pairwise interactions 
        :param True zscored: computes the z-score of the log10(data)
        :param True normalized: use the weights to normalize the data
        :param None cutoff: if defined, only the zscores above the cutoff will
           be writen to the file
        :param False uniq: only writes one representent per interacting pair
           
        """
        if not self._zscores and zscored:
            for i in xrange(self.size):
                for j in xrange(self.size):
                    self._zscores.setdefault(i, {})
                    self._zscores[i][j] = self.hic_data[0][i * self.size + j]
        if not self.wght:
            raise Exception('Experiment not normalized.')
        # write to file
        out = open(fname, 'w')
        if header:
            out.write('elt1\telt2\t%s\n' % ('zscore' if zscored else 
                                            'normalized hi-c' if normalized 
                                            else 'raw hi-c'))
        for i in xrange(self.size):
            if i in self._zeros:
                continue
            start = i if uniq else 0
            for j in xrange(start, self.size):
                if j in self._zeros:
                    continue
                if not diagonal and i == j:
                    continue
                if zscored:
                    try:
                        if self._zscores[i][j] < cutoff:
                            continue
                        if self._zscores[i][j] == -99:
                            continue
                    except KeyError:
                        continue
                    val = self._zscores[i][j]
                elif normalized:
                    val = self.wght[0][self.size*i+j]
                else:
                    val = self.hic_data[0][self.size*i+j]
                if remove_zeros and not val:
                    continue
                if true_position:
                    out.write('%s\t%s\t%s\n' % (self.resolution * (i + 1),
                                                self.resolution * (j + 1), val))
                else:
                    out.write('%s\t%s\t%s\n' % (i + 1, j + 1, val))
        out.close()


    def get_hic_matrix(self):
        """
        Return the Hi-C matrix.

        :returns: list of lists representing the Hi-C data matrix of the
           current experiment
        """
        siz = self.size
        hic = self.hic_data[0]
        return [[hic[i+siz * j] for i in xrange(siz)] for j in xrange(siz)]


    def print_hic_matrix(self, print_it=True):
        """
        Return the Hi-C matrix as string

        :returns: list of lists representing the Hi-C data matrix of the
           current experiment
        """
        siz = self.size
        hic = self.hic_data[0]
        out = '\n'.join(['\t'.join([str(hic[i+siz * j]) \
                                    for i in xrange(siz)]) \
                         for j in xrange(siz)])
        if print_it:
            print out
        else:
            return out + '\n'


    # def generate_densities(self):
    #     """
    #     Related to the generation of 3D models.
    #     In the case of Hi-C data, the density is equal to the number of
    #     nucleotides in a bin, which is equal to the experiment resolution.
    #     """
    #     dens = {}
    #     for i in self.size:
    #         dens[i] = self.resolution
    #     return dens
