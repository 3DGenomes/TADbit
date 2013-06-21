"""
20 Feb 2013


"""

from pytadbit.parsers.hic_parser import read_matrix
from pytadbit.utils import nicer, zscore, hic_filtering_for_modelling
from pytadbit.parsers.tad_parser import parse_tads
from warnings import warn
from math import sqrt


class Experiment(object):
    """
    Hi-C experiment.

    :param name: name of the experiment
    :param resolution: resolution of the experiment (size of a bin in bases)
    :param None xp_handler: whether a file or a list of lists corresponding to
       the hi-c data
    :param None tad_handler: a file or a dict with precomputed TADs for this
       experiment
    :param None parser: a parser function that returns a tuple of lists
       representing the data matrix, and the length of a row/column, with
       this file example.tsv:

       ::
       
         chrT_001	chrT_002	chrT_003	chrT_004
         chrT_001	629	164	88	105
         chrT_002	164	612	175	110
         chrT_003	88	175	437	100
         chrT_004	105	110	100	278

       the output of parser('example.tsv') might be:
       ``[([629, 164, 88, 105, 164, 612, 175, 110, 88, 175, 437, 100, 105,
       110, 100, 278]), 4]``
    :param None conditions: :py:func:`list` of experimental conditions, can be 
       the cell type, the enzyme... (i.e.: ['HindIII', 'cortex', 'treatment']).
       This parameter may be used to compare the effect of this conditions on
       the TADs.
       
    """


    def __init__(self, name, resolution, xp_handler=None, tad_handler=None,
                 parser=None, no_warn=False, weights=None,
                 conditions=None):
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
        self._zeros          = None
        self._zscores        = {}
        if xp_handler:
            self.load_experiment(xp_handler, parser)
        if tad_handler:
            self.load_tad_def(tad_handler, weights=weights)
        elif not xp_handler and not no_warn:
            warn('WARNING: this is an empty shell, no data here.\n')


    def __repr__(self):
        return 'Experiment {} (resolution: {}, TADs: {}, Hi-C rows: {})'.format(
            self.name, nicer(self.resolution), len(self.tads) or None,
            self.size)


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
            
        xpr = Experiment(name='{}+{}'.format(self.name, other.name),
                         resolution=resolution,
                         xp_handler=tuple([i + j for i, j in zip(
                             self.hic_data[0], other.hic_data[0])]))
        self.set_resolution(reso1)
        other.set_resolution(reso2)
        xpr.crm = self.crm
        return xpr


    def set_resolution(self, resolution, keep_original=True):
        """
        Set a new value for resolution. copy original data into
        Experiment._ori_hic and replaces the Experiment.hic_data
        with the data corresponding to new data 
        (:func:`pytadbit.Chromosome.compare_condition`).

        :param resolution: an integer, representing resolution. This number
           must be a multiple of the original resolution, and higher than it.
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
                    if i + k >= size: break
                    for l in  xrange(fact):
                        if j + l >= size: break
                        val += self._ori_hic[0][(i + k) * size + j + l]
                self.hic_data[0].append(val)
        # hic_data needs always to be stored as tuple
        self.hic_data[0] = tuple(self.hic_data[0])
        if not keep_original:
            del(self._ori_hic)


    def load_experiment(self, handler, parser=None, resolution=None):
        """
        Add Hi-C experiment to Chromosome
        
        :param f_name: path to tsv file
        :param name: name of the experiment
        :param False force: overwrite experiments loaded under the same name
        :param None parser: a parser function that returns a tuple of lists
           representing the data matrix, and the length of a row/column, with
           this file example.tsv:

           ::
           
             chrT_001	chrT_002	chrT_003	chrT_004
             chrT_001	629	164	88	105
             chrT_002	86	612	175	110
             chrT_003	159	216	437	105
             chrT_004	100	111	146	278
           
           the output of parser('example.tsv') might be:
           ``[([629, 86, 159, 100, 164, 612, 216, 111, 88, 175, 437, 146, 105,
           110, 105, 278]), 4]``
        :param None resolution: resolution of the experiment in the file, it
           will be adjusted to the resolution of the experiment. By default the
           file is expected to contain an hi-c experiment at the same resolution
           as the :class:`pytadbit.Experiment` created, and no change is made.
        
        """
        nums, size = read_matrix(handler, parser=parser)
        self.hic_data = nums
        self.size     = size
        resolution = resolution or self.resolution
        self.set_resolution(resolution, keep_original=False)
        # self._zeros   = [int(pos) for pos, raw in enumerate(
        #     xrange(0, self.size**2, self.size))
        #                  if sum(self.hic_data[0][raw:raw + self.size]) <= 100]
        self._zeros = hic_filtering_for_modelling(self.get_hic_matrix())
        

    def load_tad_def(self, handler, weights=None):
        """
         Add Topologically Associated Domains definition detection to Slice
        
        :param f_name: path to file
        :param None name: name of the experiment, if None f_name will be used
        :param None weights: Store information about the weights, corresponding
           to the normalization of the Hi-C data (see tadbit function
           documentation)
        
        """
        tads, wght = parse_tads(handler)
        self.tads = tads
        self.wght  = weights or wght
        

    def normalize_hic(self, method='sqrt'):
        """
        Normalize Hi-C data. This normalize step is an exact replicate of what
        is done inside :func:`pytadbit.tadbit.tadbit` (default parameters),

        It fills the Experiment.wght variable.

        the weight of a given cell in column i and row j corresponds to the
        square root of the product of the sum of the column i by the sum of row
        j.

        :param sqrt method: either 'sqrt' or 'over_tot'. Depending on this param
           the weight of the Hi-C count in row I, column J of the Hi-C matrix
           would be, under 'sqrt':
           ::
              
                                   _________________________________________
                         \        / N                    N                  |
                          \      / ___                  ___             
            weight(I,J) =  \    /  \                    \           
                            \  /   /__ (matrix(J, i)) * /__  (matrix(j, I))
                             \/    i=0                  j=0

           and under 'over_tot': 
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
   
   
           
           N being the number or rows/columns of the Hi-C matrix in both cases.

           Note that the default behavior (also used in
           :func:`pytadbit.tadbit.tadbit`)
           corresponds to method='sqrt'.
        """
        if not self.hic_data:
            raise Exception('ERROR: No Hi-C data loaded\n')
        if not method in ['sqrt', 'bytot']:
            raise LookupError('Only "sqrt" and "bytot" methods are implemented')
        if self.wght:
            warn('WARNING: removing previous weights\n')
        forbidden = [i for i in xrange(self.size) if not self.hic_data[0][i*self.size+i]]
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
        if method == 'bytot':
            total = sum(rowsums)
            func = lambda x, y: float(rowsums[x] * rowsums[y]) / total
        elif method == 'sqrt':
            func = lambda x, y: sqrt(rowsums[x] * rowsums[y])
        for i in xrange(self.size):
            for j in xrange(self.size):
                if i in forbidden or j in forbidden:
                    self.wght[0][i * self.size + j] = 0.0
                else:
                    self.wght[0][i * self.size + j] = func(i, j)


    def get_hic_zscores(self, normalized=True, zscored=True):
        """
        Computes a normalization of Hi-C raw data. Result will be stored into
        the private Experiment._zscore list

        :param True normalized: whether to normalize the result using the
           weights (see :func:`normalize_hic`)
        :param True zscored: apply a z-score transform over the data.
        
        """
        values = []
        if normalized:
            for i in xrange(self.size):
                if i in self._zeros:
                    continue
                for j in xrange(self.size):
                    if j in self._zeros:
                        continue
                    try:
                        values.append(
                            self.hic_data[0][i * self.size + j] /\
                            self.wght[0][i * self.size + j])
                    except ZeroDivisionError:
                        values.append(0.0)
        else:
            for i in xrange(self.size):
                if i in self._zeros:
                    continue
                for j in xrange(self.size):
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
            for j in xrange(self.size):
                if j in self._zeros:
                    continue
                zsc = iterval.next()
                self._zscores.setdefault(i, {})
                self._zscores[i][j] = zsc
                self._zscores.setdefault(j, {})
                self._zscores[j][i] = zsc


    def write_interaction_pairs(self, fname, normalized=True, zscored=True,
                                diagonal=False, cutoff=None,
                                true_position=False, uniq=False):
        """
        Creates a tab separated file with all interactions
        
        :param fname: file name to write the interactions pairs 
        :param True zscored: computes the z-score of the log10(data)
        :param True normalized: use weights to normalize data
        :param None cutoff: if defined, only zscores above the cutoff will be
           writen to file
        :param False uniq: only writes on representent of an interacting pair
           
        """
        if not self._zscores:
            for i in xrange(self.size):
                for j in xrange(self.size):
                    self._zscores.setdefault(i, {})
                    self._zscores[i][j] = self.hic_data[0][i * self.size + j]
        # write to file
        out = open(fname, 'w')
        out.write('elt1\telt2\t{}\n'.format('zscore' if zscored else \
                                            'normalized hi-c' if normalized \
                                            else 'raw hi-c'))
        for i in xrange(self.size):
            if i in self._zeros:
                continue
            start = i if uniq else 0
            for j in xrange(start, self.size):
                if j in self._zeros:
                    continue
                if self._zscores[i][j] < cutoff:
                    continue
                if self._zscores[i][j] == -99:
                    continue
                if diagonal and i==j:
                    continue
                if true_position:
                    out.write('{}\t{}\t{}\n'.format(self.resolution * (i + 1),
                    self.resolution * (j + 1),
                    self._zscores[i][j]))
                else:
                    out.write('{}\t{}\t{}\n'.format(i + 1, j + 1,
                                                    self._zscores[i][j]))
        out.close()


    def get_hic_matrix(self):
        """
        Returns the Hi-C matrix

        :returns: list of lists representing Hi-C data matrix of current
           experiment
        """
        siz = self.size
        hic = self.hic_data[0]
        return [[hic[i+siz * j] for i in xrange(siz)] for j in xrange(siz)]


    def print_hic_matrix(self, print_it=True):
        """
        Returns the Hi-C matrix as string

        :returns: list of lists representing Hi-C data matrix of current
           experiment
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


    def generate_densities(self):
        """
        Related to the generation of 3D models.
        In the case of Hi-C data, the density is equal to the number of
        nucleotides in a bin, that is equal to the resolution
        """
        dens = {}
        for i in self.size:
            dens[i] = self.resolution
        return dens



