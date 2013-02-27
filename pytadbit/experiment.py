"""
20 Feb 2013


"""

from pytadbit.parsers.hic_parser import read_matrix
from pytadbit.utils import nicer
from pytadbit.parsers.tad_parser import parse_tads
from warnings import warn
from math import sqrt, log10
from scipy.stats.mstats import zscore


class Experiment(object):
    """
    Hi-C experiment.

    :param name: name of the experiment
    :param resolution: resolution of the experiment (size of a bin in bases)
    :param None xp_handler: whether a file or a list of lists corresponding to
       the hi-c data
    :param None tad_handler: a file or a dict with pre-calculated TADs for this
       experiment
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
    :param None max_tad_size: filter TADs longer than this value
       
    """


    def __init__(self, name, resolution, xp_handler=None, tad_handler=None,
                 parser=None, max_tad_size=None, no_warn=False):
        self.name       = name
        self.resolution = resolution
        self.hic_data   = None
        self.size       = None
        self.tads       = {}
        self.brks       = []
        self.wght       = None
        if xp_handler:
            self.load_experiment(xp_handler, parser)
        if tad_handler:
            self.load_tad_def(tad_handler, max_tad_size=max_tad_size)
        elif not xp_handler and not no_warn:
            warn('WARNING: this is an empty shell, no data here.\n')


    def __repr__(self):
        return 'Experiment {} (resolution: {}, TADs: {}, Hi-C rows: {})'.format(
            self.name, nicer(self.resolution), len(self.tads) or None,
            self.size)


    def load_experiment(self, handler, parser=None):
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
        
        """
        nums, size = read_matrix(handler, parser=parser)
        self.hic_data = nums
        self.size     = size
        

    def load_tad_def(self, handler, weights=None, max_tad_size=None):
        """
         Add Topologically Associated Domains definition detection to Slice
        
        :param f_name: path to file
        :param None name: name of the experiment, if None f_name will be used
        :param None weights: Store information about the weights, corresponding
           to the normalization of the Hi-C data (see tadbit function
           documentation)
        :param None max_tad_size: filter TADs longer than this value
        
        """
        tads = parse_tads(handler, max_size=max_tad_size,
                                     bin_size=self.resolution)
        self.tads = tads
        self.brks = [t['brk'] for t in tads.values() if t['brk']]
        self.wght  = weights or None
        

    def normalize_hic(self):
        """
        Normalize Hi-C data. This normalize step is an exact replicate of what
        is done inside :func:`pytadbit.tadbit.tadbit`,

        It fills the Experiment.wght variable.

        the weight of a given cell in column i and row j corresponds to the
        square root of the product of the sum of the column i by the sum of row
        j.

        the weight of the Hi-C count in row I, column J of the Hi-C matrix
        would be:
        ::
           
                                _________________________________________
                      \        / N                    N                  |
                       \      / ___                  ___             
         weight(I,J) =  \    /  \                    \           
                         \  /   /__ (matrix(J, i)) * /__  (matrix(j, I))
                          \/    i=0                  j=0

        
        N being the number or rows/columns of the Hi-C matrix
        """
        if not self.hic_data:
            raise Exception('ERROR: No Hi-C data loaded\n')
        if self.wght:
            warn('WARNING: removing previous weights\n')
        rowsums = []
        for i in xrange(self.size):
            i *= self.size
            rowsums.append(0)
            for j in xrange(self.size):
                rowsums[-1] += self.hic_data[0][i + j]
        self.wght = [[0 for _ in xrange(self.size * self.size)]]
        for i in xrange(self.size):
            for j in xrange(self.size):
                self.wght[0][i * self.size + j] = sqrt(rowsums[i] * rowsums[j])


    def write_interaction_pairs(self, fname, normalized=True, zscored=True):
        """
        Creates a tab separated file with all interactions
        
        :param fname: file name to write the interactions pairs 
        :param True zscored: computes the zscore of the log10(data)
        :param True normalized: use weights to normalize data
        """
        zeros =  [float(pos) for pos, raw in enumerate(xrange(0, self.size**2,
                                                              self.size))
                  if sum(self.hic_data[0][raw:raw + self.size]) == 0]
        values = []
        if normalized:
            for i in xrange(self.size):
                if i in zeros:
                    continue
                for j in xrange(self.size):
                    values.append(self.hic_data[0][i * self.size + j] *\
                                  self.wght[0][i * self.size + j])
        else:
            for i in xrange(self.size):
                if i in zeros:
                    continue
                for j in xrange(self.size):
                    values.append(self.hic_data[0][i * self.size + j])
        # remove negative values from list
        notin = []
        yesin = []
        for i, val in enumerate(values):
            if val < 0:
                notin.append(i)
            else:
                yesin.append(val)
        # compute Z-score
        if zscored:
            values = list(zscore([log10(v) if v > 0 else -99.0 if v else 0.0 \
                                  for v in yesin]))
        # put back negative values
        for i in notin:
            values.insert(i, -99)
        # write to file
        out = open(fname, 'w')
        if normalized:
            out.write('elt1\telt2\tnorm_count\n')
        else:
            out.write('elt1\telt2\tcount\n')
        iterval = values.__iter__()
        for i in xrange(self.size):
            if i in zeros:
                continue
            for j in xrange(self.size):
                out.write('{}\t{}\t{}\n'.format(
                    i + 1, j + 1, iterval.next()))
        out.close()
