"""
20 Feb 2013


"""

from pytadbit.parsers.hic_parser import read_matrix
from pytadbit.utils import nicer
from pytadbit.parsers.tad_parser import parse_tads
from warnings import warn


class Experiment(object):
    """
    Hi-C experiment.

    :param name: name of the experiment
    :param resolution: resolution of the experiment (size of a bin in bases)
    :param None xp_handler: wether a file or a list of lists corresponding to
       the hi-c data
    :param None tad_handler: a file or a dict with precalculated TADs for this
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
            warn('WARNING: thi is an empty shell, no data here.\n')


    def __repr__(self):
        return 'Experiment {} (resolution: {}, TADs: {}, Hi-C rows: {})'.format(
            self.name, nicer(self.resolution), len(self.tads) or None, self.size)


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
         Add Topologically Associated Domains defintinion detection to Slice
        
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
        
        
        
