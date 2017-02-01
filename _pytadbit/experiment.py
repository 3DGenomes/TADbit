"""
20 Feb 2013


"""

from pytadbit                     import HiC_data
from pytadbit.parsers.hic_parser  import read_matrix
from pytadbit.utils.extraviews    import nicer
from pytadbit.utils.extraviews    import tadbit_savefig
from pytadbit.utils.tadmaths      import zscore, nozero_log_matrix
from pytadbit.utils.normalize_hic import iterative
from pytadbit.utils.hic_filtering import hic_filtering_for_modelling
from pytadbit.parsers.tad_parser  import parse_tads
from math                         import isnan
from numpy                        import log2, array
from pytadbit.modelling.IMP_CONFIG import CONFIG
from copy                         import deepcopy as copy
from sys                          import stderr
from warnings                     import warn

try:
    from pytadbit.modelling.impoptimizer  import IMPoptimizer
    from pytadbit.modelling.imp_modelling import generate_3d_models
except ImportError:
    stderr.write('IMP not found, check PYTHONPATH\n')

try:
    import matplotlib.pyplot as plt
    from matplotlib.cm import jet
except ImportError:
    stderr.write('matplotlib not found\n')


def load_experiment_from_reads(name, fnam, genome_seq, resolution, 
                               conditions=None, identifier=None, cell_type=None,
                               enzyme=None, exp_type='Hi-C', **kw_descr):
    """
    Loads an experiment object from TADbit-generated read files, that are lists
    of pairs of reads mapped to a reference genome.
    
    :param fnam: tsv file with reads1 and reads2
    :param name: name of the experiment
    :param resolution: the resolution of the experiment (size of a bin in
       bases)
    :param None identifier: some identifier relative to the Hi-C data
    :param None cell_type: cell type on which the experiment was done
    :param None enzyme: restriction enzyme used in  the Hi-C experiment
    :param Hi-C exp_type: name of the experiment used (currently only Hi-C is
       supported)
    :param None conditions: :py:func:`list` of experimental conditions, e.g. 
       the cell type, the enzyme... (i.e.: ['HindIII', 'cortex', 'treatment']).
       This parameter may be used to compare the effect of this conditions on
       the TADs
    :param None kw_descr: any other argument passed would be stored as
       complementary descriptive field. For example::
       
           exp  = Experiment('k562_rep2', resolution=100000,
                             identifier='SRX015263', cell_type='K562',
                             enzyme='HindIII', cylce='synchronized')
           print exp

           # Experiment k562_rep2:
           #    resolution        : 100Kb
           #    TADs              : None
           #    Hi-C rows         : None
           #    normalized        : None
           #    identifier        : SRX015263
           #    cell type         : K562
           #    restriction enzyme: HindIII
           #    cylce             : synchronized

       *note that these fields may appear in the header of generated out files*
    """
    size = 0
    section_sizes = {}
    sections = []
    for crm in genome_seq:
        len_crm = int(float(len(genome_seq[crm])) / resolution + 1)
        section_sizes[(crm,)] = len_crm
        size += len_crm + 1
        sections.extend([(crm, '%04d' % i) for i in xrange(len_crm + 1)])
    imx = HiC_data((), size)
    dict_sec = dict([(j, i) for i, j in enumerate(sections)])
    for line in open(fnam):
        _, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.split('\t', 9)
        ps1 = dict_sec[(cr1, '%04d' % (int(ps1) / resolution))]
        ps2 = dict_sec[(cr2, '%04d' % (int(ps2) / resolution))]
        imx[ps1 + ps2 * size] += 1
        imx[ps2 + ps1 * size] += 1
    
    return Experiment(name, resolution=resolution, hic_data=imx,
                      conditions=conditions, identifier=identifier,
                      cell_type=cell_type, enzyme=enzyme, exp_type=exp_type,
                      **kw_descr)

class Experiment(object):
    """
    Hi-C experiment.

    :param name: name of the experiment
    :param resolution: the resolution of the experiment (size of a bin in
       bases)
    :param None identifier: some identifier relative to the Hi-C data
    :param None cell_type: cell type on which the experiment was done
    :param None enzyme: restriction enzyme used in  the Hi-C experiment
    :param Hi-C exp_type: name of the experiment used (currently only Hi-C is
       supported)
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
    :param None kw_descr: any other argument passed would be stored as
       complementary descriptive field. For example::
       
           exp  = Experiment('k562_rep2', resolution=100000,
                             identifier='SRX015263', cell_type='K562',
                             enzyme='HindIII', cylce='synchronized')
           print exp

           # Experiment k562_rep2:
           #    resolution        : 100Kb
           #    TADs              : None
           #    Hi-C rows         : None
           #    normalized        : None
           #    identifier        : SRX015263
           #    cell type         : K562
           #    restriction enzyme: HindIII
           #    cylce             : synchronized

       *note that these fields may appear in the header of generated out files*

    TODO: doc conditions
    TODO: normalization
    """


    def __init__(self, name, resolution, hic_data=None, norm_data=None,
                 tad_def=None, parser=None, no_warn=False, weights=None,
                 conditions=None, identifier=None,
                 cell_type=None, enzyme=None, exp_type='Hi-C', **kw_descr):
        self.name            = name
        self.resolution      = resolution
        self.identifier      = identifier
        self.cell_type       = cell_type
        self.enzyme          = enzyme
        self.description     = kw_descr
        self.exp_type        = exp_type
        self.crm             = None
        self._ori_resolution = resolution
        self.hic_data        = None
        self._ori_hic        = None
        self._ori_norm       = None
        self._ori_size       = None
        self.conditions      = sorted(conditions) if conditions else []
        self.size            = None
        self.tads            = {}
        self.norm            = None
        self.bias            = None 
        self._normalization  = None
        self._filtered_cols  = False
        self._zeros          = {}
        self._zscores        = {}
        if hic_data:
            self.load_hic_data(hic_data, parser, **kw_descr)
        if norm_data:
            self.load_norm_data(norm_data, parser, **kw_descr)
        if tad_def:
            self.load_tad_def(tad_def, weights=weights)
        elif not hic_data and not no_warn and not norm_data:
            stderr.write('WARNING: this is an empty shell, no data here.\n')


    def __repr__(self):
        return 'Experiment %s (resolution: %s, TADs: %s, Hi-C rows: %s, normalized: %s)' % (
            self.name, nicer(self.resolution), len(self.tads) or None,
            self.size, self._normalization if self._normalization else 'None')


    def __str__(self):
        outstr = 'Experiment %s:\n' % (self.name)
        outstr += '   resolution        : %s\n' % (nicer(self.resolution))
        outstr += '   TADs              : %s\n' % (len(self.tads) or None)
        outstr += '   Hi-C rows         : %s\n' % (self.size)
        outstr += '   normalized        : %s\n' % (self._normalization or None)
        ukw = 'UNKNOWN'
        try: # new in version post-CSDM13
            outstr += '   identifier        : %s\n' % (self.identifier or ukw)
            outstr += '   cell type         : %s\n' % (self.cell_type  or ukw)
            outstr += '   restriction enzyme: %s\n' % (self.enzyme     or ukw)
            for desc in self.description:
                outstr += '   %-18s: %s\n' % (desc, self.description[desc])
        except AttributeError:
            pass
        return outstr
        

    def __add__(self, other, silent=False):
        """
        sum Hi-C data of experiments into a new one.
        """
        reso1, reso2 = self.resolution, other.resolution
        if self.resolution == other.resolution:
            resolution = self.resolution
            changed_reso = False
        else:
            resolution = max(reso1, reso2)
            self.set_resolution(resolution)
            other.set_resolution(resolution)
            if not silent:
                stderr.write('WARNING: experiments of different resolution, ' +
                             'setting both resolution of %s, and normalizing ' +
                             'at this resolution\n' % (resolution))
            norm1 = copy(self.norm)
            norm2 = copy(other.norm)
            if self._normalization:
                self.normalize_hic()
            if other._normalization:
                other.normalize_hic()
            changed_reso = True
        if self.hic_data:
            new_hicdata = HiC_data([], size=self.size)
            for i in self.hic_data[0]:
                new_hicdata[i] = self.hic_data[0].get(i)
            for i in other.hic_data[0]:
                new_hicdata[i] += other.hic_data[0].get(i)
        else:
            new_hicdata = None
        xpr = Experiment(name='%s+%s' % (self.name, other.name),
                         resolution=resolution,
                         hic_data=new_hicdata, no_warn=True)
        # check if both experiments are normalized with the same method
        # and sum both normalized data
        if self._normalization != None and other._normalization != None:
            if (self._normalization.split('_factor:')[0] ==
                other._normalization.split('_factor:')[0]):
                xpr.norm = [HiC_data([], size=self.size)]
                for i in self.norm[0]:
                    xpr.norm[0][i] = self.norm[0].get(i)
                for i in other.norm[0]:
                    xpr.norm[0][i] += other.norm[0].get(i)
                # The final value of the factor should be the sum of each
                try:
                    xpr._normalization = (
                        self._normalization.split('_factor:')[0] +
                        '_factor:' +
                        str(int(self._normalization.split('_factor:')[1]) +
                            int(other._normalization.split('_factor:')[1])))
                except IndexError: # no factor there
                    xpr._normalization = (self._normalization)
        elif self.norm or other.norm:
            try:
                if (self.norm[0] or other.norm[0]) != {}:
                    if not silent:
                        raise Exception('ERROR: normalization differs between' +
                                        ' each experiment\n')
                else:
                    if not silent:
                        stderr.write('WARNING: experiments should be ' +
                                     'normalized before being summed\n')
            except TypeError:
                if not silent:
                    stderr.write('WARNING: experiments should be normalized ' +
                                 'before being summed\n')
        else:
            if not silent:
                stderr.write('WARNING: experiments should be normalized ' +
                             'before being summed\n')
        if changed_reso:
            self.set_resolution(reso1)
            self.norm = norm1
            other.set_resolution(reso2)
            other.norm = norm2
        xpr.crm = self.crm
        if not xpr.size:
            xpr.size = len(xpr.norm[0])
        
        def __merge(own, fgn):
            "internal function to merge descriptions"
            if own == fgn:
                return own
            return '%s+%s' % (own , fgn)
        
        xpr.identifier  = __merge(self.identifier , other.identifier )
        xpr.cell_type   = __merge(self.cell_type  , other.cell_type  )
        xpr.enzyme      = __merge(self.enzyme     , other.enzyme     )
        xpr.description = __merge(self.description, other.description)
        xpr.exp_type    = __merge(self.exp_type   , other.exp_type   )
        
        for des in self.description:
            if not des in other.description:
                continue
            xpr.description[des] = __merge(self.description[des],
                                           other.description[des])
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
            self.norm       = self._ori_norm
            self.size       = self._ori_size
            self.resolution = self._ori_resolution
            return
        # if current resolution is the original one
        if self.resolution == self._ori_resolution:
            if self.hic_data:
                self._ori_hic  = copy(self.hic_data)
            if self.norm:
                self._ori_norm = self.norm[:]
                # change the factor value in normalization description
                try:
                    self._normalization = (
                        self._normalization.split('_factor:')[0] +
                        '_factor:'+
                        str(int(self._normalization.split('factor:')[1])
                            * (resolution / self.resolution)))
                except IndexError: # no factor there
                    pass
        self.resolution = resolution
        fact = self.resolution / self._ori_resolution
        # super for!
        try:
            size = len(self._ori_hic[0])
        except TypeError:
            size = len(self._ori_norm[0])
        self.size     = size / fact
        rest = size % fact
        if rest:
            self.size += 1
        self.hic_data = [HiC_data([], size / fact + (1 if rest else 0))]
        self.norm     = [HiC_data([], size / fact + (1 if rest else 0))]
        def resize(mtrx, copee):
            "resize both hic_data and normalized data"
            for i in xrange(0, size, fact):
                for j in xrange(0, size, fact):
                    val = 0
                    for k in xrange(fact):
                        if i + k >= size:
                            break
                        for l in  xrange(fact):
                            if j + l >= size:
                                break
                            val += copee[(i + k) * size + j + l]
                    if val:
                        mtrx[i/fact * self.size + j/fact] = val
        try:
            resize(self.hic_data[0], self._ori_hic[0])
        except TypeError:
            pass
        try:
            resize(self.norm[0], self._ori_norm[0])
        except TypeError:
            pass
        # we need to recalculate zeros:
        if self._filtered_cols:
            stderr.write('WARNING: definition of filtered columns lost at ' +
                         'this resolution\n')
            self._filtered_cols = False
        if not keep_original:
            del(self._ori_hic)
            del(self._ori_norm)


    def filter_columns(self, silent=False, draw_hist=False, savefig=None,
                       diagonal=True, perc_zero=90, auto=True, min_count=None):
        """
        Call filtering function, to remove artefactual columns in a given Hi-C
        matrix. This function will detect columns with very low interaction
        counts; columns passing through a cell with no interaction in the
        diagonal; and columns with NaN values (in this case NaN will be replaced
        by zero in the original Hi-C data matrix). Filtered out columns will be
        stored in the dictionary Experiment._zeros.

        :param False silent: does not warn for removed columns
        :param False draw_hist: shows the distribution of mean values by column
           the polynomial fit, and the cut applied.
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param True diagonal: remove row/columns with zero in the diagonal
        :param 90 perc_zero: maximum percentage of cells with no interactions
           allowed.
        :param None min_count: minimum number of reads mapped to a bin (recommended
           value could be 2500). If set this option overrides the perc_zero
           filtering... This option is slightly slower.
        :param True auto: if False, only filters based on the given percentage
           zeros

        """
        try:
            data = self.hic_data[0]
        except:
            data = self.norm[0]
            diagonal = True
        self._zeros, has_nans = hic_filtering_for_modelling(
            data, silent=silent, draw_hist=draw_hist, savefig=savefig,
            diagonal=diagonal, perc_zero=perc_zero, auto=auto,
            min_count=min_count)
        if has_nans: # to make it simple
            for i in xrange(self.hic_data[0]._size2):
                if repr(self.hic_data[0][i]) == 'nan':
                    del(self.hic_data[0][i])
        # Also remove columns where there is no data in the diagonal
        size = self.size
        # else:
        #     self._zeros.update(dict([(i, None) for i in xrange(size)
        #                              if not self.norm[0][i * size + i]]))
        self._filtered_cols = True


    def load_hic_data(self, hic_data, parser=None, wanted_resolution=None,
                      data_resolution=None, silent=False, **kwargs):
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
        :param True filter_columns: filter the columns with unexpectedly high 
           content of low values
        :param False silent: does not warn for removed columns

        """
        self.hic_data = read_matrix(hic_data, parser=parser, one=False)
        self._ori_size       = self.size       = len(self.hic_data[0])
        self._ori_resolution = self.resolution = (data_resolution or
                                                  self._ori_resolution)
        wanted_resolution = wanted_resolution or self.resolution
        self.set_resolution(wanted_resolution, keep_original=False)
        if self.hic_data[0].bads:
            self._zeros = self.hic_data[0].bads
            self._filtered_cols = True

    def load_norm_data(self, norm_data, parser=None, resolution=None,
                       normalization='visibility', **kwargs):
        """
        Add a normalized Hi-C experiment to the Chromosome object.
        
        :param None norm_data: whether a file or a list of lists corresponding to
           the normalized Hi-C data
        :param name: name of the experiment
        :param False force: overwrite the experiments loaded under the same 
           name
        :param None parser: a parser function that returns a tuple of lists
           representing the data matrix and the length of a row/column. 
           With the file example.tsv:

           ::
           
             chrT_001	chrT_002	chrT_003	chrT_004
             chrT_001	12.5	164	8.8	0.5
             chrT_002	8.6	61.2	1.5	1.1
             chrT_003	15.9	21.6	3.7	0.5
             chrT_004	0.0	1.1	1.6	2.8
           
        :param None resolution: resolution of the experiment in the file; it
           will be adjusted to the resolution of the experiment. By default the
           file is expected to contain a Hi-C experiment with the same resolution
           as the :class:`pytadbit.Experiment` created, and no change is made
        :param True filter_columns: filter the columns with unexpectedly high 
           content of low values
        :param False silent: does not warn for removed columns
        
        """
        self.norm = read_matrix(norm_data, parser=parser, hic=False, one=False)
        self._ori_size       = self.size       = len(self.norm[0])
        self._ori_resolution = self.resolution = resolution or self._ori_resolution
        if not self._zeros: # in case we do not have original Hi-C data
            for i in xrange(self.size):
                if all([isnan(j) for j in
                        [self.norm[0][k] for k in
                         xrange(i * self.size, i * self.size + self.size)]]):
                    self._zeros[i] = None
        # remove NaNs, we do not need them as we have zeroes
        for i in self.norm[0].keys():
            if isnan(self.norm[0][i]):
                del(self.norm[0][i])
        self._normalization = normalization
        if self.norm[0].bads:
            self._zeros = self.norm[0].bads
            self._filtered_cols = True


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
        tads, norm = parse_tads(tad_def)
        last = max(tads.keys())
        if not self.size:
            self.size = tads[last]['end']
        self.tads = tads
        if not self.norm:
            self.norm  = weights or norm
            if self.norm:
                self._normalization = 'visibility'
        if self._normalization:
            norms = self.norm[0]
        elif self.hic_data:
            norms = self.hic_data[0]
        else:
            warn("WARNING: raw Hi-C data not available, " +
                 "TAD's height fixed to 1")
            norms = None
        diags = []
        siz = self.size
        sp1 = siz + 1
        zeros = self._zeros or {}
        if norms:
            for k in xrange(1, siz):
                s_k = siz * k
                diags.append(sum([norms[i * sp1 + s_k]
                                 if not (i in zeros
                                         or (i + k) in zeros) else 1. # 1 is the mean
                                  for i in xrange(siz - k)]) / (siz - k))
        for tad in tads:
            start, end = (int(tads[tad]['start']) + 1,
                          int(tads[tad]['end']) + 1)
            if norms:
                matrix = sum([norms[i + siz * j]
                             if not (i in zeros
                                     or j in zeros) else 1. 
                              for i in xrange(start - 1, end - 1)
                              for j in xrange(i + 1, end - 1)])
            try:
                if norms:
                    height = float(matrix) / sum(
                        [diags[i-1] * (end - start - i)
                         for i in xrange(1, end - start)])
                else:
                    height = tads[tad].get('height', 1.0)
            except ZeroDivisionError:
                height = 0.
            tads[tad]['height'] = height

    def normalize_hic(self, factor=1, iterations=0, max_dev=0.1, silent=False,
                      rowsums=None):
        """
        Normalize the Hi-C data. This normalization step does the same of
        the :func:`pytadbit.tadbit.tadbit` function (default parameters),

        It fills the Experiment.norm variable with the Hi-C values divided by
        the calculated weight.

        The weight of a given cell in column i and row j corresponds to the
        square root of the product of the sum of column i by the sum of row
        j.

        normalization is done according to this formula:

        .. math::

          weight_{(I,J)} = \\frac{\\sum^N_{j=0}\\sum^N_{i=0}(matrix(i,j))}
                                 {\\sum^N_{i=0}(matrix(i,J)) \\times \\sum^N_{j=0}(matrix(I,j))}
 
        with N being the number or rows/columns of the Hi-C matrix in both
        cases.
        
        :param 1 factor: final mean number of normalized interactions wanted
           per cell
        :param False silent: does not warn when overwriting weights
        :param None rowsums: input a list of rowsums calculated elsewhere
        """

        if not self.hic_data:
            raise Exception('ERROR: No Hi-C data loaded\n')
        if self.norm and not silent:
            stderr.write('WARNING: removing previous weights\n')
        size = self.size
        self.bias = iterative(self.hic_data[0], iterations=iterations,
                              max_dev=max_dev, bads=self._zeros,
                              verbose=not silent)
        self.norm = [HiC_data([(i + j * size, float(self.hic_data[0][i, j]) /
                                self.bias[i] /
                                self.bias[j] * size)
                               for i in self.bias for j in self.bias], size)]
        # no need to use lists, tuples use less memory
        if factor:
            self._normalization = 'visibility_factor:' + str(factor)
            factor = sum(self.norm[0].values()) / (self.size * self.size * factor)
            for n in self.norm[0]:
                self.norm[0][n] = self.norm[0][n] / factor
        else:
            self._normalization = 'visibility'


    def get_hic_zscores(self, normalized=True, zscored=True, remove_zeros=True):
        """
        Normalize the Hi-C raw data. The result will be stored into
        the private Experiment._zscore list.

        :param True normalized: whether to normalize the result using the
           weights (see :func:`normalize_hic`)
        :param True zscored: calculate the z-score of the data
        :param False remove_zeros: remove null interactions. Dangerous, null
           interaction are informative.

        """
        values = {}
        zeros  = {}
        self._zscores = {}
        if normalized:
            for i in xrange(self.size):
                # zeros are rows or columns having a zero in the diagonal
                if i in self._zeros:
                    continue
                for j in xrange(i + 1, self.size):
                    if j in self._zeros:
                        continue
                    if (not self.norm[0][i * self.size + j]
                        and remove_zeros):
                        zeros[(i, j)] = None
                        continue
                    values[(i, j)] = self.norm[0][i * self.size + j]
        else:
            for i in xrange(self.size):
                if i in self._zeros:
                    continue
                for j in xrange(i + 1, self.size):
                    if j in self._zeros:
                        continue
                    values[(i, j)] = self.hic_data[0][i * self.size + j]
        # compute Z-score
        if zscored:
            zscore(values)
        for i in xrange(self.size):
            if i in self._zeros:
                continue
            for j in xrange(i + 1, self.size):
                if j in self._zeros:
                    continue
                if (i, j) in zeros and remove_zeros:
                    continue
                self._zscores.setdefault(str(i), {})
                self._zscores[str(i)][str(j)] = values[(i, j)]


    def model_region(self, start=1, end=None, n_models=5000, n_keep=1000,
                     n_cpus=1, verbose=0, keep_all=False, close_bins=1,
                     outfile=None, config=CONFIG['dmel_01'],
                     container=None):
        """
        Generates of three-dimentional models using IMP, for a given segment of
        chromosome.
        
        :param 1 start: first bin to model (bin number)
        :param None end: last bin to model (bin number). By default goes to the
           last bin.
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
        :param None container: restrains particle to be within a given object. Can 
           only be a 'cylinder', which is, in fact a cylinder of a given height to 
           which are added hemispherical ends. This cylinder is defined by a radius, 
           its height (with a height of 0 the cylinder becomes a sphere) and the 
           force applied to the restraint. E.g. for modeling E. coli genome (2 
           micrometers length and 0.5 micrometer of width), these values could be 
           used: ['cylinder', 250, 1500, 50], and for a typical mammalian nuclei
           (6 micrometers diameter): ['cylinder', 3000, 0, 50]
        :param CONFIG['dmel_01'] config: a dictionary containing the standard
           parameters used to generate the models. The dictionary should
           contain the keys kforce, maxdist, upfreq and lowfreq.
           Examples can be seen by doing:
           
           ::
           
             from pytadbit.imp.CONFIG import CONFIG

           where CONFIG is a dictionarry of dictionnaries to be passed to this
           function:
           
           ::
           
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
             
                  # How much space (radius in nm) ocupies a nucleotide
                  'scale'     : 0.005
                  }
              }

        :returns: a :class:`pytadbit.imp.structuralmodels.StructuralModels` object.

        """
        if not self._normalization:
            stderr.write('WARNING: not normalized data, should run ' +
                         'Experiment.normalize_hic()\n')
        if not end:
            end = self.size
        zscores, values, zeros = self._sub_experiment_zscore(start, end)
        coords = {'crm'  : self.crm.name,
                  'start': start,
                  'end'  : end}
        zeros = tuple([i not in zeros for i in xrange(end - start + 1)])
        nloci = end - start + 1
        if verbose:
            stderr.write('Preparing to model %s particles\n' % nloci)
        return generate_3d_models(zscores, self.resolution, nloci,
                                  values=values, n_models=n_models,
                                  outfile=outfile, n_keep=n_keep, n_cpus=n_cpus,
                                  verbose=verbose, keep_all=keep_all, first=0,
                                  close_bins=close_bins, config=config, container=container,
                                  experiment=self, coords=coords, zeros=zeros)


    def optimal_imp_parameters(self, start=1, end=None, n_models=500, n_keep=100,
                               n_cpus=1, upfreq_range=(0, 1, 0.1), close_bins=1,
                               lowfreq_range=(-1, 0, 0.1),
                               scale_range=[0.01][:],
                               maxdist_range=(400, 1400, 100),
                               dcutoff_range=[2][:],
                               outfile=None, verbose=True, corr='spearman',
                               off_diag=1, savedata=None,
                               container=None):
        """
        Find the optimal set of parameters to be used for the 3D modeling in
        IMP.

        :param 1 start: first bin to model (bin number)
        :param None end: last bin to model (bin number). By default goes to the
           last bin.
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
        :param [2] dcutoff_range: upper and lower bounds used to search for
           the optimal distance cutoff parameter (distance, in number of beads,
           from which to consider 2 beads as being close). The last value of the
           input tuple is the incremental step for scale parameter values
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
        :param True verbose: print the results to the standard output

        .. note::
        
          Each of the *_range* parameters accept tuples in the form
           *(start, end, step)*, or a list with the list of values to test.

           E.g.:
             * scale_range=[0.001, 0.005, 0.006] will test these three values.
             * scale_range=(0.001, 0.005, 0.001) will test the values 0.001,
               0.002, 0.003, 0.004 and 0.005


        :returns: an :class:`pytadbit.imp.impoptimizer.IMPoptimizer` object

        """
        if not self._normalization:
            stderr.write('WARNING: not normalized data, should run ' +
                         'Experiment.normalize_hic()\n')
        if not end:
            end = self.size
        optimizer = IMPoptimizer(self, start, end, n_keep=n_keep,
                                 n_models=n_models, close_bins=close_bins,
                                 container=container)
        optimizer.run_grid_search(maxdist_range=maxdist_range,
                                  upfreq_range=upfreq_range,
                                  lowfreq_range=lowfreq_range,
                                  scale_range=scale_range,
                                  dcutoff_range=dcutoff_range, corr=corr,
                                  n_cpus=n_cpus, verbose=verbose,
                                  off_diag=off_diag, savedata=savedata)

        if outfile:
            optimizer.write_result(outfile)

        return optimizer


    def _sub_experiment_zscore(self, start, end):
        """
        Get the z-score of a sub-region of an  experiment.

        TODO: find a nicer way to do this...

        :param start: first bin to model (bin number)
        :param end: first bin to model (bin number)

        :returns: z-score, raw values and zzeros of the experiment
        """
        if not self._normalization.startswith('visibility'):
            stderr.write('WARNING: normalizing according to visibility method\n')
            self.normalize_hic()
        from pytadbit import Chromosome
        if start < 1:
            raise ValueError('ERROR: start should be higher than 0\n')
        start -= 1 # things starts at 0 for python. we keep the end coordinate
                   # at its original value because it is inclusive
        siz = self.size
        try:
            matrix = self.get_hic_matrix()
            new_matrix = [[matrix[i][j] for i in xrange(start, end)]
                          for j in xrange(start, end)]
            tmp = Chromosome('tmp')
            tmp.add_experiment('exp1', hic_data=[new_matrix],
                               resolution=self.resolution, filter_columns=False)
            exp = tmp.experiments[0]
            # We want the weights and zeros calculated in the full chromosome
            exp.norm = [[self.norm[0][i + siz * j] for i in xrange(start, end)
                         for j in xrange(start, end)]]
        except TypeError: # no Hi-C data provided
            matrix = self.get_hic_matrix(normalized=True)
            new_matrix = [[matrix[i][j] for i in xrange(start, end)]
                           for j in xrange(start, end)]
            tmp = Chromosome('tmp')
            tmp.add_experiment('exp1', norm_data=[new_matrix],
                               resolution=self.resolution, filter_columns=False)
            exp = tmp.experiments[0]
        exp._zeros = dict([(z - start, None) for z in self._zeros
                           if start <= z <= end - 1])
        if len(exp._zeros) == (end - start):
            raise Exception('ERROR: no interaction found in selected regions')
        # ... but the z-scores in this particular region
        exp.get_hic_zscores()
        values = [[float('nan') for _ in xrange(exp.size)]
                  for _ in xrange(exp.size)]
        for i in xrange(exp.size):
            # zeros are rows or columns having a zero in the diagonal
            if i in exp._zeros:
                continue
            for j in xrange(i + 1, exp.size):
                if j in exp._zeros:
                    continue
                if (not exp.norm[0][i * exp.size + j] 
                    or not exp.norm[0][i * exp.size + j]):
                    continue
                values[i][j] = exp.norm[0][i * exp.size + j]
                values[j][i] = exp.norm[0][i * exp.size + j]
        return exp._zscores, values, exp._zeros


    def write_interaction_pairs(self, fname, normalized=True, zscored=True,
                                diagonal=False, cutoff=None, header=False,
                                true_position=False, uniq=True,
                                remove_zeros=False, focus=None, format='tsv'):
        """
        Creates a tab separated file with all the pairwise interactions.
        
        :param fname: file name where to write the  pairwise interactions 
        :param True zscored: computes the z-score of the log10(data)
        :param True normalized: use the weights to normalize the data
        :param None cutoff: if defined, only the zscores above the cutoff will
           be writen to the file
        :param False uniq: only writes one representent per interacting pair
        :param False true_position: if, true writes genomic coordinates,
           otherwise, genomic bin.
        :param None focus: writes interactions between the start and stop bin
           passed to this parameter.
        :param 'tsv' format: in which to write the file, can be tab separated
           (tsv) or JSON (json)
           
        """
        if not self._zscores and zscored:
            self.get_hic_zscores()
        if not self.norm and normalized:
            raise Exception('Experiment not normalized.')
        # write to file
        if isinstance(fname, str):
            out = open(fname, 'w')
        elif isinstance(fname, file):
            out = fname
        else:
            raise Exception('Not recognize file type\n')
        if header:
            if format == 'tsv':
                out.write('elt1\telt2\t%s\n' % ('zscore' if zscored else 
                                                'normalized hi-c' if normalized 
                                                else 'raw hi-c'))
            elif format == 'json':
                out.write('''
{
    "metadata": {
                        "formatVersion" : 3,
                        %s
                        "species" : "%s",
                        "cellType" : "%s",
                        "experimentType" : "%s",
                        "identifier" : "%s",
                        "resolution" : %s,
                        "chromosome" : "%s",
                        "start" : %s,
                        "end" : %s
                },
    "interactions": [
                ''' % ('\n'.join(['"%s": "%s",' % (k, self.description[k])
                                  for k in self.description]),
                       self.description.get('species', ''),
                       self.cell_type,
                       self.exp_type,
                       self.identifier,
                       self.resolution,
                       self.crm.name,
                       focus[0] * self.resolution if focus else 1,
                       (focus[1] * self.resolution if focus else
                        self.resolution * self.size)))
        if focus:
            start, end = focus[0], focus[1] + 1
        else:
            start, end = 0, self.size
        for i in xrange(start, end):
            if i in self._zeros:
                continue
            newstart = i if uniq else 0
            for j in xrange(newstart, end):
                if j in self._zeros:
                    continue
                if not diagonal and i == j:
                    continue
                if zscored:
                    try:
                        if self._zscores[str(i)][str(j)] < cutoff:
                            continue
                        if self._zscores[str(i)][str(j)] == -99:
                            continue
                    except KeyError:
                        continue
                    val = self._zscores[str(i)][str(j)]
                elif normalized:
                    val = self.norm[0][self.size*i+j]
                else:
                    val = self.hic_data[0][self.size*i+j]
                if remove_zeros and not val:
                    continue
                if true_position:
                    if format == 'tsv':
                        out.write('%s\t%s\t%s\n' % (
                            self.resolution * (i + 1),
                            self.resolution * (j + 1), val))
                    elif format == 'json':
                        out.write('%s,%s,%s,\n' % (
                            self.resolution * (i + 1),
                            self.resolution * (j + 1), val))
                else:
                    if format == 'tsv':
                        out.write('%s\t%s\t%s\n' % (
                            i + 1 - start, j + 1 - start, val))
                    elif format == 'json':
                        out.write('%s,%s,%s\n' % (
                            i + 1 - start, j + 1 - start, val))
        if format == 'json':
            out.write(']}\n')
        out.close()


    def get_hic_matrix(self, focus=None, diagonal=True, normalized=False):
        """
        Return the Hi-C matrix.

        :param None focus: if a tuple is passed (start, end), wil return a Hi-C
           matrix starting at start, and ending at end (all inclusive).
        :param True diagonal: replace the values in the diagonal by one. Used
           for the filtering in order to smooth the distribution of mean values
        :para False normalized: returns normalized data instead of raw Hi-C

        :returns: list of lists representing the Hi-C data matrix of the
           current experiment
        """
        siz = self.size
        if normalized:
            try:
                hic = self.norm[0]
            except TypeError:
                raise Exception('ERROR: experiment not normalized yet')
        else:
            hic = self.hic_data[0]
        if focus:
            start, end = focus
            start -= 1
        else:
            start = 0
            end   = siz
        if diagonal:
            return [[hic[i + siz * j] for i in xrange(start, end)]
                    for j in xrange(start, end)]
        else:
            mtrx = [[hic[i + siz * j] for i in xrange(start, end)]
                    for j in xrange(start, end)]
            for i in xrange(start, end):
                mtrx[i][i] = 1 if mtrx[i][i] else 0
            return mtrx
            

    def print_hic_matrix(self, print_it=True, normalized=False, zeros=False):
        """
        Return the Hi-C matrix as string

        :param True print_it: Otherwise, returns the string
        :param False normalized: returns normalized data, instead of raw Hi-C
        :param False zeros: take into account filtered columns
        :returns: list of lists representing the Hi-C data matrix of the
           current experiment
        """
        siz = self.size
        if normalized:
            hic = self.norm[0]
        else:
            hic = self.hic_data[0]
        if zeros:
            out = '\n'.join(['\t'.join(
                ['nan' if (i in self._zeros or j in self._zeros) else
                 str(hic[i+siz * j]) for i in xrange(siz)])
                             for j in xrange(siz)])
        else:
            out = '\n'.join(['\t'.join([str(hic[i+siz * j]) 
                                        for i in xrange(siz)])
                             for j in xrange(siz)])
        if print_it:
            print out
        else:
            return out + '\n'


    def view(self, tad=None, focus=None, paint_tads=False, axe=None,
             show=True, logarithm=True, normalized=False, relative=True,
             decorate=True, savefig=None, where='both', clim=None,
             cmap='jet'):
        """
        Visualize the matrix of Hi-C interactions

        :param None tad: a given TAD in the form:
           ::
           
             {'start': start,
              'end'  : end,
              'brk'  : end,
              'score': score}
              
           **Alternatively** a list of the TADs can be passed (all the TADs
           between the first and last one passed will be showed. Thus, passing
           more than two TADs might be superfluous)
        :param None focus: a tuple with the start and end positions of the 
           region to visualize
        :param False paint_tads: draw a box around the TADs defined for this
           experiment
        :param None axe: an axe object from matplotlib can be passed in order
           to customize the picture
        :param True show: either to pop-up matplotlib image or not
        :param True logarithm: show the logarithm values
        :param True normalized: show the normalized data (weights might have
           been calculated previously). *Note: white rows/columns may appear in
           the matrix displayed; these rows correspond to filtered rows (see*
           :func:`pytadbit.utils.hic_filtering.hic_filtering_for_modelling` *)*
        :param True relative: color scale is relative to the whole matrix of
           data, not only to the region displayed
        :param True decorate: draws color bar, title and axes labels
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param None clim: tuple with minimum and maximum value range for color
           scale. I.e. clim=(-4, 10)
        :param 'jet' cmap: color map from matplotlib. Can also be a 
           preconfigured cmap object.
        """
        if logarithm==True:
            fun = log2
        elif logarithm:
            fun = logarithm
        else:
            fun = lambda x: x
        size = self.size
        if normalized and not self.norm:
            raise Exception('ERROR: weights not calculated for this ' +
                            'experiment. Run Experiment.normalize_hic\n')
        if tad and focus:
            raise Exception('ERROR: only one of "tad" or "focus" might be set')
        start = end = None
        if focus:
            start, end = focus
            if start == 0:
                stderr.write('WARNING: Hi-C matrix starts at 1, setting ' +
                             'starting point to 1.\n')
                start = 1
        elif isinstance(tad, dict):
            start = int(tad['start'])
            end   = int(tad['end'])
        elif isinstance(tad, list):
            if isinstance(tad[0], dict):
                start = int(sorted(tad,
                                   key=lambda x: int(x['start']))[0 ]['start'])
                end   = int(sorted(tad,
                                   key=lambda x: int(x['end'  ]))[-1]['end'  ])
        elif self.tads:
            start = self.tads[min(self.tads)]['start'] + 1
            end   = self.tads[max(self.tads)]['end'  ] + 1
        else:
            start =  1
            end   = size
        if normalized:
            norm_data = self.norm
        else:
            hic_data = self.hic_data
        if relative and not clim:
            if normalized:
                # find minimum, if value is non-zero... for logarithm
                mini = min([i for i in norm_data[0].values() if i])
                if mini == int(mini):
                    vmin = min(norm_data[0].values())
                else:
                    vmin = mini
                vmin = fun(vmin or (1 if logarithm else 0))
                vmax = fun(max(norm_data[0].values()))
            else:
                vmin = fun(min(hic_data[0].values()) or
                           (1 if logarithm else 0))
                vmax = fun(max(hic_data[0].values()))
        elif clim:
            vmin, vmax = clim
        if axe is None:
            plt.figure(figsize=(8, 6))
            axe = plt.subplot(111)
        if tad or focus:
            if start > -1:
                if normalized:
                    matrix = [
                        [norm_data[0][i+size*j]
                         if (not i in self._zeros
                             and not j in self._zeros) else float('nan')
                         for i in xrange(int(start) - 1, int(end))]
                        for j in xrange(int(start) - 1, int(end))]
                else:
                    matrix = [
                        [hic_data[0][i+size*j]
                         for i in xrange(int(start) - 1, int(end))]
                        for j in xrange(int(start) - 1, int(end))]
            elif isinstance(tad, list):
                if normalized:
                    stderr.write('WARNING: List passed, not going to be ' +
                                 'normalized.\n')
                matrix = tad
            else:
                # TODO: something... matrix not declared...
                pass
        else:
            if normalized:
                matrix = [[norm_data[0][i+size*j]
                           if (not i in self._zeros
                               and not j in self._zeros) else float('nan')
                           for i in xrange(size)]
                          for j in xrange(size)]
            else:
                matrix = [[hic_data[0][i+size*j]\
                           for i in xrange(size)] \
                          for j in xrange(size)]
        if where == 'up':
            for i in xrange(int(end - start)):
                for j in xrange(i, int(end - start)):
                    matrix[i][j] = vmin
            alphas = array([0, 0] + [1] * 256 + [0])
            jet._init()
            jet._lut[:, -1] = alphas
        elif where == 'down':
            for i in xrange(int(end - start)):
                for j in xrange(i + 1):
                    matrix[i][j] = vmin
            alphas = array([0, 0] + [1] * 256 + [0])
            jet._init()
            jet._lut[:,-1] = alphas

        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
            cmap.set_bad('darkgrey', 1)
        if relative:
            img = axe.imshow(nozero_log_matrix(matrix, fun), origin='lower', vmin=vmin, vmax=vmax,
                             interpolation="nearest", cmap=cmap,
                             extent=(int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5,
                                     int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5))
        else:
            img = axe.imshow(nozero_log_matrix(matrix, fun), origin='lower',
                             interpolation="nearest", cmap=cmap,
                             extent=(int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5,
                                     int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5))
        if decorate:
            cbar = axe.figure.colorbar(img)
            cbar.ax.set_ylabel('%sHi-C %sinteraction count' % (
                'Log2 ' * (logarithm==True), 'normalized ' * normalized), rotation=-90)
            axe.set_title(('Chromosome %s experiment %s %s') % (
                self.crm.name, self.name,
                'focus: %s-%s' % (start, end) if tad else ''))
            axe.set_xlabel('Genomic bin (resolution: %s)' % (self.resolution))
            if paint_tads:
                axe.set_ylabel('TAD number')
            else:
                axe.set_ylabel('Genomic bin (resolution: %s)' % (
                    self.resolution))
        if not paint_tads:            
            axe.set_ylim(int(start or 1) - 0.5,
                         int(start or 1) + len(matrix) - 0.5)
            axe.set_xlim(int(start or 1) - 0.5,
                         int(start or 1) + len(matrix) - 0.5)
            if show:
                plt.show()
            return img
        pwidth = 1
        tads = dict([(t, self.tads[t]) for t in self.tads
                     if  ((int(self.tads[t]['start']) + 1 >= start
                           and int(self.tads[t]['end'  ]) + 1 <= end)
                          or not start)])
        for i, tad in tads.iteritems():
            t_start = int(tad['start']) + .5
            t_end   = int(tad['end'])   + 1.5
            nwidth = float(abs(tad['score'])) / 4
            if where in ['down', 'both']:
                axe.hlines(t_start, t_start, t_end, colors='k', lw=pwidth)
            if where in ['up', 'both']:
                axe.hlines(t_end  , t_start, t_end, colors='k', lw=nwidth)
            if where in ['up', 'both']:
                axe.vlines(t_start, t_start, t_end, colors='k', lw=pwidth)
            if where in ['down', 'both']:
                axe.vlines(t_end  , t_start, t_end, colors='k', lw=nwidth)
            pwidth = nwidth
            if tad['score'] < 0:
                for j in xrange(0, int(t_end) - int(t_start), 2):
                    axe.plot((t_start    , t_start + j),
                             (t_end   - j, t_end      ), color='k')
                    axe.plot((t_end      , t_end   - j),
                             (t_start + j, t_start    ), color='k')
        axe.set_ylim(int(start or 1) - 0.5,
                     int(start or 1) + len(matrix) - 0.5)
        axe.set_xlim(int(start or 1) - 0.5,
                     int(start or 1) + len(matrix) - 0.5)
        if paint_tads:
            ticks = []
            labels = []
            for tad, tick in [(t, tads[t]['start'] + (tads[t]['end'] -
                                                      tads[t]['start'] - 1))
                              for t in tads.keys()[::(len(tads)/11 + 1)]]:
                ticks.append(tick)
                labels.append(tad + 1)
            axe.set_yticks(ticks)
            axe.set_yticklabels(labels)
        if show:
            plt.show()
        if savefig:
            tadbit_savefig(savefig)
        return img


    def write_tad_borders(self, density=False, savedata=None, normalized=False):
        """
        Print a table summarizing the TADs found by tadbit. This function outputs
        something similar to the R function.

        :param False density: if True, adds a column with the relative
           interaction frequency measured within each TAD (value of 1 means an
           interaction frequency equal to the expectation in the experiment)
        :param None savedata: path to a file where to save the density data
           generated (1 column per step + 1 for particle number). If None, print
           a table.
        :param False normalized: uses normalized data to calculate the density
        """
        if normalized and self.norm:
            norms = self.norm[0]
        elif self.hic_data:
            if normalized:
                warn("WARNING: weights not available, using raw data")
            norms = self.hic_data[0]
        else:
            warn("WARNING: raw Hi-C data not available, " +
                 "TAD's height fixed to 1")
            norms = None
        zeros = self._zeros or {}
        table = ''
        table += '%s\t%s\t%s\t%s%s\n' % ('#', 'start', 'end', 'score',
                                        '' if not density else '\tdensity')
        tads = self.tads
        sp1 = self.size + 1
        diags = []
        if norms:
            for k in xrange(1, self.size):
                s_k = self.size * k
                diags.append(sum([norms[i * sp1 + s_k]
                                 if not (i in zeros
                                         or (i + k) in zeros) else 0.
                                  for i in xrange(
                                      self.size - k)]) / (self.size - k))
        for tad in tads:
            table += '%s\t%s\t%s\t%s%s\n' % (
                tad, int(tads[tad]['start'] + 1), int(tads[tad]['end'] + 1),
                abs(tads[tad]['score']), '' if not density else
                '\t%s' % (round(float(tads[tad]['height']), 3)))
        if not savedata:
            print table
            return
        if isinstance(savedata, file):
            out = savedata
        else:
            out = open(savedata, 'w')
        out.write(table)


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

