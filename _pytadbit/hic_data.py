"""
December 12, 2014.

"""

from pytadbit.utils.extraviews      import plot_compartments
from pytadbit.utils.extraviews      import plot_compartments_summary
from pytadbit.utils.hic_filtering   import filter_by_mean, filter_by_zero_count
from pytadbit.utils.normalize_hic   import iterative, expected
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.parsers.bed_parser    import parse_bed
from pytadbit.utils.file_handling   import mkdir
from numpy.linalg                   import LinAlgError
from numpy                          import corrcoef, nansum, array, isnan
from numpy                          import nanpercentile as npperc
from scipy.cluster.hierarchy        import linkage, fcluster
from scipy.sparse.linalg            import eigsh
from pytadbit.utils.tadmaths        import calinski_harabasz
from scipy.stats                    import ttest_ind
from collections                    import OrderedDict
from warnings                       import warn
from bisect                         import bisect_right as bisect
import os

class HiC_data(dict):
    """
    This may also hold the print/write-to-file matrix functions
    """
    def __init__(self, items, size, chromosomes=None, dict_sec=None,
                 resolution=1, masked=None):
        super(HiC_data, self).__init__(items)
        self.__size = size
        self._size2 = size**2
        self.bias = None
        self.bads = masked or {}
        self.chromosomes = chromosomes
        self.sections = dict_sec
        self.section_pos = {}
        self.resolution = resolution
        self.expected = None
        self.compartments = {}
        if self.chromosomes:
            total = 0
            for crm in self.chromosomes:
                self.section_pos[crm] = (total, total + self.chromosomes[crm])
                total += self.chromosomes[crm]
        if self.sections == {}:
            self.section_pos = {None: (0, self.__size)}
            self.sections = dict([((None, i), i)
                                  for i in xrange(0, self.__size)])

    def _update_size(self, size):
        self.__size +=  size
        self._size2 = self.__size**2
        
    def __len__(self):
        return self.__size

    def __getitem__(self, row_col):
        """
        slow one... for user
        for fast item getting, use self.get()
        """
        try:
            row, col = row_col
            pos = row * self.__size + col
            if pos > self._size2:
                raise IndexError(
                    'ERROR: row or column larger than %s' % self.__size)
            return self.get(pos, 0)
        except TypeError:
            if row_col > self._size2:
                raise IndexError(
                    'ERROR: position %d larger than %s^2' % (row_col,
                                                             self.__size))
            return self.get(row_col, 0)

    def __setitem__(self, row_col, val):
        """
        slow one... for user
        for fast item getting, use self.get()
        """
        try:
            row, col = row_col
            pos = row * self.__size + col
            if pos > self._size2:
                raise IndexError(
                    'ERROR: row or column larger than %s' % self.__size)
            super(HiC_data, self).__setitem__(pos, val)
        except TypeError:
            if row_col > self._size2:
                raise IndexError(
                    'ERROR: position %d larger than %s^2' % (row_col,
                                                             self.__size))
            super(HiC_data, self).__setitem__(row_col, val)

    def add_sections_from_fasta(self, fasta):
        """
        Add genomic coordinate to HiC_data object by getting them from a fasta
        file containing chromosome sequences

        :param fasta: path to a fasta file
        """
        genome = parse_fasta(fasta, verbose=False)
        sections = []
        genome_seq = OrderedDict()
        size = 0
        for crm in  genome:
            genome_seq[crm] = int(len(genome[crm])) / self.resolution + 1
            size += genome_seq[crm]
        section_sizes = {}
        for crm in genome_seq:
            len_crm = genome_seq[crm]
            section_sizes[(crm,)] = len_crm
            sections.extend([(crm, i) for i in xrange(len_crm)])
        dict_sec = dict([(j, i) for i, j in enumerate(sections)])
        self.chromosomes = genome_seq
        self.sections = dict_sec
        if self.chromosomes:
            total = 0
            for crm in self.chromosomes:
                self.section_pos[crm] = (total, total + self.chromosomes[crm])
                total += self.chromosomes[crm]
        if size != self.__size:
            warn('WARNING: different sizes (%d, now:%d), ' % (self.__size, size)
                 + 'should adjust the resolution')
        self.__size = size
        self._size2 = size**2

    def add_sections(self, lengths, chr_names=None, binned=False):
        """
        Add genomic coordinate to HiC_data object by getting them from a fasta
        file containing chromosome sequences. Orders matters.

        :param lengths: list of chromosome lengths
        :param None chr_names: list of corresponding chromosome names.
        :param False binned: if True, leghths will not be divided by resolution
        """
        sections = []
        genome_seq = OrderedDict()
        size = 0
        resolution = 1 if binned else self.resolution
        for crm, length in  enumerate(lengths):
            cnam = 'chr' + str(crm) if not chr_names else chr_names[crm]
            genome_seq[cnam] = int(length) / resolution + 1
            size += genome_seq[cnam]
        section_sizes = {}
        for crm in genome_seq:
            len_crm = genome_seq[crm]
            section_sizes[(crm,)] = len_crm
            sections.extend([(crm, i) for i in xrange(len_crm)])
        dict_sec = dict([(j, i) for i, j in enumerate(sections)])
        self.chromosomes = genome_seq
        self.sections = dict_sec
        if self.chromosomes:
            total = 0
            for crm in self.chromosomes:
                self.section_pos[crm] = (total, total + self.chromosomes[crm])
                total += self.chromosomes[crm]
        if size != self.__size:
            warn('WARNING: different sizes (%d, now:%d), ' % (self.__size, size)
                 + 'should adjust the resolution')
        self.__size = size
        self._size2 = size**2

    def cis_trans_ratio(self, normalized=False, exclude=None, diagonal=True,
                        equals=None):
        """
        Counts the number of interactions occuring within chromosomes (cis) with
        respect to the total number of interactions

        :param False normalized: used normalized data
        :param None exclude: exclude a given list of chromosome from the
           ratio (may want to exclude translocated chromosomes)
        :param False diagonal: replace values in the diagonal by 0 or 1
        :param None equals: can pass a function that would decide if 2 chromosomes
           have to be considered as the same. e.g. lambda x, y: x[:4]==y[:4] will
           consider chr2L and chr2R as being the same chromosome. WARNING: only
           working on consecutive chromosomes.

        :returns: the ratio of cis interactions over the total number of
           interactions. This number is expected to be between at least 40-60%
           in Human classic dilution Hi-C with HindIII as restriction enzyme.
        """
        if normalized and not self.bias:
            raise Exception('ERROR: experiment not normalized yet')
        if exclude == None:
            exclude = []
        if equals == None:
            equals = lambda x, y: x == y
        intra = 0
        if not self.chromosomes:
            return float('nan')
        # define chromosomes to be merged
        to_skip = set()
        c_prev = ''
        for c in self.chromosomes:
            if equals(c, c_prev):
                to_skip.add(c_prev)
            c_prev = c
        sections = sorted([-1] + [self.section_pos[c][1]
                                  for c in self.section_pos
                                  if not c in to_skip])
        # defines columns to be skipped
        bads = set(self.bads.keys())
        for c in exclude:
            bads.update(i for i in xrange(*self.section_pos[c]))
        # diagonal
        if diagonal:
            valid = lambda x, y: True
        else:
            valid = lambda x, y: x != y
        # normalization
        if normalized:
            transform = lambda x, y, z: x / self.bias[y] / self.bias[z]
        else:
            transform = lambda x, y, z: x
        # compute ratio
        for k, v in self.iteritems():
            i, j = divmod(k, self.__size)
            if bisect(sections, i) != bisect(sections, j):
                continue
            if i in bads or j in bads:
                continue
            if valid(i, j): # diagonal thing
                intra += transform(v, i, j)
        return float(intra) / self.sum(bias=self.bias if normalized else None, bads=bads)

    def filter_columns(self, draw_hist=False, savefig=None, perc_zero=75,
                       by_mean=True, silent=False):
        """
        Call filtering function, to remove artefactual columns in a given Hi-C
        matrix. This function will detect columns with very low interaction
        counts; columns passing through a cell with no interaction in the
        diagonal; and columns with NaN values (in this case NaN will be replaced
        by zero in the original Hi-C data matrix). Filtered out columns will be
        stored in the dictionary Experiment._zeros.

        :param False draw_hist: shows the distribution of mean values by column
           the polynomial fit, and the cut applied.
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param 75 perc_zero: maximum percentage of cells with no interactions
           allowed.
        :param True by_mean: filter columns by mean column value using
           :func:`pytadbit.utils.hic_filtering.filter_by_mean` function

        """
        self.bads = filter_by_zero_count(self, perc_zero, silent=silent)
        if by_mean:
            self.bads.update(filter_by_mean(
                self, draw_hist=draw_hist, silent=silent,
                savefig=savefig, bads=self.bads))
        if not silent:
            print 'Found %d of %d columnswith poor signal' % (len(self.bads),
                                                              len(self))

    def sum(self, bias=None, bads=None):
        """
        Sum Hi-C data matrix
        WARNING: parameters are not meant to be used by external users

        :params None bias: expects a dictionary of biases to use normalized matrix
        :params None bads: extends computed bad columns
        
        :returns: the sum of the Hi-C matrix skipping bad columns
        """
        N = self.__size
        norm_sum = 0
        bads = bads or self.bads
        if bias:
            for k, v in self.iteritems():
                i, j = divmod(k, N)
                if i in bads or j in bads:
                    continue
                norm_sum += v / (bias[i] * bias[j])
        else:
            for k, v in self.iteritems():
                i, j = divmod(k, N)
                if i in bads or j in bads:
                    continue
                norm_sum += v
        return norm_sum

    def normalize_hic(self, iterations=0, max_dev=0.1, silent=False, factor=1):
        """
        Normalize the Hi-C data.

        It fills the Experiment.norm variable with the Hi-C values divided by
        the calculated weight.

        :param 0 iteration: number of iterations
        :param 0.1 max_dev: iterative process stops when the maximum deviation
           between the sum of row is equal to this number (0.1 means 10%)
        :param False silent: does not warn when overwriting weights
        :param 1 factor: final mean number of normalized interactions wanted
           per cell (excludes filtered, or bad, out columns)
        """
        bias = iterative(self, iterations=iterations,
                         max_dev=max_dev, bads=self.bads,
                         verbose=not silent)
        if factor:
            if not silent:
                print 'rescaling to factor %d' % factor
            # get the sum on half of the matrix
            norm_sum = self.sum(bias)
            # divide biases
            target = (norm_sum / float(len(self) * len(self) * factor))**0.5
            bias = dict([(b, bias[b] * target) for b in bias])
        self.bias = bias

    def get_as_tuple(self):
        return tuple([self[i, j]
                      for j in xrange(len(self))
                      for i in xrange(len(self))])

    def write_matrix(self, fname, focus=None, diagonal=True, normalized=False):
        """
        writes the matrix to a file.
        
        :param None focus: a tuple with the (start, end) position of the desired
           window of data (start, starting at 1, and both start and end are
           inclusive). Alternatively a chromosome name can be input or a tuple
           of chromosome name, in order to retrieve a specific inter-chromosomal
           region
        :param True diagonal: if False, diagonal is replaced by zeroes
        :param False normalized: get normalized data
        """
        if focus:
            if isinstance(focus, tuple) and isinstance(focus[0], int):
                if len(focus) == 2:
                    start1, end1 = focus
                    start2, end2 = focus
                    start1 -= 1
                    start2 -= 1
                else:
                    start1, end1, start2, end2 = focus
                    start1 -= 1
                    start2 -= 1
            elif isinstance(focus, tuple) and isinstance(focus[0], str):
                start1, end1 = self.section_pos[focus[0]]
                start2, end2 = self.section_pos[focus[1]]
            else:
                start1, end1 = self.section_pos[focus]
                start2, end2 = self.section_pos[focus]
        else:
            start1 = start2 = 0
            end1   = end2   = len(self)
        out = open(fname, 'w')
        out.write('# MASKED %s\n' % (' '.join([str(k - start1)
                                               for k in self.bads.keys()
                                               if start1 <= k <= end1])))
        rownam = ['%s\t%d-%d' % (k[0],
                                 k[1] * self.resolution,
                                 (k[1] + 1) * self.resolution)
                  for k in sorted(self.sections,
                                  key=lambda x: self.sections[x])
                  if start2 <= self.sections[k] < end2]
        if rownam:
            for line in self.yield_matrix(focus=focus, diagonal=diagonal,
                                          normalized=normalized):
                out.write(rownam.pop(0) + '\t' +
                          '\t'.join([str(i) for i in line]) + '\n')
        else:
            for line in self.yield_matrix(focus=focus, diagonal=diagonal,
                                          normalized=normalized):
                out.write('\t'.join([str(i) for i in line]) + '\n')
        out.close()

    def get_matrix(self, focus=None, diagonal=True, normalized=False):
        """
        returns a matrix.
        
        :param None focus: a tuple with the (start, end) position of the desired
           window of data (start, starting at 1, and both start and end are
           inclusive). Alternatively a chromosome name can be input or a tuple
           of chromosome name, in order to retrieve a specific inter-chromosomal
           region
        :param True diagonal: if False, diagonal is replaced by ones, or zeroes
           if normalized
        :param False normalized: get normalized data

        :returns: matrix (a list of lists of values)
        """
        siz = len(self)
        if normalized and not self.bias:
            raise Exception('ERROR: experiment not normalized yet')
        if focus:
            if isinstance(focus, tuple) and isinstance(focus[0], int):
                if len(focus) == 2:
                    start1, end1 = focus
                    start2, end2 = focus
                    start1 -= 1
                    start2 -= 1
                else:
                    start1, end1, start2, end2 = focus
                    start1 -= 1
                    start2 -= 1
            elif isinstance(focus, tuple) and isinstance(focus[0], str):
                start1, end1 = self.section_pos[focus[0]]
                start2, end2 = self.section_pos[focus[1]]
            else:
                start1, end1 = self.section_pos[focus]
                start2, end2 = self.section_pos[focus]
        else:
            start1 = start2 = 0
            end1   = end2   = siz
        if normalized:
            if diagonal:
                return [[self[i, j] / self.bias[i] / self.bias[j]
                         for i in xrange(start2, end2)]
                        for j in xrange(start1, end1)]
            else:
                mtrx = [[self[i, j] / self.bias[i] / self.bias[j]
                         for i in xrange(start2, end2)]
                        for j in xrange(start1, end1)]
                if start1 == start2:
                    for i in xrange(len(mtrx)):
                        mtrx[i][i] = 0
                return mtrx
        else:
            if diagonal:
                return [[self[i, j] for i in xrange(start2, end2)]
                        for j in xrange(start1, end1)]
            else:
                mtrx = [[self[i, j] for i in xrange(start2, end2)]
                        for j in xrange(start1, end1)]
                if start1 == start2:
                    for i in xrange(len(mtrx)):
                        mtrx[i][i] = 1 if mtrx[i][i] else 0
                return mtrx

    def find_compartments(self, crms=None, savefig=None, savedata=None,
                          savecorr=None, show=False, suffix='',
                          label_compartments=True, log=None, max_mean_size=10000,
                          ev_index=None, rich_in_A=None, **kwargs):
        """
        Search for A/B copartments in each chromsome of the Hi-C matrix.
        Hi-C matrix is normalized by the number interaction expected at a given
        distance, and by visibility (one iteration of ICE). A correlation matrix
        is then calculated from this normalized matrix, and its first
        eigenvector is used to identify compartments. Changes in sign marking
        boundaries between compartments.
        Result is stored as a dictionary of compartment boundaries, keys being
        chromsome names.
        
        :param 99 perc_zero: to filter bad columns
        :param 0.05 signal_to_noise: to calculate expected interaction counts,
           if not enough reads are observed at a given distance the observations
           of the distance+1 are summed. a signal to noise ratio of < 0.05
           corresponds to > 400 reads.
        :param None crms: only runs these given list of chromosomes
        :param None savefig: path to a directory to store matrices with
           compartment predictions, one image per chromosome, stored under
           'chromosome-name.png'.
        :param False show: show the plot
        :param None savedata: path to a new file to store compartment
           predictions, one file only.
        :param None savecorr: path to a directory where to save correlation
           matrices of each chromosome
        :param -1 vmin: for the color scale of the plotted map (use vmin='auto',
           and vmax='auto' to color according to the absolute maximum found).
        :param 1 vmax: for the color scale of the plotted map (use vmin='auto',
           and vmax='auto' to color according to the absolute maximum found).
        :param False yield_ev1: if True yields one list per chromosome with the
           first eigenvector used to compute compartments.
        :param '' suffix: to be placed after file names of compartment images
        :param None ev_index: a list of number refering to the index of the
           eigenvector to be used. By default the first eigenvector is used.
           WARNING: index starts at 1, default is thus a list of ones.
        :param None rich_in_A: by default compartments are identified using mean
           number of intra-interactions (A compartments are expected to have
           less). However this measure is not very accurate. Using this
           parameter a path to a BED or BED-Graph file with a list of genes or
           active epigenetic marks can be passed, and used instead of the mean
           interactions.
        :param None log: path to a folder where to save log of the assignment
           of A/B compartments
        :param True label_compartments: label compartments into A/B categories,
           otherwise just find borders (faster).

        TODO: this is really slow...

        Notes: building the distance matrix using the amount of interactions
               instead of the mean correlation, gives generally worse results.

        :returns: a dictionary with the two first eigen vectors used to define
           compartment borders for each chromosome (keys are chromosome names)
        """
        if not self.bads:
            if kwargs.get('verbose', False):
                print 'Filtering bad columns %d' % 99
            self.filter_columns(perc_zero=kwargs.get('perc_zero', 99),
                                by_mean=False, silent=True)
            if len(self.bads) == len(self):
                self.bads = {}
                warn('WARNING: all columns would have been filtered out, '
                     'filtering disabled')
        if not self.expected:
            if kwargs.get('verbose', False):
                print 'Normalizing by expected values'
            self.expected = expected(self, bads=self.bads, **kwargs)
        if not self.bias:
            if kwargs.get('verbose', False):
                print 'Normalizing by ICE (1 round)'
            self.normalize_hic(iterations=0,
                               silent=not kwargs.get('verbose', False))
        if savefig:
            mkdir(savefig)
        if savecorr:
            mkdir(savecorr)
        if suffix != '':
            suffix = '_' + suffix
        # parse bed file
        if rich_in_A:
            rich_in_A = parse_bed(rich_in_A, resolution=self.resolution)
        cmprts = {}
        firsts = {}
        ev_nums = {}
        count = 0
        for sec in self.section_pos:
            if crms and sec not in crms:
                continue
            if kwargs.get('verbose', False):
                print 'Processing chromosome', sec
                warn('Processing chromosome %s' % (sec))
            matrix = [[(float(self[i,j]) / self.expected[abs(j-i)]
                       / self.bias[i] / self.bias[j])
                      for i in xrange(*self.section_pos[sec])
                       if not i in self.bads]
                     for j in xrange(*self.section_pos[sec])
                      if not j in self.bads]
            if not matrix: # MT chromosome will fall there
                warn('Chromosome %s is probably MT :)' % (sec))
                cmprts[sec] = []
                count += 1
                continue
            for i in xrange(len(matrix)):
                for j in xrange(i+1, len(matrix)):
                    matrix[i][j] = matrix[j][i]
            try:
                matrix = [list(m) for m in corrcoef(matrix)]
            except TypeError:
                # very small chromosome?
                warn('Chromosome %s is probably MT :)' % (sec))
                cmprts[sec] = []
                count += 1
                continue
            # write correlation matrix to file. replaces filtered row/columns by NaN
            if savecorr:
                out = open(os.path.join(savecorr, '%s_corr-matrix.tsv' % (sec)),
                           'w')
                start1, end1 = self.section_pos[sec]
                out.write('# MASKED %s\n' % (' '.join([str(k - start1)
                                                       for k in self.bads.keys()
                                                       if start1 <= k <= end1])))
                rownam = ['%s\t%d-%d' % (k[0],
                                         k[1] * self.resolution,
                                         (k[1] + 1) * self.resolution)
                          for k in sorted(self.sections,
                                          key=lambda x: self.sections[x])
                          if k[0] == sec]
                length = self.section_pos[sec][1] - self.section_pos[sec][0]
                empty = 'NaN\t' * (length - 1) + 'NaN\n'
                badrows = 0
                for row in xrange(length):
                    if row in self.bads:
                        out.write(rownam.pop(0) + '\t' +empty)
                        badrows += 1
                        continue
                    vals = []
                    badcols = 0
                    for col in xrange(length):
                        if col in self.bads:
                            vals.append('NaN')
                            badcols += 1
                            continue
                        vals.append(str(matrix[row-badrows][col-badcols]))
                    out.write(rownam.pop(0) + '\t' +'\t'.join(vals) + '\n')
                out.close()
            try:
                # This eighs is very very fast, only ask for one eigvector
                _, evect = eigsh(array(matrix), k=ev_index[count] if ev_index else 2)
            except LinAlgError:
                warn('Chromosome %s too small to compute PC1' % (sec))
                cmprts[sec] = [] # Y chromosome, or so...
                count += 1
                continue
            index = ev_index[count] if ev_index else 1
            two_first = [evect[:, -1], evect[:, -2]]
            for ev_num in range(index, 3):
                first = list(evect[:, -ev_num])
                breaks = [i for i, (a, b) in
                          enumerate(zip(first[1:], first[:-1]))
                          if a * b < 0] + [len(first) - 1]
                breaks = [{'start': breaks[i-1] + 1 if i else 0, 'end': b}
                          for i, b in enumerate(breaks)]
                if self.resolution * float(len(breaks)) / len(matrix) > max_mean_size:
                    warn('WARNING: number of compartments found with the '
                         'EigenVector number %d is too low (%d compartments)'
                         % (ev_num, len(breaks)))
                else:
                    break
            ev_nums[sec] = ev_num
            beg, end = self.section_pos[sec]
            bads = [k - beg for k in self.bads if beg <= k <= end]
            _ = [first.insert(b, 0) for b in bads]
            _ = [matrix.insert(b, [float('nan')] * len(matrix[0]))
                 for b in bads]
            _ = [matrix[i].insert(b, float('nan'))
                 for b in bads for i in xrange(len(first))]
            breaks = [i for i, (a, b) in
                      enumerate(zip(first[1:], first[:-1]))
                      if a * b < 0] + [len(first) - 1]
            breaks = [{'start': breaks[i-1] + 1 if i else 0, 'end': b}
                      for i, b in enumerate(breaks)]
            cmprts[sec] = breaks
            
            firsts[sec] = two_first
            self.__apply_metric(cmprts, sec, rich_in_A)
            if label_compartments:
                if log:
                    logf = os.path.join(log, sec + suffix + '.log')
                else:
                    logf = None
                gammas = {}
                for gamma in range(101):
                    gammas[gamma] = _find_ab_compartments(float(gamma)/100, matrix,
                                                          breaks, cmprts[sec],
                                                          rich_in_A, ev_num=ev_num,
                                                          log=logf, save=False)
                    if kwargs.get('verbose', False):
                        print gamma, gammas[gamma]
                gamma = min(gammas.keys(), key=lambda k: gammas[k][0])
                if kwargs.get('verbose', False):
                    print '   ====>  minimum:', gamma
                _ = _find_ab_compartments(float(gamma)/100, matrix, breaks,
                                          cmprts[sec], rich_in_A, save=True,
                                          log=logf, ev_num=ev_num)
            if savefig or show:
                vmin = kwargs.get('vmin', -1)
                vmax = kwargs.get('vmax',  1)
                if vmin == 'auto' == vmax:
                    vmax = max([abs(npperc(matrix, 99.5)),
                                abs(npperc(matrix, 0.5))])
                    vmin = -vmax
                plot_compartments(
                    sec, first, cmprts, matrix, show,
                    savefig + '/chr' + sec + suffix + '.pdf' if savefig else None,
                    vmin=vmin, vmax=vmax, whichpc=ev_num)
                plot_compartments_summary(
                    sec, cmprts, show,
                    savefig + '/chr' + sec + suffix + '_summ.pdf' if savefig else None)
            count += 1
        self.compartments = cmprts
        if savedata:
            self.write_compartments(savedata, chroms=self.compartments.keys(),
                                    ev_nums=ev_nums)
        return firsts

    def __apply_metric(self, cmprts, sec, rich_in_A):
        """
        calculate compartment internal density if no rich_in_A, otherwise
        sum this list
        """
        for cmprt in cmprts[sec]:
            if rich_in_A:
                beg1, end1 = cmprt['start'], cmprt['end'] + 1
                sec_matrix = [rich_in_A.get(sec, {None: 0}).get(i, 0)
                              for i in xrange(beg1, end1)
                              if not i in self.bads]
                try:
                    cmprt['dens'] = float(sum(sec_matrix)) / len(sec_matrix)
                except ZeroDivisionError:
                    cmprt['dens'] = 0.
            else:
                beg = self.section_pos[sec][0]
                beg1, end1 = cmprt['start'] + beg, cmprt['end'] + beg + 1
                sec_matrix = [(self[i,j] / self.expected[abs(j-i)]
                               / self.bias[i] / self.bias[j])
                              for i in xrange(beg1, end1) if not i in self.bads
                              for j in xrange(i, end1) if not j in self.bads]
                try:
                    cmprt['dens'] = sum(sec_matrix) / len(sec_matrix)
                except ZeroDivisionError:
                    cmprt['dens'] = 0.
        try:
            meanh = sum([cmprt['dens'] for cmprt in cmprts[sec]]) / len(cmprts[sec])
        except ZeroDivisionError:
            meanh = 1.
        for cmprt in cmprts[sec]:
            try:
                cmprt['dens'] /= meanh
            except ZeroDivisionError:
                cmprt['dens'] = 1.


    def write_compartments(self, savedata, chroms=None, ev_nums=None):
        """
        Write compartments to a file.

        :param savedata: path to a file.
        :param None chroms: write only the given list of chromosomes (default
           all chromosomes are written, note that the first column corresponding
           to chromosome name will disappear in non default case)
        """
        out = open(savedata, 'w')
        sections = chroms if chroms else self.compartments.keys()
        if ev_nums:
            for sec in sections:
                out.write('## CHR %s\tEigenvector: %d\n' % (sec, ev_nums[sec]))
        out.write('#%sstart\tend\tdensity\ttype\n'% (
            'CHR\t' if len(sections) > 1 else ''))
        try:
            out.write('\n'.join(['\n'.join(['%s%d\t%d\t%.2f\t%s' % (
                (sec + '\t') if sections else '',
                c['start'] + 1, c['end'] + 1,
                c['dens'], c['type'])
                                            for c in self.compartments[sec]])
                                 for sec in sections]) + '\n')
        except KeyError:
            out.write('\n'.join(['\n'.join(['%s%d\t%d\t%.2f\t%s' % (
                (sec + '\t') if sections else '',
                c['start'], c['end'], c['dens'], '')
                                            for c in self.compartments[sec]])
                                 for sec in sections]) + '\n')
        out.close()
        

    def yield_matrix(self, focus=None, diagonal=True, normalized=False):
        """
        Yields a matrix line by line.
        Bad row/columns are returned as null row/columns.
        
        :param None focus: a tuple with the (start, end) position of the desired
           window of data (start, starting at 1, and both start and end are
           inclusive). Alternatively a chromosome name can be input or a tuple
           of chromosome name, in order to retrieve a specific inter-chromosomal
           region
        :param True diagonal: if False, diagonal is replaced by zeroes
        :param False normalized: get normalized data

        :yields: matrix line by line (a line being a list of values)
        """
        siz = len(self)
        if normalized and not self.bias:
            raise Exception('ERROR: experiment not normalized yet')
        if focus:
            if isinstance(focus, tuple) and isinstance(focus[0], int):
                if len(focus) == 2:
                    start1, end1 = focus
                    start2, end2 = focus
                    start1 -= 1
                    start2 -= 1
                else:
                    start1, end1, start2, end2 = focus
                    start1 -= 1
                    start2 -= 1
            elif isinstance(focus, tuple) and isinstance(focus[0], str):
                start1, end1 = self.section_pos[focus[0]]
                start2, end2 = self.section_pos[focus[1]]
            else:
                start1, end1 = self.section_pos[focus]
                start2, end2 = self.section_pos[focus]
        else:
            start1 = start2 = 0
            end1   = end2   = siz
        if normalized:
            for i in xrange(start2, end2):
                # if bad column:
                if i in self.bads:
                    yield [0.0 for j in xrange(start1, end1)]
                # if we want the diagonal, or we don't but are looking at a
                # region that is not symmetric
                elif diagonal or start1 != start2:
                    yield [self[i, j] / self.bias[i] / self.bias[j]
                           for j in xrange(start1, end1)]
                # diagonal replaced by zeroes
                else:
                    yield ([self[i, j] / self.bias[i] / self.bias[j]
                            for j in xrange(start1, i)] +
                           [0.0] + 
                           [self[i, j] / self.bias[i] / self.bias[j]
                            for j in xrange(i + 1, end1)])
        else:
            for i in xrange(start2, end2):
                # if bad column:
                if i in self.bads:
                    yield [0 for j in xrange(start1, end1)]
                # if we want the diagonal, or we don't but are looking at a
                # region that is not symmetric
                elif diagonal or start1 != start2:
                    yield [self[i, j] for j in xrange(start1, end1)]
                # diagonal replaced by zeroes
                else:
                    yield ([self[i, j] for j in xrange(start1, i)] +
                           [0] + 
                           [self[i, j] for j in xrange(i + 1, end1)])

def _find_ab_compartments(gamma, matrix, breaks, cmprtsec, rich_in_A, save=True,
                          ev_num=1, log=None, verbose=False):
    # function to convert correlation into distances

    gamma += 1
    func = lambda x: -abs(x)**gamma / x
    funczero = lambda x: 0.0
    # calculate distance_matrix
    dist_matrix = [[0 for _ in xrange(len(breaks))]
                   for _ in xrange(len(breaks))]
    scores = {}
    for k, cmprt in enumerate(cmprtsec):
        beg1, end1 = cmprt['start'], cmprt['end'] + 1
        diff1 = end1 - beg1
        scores[(k,k)] = dist_matrix[k][k] = -1
        for l in xrange(k + 1, len(cmprtsec)):
            beg2, end2 = cmprtsec[l]['start'], cmprtsec[l]['end'] + 1
            val = nansum([matrix[i][j] for i in xrange(beg1, end1)
                          for j in xrange(beg2, end2)]) / (end2 - beg2) / diff1
            try:
                scores[(k,l)] = dist_matrix[k][l] = scores[(l,k)] = dist_matrix[l][k] = func(val)
            except ZeroDivisionError:
                scores[(k,l)] = dist_matrix[k][l] = scores[(l,k)] = dist_matrix[l][k] = funczero(val)
            if isnan(scores[(k,l)]):
                scores[(k,l)] = dist_matrix[k][l] = scores[(l,k)] = dist_matrix[l][k] = funczero(0)
    # cluster compartments according to their correlation score
    try:
        clust = linkage(dist_matrix, method='ward')
    except UnboundLocalError:
        print('WARNING: Chromosome probably too small. Skipping')
        warn('WARNING: Chromosome probably too small. Skipping')
        return (float('inf'), float('inf'), float('inf'))
    # find best place to divide dendrogram (only check 1, 2, 3 or 4 clusters)
    solutions = {}
    for k in clust[:,2][-3:]:
        clusters = {}
        _ = [clusters.setdefault(j, []).append(i) for i, j in
             enumerate(fcluster(clust, k, criterion='distance'))]
        solutions[k] = {'out': clusters}
        solutions[k]['score'] = calinski_harabasz(scores, clusters)
    try:
        # take best cluster according to calinski_harabasz score
        clusters = [solutions[s] for s in sorted(
            solutions, key=lambda x: solutions[x]['score'])
                    if solutions[s]['score']>0][-1]['out']
    except IndexError:
        #warn('WARNING: compartment clustering is not clear. Skipping')
        return (float('inf'), float('inf'), float('inf'))
    if len(clusters) != 2:
        #warn('WARNING: compartment clustering is too clear. Skipping')
        return (float('inf'), float('inf'), float('inf'))
    # labelling compartments. A compartments shall have lower
    # mean intra-interactions
    dens = {}
    if rich_in_A:
        test = lambda x: x >= 1
    else:
        test = lambda x: x < 1
    for k in clusters:
        val = sum([cmprtsec[c]['dens']
                   for c in clusters[k]]) / len(clusters[k])
        dens['A' if test(val) else 'B'] = [
            cmprtsec[c]['dens'] for c in clusters[k]
            if cmprtsec[c]['end'] + 1 - cmprtsec[c]['start'] > 2]
        if save:
            for c in clusters[k]:
                cmprtsec[c]['type'] = 'A' if test(val) else 'B'
    try:
        tt, pval = ttest_ind(dens['A'], dens['B'])
    except ZeroDivisionError:
        return (float('inf'), float('inf'), float('inf'))
    prop = float(len(dens['A'])) / (len(dens['A']) + len(dens['B']))
    # to avoid having all A or all B
    score1 = 5000 * (prop - 0.5)**4 - 2
    # to avoid having  consecutive As or Bs
    score2 = 0
    prev = None
    for cmprt in cmprtsec:
        if cmprt.get('type', None) == prev:
            score2 += 1
        prev = cmprt.get('type', prev)
    score2 /= float(len(cmprtsec))
    score = score1 + score2
    if verbose:
        print ('[EV%d] g:%5s %5s%% tt:%7s '
               'score-proportion:%7s score-interleave:%7s '
               'final: %7s pv:%s' % (
                   ev_num, gamma - 1, round(prop * 100, 1),
                   round(tt, 3), round(score1, 3), round(score2, 3), 
                   round(score + tt, 3), pval))
    if log:
        log = open(log, 'a')
        log.write('[EV%d] g:%5s %5s%% tt:%7s '
                  'score-proportion:%7s score-interleave:%7s '
                  'final: %7s pv:%s\n' % (
                      ev_num, gamma - 1, round(prop * 100, 1),
                      round(tt, 3), round(score1, 3), round(score2, 3), 
                      round(score + tt, 3), pval))
        log.close()
    return score + tt, tt, prop

