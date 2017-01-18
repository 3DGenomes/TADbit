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
from pytadbit.utils.hmm             import gaussian_prob, best_path, train
from numpy.linalg                   import LinAlgError
from numpy                          import corrcoef, nansum, array, isnan, mean
from numpy                          import meshgrid, asarray, exp, linspace, std
from numpy                          import nanpercentile as npperc, log as nplog
from numpy                          import nanmax
from scipy.special                  import gammaincc
from scipy.cluster.hierarchy        import linkage, fcluster, dendrogram
from scipy.sparse.linalg            import eigsh
from pytadbit.utils.tadmaths        import calinski_harabasz
from scipy.stats                    import ttest_ind
from collections                    import OrderedDict
from warnings                       import warn
from bisect                         import bisect_right as bisect
from scipy.sparse                   import csr_matrix
import os

class HiC_data(dict):
    """
    This may also hold the print/write-to-file matrix functions
    """
    def __init__(self, items, size, chromosomes=None, dict_sec=None,
                 resolution=1, masked=None, symmetricized=False):
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
        self.symmetricized = symmetricized
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
    
    def get_hic_data_as_csr(self):
        """
        Returns a scipy sparse matrix in Compressed Sparse Row format of the HiC data in the dictionary

        :returns: scipy sparse matrix in Compressed Sparse Row format
        """
        values = []
        cols = []
        rows = []
        for key, value in self.iteritems():
            row, col = round(key / self.__size), key % self.__size
            values.append(float(value))
            cols.append(col)
            rows.append(row)
                
        return csr_matrix((values, (rows, cols)), shape=(self.__size,self.__size))
        
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
        try:
            return float(intra) / self.sum(bias=self.bias if normalized else None, bads=bads)
        except ZeroDivisionError:
            return 0.

    def filter_columns(self, draw_hist=False, savefig=None, perc_zero=75,
                       by_mean=True, min_count=None, silent=False):
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
        :param None min_count: minimum number of reads mapped to a bin (recommended
           value could be 2500). If set this option overrides the perc_zero
           filtering... This option is slightly slower.
        :param True by_mean: filter columns by mean column value using
           :func:`pytadbit.utils.hic_filtering.filter_by_mean` function

        """
        self.bads = filter_by_zero_count(self, perc_zero, min_count=min_count,
                                         silent=silent)
        if by_mean:
            self.bads.update(filter_by_mean(
                self, draw_hist=draw_hist, silent=silent,
                savefig=savefig, bads=self.bads))
        if not silent:
            print 'Found %d of %d columns with poor signal' % (len(self.bads),
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
                print '  - getting the sum of the matrix'
            # get the sum on half of the matrix
            norm_sum = self.sum(bias)
            if not silent:
                print '    => %.3f' % norm_sum
                print '  - rescaling biases'
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
        if normalized and not self.bias:
            raise Exception('ERROR: experiment not normalized yet')
        start1, start2, end1, end2 = self._focus_coords(focus)
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

    def _focus_coords(self, focus):
        siz = len(self)
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
        return start1, start2, end1, end2
            
    def find_compartments(self, crms=None, savefig=None, savedata=None,
                          savecorr=None, show=False, suffix='', how='',
                          label_compartments='hmm', log=None, max_mean_size=10000,
                          ev_index=None, rich_in_A=None, max_ev=3, **kwargs):
        """
        Search for A/B copartments in each chromosome of the Hi-C matrix.
        Hi-C matrix is normalized by the number interaction expected at a given
        distance, and by visibility (one iteration of ICE). A correlation matrix
        is then calculated from this normalized matrix, and its first
        eigenvector is used to identify compartments. Changes in sign marking
        boundaries between compartments.
        Result is stored as a dictionary of compartment boundaries, keys being
        chromosome names.
        
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
        :param 3 max_ev: maximum number of EV to try
        :param None ev_index: a list of number refering to the index of the
           eigenvector to be used. By default the first eigenvector is used.
           WARNING: index starts at 1, default is thus a list of ones. Note:
           if asking for only one chromosome the list should be only of one 
           element.
        :param None rich_in_A: by default compartments are identified using mean
           number of intra-interactions (A compartments are expected to have
           less). However this measure is not very accurate. Using this
           parameter a path to a BED or BED-Graph file with a list of genes or
           active epigenetic marks can be passed, and used instead of the mean
           interactions.
        :param None log: path to a folder where to save log of the assignment
           of A/B compartments
        :param hmm label_compartments: label compartments into A/B categories,
           otherwise just find borders (faster). Can be either hmm (default), or
           cluster.
        :param 'ratio' how: ratio divide by column, subratio divide by
           compartment, diagonal only uses diagonal
           

        TODO: this is really slow...

        Notes: building the distance matrix using the amount of interactions
               instead of the mean correlation, gives generally worse results.

        :returns: a dictionary with the N (max_ev) first eigen vectors used to define
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
                for row, posx in enumerate(xrange(self.section_pos[sec][0],
                                                  self.section_pos[sec][1])):
                    if posx in self.bads:
                        out.write(rownam.pop(0) + '\t' + empty)
                        badrows += 1
                        continue
                    vals = []
                    badcols = 0
                    for col, posy in enumerate(xrange(self.section_pos[sec][0],
                                                      self.section_pos[sec][1])):
                        if posy in self.bads:
                            vals.append('NaN')
                            badcols += 1
                            continue
                        vals.append(str(matrix[row-badrows][col-badcols]))
                    out.write(rownam.pop(0) + '\t' +'\t'.join(vals) + '\n')
                out.close()

            try:
                # This eighs is very very fast, only ask for one eigvector
                _, evect = eigsh(array(matrix), k=max_ev)
            except (LinAlgError, ValueError):
                warn('Chromosome %s too small to compute PC1' % (sec))
                cmprts[sec] = [] # Y chromosome, or so...
                count += 1
                continue
            index = ev_index[count] if ev_index else 1
            n_first = [list(evect[:, -i]) for i in xrange(1, max_ev + 1)]
            for ev_num in range(index, max_ev + 1):
                first = list(evect[:, -ev_num])
                breaks = [i for i, (a, b) in
                          enumerate(zip(first[1:], first[:-1]))
                          if a * b < 0] + [len(first) - 1]
                breaks = [{'start': breaks[i-1] + 1 if i else 0, 'end': b}
                          for i, b in enumerate(breaks)]
                if (self.resolution * (len(breaks) - 1.0) / len(matrix)
                    > max_mean_size):
                    warn('WARNING: number of compartments found with the '
                         'EigenVector number %d is too low (%d compartments '
                         'in %d rows), for chromosome %s' % (
                             ev_num, len(breaks), len(matrix), sec))
                else:
                    break
            if (self.resolution * (len(breaks) - 1.0) / len(matrix)
                > max_mean_size):
                warn('WARNING: keeping first eigenvector, for chromosome %s' % (
                    sec))
                ev_num = 1
            if ev_index:
                ev_num = ev_index[count]
            first = list(evect[:, -ev_num])
            breaks = [i for i, (a, b) in
                      enumerate(zip(first[1:], first[:-1]))
                      if a * b < 0] + [len(first) - 1]
            breaks = [{'start': breaks[i-1] + 1 if i else 0, 'end': b}
                      for i, b in enumerate(breaks)]
            ev_nums[sec] = ev_num
            beg, end = self.section_pos[sec]
            bads = [k - beg for k in self.bads if beg <= k <= end]
            for evect in n_first:
                _ = [evect.insert(b, float('nan')) for b in bads]
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
            firsts[sec] = n_first
            # needed for the plotting
            self._apply_metric(cmprts, sec, rich_in_A, how=how)
            
            if label_compartments == 'cluster':
                if log:
                    logf = os.path.join(log, sec + suffix + '.log')
                else:
                    logf = None

                gammas = {}
                for n_clust in range(2, 4):
                    for gamma in range(0, 101, 1):
                        scorett, tt, prop = _cluster_ab_compartments(
                            float(gamma)/100, matrix, breaks, cmprts[sec],
                            rich_in_A, ev_num=ev_num, log=logf, save=False,
                            verbose=kwargs.get('verbose', False),
                            n_clust=n_clust)
                        gammas[gamma] = scorett, tt, prop
                    gamma = min(gammas.keys(), key=lambda k: gammas[k][0])
                    if gammas[gamma][0] - gammas[gamma][1] > 7:
                        print (' WARNING: minimum showing very low '
                               'intermeagling of A/B compartments, trying '
                               'with 3 clusters, for chromosome %s', sec)
                        gammas = {}
                        continue
                    if kwargs.get('verbose', False):
                        print '   ====>  minimum:', gamma
                    break
                _ = _cluster_ab_compartments(float(gamma)/100, matrix, breaks,
                                          cmprts[sec], rich_in_A, save=True,
                                          log=logf, ev_num=ev_num, n_clust=n_clust)

            if savefig or show:
                vmin = kwargs.get('vmin', -1)
                vmax = kwargs.get('vmax',  1)
                if vmin == 'auto' == vmax:
                    vmax = max([abs(npperc(matrix, 99.5)),
                                abs(npperc(matrix, 0.5))])
                    vmin = -vmax
                plot_compartments(
                    sec, first, cmprts, matrix, show,
                    savefig + '/chr' + str(sec) + suffix + '.pdf' if savefig else None,
                    vmin=vmin, vmax=vmax, whichpc=ev_num)
                plot_compartments_summary(
                    sec, cmprts, show,
                    savefig + '/chr' + str(sec) + suffix + '_summ.pdf' if savefig else None)
            count += 1

        if label_compartments == 'hmm':
            x = {}
            for sec in self.section_pos:
                try:
                    x[sec] = firsts[sec][ev_nums[sec] - 1]
                except KeyError:
                    continue

            # train two HMMs on the genomic data:
            #  - one with 2 states A B
            #  - one with 3 states A B I 
            #  - one with 4 states A a B b
            #  - one with 5 states A a B b I 
            models = {}
            for n in range(2, 6):
                if kwargs.get('verbose', False):
                    print ('Training HMM for %d categories of '
                           'compartments' % n)
                models[n] = _training(x, n, kwargs.get('verbose', False))

            # apply HMM models on each chromosome
            results = {}
            for sec in self.section_pos:
                if not sec in x:
                    continue
                if kwargs.get('verbose', False):
                    print 'Chromosome', sec
                beg, end = self.section_pos[sec]
                bads = [k - beg for k in self.bads if beg <= k <= end]
                # print 'CMPRTS before   ', sec, cmprts[sec]
                n_states, breaks = _hmm_refine_compartments(
                    x, sec, models, bads, kwargs.get('verbose', False))
                results[sec] = n_states, breaks
                cmprts[sec] = breaks
                # print 'CMPRTS after hmm', sec, cmprts[sec]
                self._apply_metric(cmprts, sec, rich_in_A, how=how)
                
                if rich_in_A:
                    test = lambda x: x >= 1
                else:
                    test = lambda x: x < 1
                max_type = nanmax([c['type'] for c in cmprts[sec]])

                # find which category of compartment has the highest "density"
                atyp = 0.
                alen = 0.
                btyp = 0.
                blen = 0.
                max_type = nanmax([c['type'] for c in cmprts[sec]])
                for typ in range(5):
                    subset = set([i for i, c in enumerate(cmprts[sec])
                                 if c['type'] == typ])
                    dens = sum(cmprts[sec][c]['dens'] * (cmprts[sec][c]['end'] - cmprts[sec][c]['start']) for c in subset)
                    leng = sum((cmprts[sec][c]['end'] - cmprts[sec][c]['start'])**2 / 2. for c in subset)
                    # leng = sum(1 for c in subset)
                    val = float(dens) / leng if leng else 0.
                    #if typ == 0:
                    if typ < max_type / 2.:
                        alen += leng
                        atyp += val * leng
                    # elif typ == max_type:
                    elif typ > max_type / 2.:
                        blen += leng
                        btyp += val * leng

                for i, comp in enumerate(cmprts[sec]):
                    if comp['type'] < max_type / 2.:
                        # if mean density of compartments of type 0 is higher than 1
                        # than label them as 'B', otherwise, as 'A'
                        if comp['type'] == 0:
                            comp['type'] = 'A' if test(val) else 'B'
                        else:
                            comp['type'] = 'a' if test(val) else 'b'
                    elif comp['type'] > max_type / 2.:
                        if comp['type'] == max_type:
                            comp['type'] = 'B' if test(val) else 'A'
                        else:
                            comp['type'] = 'b' if test(val) else 'a'
                    elif isnan(comp['type']):
                        comp['type'] = 'NA'
                    else:
                        comp['type'] = 'I'
        self.compartments = cmprts
        if savedata:
            self.write_compartments(savedata, chroms=self.compartments.keys(),
                                    ev_nums=ev_nums)
        return firsts

    def _apply_metric(self, cmprts, sec, rich_in_A, how='ratio'):
        """
        calculate compartment internal density if no rich_in_A, otherwise
        sum this list
        """
        # print 'SEGMENTS'
        # print sec, self.section_pos[sec]
        # for i in range(0, len(cmprts[sec]), 20):
        #     print '  ' + ''.join(['%5d/%-5d'% (s['start'], s['end']) for s in cmprts[sec][i:i+20]])
        # print 'CHROMOSOME', sec
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
                beg, end = self.section_pos[sec]
                beg1, end1 = cmprt['start'] + beg, cmprt['end'] + beg + 1
                # print 'BEG:%7d, END:%7d, LEN bias:%7d, LEN self:%7d, LEN expected:%7d' % (beg1, end1, len(self.bias),
                                                                                          # len(self), len(self.expected))
                if 'diagonal' in how:
                    sec_matrix = [(self[i,i] / self.expected[0] / self.bias[i]**2)
                                  for i in xrange(beg1, end1) if not i in self.bads]
                else: #if 'compartment' in how:
                    sec_matrix = [(self[i,j] / self.expected[abs(j-i)]
                                   / self.bias[i] / self.bias[j])
                                  for i in xrange(beg1, end1) if not i in self.bads
                                  for j in xrange(beg1, end1) if not j in self.bads]
                if '/compartment' in how: # diagonal / compartment
                    sec_column = [(self[i,j] / self.expected[abs(j-i)]
                                   / self.bias[i] / self.bias[j])
                                  for i in xrange(beg1, end1) if not i in self.bads
                                  for j in xrange(beg1, end1) if not j in self.bads]
                elif '/column' in how:
                    sec_column = [(self[i,j] / self.expected[abs(j-i)]
                                   / self.bias[i] / self.bias[j])
                                  for i in xrange(beg1, end1) if not i in self.bads
                                  for j in range(beg, end)
                                  if not j in self.bads]
                else:
                    sec_column = [1.]
                try:
                    if 'type' in cmprt and isnan(cmprt['type']):
                        cmprt['dens'] = 1.
                    else:
                        cmprt['dens'] = float(sum(sec_matrix)) / sum(sec_column)
                except ZeroDivisionError:
                    cmprt['dens'] = 1.
        # normalize to 1.0
        try:
            if 'type' in cmprt: # hmm already run and have the types definded
                meanh = (sum(cmprt['dens'] for cmprt in cmprts[sec]
                             if not isnan(cmprt['type'])) /
                         sum(1 for cmprt in cmprts[sec]
                             if not isnan(cmprt['type'])))
            else:
                meanh = (sum(cmprt['dens'] for cmprt in cmprts[sec]) /
                         sum(1 for cmprt in cmprts[sec]))
        except ZeroDivisionError:
            meanh = 1.
        for cmprt in cmprts[sec]:
            try:
                if 'type' in cmprt and isnan(cmprt['type']):
                    cmprt['dens'] = 1.0
                else:
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
                try:
                    out.write('## CHR %s\tEigenvector: %d\n' % (sec, ev_nums[sec]))
                except KeyError:
                    continue
        out.write('#%sstart\tend\tdensity\ttype\n'% (
            'CHR\t' if len(sections) > 1 else ''))
        for sec in sections:
            for c in self.compartments[sec]:
                try:
                    out.write('%s%d\t%d\t%.2f\t%s\n' % (
                        (sec + '\t') if sections else '',
                        c['start'] + 1, c['end'] + 1,
                        c['dens'], c['type']))
                except KeyError:
                    out.write('%s%d\t%d\t%.2f\t%s\n' % (
                        (sec + '\t') if sections else '',
                        c['start'], c['end'], c['dens'], ''))
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

def _hmm_refine_compartments(x, sec, models, bads, verbose):
    prevll = float('-inf')
    prevdf = 0
    results = {}
    for n in range(2, 6):
        E, pi, T = models[n]
        probs = gaussian_prob(x[sec], E)
        pathm, llm = best_path(probs, pi, T)
        pathm = asarray(map(float, pathm))
        df = n**2 - n + n * 2 + n - 1
        len_seq = len(pathm)
        lrt = gammaincc((df - prevdf) / 2., (llm - prevll) / 2.)
        bic = -2 * llm + df * nplog(len_seq)
        aic = 2 * df - 2 * llm
        if verbose:
            print 'Ll for %d states (%d df): %4.0f AIC: %4.0f BIC: %4.0f LRT=%f'% (
                n, df, llm, aic, bic, lrt)
        prevdf = df
        prevll = llm
        results[n] = {'AIC': aic,
                      'BIC': bic,
                      'LRT': lrt,
                      'PATH': pathm
                      }
    n_states = min(results, key=lambda x: results[x]['BIC'])
    results = list(results[n_states]['PATH'])
    # print 'RESULTS', results
    _ = [results.insert(b, float('nan')) for b in sorted(bads)]
    # print 'RESULTS', results
    breaks = [(i, b) for i, (a, b) in
              enumerate(zip(results[1:], results[:-1]))
              if str(a) != str(b)] + [len(results) - 1]
    # print 'BREAKS', breaks
    breaks[-1] = (breaks[-1], results[-1])
    # print 'BREAKS', breaks
    breaks = [{'start': breaks[i-1][0] + 1 if i else 0, 'end': b,
               'type': a}
              for i, (b, a) in enumerate(breaks)]
    # print 'BREAKS', breaks
    return n_states, breaks

def _training(x, n, verbose):
    """
    define default emision transition and initial states, and train the hmm
    """
    pi = [0.5 - ((n - 2) * 0.05)**2 if i == 0 or i == n - 1 else ((n - 2)*0.05)**2*2 / (n - 2) for i in range(n)]
    T = [[0.9 if i==j else 0.1/(n-1) for i in xrange(n)] for j in xrange(n)]
    E =  asarray(zip(linspace(-1, 1, n), [1./n for _ in range(n)]))

    # normalize values of the first eigenvector
    for c in x:
        this_mean = mean(x[c])
        this_std  = std (x[c])
        x[c] = [v - this_mean for v in x[c]]
        x[c] = [v / this_std  for v in x[c]]

    train(pi, T, E, x.values(), verbose=verbose, threshold=1e-6, n_iter=1000)
    return E, pi, T
    
def _cluster_ab_compartments(gamma, matrix, breaks, cmprtsec, rich_in_A, save=True,
                             ev_num=1, log=None, verbose=False, savefig=None,
                             n_clust=2):
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
    # plot
    if savefig:
        xedges = [b['start'] for b in breaks]
        yedges = [b['start'] for b in breaks]
        xedges += [breaks[-1]['end']]
        yedges += [breaks[-1]['end']]
        X, Y = meshgrid(xedges, yedges)
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10,10))
        _ = fig.add_axes([0.09,0.1,0.2,0.6])
        Z1 = dendrogram(clust, orientation='left')
        idx1 = Z1['leaves']
        idx2 = Z1['leaves']
        D = asarray(dist_matrix)[idx1,:]
        D = D[:,idx2]
        axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
        m = axmatrix.pcolormesh(X, Y, D)
        axmatrix.set_aspect('equal')
        axmatrix.set_yticks([])
        axmatrix.set_xlim((0, breaks[-1]['end']))
        axmatrix.set_ylim((0, breaks[-1]['end']))
        plt.colorbar(m)
        plt.savefig(savefig)
    try:
        # take best cluster according to calinski_harabasz score
        clusters = [solutions[s] for s in sorted(
            solutions, key=lambda x: solutions[x]['score'])
                    if solutions[s]['score']>0][1 - n_clust]['out']
    except IndexError:
        # warn('WARNING1: compartment clustering is not clear. Skipping')
        return (float('inf'), float('inf'), float('inf'))
    if len(clusters) != n_clust:
        # warn('WARNING2: compartment clustering is too clear. Skipping')
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
        for c in clusters[k]:
            cmprtsec[c]['type'] = 'A' if test(val) else 'B'
    try:
        tt, pval = ttest_ind(dens['A'], dens['B'])
    except ZeroDivisionError:
        return (float('inf'), float('inf'), float('inf'))
    prop = float(len(dens['A'])) / (len(dens['A']) + len(dens['B']))
    # to avoid having all A or all B
    # score = 5000 * (prop - 0.5)**4 - 2
    # to avoid having  consecutive As or Bs
    score = 0.
    prev = None
    for cmprt in cmprtsec:
        if cmprt.get('type', None) == prev:
            score += 1.
        prev = cmprt.get('type', prev)
    score /= len(cmprtsec)
    score = exp(10 * (score - 0.4)) # 5000 * (score - 0.5)**4 - 2
    # score = score1 + score2
    if verbose:
        print ('[EV%d CL%s] g:%5s prop:%5s%% tt:%7s '
               'score-interleave:%5s ' # score-proportion:%7s 
               'final: %7s pv:%7s' % (
                   ev_num, n_clust, gamma - 1, round(prop * 100, 1),
                   round(tt, 3), round(score, 3), #round(score2, 3), 
                   round(score + tt, 3), round(pval, 5)))
    if log:
        log = open(log, 'a')
        log.write('[EV%d CL%s] g:%5s prop:%5s%% tt:%6s '
                  'score-interleave:%6s ' # score-proportion:%7s 
                  'final: %7s pv:%s\n' % (
                      ev_num, n_clust, gamma - 1, round(prop * 100, 1),
                      round(tt, 3), round(score, 3), # round(score2, 3), 
                      round(score + tt, 3), round(pval, 4)))
        log.close()
    if not save:
        for cmprt in cmprtsec:
            if 'type' in cmprt:
                cmprt['type'] = None 
    return score + tt, tt, prop


