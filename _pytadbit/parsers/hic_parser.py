"""
November 7, 2013.

"""

from warnings import warn
from math import sqrt, isnan
from pytadbit.parsers.gzopen import gzopen
from pytadbit.utils.hic_filtering   import filter_by_mean, filter_by_zero_count
from collections import OrderedDict
from pytadbit.utils.normalize_hic  import iterative
from pytadbit.parsers.genome_parser import parse_fasta

HIC_DATA = True

# Exception to handle failed autoread.
class AutoReadFail(Exception):
    pass

# Helper functions for the autoreader.
def is_asymmetric(matrix):
    maxn = len(matrix)
    for i in range(maxn):
        maxi = matrix[i] # slightly more efficient
        for j in range(i+1, maxn):
            if maxi[j] != matrix[j][i]:
                if isnan(maxi[j]) and isnan(matrix[j][i]):
                    continue
                return True
    return False

def symmetrize(matrix):
    maxn = len(matrix)
    for i in range(maxn):
        for j in range(i+1, maxn):
            matrix[i][j] = matrix[j][i] = matrix[i][j] + matrix[j][i]


def autoreader(f):
    """
    Auto-detect matrix format of HiC data file.
    
    :param f: an iterable (typically an open file).
    
    :returns: A tuple with integer values and the dimension of
       the matrix.
    """

    # Skip initial comment lines and read in the whole file
    # as a list of lists.
    masked = {}
    for line in f:
        if line[0] != '#':
            break
        if line.startswith('# MASKED'):
            masked = dict([(int(n), True) for n in line.split()[2:]])
    items = [line.split()] + [line.split() for line in f]

    # Count the number of elements per line after the first.
    # Wrapping in a set is a trick to make sure that every line
    # has the same number of elements.
    S = set([len(line) for line in items[1:]])
    ncol = S.pop()
    # If the set 'S' is not empty, at least two lines have a
    # different number of items.
    if S:
        raise AutoReadFail('ERROR: unequal column number')

    nrow = len(items)
    # Auto-detect the format, there are only 4 cases.
    if ncol == nrow:
        try:
            _ = [float(item) for item in items[0]
                 if not item.lower() in ['na', 'nan']]
            # Case 1: pure number matrix.
            header = False
            trim = 0
        except ValueError:
            # Case 2: matrix with row and column names.
            header = True
            trim = 1
            warn('WARNING: found header')
    else:
        if len(items[0]) == len(items[1]):
            # Case 3: matrix with row information.
            header = False
            trim = ncol - nrow
            warn('WARNING: found %d colum(s) of row names' % trim)
        else:
            # Case 4: matrix with header and row information.
            header = True
            trim = ncol - nrow + 1
            warn('WARNING: found header and %d colum(s) of row names' % trim)
    # Remove header line if needed.
    if header and not trim:
        header = items.pop(0)
        nrow -= 1
    elif not trim:
        header = range(1, nrow + 1)
    elif not header:
        header = [tuple([a for a in line[:trim]]) for line in items]
    else:
        del(items[0])
        nrow -= 1
        header = [tuple([a for a in line[:trim]]) for line in items]
    # Get the numeric values and remove extra columns
    what = int if HIC_DATA else float
    try:
        items = [[what(a) for a in line[trim:]] for line in items]
    except ValueError:
        if not HIC_DATA:
            raise AutoReadFail('ERROR: non numeric values')
        try:
            # Dekker data 2009, uses integer but puts a comma... 
            items = [[int(float(a)+.5) for a in line[trim:]] for line in items]
            warn('WARNING: non integer values')
        except ValueError:
            try:
                # Some data may contain 'NaN' or 'NA'
                items = [
                    [0 if a.lower() in ['na', 'nan']
                     else int(float(a)+.5) for a in line[trim:]]
                for line in items]
                warn('WARNING: NA or NaN founds, set to zero')
            except ValueError:
                raise AutoReadFail('ERROR: non numeric values')

    # Check that the matrix is square.
    ncol -= trim
    if ncol != nrow:
        raise AutoReadFail('ERROR: non square matrix')

    if is_asymmetric(items):
        warn('WARNING: input matrix not symmetric: symmetrizing')
        symmetrize(items)
    return tuple([a for line in items for a in line]), ncol, header, masked

def _header_to_section(header, resolution):
    """
    converts row-names of the form 'chr12\t1000-2000' into sections, suitable
    to create HiC_data objects. Also creates chromosomes, from the reads
    """
    chromosomes = OrderedDict()
    sections = {}
    sections = {}
    chromosomes = None
    if (isinstance(header, list)
        and isinstance(header[0], tuple)
        and len(header[0]) > 1):
        chromosomes = OrderedDict()
        for i, h in enumerate(header):
            if '-' in h[1]:
                a, b = map(int, h[1].split('-'))
                if resolution==1:
                    resolution = abs(b - a)
                elif resolution != abs(b - a):
                    raise Exception('ERROR: found different resolution, ' +
                                    'check headers')
            else:
                a = int(h[1])
                if resolution==1 and i:
                    resolution = abs(a - b)
                elif resolution == 1:
                    b = a
            sections[(h[0], a / resolution)] = i
            chromosomes.setdefault(h[0], 0)
            chromosomes[h[0]] += 1
    return chromosomes, sections, resolution

def read_matrix(things, parser=None, hic=True, resolution=1, **kwargs):
    """
    Read and checks a matrix from a file (using
    :func:`pytadbit.parser.hic_parser.autoreader`) or a list.

    :param things: might be either a file name, a file handler or a list of
        list (all with same length)
    :param None parser: a parser function that returns a tuple of lists
       representing the data matrix,
       with this file example.tsv:
       ::
       
         chrT_001	chrT_002	chrT_003	chrT_004
         chrT_001	629	164	88	105
         chrT_002	86	612	175	110
         chrT_003	159	216	437	105
         chrT_004	100	111	146	278

       the output of parser('example.tsv') might be:
       ``([629, 86, 159, 100, 164, 612, 216, 111, 88, 175, 437, 146, 105, 110,
       105, 278])``

    :param 1 resolution: resolution of the matrix
    :param True hic: if False, TADbit assumes that files contains normalized
       data
    :returns: the corresponding matrix concatenated into a huge list, also
       returns number or rows

    """
    one = kwargs.get('one', True)
    global HIC_DATA
    HIC_DATA = hic
    parser = parser or autoreader
    if not isinstance(things, list):
        things = [things]
    matrices = []
    for thing in things:
        if isinstance(thing, HiC_data):
            matrices.append(thing)
        elif isinstance(thing, file):
            matrix, size, header = parser(thing)
            thing.close()
            chromosomes, sections, resolution = _header_to_section(header,
                                                                   resolution)
            matrices.append(HiC_data([(i, matrix[i]) for i in xrange(size**2)
                                      if matrix[i]], size, dict_sec=sections,
                                     chromosomes=chromosomes,
                                     resolution=resolution))
        elif isinstance(thing, str):
            try:
                matrix, size, header, masked = parser(gzopen(thing))
            except IOError:
                if len(thing.split('\n')) > 1:
                    matrix, size, header, masked = parser(thing.split('\n'))
                else:
                    raise IOError('\n   ERROR: file %s not found\n' % thing)
            sections = dict([(h, i) for i, h in enumerate(header)])
            chromosomes, sections, resolution = _header_to_section(header,
                                                                   resolution)
            matrices.append(HiC_data([(i, matrix[i]) for i in xrange(size**2)
                                      if matrix[i]], size, dict_sec=sections,
                                     chromosomes=chromosomes,
                                     resolution=resolution))
        elif isinstance(thing, list):
            if all([len(thing)==len(l) for l in thing]):
                matrix  = reduce(lambda x, y: x+y, thing)
                size = len(thing)
            else:
                raise Exception('must be list of lists, all with same length.')
            matrices.append(HiC_data([(i, matrix[i]) for i in xrange(size**2)
                                      if matrix[i]], size))
        elif isinstance(thing, tuple):
            # case we know what we are doing and passing directly list of tuples
            matrix = thing
            siz = sqrt(len(thing))
            if int(siz) != siz:
                raise AttributeError('ERROR: matrix should be square.\n')
            size = int(siz)
            matrices.append(HiC_data([(i, matrix[i]) for i in xrange(size**2)
                                      if matrix[i]], size))
        elif 'matrix' in str(type(thing)):
            try:
                row, col = thing.shape
                if row != col:
                    raise Exception('matrix needs to be square.')
                matrix  = thing.reshape(-1).tolist()[0]
                size = row
            except Exception as exc:
                print 'Error found:', exc
            matrices.append(HiC_data([(i, matrix[i]) for i in xrange(size**2)
                                      if matrix[i]], size))
        else:
            raise Exception('Unable to read this file or whatever it is :)')
    if one:
        return matrices[0]
    else:
        return matrices

def load_hic_data_from_reads(fnam, resolution, **kwargs):
    """
    :param fnam: tsv file with reads1 and reads2
    :param resolution: the resolution of the experiment (size of a bin in
       bases)
    :param genome_seq: a dictionary containing the genomic sequence by
       chromosome
    :param False get_sections: for very very high resolution, when the column
       index does not fit in memory
    """
    sections = []
    genome_seq = OrderedDict()
    fhandler = open(fnam)
    line = fhandler.next()
    size = 0
    while line.startswith('#'):
        if line.startswith('# CRM '):
            crm, clen = line[6:].split()
            genome_seq[crm] = int(clen) / resolution + 1
            size += genome_seq[crm]
        line = fhandler.next()
    section_sizes = {}
    if kwargs.get('get_sections', True):
        for crm in genome_seq:
            len_crm = genome_seq[crm]
            section_sizes[(crm,)] = len_crm
            sections.extend([(crm, i) for i in xrange(len_crm)])
    dict_sec = dict([(j, i) for i, j in enumerate(sections)])
    imx = HiC_data((), size, genome_seq, dict_sec, resolution=resolution)
    while True:
        _, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.split('\t', 9)
        try:
            ps1 = dict_sec[(cr1, int(ps1) / resolution)]
            ps2 = dict_sec[(cr2, int(ps2) / resolution)]
        except KeyError:
            ps1 = int(ps1) / resolution
            ps2 = int(ps2) / resolution
        imx[ps1, ps2] += 1
        imx[ps2, ps1] += 1
        try:
            line = fhandler.next()
        except StopIteration:
            break
    return imx


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
        if self.chromosomes:
            total = 0
            for crm in self.chromosomes:
                self.section_pos[crm] = (total, total + self.chromosomes[crm])
                total += self.chromosomes[crm]

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
                        equals=None, verbose=False):
        """
        Counts the number of interactions occuring within chromosomes (cis) with
        respect to the total number of interactions

        :param False normalized: used normalized data
        :param None exclude: exclude a given list of chromosome from the
           ratio (may want to exclude translocated chromosomes)
        :param False diagonal: replace values in the diagonal by 0 or 1
        :param None equals: can pass a function that would decide if 2 chromosomes
           have to be considered as the same. e.g. lambda x, y: x[:4]==y[:4] will
           consider chr2L and chr2R as being the same chromosome

        :returns: the ratio of cis interactions over the total number of
           interactions. This number is expected to be between at least 40-60%
           in Human classic dilution Hi-C with HindIII as restriction enzyme.
        """
        if exclude == None:
            exclude = []
        if equals == None:
            equals = lambda x, y: x == y
        intra = inter = 0
        if not self.chromosomes:
            return float('nan')
        for crm1 in self.chromosomes:
            for crm2 in self.chromosomes:
                if crm1 in exclude or crm2 in exclude:
                    continue
                if crm1 == crm2:
                    mtrx = self.get_matrix(
                        focus=(crm1, crm2), normalized=normalized, diagonal=diagonal)
                    val = sum([sum([mtrx[i][j] for j in xrange(len(mtrx))])
                               for i in xrange(len(mtrx))])
                    if verbose:
                        print 'INTRA', crm1, crm2, val
                    intra += val
                else:
                    val = sum([sum(d) for d in self.get_matrix(
                        focus=(crm1, crm2), normalized=normalized, diagonal=diagonal)])
                    if equals(crm1, crm2):
                        if verbose:
                            print '  INTRA', crm1, crm2, val
                        intra += val
                    else:
                        if verbose:
                            print '  INTER', crm1, crm2, val
                        inter += val
        return float(intra) / (intra + inter)
    

    def filter_columns(self, draw_hist=False, savefig=None, perc_zero=75):
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

        """
        self.bads = filter_by_zero_count(self, perc_zero, silent=False)
        self.bads.update(filter_by_mean(self, draw_hist=draw_hist,
                                        savefig=savefig, bads=self.bads))

    def normalize_hic(self, iterations=0, max_dev=0.1, silent=False):
        """
        Normalize the Hi-C data.

        It fills the Experiment.norm variable with the Hi-C values divided by
        the calculated weight.

        :param 0 iteration: number of iterations
        :param 0.1 max_dev: iterative process stops when the maximum deviation
           between the sum of row is equal to this number (0.1 means 10%)
        :param False silent: does not warn when overwriting weights
        """
        self.bias = iterative(self, iterations=iterations,
                              max_dev=max_dev, bads=self.bads,
                              verbose=not silent)

    def get_as_tuple(self):
        return tuple([self[i, j]
                      for j in xrange(len(self))
                      for i in xrange(len(self))])

    def write_matrix(self, fname, focus=None, diagonal=True, normalized=False):
        """
        writes the matrix to a file
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
        out.write('# MASKED %s\n' % (' '.join(self.bads.keys())))
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
        :param True diagonal: if False, diagonal is replaced by zeroes
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
                        mtrx[i][i] = 1 if mtrx[i][i] else 0
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
