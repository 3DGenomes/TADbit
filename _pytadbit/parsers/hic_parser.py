"""
November 7, 2013.

"""

from warnings import warn
from math import sqrt, isnan
from pytadbit.parsers.gzopen import gzopen
from pytadbit.utils.hic_filtering   import filter_by_mean
from collections import OrderedDict
from pytadbit.utils.normalize_hic  import iterative

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
    for line in f:
        if line[0] != '#':
            break
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
            _ = [float(item) for item in items[0]]
            # Case 1: pure number matrix.
            header = False
            trim = 0
        except ValueError:
            # Case 2: matrix with row and column names.
            header = True
            trim = 1
    else:
        if len(items[0]) == len(items[1]):
            # Case 3: matrix with row information.
            header = False
            trim = ncol - nrow
        else:
            # Case 4: matrix with header and row information.
            header = True
            trim = ncol - nrow + 1
    # Remove header line if needed.
    if header:
        del(items[0])
        nrow -= 1

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
    if ncol != nrow: raise AutoReadFail('ERROR: non square matrix')

    if is_asymmetric(items):
        warn('WARNING: input matrix not symmetric: symmetrizing')
        symmetrize(items)

    return tuple([a for line in items for a in line]), ncol


def read_matrix(things, parser=None, hic=True):
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


    :param True hic: if False, TADbit assumes that files contains normalized
       data
    :returns: the corresponding matrix concatenated into a huge list, also
       returns number or rows

    """
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
            matrix, size = parser(thing)
            thing.close()
            matrices.append(HiC_data([(i, matrix[i]) for i in xrange(size**2)
                                      if matrix[i]], size))
        elif isinstance(thing, str):
            try:
                matrix, size = parser(gzopen(thing))
            except IOError:
                if len(thing.split('\n')) > 1:
                    matrix, size = parser(thing.split('\n'))
                else:
                    raise IOError('\n   ERROR: file %s not found\n' % thing)
            matrices.append(HiC_data([(i, matrix[i]) for i in xrange(size**2)
                                      if matrix[i]], size))
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
        
    return matrices


def load_hic_data_from_reads(fnam, resolution):
    """
    :param fnam: tsv file with reads1 and reads2
    :param resolution: the resolution of the experiment (size of a bin in
       bases)
    :param genome_seq: a dictionary containing the genomic sequence by
       chromosome
    """
    sections = []
    genome_seq = OrderedDict()
    fhandler = open(fnam)
    line = fhandler.next()
    while line.startswith('#'):
        if line.startswith('# CRM '):
            crm, clen = line[6:].split()
            genome_seq[crm] = int(clen) / resolution + 1
        line = fhandler.next()
    size = 0
    section_sizes = {}
    for crm in genome_seq:
        len_crm = genome_seq[crm]
        section_sizes[(crm,)] = len_crm
        size += len_crm
        sections.extend([(crm, i) for i in xrange(len_crm)])
    dict_sec = dict([(j, i) for i, j in enumerate(sections)])
    imx = HiC_data((), size, genome_seq, dict_sec)
    while True:
        _, cr1, ps1, _, _, _, _, cr2, ps2, _ = line.split('\t', 9)
        ps1 = dict_sec[(cr1, int(ps1) / resolution)]
        ps2 = dict_sec[(cr2, int(ps2) / resolution)]
        imx[ps1 + ps2 * size] += 1
        imx[ps2 + ps1 * size] += 1
        try:
            line = fhandler.next()
        except StopIteration:
            break
    return imx


class HiC_data(dict):
    """
    This may also hold the print/write-to-file matrix functions
    """
    def __init__(self, items, size, chromosomes=None, dict_sec=None):
        super(HiC_data, self).__init__(items)
        self.__size = size
        self._size2 = size**2
        self.bias = None
        self.bads = None
        self.chromosomes = chromosomes
        self.sections = dict_sec
        self.section_pos = {}
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


    def filter_columns(self, draw_hist=False, savefig=None):
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

        """
        self.bads = filter_by_mean(self, draw_hist=draw_hist, savefig=savefig)

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

    def get_matrix(self, focus=None, diagonal=True, normalized=False):
        """
        get the matrix
        """
        siz = len(self)
        if normalized and not self.bias:
            raise Exception('ERROR: experiment not normalized yet')
        if focus:
            if isinstance(focus, tuple) and isinstance(focus[0], int):
                start1, end1 = focus
                start2, end2 = focus
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
                    for i in xrange(start1, end1):
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
                    for i in xrange(start1, end1):
                        mtrx[i][i] = 1 if mtrx[i][i] else 0
                return mtrx
