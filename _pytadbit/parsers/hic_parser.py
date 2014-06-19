"""
November 7, 2013.

"""

from warnings import warn
from math import sqrt, isnan
from pytadbit.parsers.gzopen import gzopen

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
        if all([item.isdigit() or item == 'nan' for item in items[0]]):
            # Case 1: pure integer matrix.
            header = False
            trim = 0
        else:
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
        try:
            # Dekker data 2009, uses integer but puts a comma... 
            items = [[int(float(a)+.5) for a in line[trim:]] for line in items]
            warn('WARNING: non integer values')
        except ValueError:
            try:
                # Some data may contain 'NaN' or 'NA'
                items = [
                    [float('nan') if a.lower() in ['na', 'nan']
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
                print thing
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


class HiC_data(dict):
    """
    This may also hold the print/write-to-file matrix functions
    """
    def __init__(self, items, size):
        super(HiC_data, self).__init__(items)
        self.__size = size
        self._size2 = size**2
        self.bias = None

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
                    'ERROR: position larger than %s^2' % self.__size)
            return self.get(row_col, 0)


    def get_as_tuple(self):
        return tuple([self[i, j] for j  in xrange(len(self)) for i in xrange(len(self))])

    def get_matrix(self, focus=None, diagonal=True, normalized=False):
        """
        get the matrix
        """
        siz = len(self)
        if normalized and not self.bias:
            raise Exception('ERROR: experiment not normalized yet')
        if focus:
            start, end = focus
            start -= 1
        else:
            start = 0
            end   = siz
        if normalized:
            if diagonal:
                return [[self[i, j] / self.bias[i] / self.bias[j]
                         for i in xrange(start, end)]
                        for j in xrange(start, end)]
            else:
                mtrx = [[self[i, j] / self.bias[i] / self.bias[j]
                         for i in xrange(start, end)]
                        for j in xrange(start, end)]
                for i in xrange(start, end):
                    mtrx[i][i] = 1 if mtrx[i][i] else 0
                return mtrx
        else:
            if diagonal:
                return [[self[i, j] for i in xrange(start, end)]
                        for j in xrange(start, end)]
            else:
                mtrx = [[self[i, j] for i in xrange(start, end)]
                        for j in xrange(start, end)]
                for i in xrange(start, end):
                    mtrx[i][i] = 1 if mtrx[i][i] else 0
                return mtrx
