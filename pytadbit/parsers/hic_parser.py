"""
November 7, 2013.

"""

from warnings import warn
from math import sqrt

# Exception to handle failed autoread.
class AutoReadFail(Exception):
    pass

# Helper functions for the autoreader.
def is_asymmetric(matrix):
    maxn = len(matrix)
    for i in range(maxn):
        maxi = matrix[i] # slightly more efficient
        for j in range(i+1, maxn):
            if maxi[j] != matrix[j][i]: return True
    return False

def symmetrize(matrix):
    maxn = len(matrix)
    for i in range(maxn):
        for j in range(i+1, maxn):
            matrix[i][j] = matrix[j][i] = matrix[i][j] + matrix[j][i]


def autoreader(f):
    """Auto-detect matrix format of HiC data file.
    
    ARGUMENTS:
       f: an iterable (typically an open file).
    
    RETURN:
       A tuple with integer values and the dimension of
       the matrix."""

    # Skip initial comment lines and read in the whole file
    # as a list of lists.
    for line in f:
        if line[0] != '#': break
    items = [line.split()]
    items += [line.split() for line in f]

    # Count the number of elements per line after the first.
    # Wrapping in a set is a trick to make sure that every line
    # has the same number of elements.
    S = set([len(line) for line in items[1:]])
    ncol = S.pop()
    # If the set 'S' is not empty, at least two lines have a
    # different number of items.
    if S: raise AutoReadFail('unequal column number')

    nrow = len(items)
    # Auto-detect the format, there are only 4 cases.
    if ncol == nrow:
        if all([item.isdigit() for item in items[0]]):
            # Case 1: pure integer matrix.
            header = False
            startat = 0
        else:
            # Case 2: matrix with row and column names.
            header = True
            startat = 1
    else:
        if len(items[0]) == len(items[1]):
            # Case 3: matrix with row information.
            header = False
            startat = ncol - nrow
        else:
            # Case 4: matrix with header and row information.
            header = True
            startat = ncol - nrow + 1

    # Trim the matrix.
    if header:
        # faster to do del
        del(items[0])
        nrow -= 1
    if startat > 0:
        items = [line[startat:] for line in items]
        ncol -= startat

    if ncol != nrow: raise AutoReadFail('non square matrix')

    # Time to get the numeric values.
    try:
        items = [[int(a) for a in line] for line in items]
    except ValueError:
        try:
            # dekker data 2009, uses integer but puts a comma... 
            items = [[int(float(a)) for a in line] for line in items]
            warn('WARNING: non integer values')
        except ValueError:
            raise AutoReadFail('ERROR: non numerical values')

    if is_asymmetric(items):
        warn('input matrix not symmetric: symmetrizing')
        symmetrize(items)

    # it is much faster to transform a list to a tuple than an iterator to a tuple..
    return tuple([a for line in items for a in line]), ncol


def read_matrix(things, parser=None):
    """
    Read and checks a matrix from a file or a list.

    :param things: might be either a file name, a file handler, a list of them
        or a list of list (all with same length)
    :param None parser: a parser function that returns a tuple of lists representing the data matrix,
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


    :returns: the corresponding matrix concatenated into a huge list, also
        returns number or rows

    """
    parser = parser or autoreader
    if type(things) is not list:
        things = [things]
    matrices = []
    sizes    = []
    for thing in things:
        if type(thing) is file:
            matrix, size = parser(thing)
            thing.close()
            matrices.append(matrix)
            sizes.append(size)
        elif type(thing) is str:
            try:
                matrix, size = parser(open(thing))
            except IOError:
                if len(thing.split('\n')) > 1:
                    matrix, size = parser(thing.split('\n'))
                else:
                    raise Exception('\n   ERROR: file %s not found\n' % thing)
            matrices.append(matrix)
            sizes.append(size)
        elif type(thing) is list:
            if all([len(thing)==len(l) for l in thing]):
                matrices.append(reduce(lambda x, y: x+y, thing))
                sizes.append(len(thing))
            else:
                raise Exception('must be list of lists, all with same length.')
        elif type(thing) is tuple:
            # case we know what we are doing and passing directly list of tuples
            matrices.append(thing)
            siz = sqrt(len(thing))
            if int(siz) != siz:
                raise AttributeError('ERROR: matrix should be square.\n')
            sizes.append(int(siz))
        elif 'matrix' in str(type(thing)):
            try:
                row, col = thing.shape
                if row != col:
                    raise Exception('matrix needs to be square.')
                matrices.append(thing.reshape(-1).tolist()[0])
                sizes.append(row)
            except Exception as exc:
                print 'Error found:', exc
        else:
            raise Exception('Unable to read this file or whatever it is :)')
    if all([s==sizes[0] for s in sizes]):
        return matrices, sizes[0]
    raise Exception('All matrices must have the same size ' +
                    '(same chromosome and same bins).')

