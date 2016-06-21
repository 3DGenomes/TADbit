"""
November 7, 2013.

"""

from warnings                import warn
from math                    import sqrt, isnan
from pytadbit.parsers.gzopen import gzopen
from collections             import OrderedDict
from pytadbit                import HiC_data

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
        for j in range(i, maxn):
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
        warn('WARNING: matrix not symmetric: summing cell_ij with cell_ji')
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
            matrix, size, header, masked = parser(thing)
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
                                     chromosomes=chromosomes, masked=masked,
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

