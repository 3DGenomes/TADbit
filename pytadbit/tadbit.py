"""
24 Oct 2012


"""

from os import path, listdir
from math import sqrt
from pytadbit.tadbit_py import _tadbit_wrapper



def _read_matrix(f_h):
    """
    reads from file
    """
    nums = []
    while True:
        values = f_h.next().split()
        if values[0].startswith('#'):
            # skip comments
            continue
        break
    # check if we have headers/row-names in the file
    start = 1
    if values[0].isdigit():
        try:
            nums.append([int(v) for v in values])
        except ValueError:
            nums.append([int(float(v)) for v in values])
        start = 0
    # parse the rest of the file
    for line in f_h:
        values = line.split()[start:]
        try:
            nums.append([int(v) for v in values])
        except ValueError:
            nums.append([int(float(v)) for v in values])
    f_h.close()
    size = len(nums)
    return tuple([nums[j][i] for i in xrange(size) for j in xrange(size)]), size


def read_matrix(things):
    """
    reads and checks matrix from file or list

    :argument things: might be either a file name, a file handler, a list of
    them or a list of list (all with same length)

    :returns: the corresponding matrix concatenated into a huge list, also
    returns number or rows

    """
    if type(things) is not list:
        things = [things]
    matrices = []
    sizes    = []
    for thing in things:
        if type(thing) is file:
            matrix, size = _read_matrix(thing)
            matrices.append(matrix)
            sizes.append(size)
        elif type(thing) is str:
            matrix, size = _read_matrix(open(thing))
            matrices.append(matrix)
            sizes.append(size)
        elif type(thing) is list:
            if all([len(thing)==len(l) for l in thing]):
                matrices.append(reduce(lambda x, y: x+y, thing))
                sizes.append(len(thing))
            raise Exception('must be list of lists, all with same length.')
        elif type(thing) is tuple:
            # case we know what we are doing and passing directly list of tuples
            matrices.append(thing)
            siz = sqrt(len(thing))
            if int(siz) != siz:
                raise Exception('Error: matrix should be square.')
            sizes.append(int(siz))
        elif 'matrix' in str(type(thing)):
            try:
                r, c = thing.shape
                if r!=c:
                    raise Exception('matrix needs to be square.')
                matrices.append(thing.reshape(-1).tolist()[0])
                sizes.append(r)
            except Exception as e:
                print 'Error found:', e
                pass
        else:
            raise Exception('Unable to read this file or whatever it is :)')
    if all([s==sizes[0] for s in sizes]):
        return matrices, sizes[0]
    raise Exception('All matrices must have the same size (same chromosome and same bins).')

def tadbit(x, n_cpus=None, verbose=True, max_tad_size="auto", no_heuristic=False):
    """
    The tadbit algorithm works on raw chromosome interaction count data.
    Not only is normalization not necessary, it is also not recommended
    since the data is assumed to be discrete counts.
    
    Tadbit is a breakpoint detection algorithm that returns the optimal
    segmentation of the chromosome under BIC-penalized likelihood. The
    model assumes that counts have a Poisson distribution and that the
    expected value of the counts decreases like a power-law with the
    linear distance on the chromosome. This expected value of the counts
    at position (i,j) is corrected by the counts at diagonal positions
    (i,i) and (j,j). This normalizes for different restriction enzynme
    site densities and 'mappability' of the reads in case a bin contains
    repeated regions.

    :argument x: A square matrix of interaction counts in hi-C data or a list of
    such matrices for replicated experiments. The counts must be evenly sampled
    and not normalized.
    x might be either a list of list, a path to a file or a file handler
    :argument None n_cpus: The number of CPUs to allocate to tadbit. The value default
    is the total number of CPUs minus 1.
    :argument auto max_tad_size: an integer defining maximum size of TAD.
    Default (auto) defines it to the number of rows/columns.
    :argument False no_heuristic: whether to use or not some heuristics

    :returns: the list of topologically associated domains' boundaries, and the
    corresponding list associated log likelihoods.
    """
    nums, size = read_matrix(x)
    del(x)
    max_tad_size = size if max_tad_size is "auto" else max_tad_size
    _, nbks, passages, _, _, bkpts = _tadbit_wrapper(nums,             # list of big lists representing the matrices
                                                     size,             # size of one row/column
                                                     len(nums),        # number of matrices
                                                     n_cpus or 0,      # number of threads
                                                     int(verbose),     # verbose 0/1
                                                     max_tad_size,     # max_tad_size
                                                     int(no_heuristic) # heuristic 0/1
                                                     )

    dim     = nbks*size
    breaks  = [i for i in xrange(size) if bkpts[i+dim]==1]
    scores = [p for p in passages if p > 0]

    result = {'start': [], 'end'  : [], 'score': []}
    for i in xrange(len(breaks)+1):
        result['start'].append((breaks[i-1] + 1) if i > 0 else 0)
        result['end'  ].append(breaks[i] if i < len(breaks) else size - 1)
        result['score'].append(scores[i] if i < len(breaks) else None)
        
    return result


def batch_tadbit(directory, sep='_', parser=None, **kwargs):
    """
    Use tadbit on directories of data files
    All files in the specified directory will be considered data file. The
    presence of non data files will cause the function to either crash or
    produce aberrant results.
    
    Each file has to contain the data for a single unit/chromosome. The
    files can be separated in sub-directories corresponding to single
    experiments or any other organization. Data files that should be
    considered replicates have to start with the same characters, until
    the character sep. For instance, all replicates of the unit
    'chr1' should start with 'chr1_', using the default value of sep.
    
    The data files are read through read.delim. You can pass options
    to read.delim through the list read_options. For instance
    if the files have no header, use read_options=list(header=FALSE) and if
    they also have row names, \code{read_options=list(header=FALSE, row.names=1)}.
    
    Other arguments such as \code{max_size}, \code{n_CPU} and \code{verbose} are
    passed to \code{tadbit}.
  
    :argument directory: The directory containing the data files.
    :argument _ sep: A character specifying how to identify unit/chormosome names
    (see below).
    :argument kwargs: arguments passed to :py:func:pytadbit:`tadbit` function.
    :argument None parser: a parser function that takes file name as input and
    returns a tuple representing the matrix of data. Tuple is a concatenation of
    column1 + column2 + column3 + ...

    :returns: A \code{list} where each element has the name of the unit/chromosome,
    and is the output of \code{tadbit} run on the corresponding files
    assumed to be replicates.

    """

    matrix = []
    for f_name in listdir(directory):
        if f_name.startswith('.'): continue
        f_name = path.join(directory, f_name)
        if parser:
            print 'loading file:', f_name
            matrix.append(parser(f_name))
            continue
        elif not path.isfile(f_name):
            continue
        print 'loading file:', f_name
        matrix.append(f_name)
    return tadbit(matrix, **kwargs)


def print_result_R(result, sep=' '*6, write=True):
    """
    """
    table = ''
    table += '{:<6}{:>6}{:>6}{:>6}\n'.format('#','start','end','score')
    for i in xrange(len(result['end'])):
        table += '{:<6}{:>6}{:>6}{:>6}\n'.format(i+1,
                                                 result['start'][i]+1,
                                                 result['end'][i]+1,
                                                 result['score'][i])
    if write:
        print table
    else:
        return table
