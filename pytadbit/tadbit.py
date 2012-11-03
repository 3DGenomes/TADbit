"""
24 Oct 2012


"""

from os import path, listdir
from pytadbit.tadbit_py import _tadbit_wrapper


def _read_matrix(f_h):
    """
    reads from file
    """
    nums = []
    for line in f_h:
        values = line.split()
        try:
            nums.append([float(v) for v in values[1:]])
        except ValueError:
            size = len(values)
            continue
    return reduce(lambda x,y: x+y, zip(*nums)), size


def read_matrix(things):
    """
    reads and checks matrix from file or list

    :argument things: might be either a file name, a file handler, a list of
    them or a list of list (all with same length)

    :returns: the corresponding matrix concatenated into a huge list, also
    returns number or rows

    TODO: find something better without using numpy matrices
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

def tadbit(x, n_cpus=None, verbose=True, speed=0, heuristic=True):
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
    x might be either a list of list, a file or a file handler
    :argument None n_cpus: The number of CPUs to allocate to tadbit. The value default
    is the total number of CPUs minus 1.
    :argument 0 speed: can be 0, 2, 3 or 4. Divide the calculated area by 2**(speed-1)
    :argument True heuristic: whether to use or not some heuristics

    :returns: the list of topologically associated domains' boundaries, and the
    corresponding list associated log likelihoods.
    """
    nums, size = read_matrix(x)
    del(x)
    print 'running {0} matrices of length {1}'.format(len(nums), size)
    _, nbks, likmat, _, bkpts = _tadbit_wrapper(nums,          # list of big lists representing the matrices
                                                size,          # size of one row/column
                                                len(nums),     # number of matrices
                                                n_cpus or 0,   # number of threads
                                                int(verbose),  # verbose 0/1
                                                speed,         # speed
                                                int(heuristic) # heuristic 0/1
                                                )

    dim     = nbks*size
    breaks  = [i for i in xrange(size) if bkpts[i+dim]==1]
    llikmat = [likmat[i*size:i*size+size] for i in xrange(size)]
    start   = [0]+[b+1 for b in breaks]
    end     = breaks[:]+[size-1]
    scores  = [llikmat[end[i  ]][start[i  ]] + \
               llikmat[end[i+1]][start[i+1]] - \
               llikmat[end[i+1]][start[i  ]] for i in xrange(len(breaks))]
    
    return breaks, scores


def batch_tadbit(directory, sep='_', **kwargs):
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

    :returns: A \code{list} where each element has the name of the unit/chromosome,
    and is the output of \code{tadbit} run on the corresponding files
    assumed to be replicates.

    TODO: understand what is 'sep' for...
    
    """

    matrix = []
    for f_name in listdir(directory):
        if f_name.startswith('.'): continue
        f_name = path.join(directory, f_name)
        if not path.isfile(f_name):
            continue
        print 'loading file:', f_name
        matrix.append(f_name)
    return tadbit(matrix, **kwargs)
