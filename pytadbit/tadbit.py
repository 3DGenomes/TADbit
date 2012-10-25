"""
24 Oct 2012


"""

import ctypes as ct

# TODO: fins a better way to do this:
from pytadbit import __path__
_pytadbit = ct.CDLL(__path__[0]+'/_pytadbit.so')


def c_matrix(c_type, mat):
    """Make a C matrix from a list of lists (mat)
    picked from somewhere in the web...
    """
    row_type = c_type * len(mat[0])
    mat_type = ct.POINTER(c_type) * len(mat)
    mat = mat_type(*[row_type(*row) for row in mat])
    return ct.cast(mat, ct.POINTER(ct.POINTER(c_type)))


class tadbit_output(ct.Structure):
    """
    this corresponds to the type definition of tadbit_output structure (tadbit.h)
    this object can store tadbit output
    """
    _fields_ = [("maxbreaks", ct.c_int),
                ("nbreaks_opt", ct.c_int),
                ("llikmat", ct.POINTER(ct.c_double)),
                ("mllik", ct.POINTER(ct.c_double)),
                ("bkpts", ct.POINTER(ct.c_int))]


def _read_matrix(f_h):
    """
    reads from file
    """
    nums = [[]]
    for line in f_h:
        values = line.split()
        try:
            nums[0] += [float(v) for v in values[1:]]
        except ValueError:
            size = len(values)
            continue
    return nums, size


def read_matrix(f_name):
    """
    reads and checks matrix from file or list

    :argument f_name: might be either a file name, a file handler or a list of list (all with same length)

    :returns: the corresponding matrix concatenated into a huge list, also returns number or rows
    """
    if type(f_name) is file:
        return _read_matrix(f_name)
    elif type(f_name) is str:
        return _read_matrix(open(f_name))
    elif type(f_name) is list:
        if all([len(f_name)==len(l) for l in f_name]):
            return reduce(lambda x, y: x+y, f_name), len(f_name)
        else:
            raise Exception('must be list of lists, all with same length.')
    elif 'matrix' in str(type(f_name)):
        try:
            r, c = f_name.shape
            if r!=c:
                raise Exception('matrix needs to be square.')
            return f_name.reshape(-1).tolist(), r
        except Exception as e:
            print 'Error found:', e
            pass
    raise Exception('Unable to read this file or whatever it is :)')

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
    :argument 0 speed: TBA :)
    :argument True heuristic: whether to use or not some heuristics

    :returns: the list of topologically associated domains' boundaries, and the
    corresponding list associated log likelihoods.
    """
    nums, size = read_matrix(x)
    tbo = tadbit_output()
    _pytadbit.tadbit(c_matrix(ct.c_double, nums),   # list of big lists representing the matrices
                     size,                          # size of one row/column
                     1,                             # number of matrices
                     n_cpus or 1,                   # number of threads
                     verbose,                       # verbose 0/1
                     speed,                         # speed
                     int(heuristic),                # heuristic 0/1
                     ct.POINTER(tadbit_output)(tbo) # pointer to tbo object used to store output
                     )

    dim     = tbo.nbreaks_opt*size
    breaks  = [i for i in xrange(size) if tbo.bkpts[i+dim]==1]
    llikmat = [tbo.llikmat[i*size:i*size+size] for i in xrange(size)]
    start   = [0]+[b+1 for b in breaks]
    end     = breaks[:]+[size-1]
    scores  = [llikmat[end[i  ]][start[i  ]] + \
               llikmat[end[i+1]][start[i+1]] - \
               llikmat[end[i+1]][start[i  ]] for i in xrange(len(breaks))]
    
    return breaks, scores
