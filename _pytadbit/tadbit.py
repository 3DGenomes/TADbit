"""
24 Oct 2012


"""

from os                           import path, listdir
from pytadbit.parsers.hic_parser  import read_matrix, HiC_data
from pytadbit.tadbit_py           import _tadbit_wrapper
from pytadbit.tadbitalone_py      import _tadbitalone_wrapper
from pytadbit.utils.normalize_hic import iterative
from sys                          import stderr

def tadbit(x, weights=None, remove=None, n_cpus=1, verbose=True,
           max_tad_size="max", no_heuristic=0, **kwargs):
    """
    The TADbit algorithm works on raw chromosome interaction count data.
    The normalization is neither necessary nor recommended,
    since the data is assumed to be discrete counts.
    
    TADbit is a breakpoint detection algorithm that returns the optimal
    segmentation of the chromosome under BIC-penalized likelihood. The
    model assumes that counts have a Poisson distribution and that the
    expected value of the counts decreases like a power-law with the
    linear distance on the chromosome. This expected value of the counts
    at position (i,j) is corrected by the counts at diagonal positions
    (i,i) and (j,j). This normalizes for different restriction enzynme
    site densities and 'mappability' of the reads in case a bin contains
    repeated regions.

    :param x: a square matrix of interaction counts in the HI-C data or a list
       of such matrices for replicated experiments. The counts must be evenly
       sampled and not normalized. x might be either a list of list, a path to
       a file or a file handler
    :argument 'visibility' norm: kind of normalization to use. Choose between
       'visibility' of 'Imakaev'
    :argument None remove: a python list of lists of booleans mapping positively
       columns to remove (if None only columns with a 0 in the diagonal will be
       removed)
    :param 1 n_cpus: The number of CPUs to allocate to TADbit. If
       n_cpus='max' the total number of CPUs will be used
    :param auto max_tad_size: an integer defining maximum size of TAD. Default
       (auto) defines it as the number of rows/columns
    :param False no_heuristic: whether to use or not some heuristics
    :param False get_weights: either to return the weights corresponding to the
       Hi-C count (weights are a normalization dependent of the count of each
       columns)

    :returns: the :py:func:`list` of topologically associated domains'
       boundaries, and the corresponding list associated log likelihoods.
       If no weights are given, it may also return calculated weights.
    """
    nums = [hic_data for hic_data in read_matrix(x)]
    size = len(nums[0])
    if not remove:
        remove = tuple([int(nums[0][i+i*size]==0) for i in xrange(size)])
    if not weights:
        weights = []
        for num in nums:
            B = iterative(num, remove)
            weights.append(tuple([B[i]*B[j] for i in xrange(size)
                                  for j in xrange(size)]))
    nums = [num.get_as_tuple() for num in nums]
    print nums
    n_cpus = n_cpus if n_cpus != 'max' else 0
    max_tad_size = size if max_tad_size in ["max",
                                            "auto"] else max_tad_size
    _, nbks, passages, _, _, bkpts = \
       _tadbit_wrapper(nums,             # list of lists of Hi-C data
                       remove,           # list of columns marking filtered
                       weights,
                       size,             # size of one row/column
                       len(nums),        # number of matrices
                       n_cpus,           # number of threads
                       int(verbose),     # verbose 0/1
                       max_tad_size,     # max_tad_size
                       kwargs.get('ntads', 0),
                       int(no_heuristic),# heuristic 0/1
                       )

    breaks = [i for i in xrange(size) if bkpts[i + nbks * size] == 1]
    scores = [p for p in passages if p > 0]
    print scores

    result = {'start': [], 'end'  : [], 'score': []}
    for brk in xrange(len(breaks)+1):
        result['start'].append((breaks[brk-1] + 1) if brk > 0 else 0)
        result['end'  ].append(breaks[brk] if brk < len(breaks) else size - 1)
        result['score'].append(scores[brk] if brk < len(breaks) else None)

    return result


def batch_tadbit(directory, parser=None, **kwargs):
    """
    Use tadbit on directories of data files.
    All files in the specified directory will be considered data file. The
    presence of non data files will cause the function to either crash or
    produce aberrant results.
    
    Each file has to contain the data for a single unit/chromosome. The
    files can be separated in sub-directories corresponding to single
    experiments or any other organization. Data files that should be
    considered replicates have to start with the same characters, until
    the character sep. For instance, all replicates of the unit
    'chr1' should start with 'chr1\_', using the default value of sep.
    
    The data files are read through read.delim. You can pass options
    to read.delim through the list read_options. For instance
    if the files have no header, use read_options=list(header=FALSE) and if
    they also have row names, read_options=list(header=FALSE, row.names=1).
    
    Other arguments such as max_size, n_CPU and verbose are passed to
    :func:`tadbit`.
  
    NOTE: only used externally, not from Chromosome
    
    :param directory: the directory containing the data files
    :param kwargs: arguments passed to :func:`tadbit` function
    :param None parser: a parser function that takes file name as input and
        returns a tuple representing the matrix of data. Tuple is a
        concatenation of column1 + column2 + column3 + ...

    :returns: A :py:func:`list` where each element has the name of the
        unit/chromosome, and is the output of :func:`tadbit` run on the
        corresponding files assumed to be replicates

    """

    matrix = []
    for f_name in sorted(listdir(directory)):
        if f_name.startswith('.'):
            continue
        f_name = path.join(directory, f_name)
        if parser:
            matrix.append(parser(f_name))
            continue
        elif not path.isfile(f_name):
            continue
        matrix.append(f_name)
    return tadbit(matrix, **kwargs)


def print_result_r(result, write=True):
    """
    Print a table summarizing the TADs found by tadbit. This function outputs
    something similar to the R function.

    :param result: the :py:class:`dict` that returns :func:`tadbit`
    :param True write: print table. If False, returns the string

    :returns: if write is False, returns a string corresponding to the table of
       results
    """
    table = ''
    table += '%-6s%6s%6s%6s\n' % ('#', 'start', 'end', 'score')
    for i in xrange(len(result['end'])):
        table += '%-6s%6s%6s%6s\n' % (i+1, result['start'][i]+1,
                                      result['end'][i]+1,
                                      result['score'][i])
    if write:
        print table
    else:
        return table
