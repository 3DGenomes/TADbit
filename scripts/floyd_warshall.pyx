"""
An efficient Cython implementation of the Floyd-Warshall algorithm for finding the shortest path distances between all nodes of a weighted graph.
See http://en.wikipedia.org/wiki/Floyd-Warshall_algorithm

Amit Moscovich Eiger, 2014
"""
import cython
from cython.parallel cimport prange, parallel
cimport numpy
import numpy

def floyd_warshall_single_core(adjacency_matrix):
    '''floyd_warshall_single_core(adjacency_matrix) -> shortest_path_distance_matrix

    Input
        An NxN NumPy array describing the directed distances between N nodes.

        adjacency_matrix[i,j] = distance to travel directly from node i to node j (without passing through other nodes)

        Notes:
        * If there is no edge connecting i->j then adjacency_matrix[i,j] should be equal to numpy.inf.
        * The diagonal of adjacency_matrix should be zero.

    Output
        An NxN NumPy array such that result[i,j] is the shortest distance to travel between node i and node j. If no such path exists then result[i,j] == numpy.inf
    '''
    (nrows, ncols) = adjacency_matrix.shape
    assert nrows == ncols
    cdef unsigned int n = nrows

    adj_mat_copy = adjacency_matrix.astype(numpy.double, order='C', casting='safe', copy=True)
    assert adj_mat_copy.flags['C_CONTIGUOUS']
    cdef numpy.ndarray[numpy.double_t, ndim=2, mode='c'] M = adj_mat_copy
    assert (numpy.diagonal(M) == 0.0).all()

    cdef unsigned int i, j, k
    cdef double M_ij, M_ik, cost_ikkj
    cdef double* M_ptr = &M[0,0]
    cdef double* M_i_ptr
    cdef double* M_k_ptr

    for k in range(n):
        M_k_ptr = M_ptr + n*k
        for i in range(n):
            M_i_ptr = M_ptr + n*i
            M_ik = M_i_ptr[k]
            for j in range(n):
                cost_ikkj = M_ik + M_k_ptr[j]
                M_ij = M_i_ptr[j]
                if M_ij > cost_ikkj:
                    M_i_ptr[j] = cost_ikkj

    return M


# Fastest for large matrices!
def floyd_warshall_parallelized(adjacency_matrix):
    """floyd_warshall_parallelized(adjacency_matrix) -> shortest_path_distance_matrix

    Input
        An NxN NumPy array describing the directed distances between N nodes.

        adjacency_matrix[i,j] = distance to travel directly from node i to node j (without passing through other nodes)

        Notes:
        * If there is no edge connecting i->j then adjacency_matrix[i,j] should be equal to numpy.inf.
        * The diagonal of adjacency_matrix should be zero.

    Output
        An NxN NumPy array such that result[i,j] is the shortest distance to travel between node i and node j. If no such path exists then result[i,j] == numpy.inf

    This function uses Cython's integrated OpenMP multithreaded capabilities via the prange() function.
    """
    (nrows, ncols) = adjacency_matrix.shape
    assert nrows == ncols
    cdef unsigned int n = nrows

    adj_mat_copy = adjacency_matrix.astype(numpy.double, order='C', casting='safe', copy=True)
    assert adj_mat_copy.flags['C_CONTIGUOUS']
    cdef numpy.ndarray[numpy.double_t, ndim=2, mode='c'] M = adj_mat_copy
    assert (numpy.diagonal(M) == 0.0).all()

    cdef double cost_ik, cost_ikkj
    cdef unsigned int i, j, k, w

    cdef double* M_ptr = &M[0,0]
    cdef double* M_i_ptr
    cdef double* M_k_ptr
    cdef double M_ij

    with nogil, parallel():
        for k in range(n):
            M_k_ptr = M_ptr + n*k
            for i in prange(n):
                M_i_ptr = M_ptr + n*i
                cost_ik = M_i_ptr[k]
                for j in range(n):
                    cost_ikkj = cost_ik + M_k_ptr[j]
                    M_ij = M_i_ptr[j]
                    if cost_ikkj < M_ij:
                        M_i_ptr[j] = cost_ikkj
    return M

