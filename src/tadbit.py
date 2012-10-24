"""
24 Oct 2012


"""

import ctypes as ct
from random import random

def rnd_name():
    let = 'abcdefghijklmnopqrst'
    nam = ''
    for _ in range(5):
        nam += let[int(random()*len(let))]
    return nam

def c_matrix(c_type, mat):
    """Make a C matrix from a list of lists (mat)"""

    row_type = c_type * len(mat[0])
    mat_type = ct.POINTER(c_type) * len(mat)
    mat = mat_type(* [row_type(* row) for row in mat])
    return ct.cast(mat, ct.POINTER(ct.POINTER(c_type)))


class tadbit_output(ct.Structure):
    _fields_ = [("maxbreaks", ct.c_int),
                ("nbreaks_opt", ct.c_int),
                ("llikmat", ct.POINTER(ct.c_double)),
                ("mllik", ct.POINTER(ct.c_double)),
                ("bkpts", ct.POINTER(ct.c_int))]

def read_matrix(f_name):
    """
    returns a big list just as tadbit likes
    """
    nums = [[]]
    for line in open(f_name):
        values = line.split()
        try:
            nums[0] += [float(v) for v in values[1:]]
        except ValueError:
            size = len(values)
            continue
    return nums, size

tadbit = ct.CDLL('/home/fransua/Tools/tadbit/src/tadbit_py.so')

#size = 1000
#out = open('/home/fransua/Projects/tad-model_guillaume/db/chr21/test.tsv', 'w')
#inp = [[0 for j in range(0,size)] for i in range(0,size)]
#for i in range(size):
#    for j in range(size):
#        if j<i: continue
#        rnd = int(random()*30)
#        rnd = 0 if random() < 0.9 else rnd
#        inp[i][j] = float(rnd)
#        inp[j][i] = float(rnd)
#out.write('\t'.join([rnd_name() for j in range(len(inp))]) + '\n')
#for i in range(size):
#    out.write(rnd_name()+'\t'+'\t'.join([str(j) for j in inp[i]]) + '\n')
#out.close()
#

def main():
    """
    playing around
    TODO: need to convert tadbit output to something readable....
    """
        
    
    chr21 = '/home/fransua/Projects/tad-model_guillaume/db/chr21/chr21_t0_ncoi_40kb_matrix.txt'
    nums, size = read_matrix(chr21)
    tbo = tadbit_output()
    tadbit.tadbit(c_matrix(ct.c_double, nums),   # list of big lists representing the matrices
                  size,                          # size of one row/column
                  1,                             # number of matrices
                  0,                             # number of threads
                  1,                             # verbose 0/1
                  4,                             # speed
                  1,                             # heuristic 0/1
                  ct.POINTER(tadbit_output)(tbo) # 
                  )
