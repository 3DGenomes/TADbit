"""
26 mai 2014


"""



from tadbit_py import _tadbit_wrapper
from random import random


hic = []
for line in open('test/data/hESC_chr19-rep1.txt').readlines()[1000:1100]:
    hic.extend([int(i) for i in line.split()] [1000:1100])

print _tadbit_wrapper([tuple(hic)], [tuple([random() for _ in xrange(100**2)])],
                      1, 100, 1, 1, 1, 100, 0, 0)
