"""
25 Oct 2012


"""

from pytadbit import tadbit, batch_tadbit
from sys import argv
from matplotlib import pyplot as plt
from numpy import log2, matrix


def get_matrix(f_name):
    nums = []
    for line in open(f_name):
        values = line.split()
        try:
            nums.append([float(v) for v in values[1:]])
        except ValueError:
            continue
    return matrix(nums)

def main():
    """
    main function
    """
    chrom =  argv[0]
    chrom = '/home/fransua/Tools/tadbit/pytadbit/test/test_a.tsv'
    chrom = get_matrix(chrom)

    out = tadbit(chrom, speed=4, verbose=True, heuristic=False)

    plt.imshow(log2(chrom.T), origin='lower')
    plt.vlines(out[0], 0, chrom.shape[0])
    plt.hlines(out[0], 0, chrom.shape[0])
    plt.show()

    chrom_path = '../../db/chr21/'
    out_batch = batch_tadbit(chrom_path, speed=4, heuristic=True)
    plt.imshow(log2(chrom.T), origin='lower')
    plt.vlines(out[0], 0, chrom.shape[0])
    plt.hlines(out[0], 0, chrom.shape[0])
    plt.vlines(out_batch[0], 0, chrom.shape[0], color='red')
    plt.hlines(out_batch[0], 0, chrom.shape[0], color='red')
    plt.show()
    

if __name__ == "__main__":
    exit(main())
