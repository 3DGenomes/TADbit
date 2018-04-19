
from pytadbit import Chromosome

import sys

def main():
    matrix_path   = sys.argv[1]
    config_string = sys.argv[2]
    compute_keep = sys.argv[3]

    uf, lf, md = config_string.split(':')
    lf = float(lf)
    uf = float(uf)
    md = int  (md)
    config = {'reference' : '', 'kforce'    : 5,
              'maxdist'   : md,
              'upfreq'    : uf,
              'lowfreq'   : lf,
              'scale'     : 0.01,
              'kbending'  : 0.0,
              }

    compute, keep = map(int, compute_keep.split(':'))

    chrom = Chromosome('chr')
    chrom.add_experiment('sample', norm_data=matrix_path, resolution=15000)
    exp = chrom.experiments[0]

    models = exp.model_region(n_models=compute, n_keep=keep, n_cpus=8, config=config)

    models.save_models('models_%s.pickle' % (config_string))

if __name__ == '__main__':
    exit(main())


