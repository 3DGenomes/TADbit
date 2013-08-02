"""
04 Jul 2013
"""
from pytadbit import Chromosome
from pytadbit.imp.imp_model import generate_3d_models
from optparse   import OptionParser

def parse_zscores(zscore_f):
    """
    """
    zscores = {}
    for line in open(zscore_f):
        x, y, z = line.split()
        z = float(z)
        x = x
        y = y
        zscores.setdefault(x, {})
        zscores[x][y] = z
        # 
        zscores.setdefault(y, {})
        zscores[y][x] = z
    return zscores
    

def main():

    opts = get_options()
    if opts.inabc:
        zscores = parse_zscores(opts.inabc)

    models = generate_3d_models(zscores, start=1, n_models=10, n_keep=10,
                                close_bins=1, n_cpus=8, keep_all=False, verbose=False,
                                outfile=None, resolution=10000)
    
    for i in range(10):
        models.write_cmm(i, '.')

    exit()
    
    crm  = '2R'
    xnam = 'TR1'
    crmbit=Chromosome('2R')
    # crmbit.add_experiment(xnam, resolution=10000, xp_handler='/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
    crmbit.add_experiment(xnam, resolution=10000, xp_handler=opts.incrm)
    # crmbit.experiments[xnam].load_experiment('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
    exp = crmbit.experiments[xnam]

    # for xnam in ['TR2', 'TR1', 'BR']:
    #     crmbit.experiments[xnam].load_experiment('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
    # exp = crmbit.experiments['TR1'] + crmbit.experiments['TR2'] + crmbit.experiments['BR']


    models = exp.model_region(start=190, end=295, n_models=5000, n_keep=1000,
                              n_cpus=8, verbose=False, keep_all=True)


    self = exp.model_region(start=190, end=295, n_models=10, n_keep=10,
                            n_cpus=8, verbose=True, keep_all=True)

    for i in range(10):
        models.write_cmm(i, '.')

    exit()
    
    models.cluster_models(dcutoff=200)

    models.cluster_analysis_dendrogram(n_best_clusters=10)

    models.model_consistency()


    md1 = models.models[0]
    md2 = models.models[2]

    for m in range(1, 100):
        m = models.fetch_model_by_rand_init(m, all_models=True)
        models.write_xyz(m, '.')


    models.objective_function(5)
    models.write_xyz(5, '.')

    models.write_cmm(589, '.')

    

def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        usage=("%prog [options] file [options] file [options] " +
               "file [options [file ...]]"))
    parser.add_option('--incrm', dest='incrm', metavar="PATH",
                      default=None,
                      help='''path to a hi-c matrix file.''')
    parser.add_option('--inabc', dest='inabc', metavar="PATH",
                      default=None,
                      help='''path to a zscore file.''')
    opts = parser.parse_args()[0]
    # if not opts.incrm:
    #     exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
