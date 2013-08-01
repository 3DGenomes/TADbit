"""
04 Jul 2013
"""
from pytadbit import Chromosome
from optparse   import OptionParser


def main():

    opts = get_options()
    
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

    for i in range(16):
        models.write_cmm(i, '.')

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
    opts = parser.parse_args()[0]
    if not opts.incrm:
        exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
