"""
04 Jul 2013
"""
from pytadbit import load_chromosome, Chromosome
from matplotlib import pyplot as plt
crm  = '2R'
xnam = 'TR1'
crmbit=Chromosome('2R')
crmbit.add_experiment(xnam, resolution=10000, xp_handler='/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
crmbit = load_chromosome('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/chr{0}.tdb'.format(crm))

# crmbit.experiments[xnam].load_experiment('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
# exp = crmbit.experiments[xnam]

for xnam in ['TR2', 'TR1', 'BR']:
    crmbit.experiments[xnam].load_experiment('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))
exp = crmbit.experiments['TR1'] + crmbit.experiments['TR2'] + crmbit.experiments['BR']


models = exp.model_region(start=190, end=295, n_models=500, n_keep=100,
                          n_cpus=8, verbose=False, keep_all=True)


models.write_cmm(0, '.')

models.cluster_models(dcutoff=200)

models.cluster_analysis_dendrogram(n_best_clusters=10)

models.model_consistency()


md1 = models.models[0]
md2 = models.models[2]

for m in range(1, 100):
    m = models.fetch_model_by_rand_init(m, all_models=True)
    # models.write_xyz(m, '.')
    write_xyz(models,m, '/home/fransua/scratch/dbau_test/lalal/')


models.objective_function(5)
models.write_xyz(5, '.')

models.write_cmm(589, '.')





def main():
    """
    main function
    """

    crm  = '2R'
    xnam = 'TR2'
    crmbit = load_chromosome('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/chr{0}.tdb'.format(crm))
    exp = crmbit.experiments[xnam]
    exp.load_experiment('/home/fransua/db/hi-c/corces_dmel/10Kb/{0}/{0}_{1}_10Kb.txt'.format(crm, xnam))

    matrix = exp.get_hic_matrix()

    new_matrix = [[] for _ in range(105)]

    for i in xrange(190, 295):
        for j in xrange(190, 295):
            new_matrix[i-190].append(matrix[i][j])

    crmbit = Chromosome('lala')
    crmbit.add_experiment('exp1', xp_handler=[new_matrix], resolution=10000)

    # crmbit.visualize('exp1', show=True)

    exp = crmbit.experiments[0]

    matrix = exp.get_hic_matrix()

    exp.normalize_hic(method='bytot')
    exp.get_hic_zscores()

    # plt.hist(reduce(lambda x, y: x+ y,
    #                 [v.values() for v in exp._zscores.values()]), bins=50)
    # plt.show()


    # here we are.
    # TODO: all this until here should be done by one simple matrix = exp.get_region()


    from pytadbit.imp.imp_model import IMPmodel
    
    model = IMPmodel('lala', exp._zscores, 10000)
    

    out = open('today.xyz', 'w')
    for i in exp._zscores:
        for j in exp._zscores[i]:
            out.write('{}\t{}\t{}\n'.format(i, j, exp._zscores[i][j]))
    out.close()


    

if __name__ == "__main__":
    exit(main())
