from pytadbit import Chromosome, Experiment
from pytadbit.modelling.impoptimizer import IMPoptimizer

# 1 - Load the Hi-C maps at each time-point
PATH = "./input_matrices/"

nrm_files = [PATH+"nrm_Sox2_B.txt" ,PATH+"nrm_Sox2_Ba.txt",
             PATH+"nrm_Sox2_D2.txt",PATH+"nrm_Sox2_D4.txt",
             PATH+"nrm_Sox2_D6.txt",PATH+"nrm_Sox2_D8.txt",
             PATH+"nrm_Sox2_PSC.txt"]
raw_files = [PATH+"raw_Sox2_B.txt" ,PATH+"raw_Sox2_Ba.txt",
             PATH+"raw_Sox2_D2.txt",PATH+"raw_Sox2_D4.txt",
             PATH+"raw_Sox2_D6.txt",PATH+"raw_Sox2_D8.txt",
             PATH+"raw_Sox2_PSC.txt"]


test_Sox2 = Chromosome(name="Test Sox2")
test_Sox2.add_experiment("Sox2", 5000, 
                         norm_data=nrm_files,
                         hic_data=raw_files,
                         silent=True)
exp = test_Sox2.experiments[0]
for timepoint in xrange(len(nrm_files)):
    exp.filter_columns(silent=True, index=timepoint)

# 2 - Produce the TADbit models at stage 0 (B cells) using the 
#optimizer = IMPoptimizer(exp, 1, 300, n_models=500, n_keep=100,  tool='imp', index=0, tmp_folder='./tmp/')
#optimizer.run_grid_search(n_cpus=16, lowfreq_range=(-3.0, 0.0, 1.0), upfreq_range=(0.0, 3.0, 1.0), maxdist_range=(150,400,50), verbose=1)
#optimizer.write_result('results.log')

# 3 - Produce the TADdyn trajectory
# Optimal parameters:
# lowfreq=-1.0
# upfreq=1.0
# maxdist=6.0

for nparticles in [50]:

    models = exp.model_region(1, nparticles, n_models=500, n_keep=100,
                              n_cpus=16, cleanup=False, hide_log=False,
                              initial_conformation='tadbit',                          
                              timesteps_per_k=10000, stages=[0,1,2,3,4,5,6],
                              config={'scale'  : 0.01 , 'kbending': 0.0,
                                      'maxdist': 300  , 'upfreq'  : 1.0, 
                                      'lowfreq': -1.0},
                              tool='lammps',
                              tmp_folder='./TADdyn_on_Sox2_test_%sparticles/' % nparticles,
                              timeout_job=12000)

    models.save_models("TADdyn_on_Sox2_test_%sparticles.pickle" % nparticles)
