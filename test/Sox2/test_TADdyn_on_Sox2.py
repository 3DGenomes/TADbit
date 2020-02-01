from pytadbit import Chromosome, Experiment
from pytadbit.modelling.impoptimizer import IMPoptimizer
from pytadbit.modelling.structuralmodels import load_structuralmodels
from pytadbit.utils.three_dim_stats import calc_eqv_rmsd
from pytadbit.utils.file_handling import mkdir

from numpy import median, zeros, mean, dot

import matplotlib.pyplot as plt
plt.switch_backend('agg')

import mpl_toolkits.mplot3d as a3
import scipy as sp
from scipy.optimize import linprog

from scipy.spatial import ConvexHull

### Convex-hull calculation ###
def pnt_in_cvex_hull(hull, point, tolerance=1e-3):
    return all(
        (dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)
### Convex-hull calculation ###

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
optimizer = IMPoptimizer(exp, 1, 300, n_models=500, n_keep=100,  tool='imp', index=0, tmp_folder='./tmp/')
optimizer.run_grid_search(n_cpus=16, lowfreq_range=(-3.0, 0.0, 1.0), upfreq_range=(0.0, 3.0, 1.0), maxdist_range=(150,400,50), verbose=1)
optimizer.write_result('results.log')

# 3 - Produce the TADdyn trajectory
# Optimal parameters:
# lowfreq=-1.0
# upfreq=1.0
# maxdist=6.0

for nparticles in [300]:

    models = exp.model_region(1, nparticles, n_models=500, n_keep=100,
                              n_cpus=16, cleanup=False, hide_log=False,
                              initial_conformation='tadbit',                          
                              timesteps_per_k=10000, stages=[0,1,2,3,4,5,6],
                              config={'scale'  : 0.01 , 'kbending': 0.0,
                                      'maxdist': 300  , 'upfreq'  : 1.0, 
                                      'lowfreq': -1.0},
                              tool='lammps',
                              tmp_folder='./TADdyn_on_Sox2_test_%sparticles/' % nparticles,
                              timeout_job=90000, useColvars=True)

    models.save_models("TADdyn_on_Sox2_test_%sparticles.pickle" % nparticles)

# 4 - Analysis of the TADdyn models
print "Models analysis:"
print "0 - Loading the models"
sm = load_structuralmodels('TADdyn_on_Sox2_test_300particles.pickle')
mkdir("Results")


print "1 - Z-score plots"
# A. TADbit z-score plot per stage (Figure 1B)
mkdir("Results/Zscore_plots/")
for stage in [0, 1, 2, 3, 4, 5, 6]:
     sm.zscore_plot(stage=stage, savefig="./Results/Zscore_plots/Zscore_plot_at_stage_%d.pdf" % stage)
print "DONE"

print "2 - Snapshots"
# B. Model conformations visualization (Figure 1C)
mkdir("Results/Snapshots/")
for stage in [0, 1, 2, 3, 4, 5, 6]:
    sm.view_models(models=[int(stage*100*100)], tool="plot", savefig="./Results/Snapshots/Snapshot_trajectory_1_timepoint_%d.png" % int(stage*100))
print "DONE"


print "3 - Contact matrices"
# C. Contact matrices at 200nm (Figure 1D)
mkdir("Results/Contact_maps/")
for timestep in xrange(0,601):
    sm.contact_map(stage=timestep, cutoff=200., savefig="./Results/Contact_maps/Contact_map_at_timestep_%d.png" % timestep, cmap="jet")
print "DONE"


print "4 - Spearman correlation"
# D. Spearman correlation matrix between Hi-C and TADdyn models contact matrices (Figure 2A)
mkdir("Results/Correlation_with_real_data/")
fp_output=open("./Results/Correlation_with_real_data/correlation_with_real_data.tab", "w")
fp_output.write("%s\t%s\t%s\n" % ("#HiCStage", "timestep", "SpearmanCorr"))
for index in xrange(6,7): # Hi-C stage
    for timestep in xrange(601): # TADdyn timestep         
        fp_output.write("%s\t%s\t%lf\n" % (index,timestep,sm.correlate_with_real_data(stage=timestep, index=index, cutoff=200, savefig="./Results/Correlation_with_real_data/correlation_plot_HiC_stage_%s_TADdyn_timestep_%s.png" % (index, timestep), cmap="jet")))
print "DONE"


print "5 - Models dRMSD"
# E. Median dRMSD of TADdyn models (Figure 2B)
mkdir("Results/Models_dRMSDs")
# Compute all-vs-all dRMSD
nsteps = 10
models = [sm[mod] for x in xrange(0,601,nsteps) for mod in sm.stages[x]]
dRMSDs = calc_eqv_rmsd(models=models, nloci=300, zeros=[True]*300, what="dRMSD", normed=False)
fp_output=open("./Results/Models_dRMSDs/all_dRMSD_matrix.tab", "w")
all_dRMSD = {}
fp_output.write("%s\t%s\t%s\t%s\t%s\n" % ("#Trajectory_i", "timestep_i", "trajectory_j", "timestep_j", "dRMSD"))
for pair in dRMSDs:
    trajectory_i  = (int(pair[0])%100)+1
    timestep_i    = (int(pair[0])/100)*nsteps
    trajectory_j  = (int(pair[1])%100)+1
    timestep_j    = (int(pair[1])/100)*nsteps
    if (timestep_i,timestep_j) not in all_dRMSD:
        all_dRMSD[(timestep_i,timestep_j)] = []
    all_dRMSD[(timestep_i,timestep_j)].append(dRMSDs[pair])
    fp_output.write("%d\t%d\t%d\t%d\t%lf\n" % (trajectory_i, timestep_i, trajectory_j, timestep_j, dRMSDs[pair]))

# Compute median dRMSD
fp_output=open("./Results/Models_dRMSDs/median_dRMSD_matrix.tab", "w")
fp_output.write("%s\t%s\t%s\n" % ("#Timestep_i", "timestep_j", "dRMSD"))
for timestep_i in xrange(0,601,nsteps):
    for timestep_j in xrange(0,601,nsteps):
        fp_output.write("%d\t%d\t%lf\n" % (timestep_i, timestep_j, median(all_dRMSD[(timestep_i,timestep_j)])))
print "DONE"


print "6 - Accessibility"
# F. 1.0 - Accessibility (Figure 3A)
radius=50
mkdir("./Results/Accessibility_for_molecules_of_radius_%dnm" % radius)
for timestep in xrange(0,601,1):
    #Note: For the final Figures in Ref. (1) nump=100 was used. 
    sm.accessibility(radius=radius, nump=10, models=sm.stages[timestep], 
                     savedata="./Results/Accessibility_for_molecules_of_radius_%dnm/accessibility_timestep_%d.txt" % (radius, timestep),
                     savefig ="./Results/Accessibility_for_molecules_of_radius_%dnm/accessibility_timestep_%d.png" % (radius, timestep))
    print "DONE"

TSSparticle=150
print "7 - TSS particles spanned volume"
# G. Convexhull calculation volume for TSS particles spanned volume (Figure 3B)
mkdir("Results/Convexhull/")
ntimepoints=50
fp_output=open("./Results/Convexhull/convexhull_TSS_particles.txt", "w")
fp_output.write("%s\t%s\t%s\t%s\n" % ("#Trajectory","particle","timesteps","ConvexhullVolume"))

for trajectory in xrange(100):
    models = []
    for timestep in xrange(0,601,1):
        models.append(sm[sm.stages[timestep][trajectory]])

    for timestep_i in xrange(0,600,ntimepoints):        
        points = zeros((ntimepoints,3),dtype=float)
        for time in xrange(timestep_i,timestep_i+ntimepoints):
            points[time-timestep_i][0] = models[time]['x'][TSSparticle-1]
            points[time-timestep_i][1] = models[time]['y'][TSSparticle-1]
            points[time-timestep_i][2] = models[time]['z'][TSSparticle-1]
        hull = ConvexHull(points, incremental=True)
        fp_output.write("%d\t%d\t%d-%d\t%lf\n" % (trajectory+1,TSSparticle,timestep_i,timestep_i+ntimepoints,hull.volume))
print "DONE"


print "8 - Input for IS analysis"
# H. Insulation score analysis (Figure 4A)
for timestep in xrange(0,601):
    fp_output=open("./Results/Contact_maps/Contact_map_at_timestep_%d_for_IS_analysis.mat" % timestep, "w")
    matrix = sm.get_contact_matrix(cutoff=200,distance=False,models=sm.stages[timestep])
    
    # Write header
    for i in xrange(len(matrix)):
        fp_output.write("%d|mm9|chr3:%d-%d\t" % (i+1,5000*i,5000*(i+1)))
    fp_output.write("\n")                                                                                             
    for i in xrange(len(matrix)):
        fp_output.write("%d|mm9|chr3:%d-%d\t" % (i+1,5000*i,5000*(i+1)))
        for j in xrange(len(matrix[i])):
            fp_output.write("%f\t" % ((matrix[i][j]+matrix[i][j])*0.5))
        fp_output.write("\n")
#Perform the IS analysis using the script matrix2insulation.ol from XXX repository on the files Contact_map_at_timestep_%d.mat
#using the parameters perl matrix2insulation.pl -i ./Results/Contact_maps/Contact_map_at_timestep_${timestep}.mat -v -o timestep_${timestep} --is 100000 --ids 50000 --ez --im mean --nt 0.1 --bmoe 3 > /dev/null 2> /dev/null
print "DONE"


print "9 - Distances to TSS"
# I. TSS distance matrix (Figure 4B)
# Compute all TSS distances
mkdir("./Results/TSS_distance")
fp_output=open("./Results/TSS_distance/all_TSS_distances.txt", "w")
fp_output1=open("./Results/TSS_distance/mean_TSS_distances.txt", "w")
fp_output.write("%s\t%s\t%s\t%s\n" % ("#Trajectory", "timestep", "particle", "TSS_distance"))
fp_output1.write("%s\t%s\t%s\n" % ("#Timestep", "particle", "TSS_distance"))

for timestep in xrange(0,601,1):
    models = [sm[mod] for mod in sm.stages[timestep]]
    values = zeros((300,100), dtype="float")
    for trajectory,model in enumerate(models):        
        for particle in xrange(1,301):                    
            distance = model.distance(TSSparticle, particle)
            fp_output.write("%s\t%s\t%s\t%s\n" % (trajectory+1, timestep, particle, distance))
            values[particle-1][trajectory] = distance

    for particle in xrange(1,301):
        fp_output1.write("%d\t%d\t%lf\n" % (timestep, particle, mean(values[particle-1])))
print "DONE"


## L. Particle category analysis (Figure 4C)

