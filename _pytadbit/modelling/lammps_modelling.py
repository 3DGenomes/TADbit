"""
07 Nov 2016


"""
from pytadbit.modelling import LAMMPS_CONFIG as CONFIG
from pytadbit.modelling.lammpsmodel import LAMMPSmodel
from pytadbit.modelling.structuralmodels import StructuralModels
from os.path import exists
from random import randint, seed
from cPickle import load, dump
import multiprocessing as mu
import random 
from math import abs
import numpy as np

from lammps import lammps

import os

def init_lammps_run(lmp, initial_conformation,
                    neighbor=CONFIG.neighbor, minimize=True):
    
    """
    Initialise the parameters for the computation in lammps job
    
    :param lmp: lammps instance object.
    :param initial_conformation: lammps input data file with the particles initial conformation.
    :param CONFIG.neighbor neighbor: see LAMMPS_CONFIG.py.
    :param True minimize: whether to apply minimize command or not.

    """
    
    #######################################################
    # Box and units  (use LJ units and period boundaries) #
    #######################################################
    lmp.command("units %s" % CONFIG.units)
    lmp.command("atom_style %s" % CONFIG.atom_style) #with stiffness
    lmp.command("boundary %s" % CONFIG.boundary)
    
    ##########################
    # READ "start" data file #
    ##########################
    lmp.command("read_data %s" % initial_conformation)
    lmp.command("mass %s" % CONFIG.mass)
        
    ##################################################################
    # Pair interactions require lists of neighbours to be calculated #
    ##################################################################
    lmp.command("neighbor %s" % neighbor)
    lmp.command("neigh_modify %s" % CONFIG.neigh_modify)
    
    ##############################################################
    # Sample thermodynamic info  (temperature, energy, pressure) #
    ##############################################################
    lmp.command("thermo %i" % CONFIG.thermo)
    lmp.command("thermo_style   custom   step temp epair emol")
    
    ###############################
    # Stiffness term              #
    # E = K * (1+cos(theta)), K>0 #
    ###############################
    lmp.command("angle_style %s" % CONFIG.angle_style) # Write function for kinks 
    lmp.command("angle_coeff * %f" % CONFIG.persistence_length)
    
    ###################################################################
    # Pair interaction between non-bonded atoms                       #
    #                                                                 #
    #  Lennard-Jones 12-6 potential with cutoff:                      #
    #  potential E=4epsilon[ (sigma/r)^12 - (sigma/r)^6]  for r<r_cut #
    #  r_cut =1.12246 = 2^(1/6) is the minimum of the potential       #
    ###################################################################
    lmp.command("pair_style lj/cut %f" % CONFIG.PurelyRepulsiveLJcutoff)
    
    ################################################################
    #  pair_modify shift yes adds a constant to the potential such #
    #  that E(r_cut)=0. Forces remains unchanged.                  #
    ################################################################
    lmp.command("pair_modify shift yes")
    
    ######################################
    #  pair_coeff for lj/cut, specify 4: #
    #    * atom type interacting with    #
    #    * atom type                     #
    #    * epsilon (energy units)        #
    #    * sigma (distance units)        #
    ######################################
    lmp.command("pair_coeff * * %f %f %f" % (CONFIG.PurelyRepulsiveLJepsilon, CONFIG.PurelyRepulsiveLJsigma, CONFIG.PurelyRepulsiveLJcutoff))
    
    #########################################################
    # Pair interaction between bonded atoms                 #
    #                                                       #
    # Fene potential + Lennard Jones 12-6:                  #
    #  E= - 0.5 K R0^2 ln[ 1- (r/R0)^2]                     #
    #     + 4epsilon[ (sigma/r)^12 - (sigma/r)^6] + epsilon #
    #########################################################
    lmp.command("bond_style fene")
    
    ########################################
    # For style fene, specify:             #
    #   * bond type                        #
    #   * K (energy/distance^2)            #
    #   * R0 (distance)                    #
    #   * epsilon (energy)  (LJ component) #
    #   * sigma (distance)  (LJ component) #
    ########################################
    lmp.command("bond_coeff * %f %f %f %f" % (CONFIG.FENEK, CONFIG.FENER0, CONFIG.FENEepsilon, CONFIG.FENEsigma))
    
    lmp.command("special_bonds fene") #<=== I M P O R T A N T (new command)
    
    ##############################
    # set timestep of integrator #
    ##############################
    lmp.command("timestep %f" % CONFIG.timestep)
    
    if minimize:
        lmp.command("minimize 1.0e-4 1.0e-6 100000 100000")
        
def lammps_simulate(initial_conformation, run_time, colvars=None,
                    initial_seed=None, n_models=500,
                    resolution=10000, description=None,
                    neighbor=CONFIG.neighbor, tethering=False, 
                    minimize=True, keep_restart_step=1000000, 
                    keep_restart_out_dir=None, outfile=None, n_cpus=1):

    """
    This function launches jobs to generate three-dimensional models in lammps
    
    :param initial_conformation: lammps input data file with the particles initial conformation. http://lammps.sandia.gov/doc/2001/data_format.html
    :param run_time: # of timesteps.
    :param None colvars: space-separated input file with particles contacts http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf.
            Should at least contain Chromosome, loci1, loci2 as 1st, 2nd and 3rd column 
    :param None initial_seed: Initial random seed. If None then computer time is taken.
    :param 500 n_models: number of models to generate.
    :param 10000 resolution: resolution to specify for the StructuralModels object.
    :param None description: description to specify for the StructuralModels object.
    :param CONFIG.neighbor neighbor: see LAMMPS_CONFIG.py.
    :param True minimize: whether to apply minimize command or not. 
    :param 1000000 keep_restart_step: step to recover stopped computation. To be implemented.
    :param None keep_restart_out_dir: recover stopped computation. To be implemented.
    :param None outfile: store result in outfile
    :param 1 n_cpus: number of CPUs to use.
    
    :returns: a StructuralModels object

    """
    
    #===========================================================================
    # ##########################################################
    # # Generate RESTART file, SPECIAL format, not a .txt file #
    # # Useful if simulation crashes                           #
    # ##########################################################
    # 
    # if keep_restart_out_dir:
    #     if not os.path.exists(keep_restart_out_dir):
    #         os.makedirs(keep_restart_out_dir)
    #     lmp.command("restart %i %s/relaxation_%i_*.restart" % (keep_restart_step, keep_restart_out_dir, kseed))
    #===========================================================================
        
    if initial_seed:
        seed(initial_seed)

    pool = mu.Pool(n_cpus)
    
    kseeds = []
    
    for k in xrange(n_models):
        kseeds.append(randint(1,100000))
        
    jobs = {}
    for k in kseeds:
 
        jobs[k] = pool.apply_async(run_lammps,
                                           args=(k, initial_conformation, run_time, colvars,
                                                neighbor, tethering, 
                                                minimize,))

    pool.close()
    pool.join()
        
    
    
    results = []
    
    for k in kseeds:
        results.append((k, jobs[k].get()))
         
    nloci = 0
    models = {}
    for i, (_, m) in enumerate(
        sorted(results, key=lambda x: x[1]['objfun'])):
        models[i] = m
        nloci = len(m['x'])
        
    if outfile:
        if exists(outfile):
            old_models = load(open(outfile))
        else:
            old_models = {}
        models.update(old_models)
        out = open(outfile, 'w')
        dump((models), out)
        out.close()
    else:
        return StructuralModels(
            nloci, models, [], resolution,description=description, zeros=tuple([1 for i in xrange(nloci)]))

def run_lammps(kseed, initial_conformation, run_time, colvars=None,
                    neighbor=CONFIG.neighbor, tethering=False, 
                    minimize=True):
    """
    Generates one lammps model
    
    :param kseed: Random number to identify the model.
    :param initial_conformation: lammps input data file with the particles initial conformation. http://lammps.sandia.gov/doc/2001/data_format.html
    :param run_time: # of timesteps.
    :param None colvars: particles contacts file from colvars fix http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf. 
    :param CONFIG.neighbor neighbor: see LAMMPS_CONFIG.py.
    :param False tethering: whether to apply tethering command or not.
    :param True minimize: whether to apply minimize command or not. 
    
    :returns: a LAMMPSModel object

    """
    
    lmp = lammps()
        
    
    init_lammps_run(lmp, initial_conformation,
                neighbor=neighbor, minimize=minimize)
    
    #######################################################
    # Set up fixes                                        #
    # use NVE ensemble                                    #
    # Langevin integrator Tstart Tstop 1/friction rndseed #
    # => sampling NVT ensamble                            #
    #######################################################
    lmp.command("fix 1 all nve")
    lmp.command("fix 2 all langevin 1.0  1.0  2.0 %i" % kseed)
    if tethering:
        lmp.command("fix 3 all spring tether 50.0 0.0 0.0 0.0 0.0")
    if colvars:
        lmp.command("fix 4 all colvars %s" % colvars)
    
    lmp.command("reset_timestep 0")
    
    lmp.command("run %i" % run_time)

    xc = np.array(lmp.gather_atoms("x",1,3))
    
    result = LAMMPSmodel({'x'          : [],
                       'y'          : [],
                       'z'          : [],
                       'cluster'    : 'Singleton',
                       'objfun'     : 0,
                       'radius'     : 1,
                       'rand_init'  : str(kseed)})
                       
    lmp.close()
    
    for i in xrange(0,len(xc),3):
        result['x'].append(xc[i])
        result['y'].append(xc[i+1])
        result['z'].append(xc[i+2])
    
    return result
 
def linecount(filename):
    """
    Count valid lines of input colvars file
    
    :param filename: input colvars file.
    
    :returns: number of valid contact lines

    """
    
    k = 0
    tfp = open(filename)
    for i, line in enumerate(tfp):   
        if line.startswith('#'):
            continue
        cols_vals = line.split(' ')
        if cols_vals[1] == cols_vals[2]:
            continue
        k += 1
        
    return k
    
def generate_colvars_list(target_pairs_file,outfile,kappa_vs_genomic_distance,
                          chrlength=0, copies=['A'],kbin=10000000, 
                          binsize=50000, percentage_enforced_contacts=10,
                          colvars_header='# collective variable: monitor distances\n\ncolvarsTrajFrequency 1000 # output every 1000 steps\ncolvarsRestartFrequency 10000000',
                          colvars_template='''

colvar {
  name %s
  # %s %s %i
  distance {
      group1 {
        atomNumbers %i
      }
      group2 {
        atomNumbers %i
      }
  }
}''',
                            colvars_tail = '''

harmonic {
  name h_pot_%i
  colvars %s
  centers %s
  forceConstant %f # %f
}'''
                            ):
                            
    """
    Generates lammps colvars file http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf
    
    :param target_pairs_file:space-separated input file with particles contacts http://lammps.sandia.gov/doc/PDF/colvars-refman-lammps.pdf.
            Should at least contain Chromosome, loci1, loci2 as 1st, 2nd and 3rd column.
    :param outfile: file where to store resulting colvars file
    :param kappa_vs_genomic_distance: space-separated file that relates genomic distance and corresponding spring constant kappa.
            Should at least contain Distance, kappa as 1st and 2nd column
    :param 0 chrlength: Chromosome lengths if more than one copy.
    :param ['A'] copies: list of chromosome copies.
    :param 10000000 kbin: bin size of kappa vs genomic distance distribution values.
    :param 50000 binsize: size of each bin in genomic distance.
    :param 10 percentage_enforced_contacts: Percentage of enforced contacts over the total existing in the target_pairs_file
    :param exisiting_template colvars_header: header template for colvars file.
    :param exisiting_template colvars_template: contact template for colvars file.
    :param exisiting_template colvars_tail: tail template for colvars file.

    """

    totcolvars = linecount(target_pairs_file)
    ncolvars = int(totcolvars*(float(percentage_enforced_contacts)/100))
    
    print "Number of enforced contacts = %i over %i" % (ncolvars,totcolvars)
    rand_positions = random.sample(list(range(totcolvars)), ncolvars)
    rand_positions = sorted(rand_positions)
    rand_lines = []

    i=0
    j=0
    tfp = open(target_pairs_file)
    with open(target_pairs_file) as f:
        for line in f:
            if j >= ncolvars:
                break
            if line.startswith('#'):
                continue
            cols_vals = line.split(' ')
            if cols_vals[1] == cols_vals[2]:
                continue
        
            if i == rand_positions[j]:
                rand_lines.append(line)
                j += 1
            i += 1
    tfp.close()
    
    seqdists = {}
    poffset=0
    outf = open(outfile,'w')
    outf.write(colvars_header)
    for copy_nbr in copies:
        i = 1
        for line in rand_lines:   
            cols_vals = line.split(' ')
            part1_start = int(cols_vals[1])*binsize
            part1_end = (int(cols_vals[1])+1)*binsize
            part2_start = int(cols_vals[2])*binsize
            part2_end = (int(cols_vals[2])+1)*binsize
            name = str(i)+copy_nbr  
            seqdist = int((int(abs((part1_start-part2_start))/kbin)+0.5)*kbin)
            region1 = cols_vals[0] + '_' + str(part1_start) + '_' + str(part1_end)
            region2 = cols_vals[0] + '_' + str(part2_start) + '_' + str(part2_end)
            particle1 = int(cols_vals[1]) + 1 + poffset
            particle2 = int(cols_vals[2]) + 1 + poffset
            seqdists[name] = seqdist
            outf.write(colvars_template % (name,region1,region2,seqdist,particle1,particle2))
            i += 1
        poffset += chrlength
            
    outf.flush()
    
    kappa_values = {}
    with open(kappa_vs_genomic_distance) as kgd:
        for line in kgd:
            line_vals = line.split(' ')
            kappa_values[int(line_vals[0])] = float(line_vals[1])
        
    for seqd in set(seqdists.values()):
        kappa = 0
        if seqd in kappa_values:
            kappa = kappa_values[seqd]*0.1
        else:
            for kappa_key in sorted(kappa_values, key=int):
                if int(kappa_key) > seqd:
                    break
                kappa = kappa_values[kappa_key]*0.1
        centres=''
        names=''
        for seq_name in seqdists:
            if seqdists[seq_name] == seqd:
                centres += ' 1.0'
                names += ' '+seq_name
      
        outf.write(colvars_tail % (seqd,names,centres,kappa,kappa))  
        
    outf.flush()
    
    outf.close()  
        
        
        
      
    
       
