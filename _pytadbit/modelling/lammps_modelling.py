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
from math import atan2
import numpy as np
import sys
import copy
from numpy import sin, cos, arccos, sqrt, fabs, asarray, pi
from itertools import combinations

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
        
### TODO Let to organize chromosomes in a cubic confinement with periodic boundary condition (PBC)
###      Add the possibility to put also spheres of different radii (e.g. nucleoli)

##########

def generate_chromosome_random_walks_conformation ( chromosome_particle_numbers ,
                                                    confining_environment=['sphere',100.] ,
                                                    particle_radius=0.5 ,
                                                    seed_of_the_random_number_generator=4242424 ,
                                                    number_of_conformations=1,
                                                    outfile="Initial_random_walk_conformation.dat"):
    """
    Generates lammps initial conformation file by random walks
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param ['sphere',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
            ['sphere',radius]
            ['ellipsoid',x-semiaxes, y-semiaxes, z-semiaxes]
            ['cylinder', basal radius, height]
    :param 0.5 particle_radius: Radius of each particle.
    :param 4242424 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    random.seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)

    for cnt in xrange(number_of_conformations):

        final_random_walks = generate_random_walks(chromosome_particle_numbers,
                                                   particle_radius,
                                                   confining_environment)

        # Writing the final_random_walks conformation
        print "Succesfully generated conformation number %d\n" % (cnt+1)
        write_initial_conformation_file(final_random_walks,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        out_file=outfile)

##########
        
def generate_chromosome_rosettes_conformation ( chromosome_particle_numbers ,
                                                fractional_radial_positions=None,
                                                confining_environment=['sphere',100.] ,
                                                rosette_radius=12.0 , particle_radius=0.5 ,
                                                seed_of_the_random_number_generator=4242424 ,
                                                number_of_conformations=1,
                                                outfile = "Initial_rosette_conformation.dat" ):
    """
    Generates lammps initial conformation file by rosettes conformation
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param None fractional_radial_positions: list with fractional radial positions for all the chromosomes.
    :param ['sphere',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
            ['sphere',radius]
            ['ellipsoid',x-semiaxes, y-semiaxes, z-semiaxes]
            ['cylinder', basal radius, height]
    :param 0.5 particle_radius: Radius of each particle.
    :param 4242424 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    random.seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)    

    initial_rosettes , rosettes_lengths = generate_rosettes(chromosome_particle_numbers,
                                                            rosette_radius,
                                                            particle_radius)
    print rosettes_lengths

    
    # Constructing the rosettes conformations
    for cnt in xrange(number_of_conformations):

        particle_inside   = 0 # 0 means a particle is outside
        particles_overlap = 0 # 0 means two particles are overlapping
        while particle_inside == 0 or particles_overlap == 0:
            particle_inside   = 1
            particles_overlap = 1
            segments_P1 = []
            segments_P0 = []
            side = 0
            init_rosettes = copy.deepcopy(initial_rosettes)
            
            # Guess of the initial segment conformation:
            # 1 - each rod is placed inside the confining evironment
            # in a random position and with random orientation
            # 2 - possible clashes between generated rods are checked
            if fractional_radial_positions:
                if len(fractional_radial_positions) != len(chromosome_particle_numbers):
                    print "Please provide the desired fractional radial positions for all the chromosomes"
                    sys.exit()
                segments_P1 , segments_P0 = generate_rods_biased_conformation(rosettes_lengths, rosette_radius,
                                                                              confining_environment,
                                                                              fractional_radial_positions,
                                                                              max_number_of_temptative=100000)
            else:
                segments_P1 , segments_P0 = generate_rods_random_conformation(rosettes_lengths, rosette_radius,
                                                                              confining_environment,
                                                                              max_number_of_temptative=100000)

            # Roto-translation of the rosettes according to the segment position and orientation 
            final_rosettes = rosettes_rototranslation(init_rosettes, segments_P1, segments_P0)
            
            # Checking that the beads are all inside the confining environment and are not overlapping
            for rosette_pair in list(combinations(final_rosettes,2)):
                
                for x0,y0,z0 in zip(rosette_pair[0]['x'],rosette_pair[0]['y'],rosette_pair[0]['z']):
                    particle_inside = check_point_inside_the_confining_environment(x0, y0, z0,
                                                                                   particle_radius,
                                                                                   confining_environment)
                    
                    if particle_inside == 0: # 0 means that the particle is outside -> PROBLEM!!!
                        print "Particle",x0,y0,z0,"is out of the confining environment\n"
                        break

                    for x1,y1,z1 in zip(rosette_pair[1]['x'],rosette_pair[1]['y'],rosette_pair[1]['z']):
                        particles_overlap = check_particles_overlap(x0,y0,z0,x1,y1,z1,particle_radius)

                        if particles_overlap == 0: # 0 means that the particles are overlapping -> PROBLEM!!!
                            print "Particle",x0,y0,z0,"and",x1,y1,z1,"overlap\n"         
                            break
                    if particle_inside == 0 or particles_overlap == 0:
                        break
                if particle_inside == 0 or particles_overlap == 0:
                    break
            
        # Writing the final_rosettes conformation
        print "Succesfully generated conformation number %d\n" % (cnt+1)
        write_initial_conformation_file(final_rosettes,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        out_file=outfile)

##########

def generate_rosettes(chromosome_particle_numbers, rosette_radius, particle_radius):
    # Genaration of the rosettes
    # XXXA. Rosa publicationXXX
    # List to contain the rosettes and the rosettes lengths
    rosettes = []
    rosettes_lengths  = []

    for number_of_particles in chromosome_particle_numbers:
        
        # Variable to build the chain
        phi = 0.0

        # Dictory of lists to contain the rosette
        rosette      = {}
        rosette['x'] = []
        rosette['y'] = []
        rosette['z'] = []        
        
        # Position of the first particle (x_0, 0.0, 0.0)
        rosette['x'].append(rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * cos(phi))
        rosette['y'].append(rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * sin(phi))
        rosette['z'].append(phi / (2.0 * pi))
        #print "First bead is in position: %f %f %f" % (rosette['x'][0], rosette['y'][0], rosette['z'][0])

        # Building the chain: The rosette is growing along the positive z-axes
        for particle in xrange(1,number_of_particles):

            distance = 0.0
            while distance < (particle_radius*2.0): 
                phi   = phi + 0.001
                x_tmp = rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * cos(phi)
                y_tmp = rosette_radius * (0.38 + (1 - 0.38) * cos(6*phi) * cos(6*phi)) * sin(phi)
                z_tmp = phi / (2.0 * pi)     
                distance  = sqrt((x_tmp - rosette['x'][-1])*(x_tmp - rosette['x'][-1]) +
                                 (y_tmp - rosette['y'][-1])*(y_tmp - rosette['y'][-1]) +
                                 (z_tmp - rosette['z'][-1])*(z_tmp - rosette['z'][-1]))

            rosette['x'].append(x_tmp)
            rosette['y'].append(y_tmp)
            rosette['z'].append(z_tmp)
            if distance > ((particle_radius*2.0)*1.2):
                print "%f %d %d %d" % (distance, particle-1, particle)
            
        rosettes.append(rosette)
        rosettes_lengths.append(rosette['z'][-1]-rosette['z'][0])
        
    return rosettes , rosettes_lengths

##########

def generate_rods_biased_conformation(rosettes_lengths, rosette_radius,
                                      confining_environment,
                                      fractional_radial_positions,
                                      max_number_of_temptative=100000):
    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []

    if confining_environment[0] != 'sphere':
        print "ERROR: Biased chromosome positioning is currently implemented"
        print "only for spherical confinement. If you need other shapes, please"
        print "contact the developers"
    
    for length , target_radial_position in zip(rosettes_lengths,fractional_radial_positions):
        tentative            = 0
        clashes              = 0 # 0 means that there is an clash -> PROBLEM
        best_radial_position = 1.0
        best_radial_distance = 1.0
        best_segment_P0      = []
        best_segment_P1      = []
        
        # Positioning the rods
        while tentative < 100000 and best_radial_distance > 0.00005:                

            print "Length = %f" % length

            print "Trying to position terminus 0"
            segment_P0_tmp = []
            segment_P0_tmp = draw_point_inside_the_confining_environment(confining_environment,
                                                                         rosette_radius)
            print "Successfully positioned terminus 0: %f %f %f" % (segment_P0_tmp[0], segment_P0_tmp[1], segment_P0_tmp[2])
            
            print "Trying to position terminus 1"
            segment_P1_tmp = []                            
            segment_P1_tmp = draw_second_extreme_of_a_segment_inside_the_confining_environment(segment_P0_tmp[0],
                                                                                               segment_P0_tmp[1],
                                                                                               segment_P0_tmp[2],
                                                                                               length,
                                                                                               rosette_radius,
                                                                                               confining_environment)
            print "Successfully positioned terminus 1: %f %f %f" % (segment_P1_tmp[0], segment_P1_tmp[1], segment_P1_tmp[2])

            # Check clashes with the previously positioned rods
            clashes = 1
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes(segment_P1,
                                                 segment_P0,
                                                 segment_P1_tmp,
                                                 segment_P0_tmp,
                                                 rosette_radius)
                if clashes == 0:
                    break                

            if clashes == 1:
                # Check whether the midpoint of the segment is close to the target radial position
                segment_midpoint = []
                segment_midpoint.append((segment_P0_tmp[0] + segment_P1_tmp[0])*0.5)
                segment_midpoint.append((segment_P0_tmp[1] + segment_P1_tmp[1])*0.5)
                segment_midpoint.append((segment_P0_tmp[2] + segment_P1_tmp[2])*0.5)

                radial_position = sqrt( ( segment_midpoint[0] * segment_midpoint[0] +
                                          segment_midpoint[1] * segment_midpoint[1] +
                                          segment_midpoint[2] * segment_midpoint[2] ) /
                                        (confining_environment[1]*confining_environment[1]))

                radial_distance = fabs(radial_position-target_radial_position)

                print radial_position , target_radial_position , radial_distance , best_radial_distance , tentative
                
                # If the midpoint of the segment is closer to the target radial position than the
                # previous guesses. Store the points as the best guesses!
                if radial_distance < best_radial_distance:                    
                    best_radial_distance = radial_distance
                    best_radial_position = radial_position
                    best_tentative       = tentative+1 # The variable tentative starts from 0

                    best_segment_P0 = []
                    best_segment_P1 = []
                    for component_P0 , component_P1 in zip(segment_P0_tmp,segment_P1_tmp):
                        best_segment_P0.append(component_P0)
                        best_segment_P1.append(component_P1)
                    
                tentative = tentative + 1
                
        if best_segment_P0 == []:
            print "Valid placement not found for chromosome rosette after %d tentatives" % tentative
            sys.exit()

        print "Successfully positioned chromosome of length %lf at tentative %d of %d tentatives" % (length, best_tentative, tentative)        
        segments_P0.append(best_segment_P0)
        segments_P1.append(best_segment_P1)

    print "Successfully generated rod conformation!"
    return segments_P1 , segments_P0
    
##########

def generate_rods_random_conformation(rosettes_lengths, rosette_radius,
                                      confining_environment,
                                      max_number_of_temptative=100000):
    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []
    
    for length in rosettes_lengths:
        tentative = 0
        clashes   = 0
        # Random positioning of the rods
        while tentative < 100000 and clashes == 0:                

            tentative += 1
            clashes    = 1
            #print "Length = %f" % length

            print "Trying to position terminus 0"
            #pick uniformly within the confining environment using the rejection method 
            first_point = []
            first_point = draw_point_inside_the_confining_environment(confining_environment,
                                                                      rosette_radius)

            print "Successfully positioned terminus 0: %f %f %f" % (first_point[0], first_point[1], first_point[2])
            
            print "Trying to position terminus 1"
            #pick from P0 another point one the sphere of radius length inside the confining environment
            last_point = []
            last_point = draw_second_extreme_of_a_segment_inside_the_confining_environment(first_point[0],
                                                                                           first_point[1],
                                                                                           first_point[2],
                                                                                           length,
                                                                                           rosette_radius,
                                                                                           confining_environment)
            
            print "Successfully positioned terminus 1: %f %f %f" % (last_point[0], last_point[1], last_point[2])
                
            # Check clashes with the previously positioned rods
            clashes = 1 
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes(segment_P1,
                                                 segment_P0,
                                                 last_point,
                                                 first_point,
                                                 rosette_radius)
                if clashes == 0:
                    break                

            #print clashes
        print "Successfully positioned chromosome of length %lf at tentative %d\n" % (length, tentative)        
        segments_P1.append(last_point)
        segments_P0.append(first_point)            

    print "Successfully generated rod conformation!"
    return segments_P1 , segments_P0

##########

def generate_random_walks(chromosome_particle_numbers,
                          particle_radius,
                          confining_environment):
    # Construction of the random walks initial conformation 
    random_walks = []
    
    for number_of_particles in chromosome_particle_numbers:
        print "Trying to position random walk"
        random_walk      = {}
        random_walk['x'] = []
        random_walk['y'] = []
        random_walk['z'] = []        


        print "Positioning first particle"            
        particle_overlap = 0
        while particle_overlap == 0:
            particle_overlap = 1
            first_particle = []
            first_particle = draw_point_inside_the_confining_environment(confining_environment,
                                                                         particle_radius)
            print first_particle
            # Check if the particle is overlapping with any other particle in the system
            for rand_walk in random_walks:
                particle_overlap = check_particle_vs_all_overlap(first_particle[0],
                                                                 first_particle[1],
                                                                 first_particle[2],
                                                                 rand_walk,
                                                                 2.0*particle_radius)            
                if particle_overlap == 0:
                    break
        random_walk['x'].append(first_particle[0])        
        random_walk['y'].append(first_particle[1])
        random_walk['z'].append(first_particle[2])

        for particle in xrange(1,number_of_particles):
            print "Positioning particle %d" % (particle+1)
            particle_overlap = 0 # 0 means that there is an overlap -> PROBLEM
            while particle_overlap == 0:
                particle_overlap = 1
                new_particle = []
                new_particle = draw_second_extreme_of_a_segment_inside_the_confining_environment(random_walk['x'][-1],
                                                                                                 random_walk['y'][-1],
                                                                                                 random_walk['z'][-1],
                                                                                                 2.0*particle_radius,
                                                                                                 2.0*particle_radius,
                                                                                                 confining_environment)

                # Check if the particle is overlapping with any other particle in the system
                for rand_walk in random_walks:
                    particle_overlap = check_particle_vs_all_overlap(new_particle[0],
                                                                     new_particle[1],
                                                                     new_particle[2],
                                                                     rand_walk,
                                                                     2.0*particle_radius)
                    if particle_overlap == 0:
                        break
                if particle_overlap == 0:
                    continue
                
                # The current random walk is not yet in the list above
                particle_overlap = check_particle_vs_all_overlap(new_particle[0],
                                                                 new_particle[1],
                                                                 new_particle[2],
                                                                 random_walk,
                                                                 2.0*particle_radius)
                if particle_overlap == 0:
                    continue
                
            random_walk['x'].append(new_particle[0])        
            random_walk['y'].append(new_particle[1])
            random_walk['z'].append(new_particle[2])
                    
        print "Successfully positioned random walk of %d particles" % number_of_particles
        random_walks.append(random_walk)

    print "Successfully generated random walk conformation!"
    return random_walks

##########

def check_particle_vs_all_overlap(x,y,z,chromosome,overlap_radius):    
    particle_overlap = 1

    for x0, y0, z0 in zip(chromosome['x'],chromosome['y'],chromosome['z']):
        particle_overlap = check_particles_overlap(x0,y0,z0,x,y,z,overlap_radius)
        if particle_overlap == 0:
            return particle_overlap
        
    return particle_overlap
            
##########

def draw_second_extreme_of_a_segment_inside_the_confining_environment(x0, y0, z0, segment_length, object_radius, confining_environment):
    inside = 0
    while inside == 0:
        particle = []
        temp_theta  = arccos(2.0*random.random()-1.0)
        temp_phi    = 2*pi*random.random()
        particle.append(x0 + segment_length * cos(temp_phi) * sin(temp_theta))
        particle.append(y0 + segment_length * sin(temp_phi) * sin(temp_theta))
        particle.append(z0 + segment_length * cos(temp_theta))
        # Check if the particle is inside the confining_environment
        inside = check_point_inside_the_confining_environment(particle[0],
                                                              particle[1],
                                                              particle[2],
                                                              object_radius,
                                                              confining_environment)

    return particle

##########

def draw_point_inside_the_confining_environment(confining_environment, object_radius):
    #pick a point uniformly within the confining environment using the rejection method 

    if confining_environment[0] == 'cube':
        dimension_x = confining_environment[1] * 0.5
        dimension_y = confining_environment[1] * 0.5
        dimension_z = confining_environment[1] * 0.5        
        if len(confining_environment) > 2:
            print "# WARNING: Defined a cubical confining environment with reduntant paramenters."
            print "# Only 2 are needed the identifier and the side"
        
    if confining_environment[0] == 'sphere':
        dimension_x = confining_environment[1]
        dimension_y = confining_environment[1]
        dimension_z = confining_environment[1]
        if len(confining_environment) > 2:
            print "# WARNING: Defined a spherical confining environment with reduntant paramenters."
            print "# Only 2 are needed the identifier and the radius"
        
    if confining_environment[0] == 'ellipsoid':
        if len(confining_environment) < 4:
            print "# ERROR: Defined an ellipsoidal confining environment without the necessary paramenters."
            print "# 4 are needed the identifier, the x-semiaxes, the y-semiaxes, and the z-semiaxes"
            sys.exit()
        dimension_x = confining_environment[1]
        dimension_y = confining_environment[2]
        dimension_z = confining_environment[3]

    if confining_environment[0] == 'cylinder':
        if len(confining_environment) < 3:
            print "# WARNING: Defined a cylindrical confining environment without the necessary paramenters."
            print "# 3 are needed the identifier, the basal radius, and the height"
            sys.exit()
        dimension_x = confining_environment[1]
        dimension_y = confining_environment[1]
        dimension_z = confining_environment[2]
            
    inside = 0
    while inside == 0:
        particle = []
        particle.append((2.0*random.random()-1.0)*(dimension_x - object_radius))
        particle.append((2.0*random.random()-1.0)*(dimension_y - object_radius))
        particle.append((2.0*random.random()-1.0)*(dimension_z - object_radius))
        # Check if the particle is inside the confining_environment
        inside = check_point_inside_the_confining_environment(particle[0],
                                                              particle[1],
                                                              particle[2],
                                                              object_radius,
                                                              confining_environment)                        
    
    return particle
    
##########        

def check_point_inside_the_confining_environment(Px, Py, Pz,
                                                 object_radius,
                                                 confining_environment):
    # The shapes are all centered in the origin
    # - sphere    : radius r
    # - cube      : side
    # - cylinder  : basal radius , height
    # - ellipsoid : semi-axes a , b , c ;

    if confining_environment[0] == 'sphere':
        radius = confining_environment[1] - object_radius
        if ((Px*Px)/(radius*radius) + (Py*Py)/(radius*radius) + (Pz*Pz)/(radius*radius)) < 1.0 : return 1

    if confining_environment[0] == 'ellipsoid':
        a = confining_environment[1] - object_radius
        b = confining_environment[2] - object_radius
        c = confining_environment[3] - object_radius
        if ((Px*Px)/(a*a) + (Py*Py)/(b*b) + (Pz*Pz)/(c*c)) < 1.0 : return 1

    if confining_environment[0] == 'cube':
        hside = confining_environment[1] * 0.5 - object_radius
        if (((Px*Px)/(hside*hside)) < 1.0) and (((Py*Py)/(hside*hside)) < 1.0) and (((Pz*Pz)/(hside*hside)) < 1.0) : return 1

    if confining_environment[0] == 'cylinder':
        radius      = confining_environment[1]     - object_radius
        half_height = confining_environment[2]*0.5 - object_radius
        if (((Px*Px)/(radius*radius) + (Py*Py)/(radius*radius)) < 1.0) and (((Pz*Pz)/(half_height*half_height)) < 1.0): return 1
            
    return 0

##########

def check_segments_clashes(s1_P1, s1_P0, s2_P1, s2_P0, rosette_radius):
    # Check steric clashes without periodic boundary conditions

    if distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0) < 2.0*rosette_radius:
        # print "Clash between segments",s1_P1,s1_P0,"and",s2_P1_tmp,s2_P0_tmp,"at distance", distance
        return 0

    return 1

##########

def distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0):

    # Inspiration: http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm 
    # Copyright 2001, softSurfer (www.softsurfer.com)
    # This code may be freely used and modified for any purpose
    # providing that this copyright notice is included with it.
    # SoftSurfer makes no warranty for this code, and cannot be held
    # liable for any real or imagined damage resulting from its use.
    # Users of this code must verify correctness for their application.

    u  = []
    v  = []
    w  = []
    dP = []

    for c_s1_P1,c_s1_P0,c_s2_P1,c_s2_P0 in zip(s1_P1, s1_P0, s2_P1, s2_P0):        
        u.append(c_s1_P1 - c_s1_P0)
        v.append(c_s2_P1 - c_s2_P0)
        w.append(c_s1_P0 - c_s2_P0)
    
    a  = scalar_product(u, u)
    b  = scalar_product(u, v)
    c  = scalar_product(v, v)
    d  = scalar_product(u, w)
    e  = scalar_product(v, w)

    D  = a*c - b*b
    sD = tD = D
        
    if D < (1.0e-7):
        # Segments almost parallel 
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else:
        # Get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if (sN < 0.0):            
            # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD: # sc > 1 => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0: # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # Recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a        
    
    elif tN > tD: # tc > 1 => the t=1 edge is visible
        tN = tD
        # Recompute sc for this edge
        if (-d + b) < 0.0:
            sN = 0
        elif (-d + b) > a:
            sN = sD;
        else:
            sN = (-d + b)
            sD = a

    # Finally do the division to get sc and tc
    if abs(sN) < (1.0e-7):
        sc = 0.0
    else:
        sc = sN / sD
        
    if abs(tN) < (1.0e-7):
        tc = 0.0
    else:
        tc = tN / tD
     
    # Get the difference of the two closest points
    for i in xrange(len(w)):    
        dP.append(w[i] + ( sc * u[i] ) - ( tc * v[i] )) # = S1(sc) - S2(tc)
    
    return norm(dP)   # return the closest distance

##########

def rosettes_rototranslation(rosettes, segments_P1, segments_P0):

    for i in xrange(len(segments_P1)):
        vector = []
        theta  = []

        for component_P1 , component_P0 in zip(segments_P1[i], segments_P0[i]):
            vector.append(component_P1-component_P0)
            
        # Rotation Angles
        theta.append(atan2(vector[1],vector[2]))
      
        x_temp_2 =  vector[0]
        y_temp_2 =  cos(theta[0]) * vector[1] - sin(theta[0]) * vector[2]                
        z_temp_2 =  sin(theta[0]) * vector[1] + cos(theta[0]) * vector[2]        
        theta.append(atan2(x_temp_2,z_temp_2))
        
        x_temp_1 =  cos(theta[1]) * x_temp_2 - sin(theta[1]) * z_temp_2
        y_temp_1 =  y_temp_2
        z_temp_1 =  sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2
        
        if(z_temp_1 < 0.0):            
            z_temp_1 = -z_temp_1
            theta.append(pi)
        else:
            theta.append(0.0)
        #print x_temp_1 , y_temp_1 , z_temp_1 
        
        # Chromosome roto-translations
        for particle in xrange(len(rosettes[i]['x'])):

            x_temp_2 =   rosettes[i]['x'][particle]
            y_temp_2 =   cos(theta[2]) * rosettes[i]['y'][particle] + sin(theta[2]) * rosettes[i]['z'][particle]
            z_temp_2 = - sin(theta[2]) * rosettes[i]['y'][particle] + cos(theta[2]) * rosettes[i]['z'][particle]
            
            x_temp_1 =   cos(theta[1]) * x_temp_2 + sin(theta[1]) * z_temp_2
            y_temp_1 =   y_temp_2
            z_temp_1 = - sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2

            x =   x_temp_1;
            y =   cos(theta[0]) * y_temp_1 + sin(theta[0]) * z_temp_1;
            z = - sin(theta[0]) * y_temp_1 + cos(theta[0]) * z_temp_1;
            
            # Chromosome translations
            rosettes[i]['x'][particle] = segments_P0[i][0] + x;
            rosettes[i]['y'][particle] = segments_P0[i][1] + y;
            rosettes[i]['z'][particle] = segments_P0[i][2] + z;
    return rosettes

##########

def scalar_product(a, b):

    scalar = 0.0
    for c_a,c_b in zip(a,b):
        scalar = scalar + c_a*c_b 

    return scalar

##########

def norm(a):

    return sqrt(scalar_product(a, a))

##########

def write_initial_conformation_file(chromosomes,
                                    chromosome_particle_numbers,
                                    confining_environment,
                                    out_file="Initial_conformation.dat"):
    # Choosing the appropriate xlo, xhi...etc...depending on the confining environment
    xlim = []
    ylim = []
    zlim = []
    if confining_environment[0] == 'sphere':
        radius = confining_environment[1] + 1.0
        xlim.append(-radius)
        xlim.append(radius)
        ylim.append(-radius)
        ylim.append(radius)
        zlim.append(-radius)
        zlim.append(radius)
        
    if confining_environment[0] == 'ellipsoid':
        a = confining_environment[1] + 1.0
        b = confining_environment[2] + 1.0
        c = confining_environment[3] + 1.0
        xlim.append(-a)
        xlim.append(a)
        ylim.append(-b)
        ylim.append(b)
        zlim.append(-c)
        zlim.append(c)

    if confining_environment[0] == 'cube':
        hside = confining_environment[1] * 0.5 + 1.0
        xlim.append(-hside)
        xlim.append(hside)
        ylim.append(-hside)
        ylim.append(hside)
        zlim.append(-hside)
        zlim.append(hside)
        
    if confining_environment[0] == 'cylinder':
        radius      = confining_environment[1]   + 1.0
        hheight = confining_environment[2] * 0.5 + 1.0
        xlim.append(-radius)
        xlim.append(radius)
        ylim.append(-radius)
        ylim.append(radius)
        zlim.append(-hheight)
        zlim.append(hheight)
    
    fileout = open(out_file,'w')
    angle_types = 1
    bond_types  = 1
    atom_types  = 1
    n_chr=len(chromosomes)
    n_atoms=0
    for n in chromosome_particle_numbers:
        n_atoms+=n    
        
    fileout.write("LAMMPS input data file \n\n")
    fileout.write("%9d atoms\n" % (n_atoms))
    fileout.write("%9d bonds\n" % (n_atoms-n_chr))
    fileout.write("%9d angles\n\n" % (n_atoms-2*n_chr))
    fileout.write("%9s atom types\n" % atom_types)
    fileout.write("%9s bond types\n" % bond_types)
    fileout.write("%9s angle types\n\n" % angle_types)
    fileout.write("%6.3lf    %6.3lf     xlo xhi\n" % (xlim[0], xlim[1]))
    fileout.write("%6.3lf    %6.3lf     ylo yhi\n" % (ylim[0], ylim[1]))
    fileout.write("%6.3lf    %6.3lf     zlo zhi\n" % (zlim[0], zlim[1]))
  
    fileout.write("\n Atoms \n\n")
    particle_number = 1
    for chromosome in chromosomes:
        for x,y,z in zip(chromosome['x'],chromosome['y'],chromosome['z']):          
            fileout.write("%-8d %s %s %7.4lf %7.4lf %7.4lf\n" % (particle_number, "1", "1", x, y, z))
            particle_number += 1
            
    # for(i = 0; i < N_NUCL; i++)
    # {
    #    k++;
    #    fileout.write("%5d %s %s %7.4lf %7.4lf %7.4lf \n", k, "1", "1", P[i][0], P[i][1], P[i][2]);
    # }
  
    fileout.write("\n Bonds \n\n")
    bond_number          = 1
    first_particle_index = 1
    for chromosome in chromosomes:
        for i in xrange(len(chromosome['x'])-1):
            fileout.write("%-4d %s %4d %4d\n" % (bond_number, "1", first_particle_index, first_particle_index+1))
            bond_number          += 1
            first_particle_index += 1
        first_particle_index += 1 # I have to go to the end of the chromosome!
          
    fileout.write("\n Angles \n\n")
    angle_number         = 1
    first_particle_index = 1
    for chromosome in chromosomes:        
        for i in xrange(len(chromosome['x'])-2):
            fileout.write("%-4d %s %5d %5d %5d\n" % (angle_number, "1", first_particle_index, first_particle_index+1, first_particle_index+2))
            angle_number         += 1    
            first_particle_index += 1
            first_particle_index += 2 # I have to go to the end of the chromosome!

    fileout.close()

##########

def distance(x0,y0,z0,x1,y1,z1):
    return sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))

##########

def check_particles_overlap(x0,y0,z0,x1,y1,z1,overlap_radius):
    if distance(x0,y0,z0,x1,y1,z1) < overlap_radius:
        #print "Particle %f %f %f and particle %f %f %f are overlapping\n" % (x0,y0,z0,x1,y1,z1)
        return 0
    return 1

    
      
    
       
