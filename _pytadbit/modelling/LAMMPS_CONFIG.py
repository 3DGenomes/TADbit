"""
25 Oct 2016


"""
###############################################################################
# Parameters to implement Kremer&Grest polymer model                          #
# Reference paper:                                                            #
# K. Kremer and G. S. Grest                                                   #
# Dynamics of entangled linear polymer melts: A molecular-dynamics simulation #
# J Chem Phys 92, 5057 (1990)                                                 #
###############################################################################


# units http://lammps.sandia.gov/doc/units.html 
units = "lj"

# atom_style http://lammps.sandia.gov/doc/atom_style.html
atom_style = "angle"

# boundary conditions http://lammps.sandia.gov/doc/boundary.html
#boundary = "p p p"
boundary = "f f f"

# mass http://lammps.sandia.gov/doc/mass.html
mass = "* 1.0"

# neighbor http://lammps.sandia.gov/doc/neighbor.html
neighbor = "3.0 nsq" # Optional for small and low density systems
#neighbor = "0.3 bin" # Standard choice for large (> 10,000 particles) systems
neigh_modify = "every 1 delay 1 check yes"

# thermo
run = 10000
thermo = 1000 #int(float(run)/100)
#thermo_style    custom   step temp epair emol press pxx pyy pzz pxy pxz pyz vol

# Excluded volume term: Purely repulsive Lennard-Jones or Truncated and Shifted Lennard-Jones
###################################################################
#  Lennard-Jones 12-6 potential with cutoff (=truncated):         #
#  potential E=4epsilon[ (sigma/r)^12 - (sigma/r)^6]  for r<r_cut #
#  r_cut =1.12246 = 2^(1/6) is the minimum of the potential       #
###################################################################
PurelyRepulsiveLJepsilon = 1.0
PurelyRepulsiveLJsigma   = 1.0
#PurelyRepulsiveLJcutoff  = PurelyRepulsiveLJsigma * 1.12246152962189
PurelyRepulsiveLJcutoff  = 2.0 * 1.12246152962189
# Chain connectivity term: FENE potential
#########################################################
# Fene potential + Lennard Jones 12-6:                  #
#  E= - 0.5 K R0^2 ln[ 1- (r/R0)^2]                     #
#     + 4epsilon[ (sigma/r)^12 - (sigma/r)^6] + epsilon #
#########################################################
#FENEK  = 30.0
FENEK  = 300.0
FENER0 = 1.5
#FENEepsilon = 1.0
FENEepsilon = 10.0
FENEsigma   = 1.0

# Bending rigidity term: 
angle_style = "cosine"
persistence_length = 0.00

##############################
# set timestep of integrator #
##############################
#timestep = 0.01
timestep = 0.001


# HiC Config
HiC = {
    # use these paramaters with the Hi-C data from:
    'reference' : 'sample reference',
    
    # Force applied to the restraints inferred to neighbor particles
    'kforce'    : 5,

    # How much space (in nm) ocupies a nucleotide
    'scale'     : 0.01, 
    
    # Strength of the bending interaction
    'kbending'     : 0.0, # OPTIMIZATION: 
    
    # Maximum experimental contact distance
    'maxdist'   : 600, # OPTIMIZATION: 500-1200
    
    # Minimum thresholds used to decide which experimental values have to be
    # included in the computation of restraints. Z-score values bigger than upfreq
    # and less that lowfreq will be include, whereas all the others will be rejected
    'lowfreq'   : -0.7, # OPTIMIZATION: min/max Z-score

    # Maximum thresholds used to decide which experimental values have to be
    # included in the computation of restraints. Z-score values bigger than upfreq
    # and less that lowfreq will be include, whereas all the others will be rejected
    'upfreq'    : 0.3 # OPTIMIZATION: min/max Z-score

    
    }
