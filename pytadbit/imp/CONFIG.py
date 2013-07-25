"""
07 Feb 2013


"""

IMPATH = '/usr/src/imp/imp-r12787-git/imp/build/lib/'

# GENERAL
#########
# Integer number used to initialize the random number generator, set for
# reproducibility concerns
RAND_INIT = 1
# Maximum experimental contact distance
CONSDIST  = 600
# Force applied to the restraints inferred to neighbor particles
KFORCE    = 5
# Minimum distance between two non-bonded particles
LOWRDIST  = 100
# Particles initial radius
RADIUS    = 50
# Maximum thresholds used to decide which experimental values have to be
# included in the computation of restraints. Z-score values bigger than upfreq
# and less that lowfreq will be include, whereas all the others will be rejected
UPFREQ    = 0.3
# Minimum thresholds used to decide which experimental values have to be
# included in the computation of restraints. Z-score values bigger than upfreq
# and less that lowfreq will be include, whereas all the others will be rejected
LOWFREQ   = -0.7

# MonteCarlo optimizer parameters
#################################
# number of iterations
NROUNDS   = 10000
# number of MonteCarlo steps per round
STEPS     = 1
# number of local steps per round
LSTEPS    = 5
