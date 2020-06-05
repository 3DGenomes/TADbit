"""
07 Feb 2013


"""

# GENERAL
#########

CONFIG = {
    # use these paramaters with the Hi-C data from:
    'reference' : 'example dataset',

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


# MonteCarlo optimizer parameters
#################################
# number of iterations
NROUNDS   = 10000
# number of MonteCarlo steps per round
STEPS     = 1
# number of local steps per round
LSTEPS    = 5
