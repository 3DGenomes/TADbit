"""
07 Feb 2013


"""

# GENERAL
#########

CONFIG = {
    'dmel_01': {
        'reference' : 'victor corces dataset 2013',
        # Maximum experimental contact distance
        'consdist'  : 600,
        # Force applied to the restraints inferred to neighbor particles
        'kforce'    : 5,
        # Minimum distance between two non-bonded particles
        'lowrdist'  : 100,
        # Maximum thresholds used to decide which experimental values have to be
        # included in the computation of restraints. Z-score values bigger than upfreq
        # and less that lowfreq will be include, whereas all the others will be rejected
        'upfreq'    : 0.3,
        # Minimum thresholds used to decide which experimental values have to be
        # included in the computation of restraints. Z-score values bigger than upfreq
        # and less that lowfreq will be include, whereas all the others will be rejected
        'lowfreq'   : -0.7}
    }


# MonteCarlo optimizer parameters
#################################
# number of iterations
NROUNDS   = 10000
# number of MonteCarlo steps per round
STEPS     = 1
# number of local steps per round
LSTEPS    = 5
