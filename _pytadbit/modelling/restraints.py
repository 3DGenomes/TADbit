from math           import fabs, pow as power
from collections    import OrderedDict

from scipy          import polyfit

class HiCBasedRestraints(object):

    """
    This class contains distance restraints based on Hi-C zscores.

    :param nloci: number of particles to model (may not all be present in
       zscores)
    :param particle_radius: radius of each particle in the model.
    :param None config: a dictionary containing the standard
       parameters used to generate the models. The dictionary should contain
       the keys kforce, lowrdist, maxdist, upfreq and lowfreq. Examples can be
       seen by doing:

       ::

         from pytadbit.modelling.CONFIG import CONFIG

         where CONFIG is a dictionary of dictionaries to be passed to this function:

       ::

         CONFIG = {
          'dmel_01': {
              # Paramaters for the Hi-C dataset from:
              'reference' : 'victor corces dataset 2013',

              # Force applied to the restraints inferred to neighbor particles
              'kforce'    : 5,

              # Space occupied by a nucleotide (nm)
              'scale'     : 0.005

              # Strength of the bending interaction
              'kbending'     : 0.0, # OPTIMIZATION:

              # Maximum experimental contact distance
              'maxdist'   : 600, # OPTIMIZATION: 500-1200

              # Minimum thresholds used to decide which experimental values have to be
              # included in the computation of restraints. Z-score values bigger than upfreq
              # and less that lowfreq will be include, whereas all the others will be rejected
              'lowfreq'   : -0.7 # OPTIMIZATION: min/max Z-score

              # Maximum threshold used to decide which experimental values have to be
              # included in the computation of restraints. Z-score values greater than upfreq
              # and less than lowfreq will be included, while all the others will be rejected
              'upfreq'    : 0.3 # OPTIMIZATION: min/max Z-score

              }
          }
    :param resolution:  number of nucleotides per Hi-C bin. This will be the
       number of nucleotides in each model's particle
    :param zscores: the dictionary of the Z-score values calculated from the
       Hi-C pairwise interactions
    :param 1 close_bins: number of particles away (i.e. the bin number
       difference) a particle pair must be in order to be considered as
       neighbors (e.g. 1 means consecutive particles)
    :param None first: particle number at which model should start (0 should be
       used inside TADbit)
    :param None remove_rstrn: list of particles which must not have restrains


    """
    def __init__(self, nloci, particle_radius,CONFIG,resolution,zscores,
                 chromosomes, close_bins=1,first=None, min_seqdist=0,
                 remove_rstrn=[]):

        self.particle_radius       = particle_radius
        self.nloci = nloci
        self.CONFIG = CONFIG
        self.resolution = resolution
        self.nnkforce = CONFIG['kforce']
        self.min_seqdist = min_seqdist
        self.chromosomes = OrderedDict()
        self.remove_rstrn = remove_rstrn
        if chromosomes:
            if isinstance(chromosomes,dict):
                self.chromosomes[chromosomes['crm']] = chromosomes['end'] - chromosomes['start'] + 1
            else:
                tot = 0
                for k in chromosomes:
                    tot += k['end'] - k['start'] + 1
                    self.chromosomes[k['crm']] = tot
        else:
            self.chromosomes['UNKNOWN'] = nloci

        self.CONFIG['lowrdist'] = self.particle_radius * 2.

        if self.CONFIG['lowrdist'] > self.CONFIG['maxdist']:
            raise TADbitModelingOutOfBound(
                ('ERROR: we must prevent you from doing this for the safe of our ' +
                 'universe...\nIn this case, maxdist must be higher than %s\n' +
                 '   -> resolution times scale -- %s*%s)') % (
                    self.CONFIG['lowrdist'], self.resolution, self.CONFIG['scale']))

        # print 'config:', self.CONFIG
        # get SLOPE and regression for all particles of the z-score data

        zsc_vals = [zscores[i][j] for i in zscores for j in zscores[i]
                    if abs(int(i) - int(j)) > 1] # condition is to avoid
                                                 # taking into account selfies
                                                 # and neighbors
        self.SLOPE, self.INTERCEPT   = polyfit([min(zsc_vals), max(zsc_vals)],
                                     [self.CONFIG['maxdist'], self.CONFIG['lowrdist']], 1)
        #print "#SLOPE = %f ; INTERCEPT = %f" % (self.SLOPE, self.INTERCEPT)
        #print "#maxdist = %f ; lowrdist = %f" % (self.CONFIG['maxdist'], self.CONFIG['lowrdist'])
        # get SLOPE and regression for neighbors of the z-score data
        xarray = [zscores[i][j] for i in zscores for j in zscores[i]
                  if abs(int(i) - int(j)) <= (close_bins + 1)]
        yarray = [self.particle_radius * 2 for _ in range(len(xarray))]
        try:
            self.NSLOPE, self.NINTERCEPT = polyfit(xarray, yarray, 1)
        except TypeError:
            self.NSLOPE, self.NINTERCEPT = 0.0, self.particle_radius * 2

        # if z-scores are generated outside TADbit they may not start at zero
        if first == None:
            first = min([int(j) for i in zscores for j in zscores[i]] +
                        [int(i) for i in zscores])
        self.LOCI  = list(range(first, nloci + first))

        # Z-scores

        self.PDIST = zscores

    def get_hicbased_restraints(self):

        # HiCbasedRestraints is a list of restraints returned by this function.
        # Each entry of the list is a list of 5 elements describing the details of the restraint:
        # 0 - particle_i
        # 1 - particle_j
        # 2 - type_of_restraint = Harmonic or HarmonicLowerBound or HarmonicUpperBound
        # 3 - the kforce of the restraint
        # 4 - the equilibrium (or maximum or minimum respectively) distance associated to the restraint

        HiCbasedRestraints = []
        nlocis = list(sorted(set(range(self.nloci)) - set(self.remove_rstrn)))
        for ni, i in enumerate(nlocis):
            chr1 = [k for k,v in list(self.chromosomes.items()) if v > i][0]
            for j in nlocis[ni+1:]:
                chr2 = [k for k,v in list(self.chromosomes.items()) if v > j][0]
                # Compute the sequence separation (in particles) depending on it the restraint changes
                seqdist = abs(j - i)

                # 1 - CASE OF TWO CONSECUTIVE LOCI (NEAREST NEIGHBOR PARTICLES)
                if seqdist == 1 and seqdist > self.min_seqdist:
                    if chr1 != chr2:
                        continue
                    RestraintType, dist = self.get_nearest_neighbors_restraint_distance(self.particle_radius, i, j)
                    kforce = self.nnkforce

                # 2 - CASE OF 2 SECOND NEAREST NEIGHBORS SEQDIST = 2
                if seqdist == 2 and seqdist > self.min_seqdist:
                    if chr1 != chr2:
                        continue
                    RestraintType, dist = self.get_second_nearest_neighbors_restraint_distance(self.particle_radius, i, j)
                    kforce = self.nnkforce

                # 3 - CASE OF TWO NON-CONSECUTIVE PARTICLES SEQDIST > 2
                if seqdist >  2 and seqdist > self.min_seqdist:

                    #CASES OF TWO NON-CONSECUTIVE PARTICLES SEQDIST > 2
                    RestraintType, kforce, dist = self.get_long_range_restraints_kforce_and_distance(i, j)
                    if RestraintType == "None":
                        #print "No HiC-based restraint between particles %d and %d" % (i,j)
                        continue

                HiCbasedRestraints.append([i, j, RestraintType, kforce, dist])

        return HiCbasedRestraints



    #Functions to add restraints: HarmonicRestraints , HarmonicUpperBoundRestraints , HarmonicLowerBoundRestraints
    #addNearestNeighborsRestraint , addSecondNearestNeighborsRestraint , addLongRangeRestraints
    def get_nearest_neighbors_restraint_distance(self, particle_radius, i, j):
        x=str(i)
        y=str(j)

        if x in self.PDIST and y in self.PDIST[x] and self.PDIST[x][y] > self.CONFIG['upfreq']:
            # When p1 and p2 have a contact propensity larger that upfreq
            # their spatial proximity and a partial overlap between them is enforced
            # The equilibrium distance of the spring is inferred from the 3C based Z-score
            RestraintType = "NeighborHarmonic"
            dist = distance(self.PDIST[x][y],self.NSLOPE,self.NINTERCEPT)
            #print "Distance = ", dist
        else:
            # When p1 and p2 have a contact propensity lower than upfreq they are simply connected to each other
            #p1 = model['particles'].get_particle(i)
            #p2 = model['particles'].get_particle(j)
            #dist = p1.get_value(model['radius']) + p2.get_value(model['radius'])
            RestraintType = "NeighborHarmonicUpperBound"

            dist = 2.0 * particle_radius
        return RestraintType , dist

    def get_second_nearest_neighbors_restraint_distance(self, particle_radius, i, j):
        # IMP COMMAND: Consider the particles i, j and the particle between i and j
        #p1      = model['particles'].get_particle(i)
        #p2      = model['particles'].get_particle(j)
        #pmiddle = model['particles'].get_particle(j-1)

        # The equilibrium distance is the sum of the radii of particles p1 and p2, and of the diameter of particle pmiddle
        RestraintType = "HarmonicUpperBound"
        dist = 4.0 * particle_radius
        #dist = p1.get_value(model['radius']) + p2.get_value(model['radius']) + 2.0 * pmiddle.get_value(model['radius'])
        #print p1.get_value(model['radius']) , p2.get_value(model['radius']) , pmiddle.get_value(model['radius'])

        #print RestraintType , dist
        return RestraintType , dist

    def get_long_range_restraints_kforce_and_distance(self, i, j):
        x = str(i)
        y = str(j)

        Zscore = float('nan')

        # For non consecutive particles the kforce is a function of the *C based Zscore
        # First we define the The kforce of the harmonic restraint. It is different for 3 scenarios...

        RestraintType = "None"
        kforce        = 0.0
        dist          = 0.0

        # 1 - If the Z-score between i and j is defined
        if x in self.PDIST and y in self.PDIST[x]:
            # Get the Zscore between particles p1 and p2
            Zscore = self.PDIST[x][y]
            kforce = k_force(Zscore)

        # 2 - If the Z-score is defined only for particle i (In the Hi-C matrix you could encounter zero values next to very high entries)
        elif x in self.PDIST:
            prevy = str(j - 1)
            posty = str(j + 1)
            # The Zscore is compute as the average of the Z-scores of p2 nearest neighbor particles with p1
            Zscore = (self.PDIST[x].get(prevy, self.PDIST[x].get(posty, float('nan'))) +
                      self.PDIST[x].get(posty, self.PDIST[x].get(prevy, float('nan')))) / 2
            kforce = 0.5 * k_force(Zscore)

        # 3 - If the Z-score is defined only for particle j
        else:
            prevx = str(i - 1)
            postx = str(i + 1)
            prevx = prevx if prevx in self.PDIST else postx
            postx = postx if postx in self.PDIST else prevx
            try:
                Zscore = (self.PDIST[prevx].get(y, self.PDIST[postx].get(y, float('nan'))) +
                          self.PDIST[postx].get(y, self.PDIST[prevx].get(y, float('nan')))) / 2
                # For non consecutive particles the kforce is a function of the *C based Zscore
            except KeyError:
                pass
            kforce = 0.5 * k_force(Zscore)


        # If the ZSCORE > UPFREQ the spatial proximity of particles p1 and p2 is favoured
        if Zscore > self.CONFIG['upfreq']:
            RestraintType = "Harmonic"
            dist = distance(Zscore, self.SLOPE, self.INTERCEPT)

        # If the ZSCORE < LOWFREQ the particles p1 and p2 are restrained to be far from each other.
        elif Zscore < self.CONFIG['lowfreq']:
            RestraintType = "HarmonicLowerBound"
            dist = distance(Zscore, self.SLOPE, self.INTERCEPT)

        #if RestraintType != "None":
        #    print i, j, Zscore, RestraintType, kforce, dist
        return RestraintType, kforce, dist

    # This is a function need for TADkit?
    def _get_restraints(self):
        """
        Same function as addAllHarmonic but just to get restraints
        """
        restraint_names = {'None'               : None,
                           'Harmonic'           : 'a',
                           'NeighborHarmonic'   : 'n',
                           'HarmonicUpperBound' : 'u',
                           'NeighborHarmonicUpperBound' : 'u',
                           'HarmonicLowerBound' : 'l'}

        #model = {'radius'     : IMP.FloatKey("radius"),
        #         'model'      : Model(),
        #         'restraints' : None, # 2.6.1 compat
        #         'particles'  : None}

        # set container
        #try:
        #    model['restraints'] = IMP.RestraintSet(model['model']) # 2.6.1 compat
        #except:
        #    pass

        #model['particles'] = ListSingletonContainer(IMP.core.create_xyzr_particles(
        #    model['model'], len(LOCI), RADIUS, 100000))

        restraints = {}
        for i, j, RestraintType, kforce, dist in self.get_hicbased_restraints():
            restraints[tuple(sorted((i, j)))] = restraint_names[RestraintType], dist, kforce
        return restraints

#Function to translate the Zscore value into distances and kforce values
def distance(Zscore, slope, intercept):
    """
    Function mapping the Z-scores into distances for neighbor and non-neighbor fragments (slope, intercept) are different
    """
    #print slope, intercept
    return (slope * Zscore) + intercept

def k_force(Zscore):
    """
    Function to assign to each restraint a force proportional to the underlying
    experimental value.
    """
    return power(fabs(Zscore), 0.5)

class TADbitModelingOutOfBound(Exception):
    pass
