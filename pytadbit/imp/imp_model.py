"""
05 Jul 2013


"""
from pytadbit.imp.CONFIG import *
from scipy import polyfit
import math

import IMP.core
import IMP.container
import IMP.algebra
import IMP.display
from IMP import Model


class IMPmodel(Model):

    def __init__(self, name, zscores, resolution, close_bins=1):
        
        zsc = set(zscores.keys() + reduce(lambda x, y: x + y,
                                          [j.keys() for j in zscores.values()]))

        zmin = min([zscores[i][j] for i in zscores for j in zscores[i]])
        zmax = max([zscores[i][j] for i in zscores for j in zscores[i]])

        xarray = [zscores[i][j] for i in zscores for j in zscores[i]
                  if abs(i - j) <= close_bins]
        yarray = [resolution * 0.01 for _ in xrange(len(xarray))]
        self.slope, self.intercept   = polyfit([zmin, zmax],
                                             [CONSDIST, LOWRDIST], 1)
        self.dens = resolution
        self.nslope, self.nintercept = polyfit(xarray, yarray, 1)
        self.loci = range(min(zsc), max(zsc))
        self.nloci = len(self.loci)
        self.pdist = zscores
        IMP.random_number_generator.seed(int(RAND_INIT))
        self.rk = IMP.FloatKey("radius")

        Model.__init__(self)


        self.ps = IMP.container.ListSingletonContainer(
            IMP.core.create_xyzr_particles(self, self.nloci, RADIUS, 1000))
        self.ps.set_name("")
        for i in range(0, self.nloci):
            p = self.ps.get_particle(i)
            p.set_name(str(self.loci[i]))
            newrk = self.dens * 0.007 # radius = diameter/2 (0.013/2, computed following the relationship with the 30nm fiber vs 40nm here)
            p.set_value(self.rk, newrk)

        # Restraints between pairs of loci proportional to the pdist
        self.pps  = IMP.ParticlePairs()

        self.addAllHarmonics()

        # Setup an excluded volume restraint between a bunch of particles with radius
        r = IMP.core.ExcludedVolumeRestraint(self.ps, KFORCE)
        self.add_restraint(r)

        print "Total number of restraints: %i" % (self.get_number_of_restraints())

        # Setup a CMM writer
        writer = IMP.display.CMMWriter("model." + str(name) + ".cmm")
        log = IMP.display.WriteOptimizerState(writer)

        # Set up optimizer
        lo = IMP.core.ConjugateGradients()
        lo.set_model(self)
        o = IMP.core.MonteCarloWithLocalOptimization(lo, LSTEPS)
        o.set_return_best(True)
        fk = IMP.core.XYZ.get_xyz_keys()
        ptmp = self.ps.get_particles()
        mov = IMP.core.NormalMover(ptmp, fk, 0.25)
        o.add_mover(mov)
        o.add_optimizer_state(log)

        # Save a conformation every 'skip' iterations, if set at command line
        if (SKIP > 0):
            log = IMP.display.WriteOptimizerState(writer)
            g = IMP.core.XYZRsGeometry(self.ps)
            log.add_geometry(g)
            log.update()

        # Optimizer's parameters
        print "nrounds: %i, steps: %i, lsteps: %i" % (NROUNDS, STEPS, LSTEPS)
        
        # Start optimization and save an VRML after 100 MC moves
        print "Start " + str(self.evaluate(False))

        #"""simulated_annealing: preform simulated annealing for at most nrounds iterations. The optimization stops if the score does not change more than
        #    a value defined by endLoopValue and for stopCount iterations. 
        #   @param endLoopCount = Counter that increments if the score of two models did not change more than a value
        #   @param stopCount = Maximum values of iteration during which the score did not change more than a specific value
        #   @paramendLoopValue = Threshold used to compute the value  that defines if the endLoopCounter should be incremented or not"""
        # IMP.fivec.simulatedannealing.partial_rounds(m, o, nrounds, steps)
        endLoopCount = 0
        stopCount = 10
        endLoopValue = 0.00001
        # alpha is a parameter that takes into account the number of particles in the model (nloci).
        # The multiplier (in this case is 1.0) is used to give a different weight to the number of particles
        alpha = 1.0 * self.nloci
        # During the firsts hightemp iterations, do not stop the optimization
        hightemp = int(0.025 * NROUNDS)
        for i in range(0, hightemp):
            temperature = alpha * (1.1 * NROUNDS - i) / NROUNDS
            o.set_kt(temperature)
            e = o.optimize(STEPS)
            prevE = e
            print str(i) + " " + str(e) + " " + str(o.get_kt())
        # After the firsts hightemp iterations, stop the optimization if the score does not change by more than a value defined by endLoopValue and
        # for stopCount iterations
        for i in range(hightemp, NROUNDS):
            temperature = alpha * (1.1 * NROUNDS - i) / NROUNDS
            o.set_kt(temperature)
            e = o.optimize(STEPS)
            print str(i) +" " + str(e) + " " + str(o.get_kt())
            # Calculate the score variation and check if the optimization can be stopped or not
            if prevE > 0:
                deltaE = math.fabs((e - prevE)/prevE)
            else:
                deltaE = e
            if (deltaE < endLoopValue and endLoopCount == stopCount):
                break
            elif (deltaE < endLoopValue and endLoopCount < stopCount):
                endLoopCount += 1
                prevE = e
            else:
                endLoopCount = 0
                prevE = e
        #"""simulated_annealing_full: preform simulated annealing for nrounds iterations"""
        # # IMP.fivec.simulatedannealing.full_rounds(m, o, nrounds, steps)
        # alpha = 1.0 * nloci
        # for i in range(0,nrounds):
        #    temperature = alpha * (1.1 * nrounds - i) / nrounds
        #    o.set_kt(temperature)
        #    e = o.optimize(steps)
        #    print str(i) + " " + str(e) + " " + str(o.get_kt())
        # Print the IMP score of the final model

        e = self.evaluate(False)
        print "Final " + str(e)

        g = IMP.core.XYZRsGeometry(self.ps)
        log.add_geometry(g)
        log.update()



    def addHarmonicPair(self, p1, p2, x, y, j, num_loci1, num_loci2, verbose):
        seqdist = num_loci2 - num_loci1

        # SHORT RANGE DISTANCE BETWEEN TWO CONSECUTIVE LOCI
        if (seqdist == 1):
            if (x in self.pdist and y in self.pdist[x]
                and float(self.pdist[x][y]) > 0):
                kforce1 = KFORCE
                self.addHarmonicNeighborsRestraints(p1, p2, kforce1, verbose)
                #print "harmo1\t%s\t%s\t%f\t%f" % ( x, y, dist1, kforce1)
            else:
                kforce1 = KFORCE
                dist1 = (p1.get_value(self.rk) + p2.get_value(self.rk))
                self.addHarmonicUpperBoundRestraints(p1, p2, dist1, kforce1, verbose)
                #print "upper1\t%s\t%s\t%f\t%f" % ( x, y, dist1, kforce1)

                # SHORT RANGE DISTANCE BETWEEN TWO SEQDIST = 2
        elif (seqdist == 2):
            if (x in self.pdist and y in self.pdist[x] and float(self.pdist[x][y]) > 0):
                kforce2 = KFORCE
                self.addHarmonicNeighborsRestraints(p1, p2, kforce2, verbose)
            else:
                p3 = self.ps.get_particle(j-1)
                kforce2 = KFORCE
                dist2 = (p1.get_value(self.rk) + p2.get_value(self.rk)) + 2.0 * p3.get_value(self.rk)
                self.addHarmonicUpperBoundRestraints(p1, p2, dist2, kforce2, verbose)
                #print "upper2\t%s\t%s\t%f\t%f" % ( x, y, dist2, kforce2)

        else:

            # LONG RANGE DISTANCE DISTANCE BETWEEN TWO NON-CONSECUTIVE LOCI
            if (x in self.pdist and y in self.pdist[x]):
                # FREQUENCY > UPFREQ
                if (float(self.pdist[x][y]) > UPFREQ):
                    kforce3 = kForce(float(self.pdist[x][y]))
                    self.addHarmonicRestraints(p1, p2, float(self.pdist[x][y]), kforce3, verbose)
                    #print "harmo3\t%s\t%s\t%f\t%f" % ( x, y, dist3, kforce3)
                # FREQUENCY > LOW THIS HAS TO BE THE THRESHOLD FOR "PHYSICAL INTERACTIONS"
                elif (float(self.pdist[x][y]) < LOWFREQ):
                    kforce3 = kForce(float(self.pdist[x][y]))
                    self.addHarmonicLowerBoundRestraints(p1, p2, float(self.pdist[x][y]), kforce3, verbose)
               #print "lower3\t%s\t%s\t%f\t%f" % ( x, y, dist3, kforce3)
                else:
                    return

            # X IN PDIST BY Y NOT IN PDIST[X]
            elif (x in self.pdist): # and y not in pdist[x]):
                if (num_loci2 > num_loci1):
                    prev_num = num_loci2 - 1
                    pnext_num = num_loci2 + 1
                else:
                    prev_num = num_loci1 - 1
                    pnext_num = num_loci1 + 1
                prev = str(prev_num)
                pnext = str(pnext_num)

                if (prev in self.pdist[x] and pnext in self.pdist[x]):
                    virt_freq = (float(self.pdist[x][prev]) + float(self.pdist[x][pnext])) / 2.0
                    if (virt_freq > UPFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "harmo4\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    elif (virt_freq < LOWFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicLowerBoundRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "lower4\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    else:
                        return

                elif (pnext in self.pdist[x]):
                    virt_freq = float(self.pdist[x][pnext])
                    if (virt_freq > UPFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "harmo5\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    elif (virt_freq < LOWFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicLowerBoundRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "lower5\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    else:
                        return

                elif (prev in self.pdist[x]):
                    virt_freq = float(self.pdist[x][prev])
                    if (virt_freq > UPFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "harmo6\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    elif (virt_freq < LOWFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicLowerBoundRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "lower6\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    else:
                        return

                else:
                    return

            # MISSING DATA (X)
            else:
                if (num_loci2 > num_loci1):
                    xprev_num = num_loci1 - 1
                    xpnext_num = num_loci1 + 1
                    prev_num = num_loci2 - 1
                    pnext_num = num_loci2 + 1
                else:
                    xprev_num = num_loci2 - 1
                    xpnext_num = num_loci2 + 1
                    prev_num = num_loci1 - 1
                    pnext_num = num_loci1 + 1
                xprev = str(xprev_num)
                xpnext = str(xpnext_num)
                prev = str(prev_num)
                pnext = str(pnext_num)

                # CASE 1
                if (xprev in self.pdist and xpnext in self.pdist):
                    if (y in self.pdist[xprev] and y in self.pdist[xpnext]):
                        virt_freq = ( float(self.pdist[xprev][y]) + float(self.pdist[xpnext][y]) ) / 2.0
                        kforce4 = 0.5 * kForce(virt_freq)
                    elif (y in self.pdist[xprev]):
                        virt_freq = float(self.pdist[xprev][y])
                        kforce4 = 0.5 * kForce(virt_freq)
                    elif (y in self.pdist[xpnext]):
                        virt_freq = float(self.pdist[xpnext][y])
                        kforce4 = 0.5 * kForce(virt_freq)
                    else:
                        return

                    if (virt_freq > UPFREQ):
                        self.addHarmonicRestraints(p1, p2, virt_freq, kforce4, verbose)
                    elif (virt_freq < LOWFREQ):
                        self.addHarmonicLowerBoundRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "lower7\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    else:
                        return

                # CASE 2
                elif (xprev in self.pdist and y in self.pdist[xprev]):
                    virt_freq = float(self.pdist[xprev][y])
                    if (virt_freq > UPFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "harmo8\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    elif (virt_freq < LOWFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicLowerBoundRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "lower8\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    else:
                        return

                # CASE 3
                elif (xpnext in self.pdist and y in self.pdist[xpnext]):
                    virt_freq = float(self.pdist[xpnext][y])
                    if (virt_freq > UPFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "harmo9\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    elif (virt_freq < LOWFREQ):
                        kforce4 = 0.5 * kForce(virt_freq)
                        self.addHarmonicLowerBoundRestraints(p1, p2, virt_freq, kforce4, verbose)
                        #print "lower9\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                    else:
                        return

                else:
                    return


    def addAllHarmonics(self, verbose=False):
        for i in range(0, self.nloci):
            p1 = self.ps.get_particle(i)
            x = p1.get_name()
            num_loci1 = int(x)

            for j in range(i+1, self.nloci):
                p2 = self.ps.get_particle(j)
                y = p2.get_name()
                num_loci2 = int(y)

                self.addHarmonicPair(p1, p2, x, y, j, num_loci1, num_loci2, verbose)

    
    def distConseq12(self, freq):
        """
        Function mapping the Z-scores into distances for neighbor fragments
        """
        return (self.nslope * freq) + self.nintercept


    def distance(self, freq):
        """
        Function mapping the Z-scores into distances for non-neighbor fragments
        """
        return (self.slope * freq) + self.intercept


    def addHarmonicNeighborsRestraints(self, p1, p2, kforce, verbose):
        dist = self.distConseq12(self.pdist[p1.get_name()][p2.get_name()])
        p = IMP.ParticlePair(p1, p2)
        self.pps.append(p)
        dr = IMP.core.DistanceRestraint(IMP.core.Harmonic(dist, kforce),p1, p2)
        self.add_restraint(dr)
        if verbose:
            print "addHn\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)


    def addHarmonicUpperBoundRestraints(self, p1, p2, dist, kforce, verbose):
        #dist = (p1.get_value(rk) + p2.get_value(rk))
        p = IMP.ParticlePair(p1, p2)
        self.pps.append(p)
        dr = IMP.core.DistanceRestraint(IMP.core.HarmonicUpperBound(dist, kforce),p1, p2)
        self.add_restraint(dr)
        if verbose:
            print "addHu\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)


    def addHarmonicRestraints(self, p1, p2, freq, kforce, verbose):
        dist = self.distance(freq)
        p = IMP.ParticlePair(p1, p2)
        self.pps.append(p)
        dr = IMP.core.DistanceRestraint(IMP.core.Harmonic(dist, kforce),p1, p2)
        self.add_restraint(dr)
        if verbose:
            print "addHa\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)

    def addHarmonicLowerBoundRestraints(self, p1, p2, freq, kforce, verbose):
        dist = self.distance(freq)
        p = IMP.ParticlePair(p1, p2)
        self.pps.append(p)
        dr = IMP.core.DistanceRestraint(IMP.core.HarmonicLowerBound(dist, kforce),p1, p2)
        self.add_restraint(dr)
        if verbose:
            print "addHl\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)


def kForce(freq):
    """
    Function to assign to each restraint a force proportional to the underlying
    experimental value.
    """
    return math.pow(math.fabs(freq), 0.5 )


