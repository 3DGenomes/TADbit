"""
08 Feb 2013


"""
import IMP.core
import IMP.container
import IMP.algebra
import IMP.display
import math
from re import findall

from CONFIG import *
from scipy import polyfit


def parse_zscores(zscore_f, dens, close_bins=1):
    """

    :param 1 close_bins: when to consider 2 bins as being close. 1 means only
       consecutive bins, while 2 means that we allow one bin in between (3 two
       bins, 4 three bins...)
    """
    zscores = {}
    loci   = []
    xarray = []
    yarray = []
    zmin   = float('inf')
    zmax   = float('-inf')
    for line in open(zscore_f):
        x, y, z = line.split()
        z = float(z)
        x = int(x)
        y = int(y)
        if z < zmin:
            zmin = z
        if z > zmax:
            zmax = z
        if x not in loci:
            loci.append(x)
        if y not in loci:
            loci.append(y)
        zscores.setdefault(x, {})
        zscores[x][y] = z
        # 
        if abs(x - y) <= close_bins:
            xarray.append(zscores[x][y])
            yarray.append(0.005 * (dens[x] + dens[y]))
        zscores.setdefault(y, {})
        zscores[y][x] = z
    loci.sort()
    nloci = len(loci)
    
    i = 0
    index = loci[0]
    while index <= loci[nloci-1]:
        if loci[i] > index:
            loci.insert( i, index )
            nloci = nloci + 1
        index = index + 1
        i = i + 1
        if i == nloci:
            break
    # Calculation of the slope and intercept of the calibration curve used to
    # derived distances for non-neighbor fragments, ie for fragments (i...i+3,n)
    global SLOPE
    global INTERCEPT
    SLOPE, INTERCEPT = polyfit([zmin, zmax], [CONSDIST, LOWRDIST], 1)

    # Calculation of the slope and intercept of the calibration curve used to
    # derived distances for neighbor fragments, ie for fragments (i...i+1,2)
    global NSLOPE
    global NINTERCEPT
    NSLOPE, NINTERCEPT = polyfit(xarray, yarray, 1)
    
    return zscores, loci, nloci



def get_loci(exp, close_bins=1):
    loci   = []
    xarray = []
    yarray = []
    zmin   = float('inf')
    zmax   = float('-inf')
    for x in xrange(exp.size):
        for y in xrange(exp.size):
            try:
                z = exp._zscores[x][y]
            except KeyError:
                continue
            if z < zmin:
                zmin = z
            if z > zmax:
                zmax = z
            if x not in loci:
                loci.append(x)
            if y not in loci:
                loci.append(y)
            # 
            if abs(x - y) <= close_bins:
                xarray.append(exp._zscores[x][y])
                yarray.append(0.005 * (exp.resolution * 2))
    loci.sort()
    nloci = len(loci)
    
    i = 0
    index = loci[0]
    while index <= loci[nloci-1]:
        if loci[i] > index:
            loci.insert( i, index )
            nloci = nloci + 1
        index = index + 1
        i = i + 1
        if i == nloci:
            break
    # Calculation of the slope and intercept of the calibration curve used to
    # derived distances for non-neighbor fragments, ie for fragments (i...i+3,n)
    global SLOPE
    global INTERCEPT
    SLOPE, INTERCEPT = polyfit([zmin, zmax], [CONSDIST, LOWRDIST], 1)

    # Calculation of the slope and intercept of the calibration curve used to
    # derived distances for neighbor fragments, ie for fragments (i...i+1,2)
    global NSLOPE
    global NINTERCEPT
    NSLOPE, NINTERCEPT = polyfit(xarray, yarray, 1)
    
    return zscores, loci, nloci
    

def parse_densities(dens_f):
    densities = {}
    for line in open(dens_f):
        line = line.split()
        densities[int(line[0])] = int(line[1])
    return densities

# Function mapping the Z-scores into distances for neighbor fragments
def distConseq12(freq):
    return (NSLOPE * freq) + NINTERCEPT

# Function mapping the Z-scores into distances for non-neighbor fragments
def distance(freq):
    return (SLOPE * freq) + INTERCEPT

# Function to assign to each restraint a force proportional to the underlying experimental value
def kForce(freq):
    return math.pow(math.fabs(freq), 0.5 )

def addHarmonicNeighborsRestraints(p1, p2, pps, freq, kforce, m, verbose):
    dist = distConseq12(freq)
    p = IMP.ParticlePair(p1, p2)
    pps.append(p)
    dr = IMP.core.DistanceRestraint(IMP.core.Harmonic(dist, kforce),p1, p2)
    dri = m.add_restraint(dr)
    if verbose:
        print "addHn\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)

def addHarmonicUpperBoundRestraints(p1, p2, pps, dist, kforce, m, verbose):
    #dist = (p1.get_value(rk) + p2.get_value(rk))
    p = IMP.ParticlePair(p1, p2)
    pps.append(p)
    dr = IMP.core.DistanceRestraint(IMP.core.HarmonicUpperBound(dist, kforce),p1, p2)
    dri = m.add_restraint(dr)
    if verbose:
        print "addHu\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)

def addHarmonicRestraints(p1, p2, pps, freq, kforce, m, verbose):
    dist = distance(freq)
    p = IMP.ParticlePair(p1, p2)
    pps.append(p)
    dr = IMP.core.DistanceRestraint(IMP.core.Harmonic(dist, kforce),p1, p2)
    dri = m.add_restraint(dr)
    if verbose:
        print "addHa\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)

def addHarmonicLowerBoundRestraints(p1, p2, pps, freq, kforce, m, verbose):
    dist = distance(freq)
    p = IMP.ParticlePair(p1, p2)
    pps.append(p)
    dr = IMP.core.DistanceRestraint(IMP.core.HarmonicLowerBound(dist, kforce),p1, p2)
    dri = m.add_restraint(dr)
    if verbose:
        print "addHl\t%s\t%s\t%f\t%f" % (p1.get_name(), p2.get_name(), dist, kforce)


def create_model(dens, loci, nloci, pdist, runname, verbose=True):
    IMP.random_number_generator.seed(int(RAND_INIT))
    rk = IMP.FloatKey("radius")
    m = IMP.Model()
    
    ps = IMP.container.ListSingletonContainer(\
        IMP.core.create_xyzr_particles(m, nloci, RADIUS, 1000))
    ps.set_name("")
    for i in range(0,nloci):
        p = ps.get_particle(i)
        p.set_name(str(loci[i]))
        if (str(loci[i]) in dens):
            newrk = float(dens[str(loci[i])]) * 0.007 # radius = diameter/2 (0.013/2, computed following the relationship with the 30nm fiber vs 40nm here)
            p.set_value(rk, newrk)
    
    # Restraints between pairs of loci proportional to the pdist
    pps  = IMP.ParticlePairs()
    for i in range(0,nloci):
        p1 = ps.get_particle(i)
        x = p1.get_name()
        num_loci1 = int(x)
    
        for j in range(i+1,nloci):
            p2 = ps.get_particle(j)
            y = p2.get_name()
            num_loci2 = int(y)
    
            seqdist =  num_loci2 - num_loci1
    
            # SHORT RANGE DISTANCE BETWEEN TWO CONSECUTIVE LOCI
            if (seqdist == 1):
                if (x in pdist and y in pdist[x] and float(pdist[x][y]) > 0):
                    kforce1 = KFORCE
                    addHarmonicNeighborsRestraints(p1, p2, pps, float(pdist[x][y]), kforce1, m, verbose)
                    #print "harmo1\t%s\t%s\t%f\t%f" % ( x, y, dist1, kforce1)
                else:
                    kforce1 = KFORCE
                    dist1 = (p1.get_value(rk) + p2.get_value(rk))
                    addHarmonicUpperBoundRestraints(p1, p2, pps, dist1, kforce1, m, verbose)
                    #print "upper1\t%s\t%s\t%f\t%f" % ( x, y, dist1, kforce1)
    
                    # SHORT RANGE DISTANCE BETWEEN TWO SEQDIST = 2
            elif (seqdist == 2):
                if (x in pdist and y in pdist[x] and float(pdist[x][y]) > 0):
                    kforce2 = KFORCE
                    addHarmonicNeighborsRestraints(p1, p2, pps, float(pdist[x][y]), kforce2, m, verbose)
                else:
                    p3 = ps.get_particle(j-1)
                    kforce2 = KFORCE
                    dist2 = (p1.get_value(rk) + p2.get_value(rk)) + 2.0 * p3.get_value(rk)
                    addHarmonicUpperBoundRestraints(p1, p2, pps, dist2, kforce2, m, verbose)
                    #print "upper2\t%s\t%s\t%f\t%f" % ( x, y, dist2, kforce2)
    
            else:
    
                # LONG RANGE DISTANCE DISTANCE BETWEEN TWO NON-CONSECUTIVE LOCI
                if (x in pdist and y in pdist[x]):
                    # FREQUENCY > UPFREQ
                    if (float(pdist[x][y]) > UPFREQ):
                        kforce3 = kForce(float(pdist[x][y]))
                        addHarmonicRestraints(p1, p2, pps, float(pdist[x][y]), kforce3, m, verbose)
                        #print "harmo3\t%s\t%s\t%f\t%f" % ( x, y, dist3, kforce3)
                    # FREQUENCY > LOW THIS HAS TO BE THE THRESHOLD FOR "PHYSICAL INTERACTIONS"
                    elif (float(pdist[x][y]) < LOWFREQ):
                        kforce3 = kForce(float(pdist[x][y]))
                        addHarmonicLowerBoundRestraints(p1, p2, pps, float(pdist[x][y]), kforce3, m, verbose)
                   #print "lower3\t%s\t%s\t%f\t%f" % ( x, y, dist3, kforce3)
                    else:
                        continue
    
                # X IN PDIST BY Y NOT IN PDIST[X]
                elif (x in pdist): # and y not in pdist[x]):
                    if (num_loci2 > num_loci1):
                        prev_num = num_loci2 - 1
                        pnext_num = num_loci2 + 1
                    else:
                        prev_num = num_loci1 - 1
                        pnext_num = num_loci1 + 1
                    prev = str(prev_num)
                    pnext = str(pnext_num)
    
                    if (prev in pdist[x] and pnext in pdist[x]):
                        virt_freq = (float(pdist[x][prev]) + float(pdist[x][pnext])) / 2.0
                        if (virt_freq > UPFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "harmo4\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        elif (virt_freq < LOWFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicLowerBoundRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "lower4\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        else:
                            continue
                    
                    elif (pnext in pdist[x]):
                        virt_freq = float(pdist[x][pnext])
                        if (virt_freq > UPFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "harmo5\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        elif (virt_freq < LOWFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicLowerBoundRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "lower5\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        else:
                            continue
                    
                    elif (prev in pdist[x]):
                        virt_freq = float(pdist[x][prev])
                        if (virt_freq > UPFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "harmo6\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        elif (virt_freq < LOWFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicLowerBoundRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "lower6\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        else:
                            continue
                    
                    else:
                        continue
    
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
                    if (xprev in pdist and xpnext in pdist):
                        if (y in pdist[xprev] and y in pdist[xpnext]):
                            virt_freq = ( float(pdist[xprev][y]) + float(pdist[xpnext][y]) ) / 2.0
                            kforce4 = 0.5 * kForce(virt_freq)
                        elif (y in pdist[xprev]):
                            virt_freq = float(pdist[xprev][y])
                            kforce4 = 0.5 * kForce(virt_freq)
                        elif (y in pdist[xpnext]):
                            virt_freq = float(pdist[xpnext][y])
                            kforce4 = 0.5 * kForce(virt_freq)
                        else:
                            continue
    
                        if (virt_freq > UPFREQ):
                            addHarmonicRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                        elif (virt_freq < LOWFREQ):
                            addHarmonicLowerBoundRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "lower7\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        else:
                            continue
    
                    # CASE 2
                    elif (xprev in pdist and y in pdist[xprev]):
                        virt_freq = float(pdist[xprev][y])
                        if (virt_freq > UPFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "harmo8\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        elif (virt_freq < LOWFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicLowerBoundRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "lower8\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        else:
                            continue
    
                    # CASE 3
                    elif (xpnext in pdist and y in pdist[xpnext]):
                        virt_freq = float(pdist[xpnext][y])
                        if (virt_freq > UPFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "harmo9\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        elif (virt_freq < LOWFREQ):
                            kforce4 = 0.5 * kForce(virt_freq)
                            addHarmonicLowerBoundRestraints(p1, p2, pps, virt_freq, kforce4, m, verbose)
                            #print "lower9\t%s\t%s\t%f\t%f" % ( x, y, dist4, kforce4)
                        else:
                            continue
                       
                    else:
                        continue
    
    # Setup an excluded volume restraint between a bunch of particles with radius
    r = IMP.core.ExcludedVolumeRestraint(ps, KFORCE)
    m.add_restraint(r)

    print "Total number of restraints: %i" % (m.get_number_of_restraints())
    
    # Setup a CMM writer
    writer = IMP.display.CMMWriter("model." + str(runname) + ".cmm")
    log = IMP.display.WriteOptimizerState(writer)
    
    # Set up optimizer
    lo = IMP.core.ConjugateGradients()
    lo.set_model(m)
    o = IMP.core.MonteCarloWithLocalOptimization(lo, LSTEPS)
    o.set_return_best(True)
    fk = IMP.core.XYZ.get_xyz_keys()
    ptmp = ps.get_particles()
    mov = IMP.core.NormalMover(ptmp, fk, 0.25)
    o.add_mover(mov)
    o.add_optimizer_state(log)
    
    # Save a conformation every 'skip' iterations, if set at command line
    if (SKIP > 0):
        log = IMP.display.WriteOptimizerState(writer)
        g = IMP.core.XYZRsGeometry(ps)
        log.add_geometry(g)
        log.update()
    
    # Optimizer's parameters
    print "nrounds: %i, steps: %i, lsteps: %i" % (NROUNDS, STEPS, LSTEPS)
    
    # Start optimization and save an VRML after 100 MC moves
    print "Start " + str(m.evaluate(False))
    
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
    alpha = 1.0 * nloci
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
    
    e = m.evaluate(False)
    print "Final " + str(e)
    
    g = IMP.core.XYZRsGeometry(ps)
    log.add_geometry(g)
    log.update()


def parse_cmm(cmm_f):
    re_line = ('<marker id="[0-9]+" x="(-?[0-9]+.[0-9]+)" ' +
               'y="(-?[0-9]+.[0-9]+)" z="(-?[0-9]+.[0-9]+)" ' +
               'radius="(-?[0-9]+)" r="0.7" g="0.7" b="0.7" note=" geometry"/>')
    residues = []
    for line in open(cmm_f):
        try:
            x, y, z, radius = findall(re_line, line)[0]
        except IndexError:
            continue
        residues.append({'x': float(x),
                         'y': float(y),
                         'z': float(z),
                         'radius': radius})
    return residues


def color_residues(residues):
    npart = len(residues)
    for n in xrange(npart):
        residues[n]['r'] = float(n+1)/npart
        residues[n]['b'] = 1 - residues[n]['r']
        residues[n]['g'] = 0


def write_residues(residues, runname):
    out = '<marker_set name=\"marker set $marker\">\n'
    form = ('<marker id=\"{0}\" x=\"{1}\" y=\"{2}\" z=\"{3}\" r=\"{4}\" ' +
            'g=\"{5}\" b=\"{6}\" radius=\"0.5\" note=\"{0}\"/>\n')
    for n in xrange(len(residues)):
        out += form.format(n + 1, residues[n]['x'], residues[n]['y'],
                           residues[n]['z'], residues[n]['r'],
                           residues[n]['g'], residues[n]['b'])
    form = ('<link id1=\"{}\" id2=\"{}\" r=\"1\" ' +
            'g=\"1\" b=\"1\" radius=\"0.1\"/>\n')
    for n in xrange(1, len(residues)):
        out += form.format(n, n + 1)
    out += '</marker_set>\n'

    out_f = open("model." + str(runname) + ".cmm", 'w')
    out_f.write(out)
    out_f.close()
    form = "{:>12}{:>12}{x:>12.3f}{y:>12.3f}{z:>12.3f}\n"
    out = ''
    for n in xrange(len(residues)):
        out += form.format('p'+str(n+1), n+1, **residues[n])
    out_f = open("model." + str(runname) + ".xyz", 'w')
    out_f.write(out)
    out_f.close()


def main():
    """
    main function
    """
    
    dens = parse_densities('/home/fransua/Tools/imp-scripts/parameters/' +
                           'Fragm_Oligos_ChrIV_dens.txt')
    pdist, loci, nloci = parse_zscores('/home/fransua/Tools/imp-scripts/' +
                                       'parameters/lane7_chrIV_input_zscore.txt',
                                       dens)

    create_model(dens, loci, nloci, pdist, 'ici', verbose=False)

    residues = parse_cmm("model.ici.cmm")

    color_residues(residues)

    write_residues(residues, 'ici')

