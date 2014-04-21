"""
03 Apr 2014

All in one script from Hi-C raw data matrix to 3D models.
"""

# MatPlotLib not asking for X11
import matplotlib as mpl
mpl.use('Agg')

import pytadbit
from argparse import ArgumentParser, HelpFormatter
from pytadbit import Chromosome
from warnings import warn
import os, sys


def main():
    """
    main function
    """
    opts = get_options()
    crm = opts.crm
    res = int(opts.res)
    ini = int(float(opts.beg) / res)
    end = int(float(opts.end) / res)
    if end - ini <= 2:
        raise Exception('"beg" and "end" parameter should be given in ' +
                        'genomic coordinates, not bin')
    datasets = opts.data
    xtra = ('_' + opts.extra) if opts.extra else ''
    nmodels_opt, nkeep_opt, ncpus = (int(opts.nmodels_opt),
                                     int(opts.nkeep_opt), int(opts.ncpus))
    nmodels_mod, nkeep_mod = int(opts.nmodels_mod), int(opts.nkeep_mod)
    if opts.xname:
        xnames = opts.xname
    else:
        xnames = [os.path.split(d)[-1] for d in datasets]

    scale   = (tuple([float(i) for i in opts.scale.split(':')  ])
               if ':' in opts.scale   else float(opts.scale)  )
    maxdist = (tuple([int(i) for i in opts.maxdist.split(':')])
               if ':' in opts.maxdist else int(opts.maxdist))
    upfreq  = (tuple([float(i) for i in opts.upfreq.split(':') ])
               if ':' in opts.upfreq  else float(opts.upfreq) )
    lowfreq = (tuple([float(i) for i in opts.lowfreq.split(':')])
               if ':' in opts.lowfreq else float(opts.lowfreq))

    dcutoff  = int(res * scale * 2)         # distance cutoff equals to 2 bins
    name = '{0}_{1}_{2}'.format(crm, ini, end)
    PATH = opts.outdir
    
    # Start the entire process
    if not os.path.exists(os.path.join(PATH, name)):
        os.makedirs(os.path.join(PATH, name))

    ############################################################################
    ############################  LOAD HI-C DATA  ##############################
    ############################################################################

    # Start reading the data
    my_chrom = Chromosome(crm, species=(
        opts.species.split('_')[0].capitalize() + opts.species.split('_')[1]
                          if '_' in opts.species else opts.species),
                          assembly=opts.assembly) # Create the chromosome object

    # Load three different experimental data sets named TR1, TR2 and BR.
    # Data obtained from Hou et al (2012) Molecular Cell.
    # doi:10.1016/j.molcel.2012.08.031
    print "\tReading input data..."
    for xnam, xpath in zip(xnames, datasets):
        my_chrom.add_experiment(
            xnam, exp_type='Hi-C', enzyme=opts.enzyme,
            cell_type=opts.cell,
            identifier=opts.identifier, # general descriptive fields
            project=opts.project, # user descriptions
            resolution=res,
            hic_data=xpath,
            silent=True)

    # Sum all experiments into a new one
    print "\tSumming experiments..."
    if len(xnames) > 1:
        exp = my_chrom.experiments[0] + my_chrom.experiments[1]
        for i in range(2, len(xnames)):
            exp += my_chrom.experiments[i]

    # Normalize the sum of the three experiments
    print "\tNormalizing HiC data..."
    exp.normalize_hic()

    # New name for the sum experiments
    xnam = '+'.join(xnames)

    # Save the chromosome
    # Chromosomes can later on be loaded to avoid re-reading the original
    # matrices. See function "load_chromosome".
    print "\tSaving the chromosome..."
    my_chrom.save_chromosome(os.path.join(PATH, name, '{0}.tdb'.format(name)),
                             force=True)


    ############################################################################
    #######################  LOAD OPTIMAL IMP PARAMETERS #######################
    ############################################################################

    if not opts.analyze_only:
        # If some optimizations have finished, we load log files into a single
        # IMPoptimizer object
        print ("\tReading optimal parameters also saved in " +
               "{0}_optimal_params{1}.log").format(
            os.path.join(PATH, name), xtra)

        from pytadbit import IMPoptimizer
        results = IMPoptimizer(exp, ini, end, n_models=nmodels_opt,
                               n_keep=nkeep_opt, cutoff=dcutoff)
        # load from log files
        if not opts.optimize_only:
            for fnam in os.listdir(os.path.join(PATH, name)):
                if fnam.endswith(xtra + '.log') and '_optimal_params' in fnam:
                    print os.path.join(PATH, name, fnam)
                    results.load_from_file(os.path.join(PATH, name, fnam))

    ############################################################################
    #########################  OPTIMIZE IMP PARAMETERS #########################
    ############################################################################

    if not opts.analyze_only:
        # Find optimal parameters
        print ("\tFinding optimal parameters for modeling " +
               "(this can take long)...")
        optname = '_{}_{}_{}'.format(opts.maxdist,
                                     opts.upfreq ,
                                     opts.lowfreq)
        logpath = os.path.join(
            PATH, name, '{}_optimal_params{}{}.log'.format(name, optname, xtra))
        results.run_grid_search(n_cpus=ncpus, off_diag=2, verbose=True,
                                lowfreq_range=lowfreq, upfreq_range=upfreq,
                                maxdist_range=maxdist, scale_range=scale)
        results.write_result(logpath)
        if opts.optimize_only:
            print 'Optimization done.\n'
            exit()

        ## get best parameters

        optpar, cc = results.get_best_parameters_dict(
            reference='Optimized for {0}'.format(name), with_corr=True)

        sc = optpar['scale']
        md = optpar['maxdist']
        uf = optpar['upfreq']
        lf = optpar['lowfreq']

        print ("\t\tOptimal values: scale:{0} maxdist:{1} upper:{2} lower:{3}" +
               " with cc: {4}").format(sc, md, uf, lf, cc)
        results.write_result(os.path.join(
            PATH, name, '{}_optimal_params{}.log'.format(name, xtra)))
        if 1 in opts.analyze:
            results.plot_2d(show_best=20,
                            savefig="{0}/{1}{2}_optimal_params.pdf".format(
                                os.path.join(PATH, name), name, xtra))
        # Optimal parameters
        kf = 5 # IMP penalty for connectivity between two consecutive particles.
               # This needs to be large enough to ensure connectivity.

        optpar['kforce'] = kf # this is already the default but it can be changed
                              # like this


    ############################################################################
    ##############################  MODEL REGION ###############################
    ############################################################################

    # if models are already calculated and we just want to load them
    if (os.path.exists(os.path.join(PATH, name, name + '.models'))
        and opts.analyze_only):
        ########################################################################
        # function for loading models
        from pytadbit.imp.structuralmodels import load_structuralmodels 
        models = load_structuralmodels(
            os.path.join(PATH, name, name + xtra + '.models'))
        ########################################################################
    else:
        # Build 3D models based on the HiC data.
        print "\tModeling (this can take long)..."
        models = exp.model_region(ini, end, n_models=nmodels_mod,
                                  n_keep=nkeep_mod, n_cpus=ncpus,
                                  keep_all=True, config=optpar)
        print models

        # Save models
        # Models can later on be loaded to avoid the CPU expensive modeling.
        # See function "load_structuralmodels".
        print "\tSaving the models..."
        models.save_models(
            os.path.join(PATH, name, name + xtra + '.models'))
    
    ############################################################################
    ##############################  ANALYZE MODELS #############################
    ############################################################################

    if 2 in opts.analyze:
        # Calculate the correlation coefficient between a set of kept models and the
        # original HiC matrix
        print "\tCorrelation with data..."
        rho, pval = models.correlate_with_real_data(
            cutoff=dcutoff,
            savefig=os.path.join(PATH, name,
                                 name + xtra + '_corre_real.pdf'),
            plot=True)
        rho, pval = models.correlate_with_real_data(
            cutoff=dcutoff,
            savefig=os.path.join(PATH, name,
                                 name + xtra + '_corre_real_bis.pdf'),
            plot=True, midplot='triple')
        print "\t\tCorrelation coefficient: {0} [p-value: {1}]".format(rho,pval)

    if 3 in opts.analyze:
        # zscore plots
        print "\tZ-score plot..."
        models.zscore_plot(
            savefig=os.path.join(PATH, name, name + xtra + '_zscores.pdf'))

    # Cluster models based on structural similarity
    print "\tClustering all models into sets of structurally similar models..."
    ffact    = 0.95 # Fraction of particles that are within the dcutoff value
    clcutoff = dcutoff - 50 # RMSD cut-off to consider two models equivalent(nm)
    for clcutoff in range(clcutoff, clcutoff * 2, 50):
        try:
            print '   cutoff =', clcutoff
            models.cluster_models(fact=ffact, dcutoff=clcutoff)
            break
        except:
            continue
    print "\tSaving again the models this time with clusters..."
    models.save_models(os.path.join(PATH, name, name + xtra + '.models'))
    # Plot the clustering
    try:
        models.cluster_analysis_dendrogram(
            color=True, savefig=os.path.join(
                PATH, name, name + xtra + '_clusters.pdf'))
    except:
        print "\t\tWARNING: plot for clusters could not be made..."

    if not (opts.not_write_xyz and opts.not_write_cmm):
        # Save the clustered models into directories for easy visualization with
        # Chimera (http://www.cgl.ucsf.edu/chimera/)
        # Move into the cluster directory and run in the prompt
        # "chimera cl_1_superimpose.cmd"            
        print "\t\tWriting models, list and chimera files..."
        for cluster in models.clusters:
            print "\t\tCluster #{0} has {1} models {2}".format(
                cluster, len(models.clusters[cluster]), models.clusters[cluster])
            if not os.path.exists(os.path.join(
                PATH, name, 'models' + xtra, 'cl_' + str(cluster))):
                os.makedirs(os.path.join(
                    PATH, name, 'models' + xtra, 'cl_' + str(cluster)))
            if not opts.not_write_cmm:
                models.write_xyz(directory=os.path.join(
                    PATH, name, 'models' + xtra, 'cl_' + str(cluster)),
                                 cluster=cluster)
            if not opts.not_write_xyz:
                models.write_cmm(directory=os.path.join(
                    PATH, name, 'models' + xtra, 'cl_' + str(cluster)),
                                 cluster=cluster)
            # Write list file
            clslstfile = os.path.join(
                PATH, name,
                'models' + xtra, 'cl_{}.lst'.format(str(cluster)))
            f = open(clslstfile,'w')
            for model_n in models.clusters[cluster]:
                f.write("model.{0}\n".format(model_n))
            f.close()
            if not opts.not_write_cmm:
                # Write chimera file
                clschmfile = os.path.join(
                    PATH, name, 'models' + xtra,
                    'cl_{}_superimpose.cmd'.format(str(cluster)))
                f = open(clschmfile, 'w')
                f.write("open " + " ".join(["cl_{0}/model.{1}.cmm".format(
                    cluster, model_n) for model_n in models.clusters[cluster]]))
                f.write("\nlabel; represent wire; ~bondcolor\n")
                for i in range(1,len(models.clusters[cluster])+1):
                    f.write("match #{0} #0\n".format(i-1))
                f.close()
    
    if 4 in opts.analyze:
        # Get the centroid model of cluster #1
        print "\tGetting centroid..."
        centroid = models.centroid_model(cluster=1)
        print ("\t\tThe model centroid (closest to the average) " +
               "for cluster 1 is: {}".format(centroid))

    if 5 in opts.analyze:
        # Calculate a consistency plot for all models in cluster #1
        print "\tGetting consistency data..."
        models.model_consistency(
            cluster=1, cutoffs=range(50, dcutoff + 50, 50),
            savefig =os.path.join(PATH, name,
                                  name + xtra + '_consistency.pdf'),
            savedata=os.path.join(PATH, name,
                                  name + xtra + '_consistency.dat'))

    if 6 in opts.analyze:
        # Calculate a DNA density plot
        print "\tGetting density data..."
        models.density_plot(
            error=True, steps=(1,3,5,7),
            savefig =os.path.join(PATH, name, name + xtra + '_density.pdf'),
            savedata=os.path.join(PATH, name, name + xtra + '_density.dat'))

    if 7 in opts.analyze:
        # Get a contact map at cut-off of 150nm for cluster #1
        print "\tGetting a contact map..."
        models.contact_map(
            cluster=1, cutoff=dcutoff,
            savedata=os.path.join(PATH, name, name + xtra + '_contact.dat'))

    if 8 in opts.analyze:
        # Get Dihedral angle plot for cluster #1
        print "\tGetting angle data..."
        models.walking_angle(
            cluster=1, steps=(1,5),
            savefig =os.path.join(PATH, name, name + xtra + '_wang.pdf'),
            savedata=os.path.join(PATH, name, name + xtra + '_wang.dat'))

    if 9 in opts.analyze:
        # Get persistence length of all models
        print "\tGetting persistence length data..."
        pltfile = os.path.join(PATH, name, name + xtra + '_pL.dat')
        f = open(pltfile,'w')
        f.write('#Model_Number\tpL\n')
        for model in models:
            try:
                f.write('%s\t%.2f\n' % (model["rand_init"], model.persistence_length()))
            except:
                warn('WARNING: failed to compute persistence length for model %s' %
                     model["rand_init"])

    if 10 in opts.analyze:
        # Get accessibility of all models
        radius = 75   # Radius of an object to calculate accessibility
        nump   = 30   # number of particles to calculate accessibility (resolution)
        print "\tGetting accessibility data (this can take long)..."
        if not os.path.exists(
            os.path.join(PATH, name, 'models' + xtra, 'asa')):
            os.makedirs(os.path.join(PATH, name, 'models' + xtra, 'asa'))
        for model in models:
            by_part = model.accessible_surface(radius, nump=nump,
                                               include_edges=False)[4]
            asafile = os.path.join(PATH, name, 'models' + xtra,
                                   'asa', 'model_{}.asa'.format(model['rand_init']))
            out = open(asafile, 'w')
            for part, acc, ina in by_part:
                try:
                    out.write('%s\t%.2f\n' % (part, 100*float(acc) / (acc + ina)))
                except ZeroDivisionError:
                    out.write('%s\t%s\n' % (part, 'nan'))
            out.close()

    if 11 in opts.analyze:
        # Get interaction data of all models at 200 nm cut-off
        print "\tGetting interaction data..."
        models.interactions(
            cutoff=dcutoff, steps=(1,3,5),
            savefig =os.path.join(PATH, name,
                                  name + xtra + '_interactions.pdf'),
            savedata=os.path.join(PATH, name,
                                  name + xtra + '_interactions.dat'),
            error=True)



def get_options():
    """
    parse option from call
   
    """
    parser = ArgumentParser(
        usage="%(prog)s [options] [--cfg CONFIG_PATH]",
        description=('Options can be passed through command line or' +
                     ' a configuration file (--cfg CONG_PATH). ' +
                     'Note: options passed through command line override the ' +
                     'ones from the configuration file.'),
        formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                   max_help_position=25))
    glopts = parser.add_argument_group('General arguments')
    optimo = parser.add_argument_group('Optimization of IMP arguments')
    modelo = parser.add_argument_group('Modeling with optimal IMP arguments')
    descro = parser.add_argument_group('Descriptive, optional arguments')
    analyz= parser.add_argument_group('Output arguments')
    
    parser.add_argument('--cfg', dest='cfg', metavar="PATH", action='store',
                      default=None, type=str,
                      help='path to a configuration file with predefined ' +
                      'parameters')
    parser.add_argument('--analyze_only', dest='analyze_only',
                        action='store_true', default=False,
                        help=('load precomputed models in outdir, ' +
                              'skip optimization, modeling.'))
    parser.add_argument('--optimize_only', dest='optimize_only', default=False,
                        action='store_true',
                        help='do the optimization of the region and exit.')
    parser.add_argument('--ncpus', dest='ncpus', metavar="INT", default='1',
                        help='[%(default)s] Number of CPUs to use')

    #########################################
    # GENERAL
    glopts.add_argument(
        '--data_path', dest='data_path', metavar="PATH", action='append',
        default='', type=str,
        help=('path to search for data files (just pass file name' +
              'in "data")'))
    glopts.add_argument(
        '-d','--data', dest='data', metavar="PATH", nargs='+',
        type=str,help=('path to file(s) with Hi-C data matrix. If' +
                       ' many, experiments will be summed up. I' +
                       '.e.: --data replicate_1.txt replicate_2.txt'))
    glopts.add_argument('-n', '--xname', dest='xname', metavar="STR", nargs='+',
                        default=[], type=str,
                        help=('experiment name(s). To be used in ' +
                              'the same order as data. By default experiment' +
                              'name is file name.'))
    glopts.add_argument('--crm', dest='crm', metavar="NAME",
                        help='chromosome name')
    glopts.add_argument('--beg', dest='beg', metavar="INT",
                        help='genomic coordinate from which to start modeling')
    glopts.add_argument('--end', dest='end', metavar="INT",
                        help='genomic coordinate where to end modeling')
    glopts.add_argument('--res', dest='res', metavar="INT",
                        help='resolution of the Hi-C experiment.')
    glopts.add_argument('--outdir', dest='outdir', metavar="PATH",
                        default=None,
                        help='out directory for results')

    #########################################
    # MODELING
    modelo.add_argument('--nmodels_mod', dest='nmodels_mod', metavar="INT",
                        default='5000',
                        help=('[%(default)s] number of models to generate for' +
                              ' modeling'))
    modelo.add_argument('--nkeep_mod', dest='nkeep_mod', metavar="INT",
                        default='1000',
                        help=('[%(default)s] number of models to keep for ' +
                        'modeling'))
    # glopts.add_argument('--scale', dest='scale', metavar="LIST",
    #                   default="0.01",
    #                   help='''[%(default)s] range of numbers to be test as optimal
    #                   scale value, i.e. 0.005:0.01:0.001 -- Can also pass only one
    #                   number.''')

    #########################################
    # OPTIMIZATION
    optimo.add_argument('--maxdist', action='store', metavar="LIST",
                        default='400', dest='maxdist',
                        help='range of numbers to be test as optimal maxdist ' +
                        'value, i.e. 400:1000:100 -- Can also pass only one ' +
                        'number.')
    optimo.add_argument('--upfreq', dest='upfreq', metavar="LIST",
                      default='0',
                        help='range of numbers to be test as optimal upfreq ' +
                        'value, i.e. 0:1.2:0.3 -- Can also pass only one ' +
                        'number.')
    optimo.add_argument('--lowfreq', dest='lowfreq', metavar="LIST",
                        default='0',
                        help='range of numbers to be test as optimal lowfreq ' +
                        'value, i.e. -1.2:0:0.3 -- Can also pass only one ' +
                        'number.')
    optimo.add_argument('--nmodels_opt', dest='nmodels_opt', metavar="INT",
                        default='500',
                        help='[%(default)s] number of models to generate for ' +
                        'optimization')
    optimo.add_argument('--nkeep_opt', dest='nkeep_opt', metavar="INT",
                        default='100',
                        help='[%(default)s] number of models to keep for ' +
                        'optimization')

    #########################################
    # DESCRIPTION
    descro.add_argument('--species', dest='species', metavar="STRING",
                        help='species name, with no spaces, i.e.: homo_sapiens')
    descro.add_argument('--cell', dest='cell', metavar="STRING",
                        help='cell type name')
    descro.add_argument('--exp_type', dest='exp_type', metavar="STRING",
                        help='experiment type name (i.e.: Hi-C)')
    descro.add_argument('--assembly', dest='assembly', metavar="STRING",
                        default=None,
                        help='''NCBI ID of the original assembly
                        (i.e.: NCBI36 for human)''')
    descro.add_argument('--enzyme', dest='enzyme', metavar="STRING",
                        default=None,
                        help='''name of the enzyme used to digest
                        chromatin (i.e. HindIII)''')
    descro.add_argument('--identifier', dest='identifier', metavar="STRING",
                        default=None,
                        help='''NCBI identifier of the experiment.''')
    descro.add_argument('--project', dest='project', metavar="STRING",
                        default=None,
                        help='''project name.''')
    descro.add_argument('--extra', dest='extra', metavar="STR",
                        default='',
                        help='extra specification to appear in the' +
                        'name of out files.')


    #########################################
    # ANALYSIS
    analyz.add_argument('--analyze', dest='analyze', nargs='*',
                        choices=range(1, 12), type=int,
                        default=range(1, 12),
                        help=('''%(default)s list of numbers representing the
                        analysis to be done. Choose between:
                        1) optimization plot
                        2) correlation real/models
                        3) z-score plot
                        4) centroid
                        5) consistency
                        6) density
                        7) contact map
                        8) walking angle
                        9) persistence length
                        10) accessibility
                        11) interaction
                        '''))    
    analyz.add_argument('--not_write_cmm', dest='not_write_cmm',
                        default=False, action='store_true',
                        help=('do not generate cmm files for each model'))    
    analyz.add_argument('--not_write_xyz', dest='not_write_xyz',
                        default=False, action='store_true',
                        help=('do not generate xyz files for each model'))    

    parser.add_argument_group(optimo)
    parser.add_argument_group(modelo)
    parser.add_argument_group(descro)
    parser.add_argument_group(analyz)
    opts = parser.parse_args()

    args = [i.strip('-') for i in sys.argv]
    new_opts = {}
    if opts.cfg:
        for line in open(opts.cfg):
            if not '=' in line:
                continue
            if line.startswith('#'):
                continue
            key, value = line.split('#')[0].strip().split('=')
            key = key.strip()
            value = value.strip()
            if key in ['data', 'xname']:
                new_opts.setdefault(key, []).append(value)
                continue
            new_opts[key] = value
        # bad key in configuration file
        for bad_k in set(new_opts.keys()) - set(opts.__dict__.keys()):
            warn('WARNING: parameter "%s" not recognized' % (bad_k))
        for key in sorted(opts.__dict__.keys()):
            if key in args:
                print '  * Command setting   %12s to %s' % (key, opts.__dict__[key])
            elif key in new_opts:
                opts.__dict__[key] = new_opts[key]
                print '  - Config. setting   %12s to %s' % (key, new_opts[key])
            else:
                print '  o Default setting   %12s to %s' % (key, opts.__dict__[key])
    
    if not opts.data:
        print 'MISSING data'
        exit(parser.print_help())
    if not opts.outdir:
        print 'MISSING outdir'
        exit(parser.print_help())
    if not opts.crm:
        print 'MISSING crm NAME'
        exit(parser.print_help())
    if not opts.beg:
        print 'MISSING beg COORDINATE'
        exit(parser.print_help())
    if not opts.end:
        print 'MISSING end COORDINATE'
        exit(parser.print_help())
    if not opts.res:
        print 'MISSING resolution'
        exit(parser.print_help())
    if not opts.analyze_only:
        if not opts.maxdist:
            print 'MISSING maxdist'
            exit(parser.print_help())
        if not opts.lowfreq:
            print 'MISSING lowfreq'
            exit(parser.print_help())
        if not opts.upfreq:
            print 'MISSING upfreq'
            exit(parser.print_help())
            
    # this options should stay as this now
    opts.scale = '0.01'

    if opts.data_path:
        for i in xrange(len(opts.data)):
            print os.path.join(opts.data_path, opts.data[i])
            opts.data[i] = os.path.join(opts.data_path, opts.data[i])

    return opts



if __name__ == "__main__":
    exit(main())
