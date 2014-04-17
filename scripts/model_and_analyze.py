"""
03 Apr 2014


"""

# MatPlotLib not asking for X11
import matplotlib as mpl
mpl.use('Agg')

import pytadbit
from optparse import OptionParser, OptionGroup
from pytadbit import Chromosome
from warnings import warn
import os


def main():
    """
    main function
    """
    opts = get_options()
    crm = opts.crm
    res = int(opts.res)
    ini = int(float(opts.beg) / res)
    end = int(float(opts.end) / res)
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
    #########################  OPTIMIZE IMP PARAMETERS #########################
    ############################################################################

    # Find optimal paramaters
    if opts.optimize or not opts.load_optimization:
        print "\tFinding optimal parameters for modeling (this can take long)..."
        optname = '_{}_{}_{}'.format(opts.maxdist,
                                     opts.upfreq ,
                                     opts.lowfreq)
        logpath = os.path.join(
            PATH, name, '{}_optimal_params{}{}.log'.format(name, optname, xtra))
        mdlpath = os.path.join(
            PATH, name, '{}_optimal_params{}{}.mdls'.format(name, optname, xtra))

        results = exp.optimal_imp_parameters(
            ini, end, n_models=nmodels_opt, n_keep=nkeep_opt, n_cpus=ncpus,
            lowfreq_range=lowfreq, upfreq_range=upfreq, maxdist_range=maxdist,
            scale_range=scale, verbose=True, outfile=logpath, 
            corr='spearman', cutoff=dcutoff, off_diag=2, savedata=mdlpath)
        
        if opts.optimize:
            print 'Optimization done.\n'
            exit()

    ############################################################################
    #######################  LOAD OPTIMAL IMP PARAMETERS #######################
    ############################################################################

    # If all optimizations have finished, we load all log file into a single
    # IMPoptimizer object
    if opts.load_optimization:
        # Define optimal paramaters
        print ("\tReading optimal parameters also saved in " +
               "{0}_optimal_params{1}.mdls").format(
            os.path.join(PATH, name), xtra)
        
        from pytadbit import IMPoptimizer
        results = IMPoptimizer(exp, ini, end, n_models=nmodels_mod,
                               n_keep=nkeep_mod, cutoff=dcutoff)
        # load from saved models
        fnams = []
        for fnam in os.listdir(os.path.join(PATH, name)):
            if fnam.endswith(xtra + '.mdls') and '_optimal_params' in fnam:
                fnams.append(os.path.join(PATH, name, fnam))
        results.load_grid_search(fnams, corr='spearman',
                                 off_diag=2,
                                 n_cpus=ncpus)
    optpar, cc = results.get_best_parameters_dict(
        reference='Optimized for {0}'.format(name), with_corr=True)

    sc = optpar['scale']
    md = optpar['maxdist']
    uf = optpar['upfreq']
    lf = optpar['lowfreq']

    print ("\t\tOptimal values: scale:{0} maxdist:{1} upper:{2} lower:{3} " +
           "with cc: {4}").format(sc, md, uf, lf, cc)
    if opts.view_optimization:
        results.write_result(os.path.join(
            PATH, name, '{}_optimal_params{}.log'.format(name, xtra)))
        # import cPickle
        # models = {}
        # for fnam in os.listdir(os.path.join(PATH, name)):
        #     if fnam.endswith(xtra + '.mdls') and '_optimal_params' in fnam:
        #         models.update(cPickle.load(open(os.path.join(PATH, name, fnam))))
                
        # results.plot_2d(show_best=20, savefig="{0}/{1}{2}_optimal_params.pdf".format(
        #    os.path.join(PATH, name), name, xtra))
        exit()
    results.plot_2d(show_best=20, savefig="{0}/{1}{2}_optimal_params.pdf".format(
        os.path.join(PATH, name), name, xtra))
    # Optimal parameters
    kf = 5 # IMP penalty for connectivity between two consequtive particles.
           # This needs to be large enough to ensure connectivity.
           
    optpar['kforce'] = kf # this is already the default but it can be changed
                          # like this


    ############################################################################
    ##############################  MODEL REGION ###############################
    ############################################################################

    load = False # TODO opts.load_models
    # if models are already calculated and we just want to load them
    if (os.path.exists(os.path.join(PATH, name, name + '.models'))
        and load):
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
        models.write_xyz(directory=os.path.join(
            PATH, name, 'models' + xtra, 'cl_' + str(cluster)),
                         cluster=cluster)
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
        # Write chimera file
        clschmfile = os.path.join(PATH, name, 'models' + xtra,
                                  'cl_{}_superimpose.cmd'.format(str(cluster)))
        f = open(clschmfile, 'w')
        f.write("open " + " ".join(["cl_{0}/model.{1}.cmm".format(cluster,
                                                                  model_n)
                                    for model_n in models.clusters[cluster]]))
        f.write("\nlabel; represent wire; ~bondcolor\n")
        for i in range(1,len(models.clusters[cluster])+1):
            f.write("match #{0} #0\n".format(i-1))
        f.close()
    
    # Get the centroid model of cluster #1
    print "\tGetting centroid..."
    centroid = models.centroid_model(cluster=1)
    print ("\t\tThe model centroid (closest to the average) " +
           "for cluster 1 is: {}".format(centroid))

    # Calculate a consistency plot for all models in cluster #1
    print "\tGetting consistency data..."
    models.model_consistency(
        cluster=1, cutoffs=range(50, dcutoff + 50, 50),
        savefig =os.path.join(PATH, name,
                              name + xtra + '_consistency.pdf'),
        savedata=os.path.join(PATH, name,
                              name + xtra + '_consistency.dat'))

    # Calculate a DNA density plot
    print "\tGetting density data..."
    models.density_plot(
        error=True, steps=(1,3,5,7),
        savefig =os.path.join(PATH, name, name + xtra + '_density.pdf'),
        savedata=os.path.join(PATH, name, name + xtra + '_density.dat'))

    # Get a contact map at cut-off of 150nm for cluster #1
    print "\tGetting a contact map..."
    models.contact_map(
        cluster=1, cutoff=dcutoff,
        savedata=os.path.join(PATH, name, name + xtra + '_contact.dat'))

    # Get Dihedral angle plot for cluster #1
    print "\tGetting angle data..."
    models.walking_angle(
        cluster=1, steps=(1,5),
        savefig =os.path.join(PATH, name, name + xtra + '_wang.pdf'),
        savedata=os.path.join(PATH, name, name + xtra + '_wang.dat'))

    # Get persistence length of all models
    print "\tGetting persistence length data..."
    pltfile = os.path.join(PATH, name, name + xtra + '_pL.dat')
    f = open(pltfile,'w')
    f.write('#Model_Number\tpL\n')
    for model in models:
        f.write('%s\t%.2f\n' % (model["rand_init"], model.persistence_length()))

    # Get accessibility of all models
    radius = 75   # Radius of an object to caluculate accessibility
    nump   = 30   # number of particles to caluculate accessibility (resolution)
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
    parser = OptionParser(
        usage=('%prog [options]'))
    modelo = OptionGroup(
        parser, 'Modeling options\n  ' + ':' * 15,
        'These options control the chromatin final modeling.')
    optimo = OptionGroup(
        parser, 'Optimization options\n  ' + ':' * 19,
        'These options control the optimization of IMP parameters.')
    descro = OptionGroup(
        parser, 'Descriptive options\n  ' + ':' * 18,
        '''These options are to reference your analysis, they all can be
        grouped into a config file.''')

    #########################################
    # GENERAL
    parser.add_option(
        '--data', dest='data', metavar="PATH", action='append', default=[],
        type=str,help=('path to a directory with a Hi-C experiment matrix. If' +
                       ' used several times, experiments will be summed up. I' +
                       '.e.: --data hic_replicate_1.txt --data hic_replicate_' +
                       '2.txt'))
    parser.add_option(
        '--data_path', dest='data_path', metavar="PATH", action='append',
        default='', type=str,
        help=('path to a directory where to search for Hi-C experiment data f' +
              'iles.'))
    parser.add_option('--xname', dest='xname', metavar="STR", action='append',
                      default=[], type=str,
                      help=('Name of the Hi-C experiments. Same order as data' +
                            '. To be used as the data option. If ommitted exp' +
                            'eriments will be named using the data file name.'))
    parser.add_option('--cfg', dest='cfg', metavar="PATH", action='store',
                      default=None, type=str,
                      help='Path to a config file where can be stored parameters.')
    parser.add_option('--crm', dest='crm', metavar="NAME", help='Chromosome name')
    parser.add_option('--beg', dest='beg', metavar="INT",
                      help='genomic coordinate from which to start modeling')
    parser.add_option('--end', dest='end', metavar="INT", 
                      help='genomic coordinate where to end modeling')
    parser.add_option('--res', dest='res', metavar="INT", 
                      help='Resolution of the Hi-C experiment.')
    parser.add_option('--outdir', dest='outdir', metavar="PATH",
                      default=None,
                      help='out directory for results')
    parser.add_option('--ncpus', dest='ncpus', metavar="INT", default='1',
                      help='[%default] Number of CPUs to use')
    parser.add_option('--extra', dest='extra', metavar="PATH",
                      default='',
                      help='''[%default] extra specification to appear in the
                      name of generated files''')
    parser.add_option('--load_models', dest='load_models',
                      action='store_true', default=False,
                      help=('Load models generated during a previous modeling' +
                            ' from outdir, to go on with the analysis.'))

    #########################################
    # MODELING
    parser.add_option('--load_optimization', dest='load_optimization',
                      action='store_true', default=False,
                      help=('Load models generated during a previous optimiza' +
                            'tion from outdir, to find the best combination o' +
                            'f parameters and go on with the modeling.'))
    parser.add_option('--view_optimization', dest='view_optimization',
                      action='store_true', default=False,
                      help=('Generates optimization plot based on loaded ' +
                            'optimizations and exit. To be used with "view_op' +
                            'timization"'))
    modelo.add_option('--nmodels_mod', dest='nmodels_mod', metavar="INT",
                      default='5000',
                      help=('[%default] Number of models to generate for mode' +
                      'ling'))
    modelo.add_option('--nkeep_mod', dest='nkeep_mod', metavar="INT",
                      default='1000',
                      help='[%default] Number of models to keep for modeling')
    # parser.add_option('--scale', dest='scale', metavar="LIST",
    #                   default="0.01",
    #                   help='''[%default] range of numbers to be test as optimal
    #                   scale value, i.e. 0.005:0.01:0.001 -- Can also pass only one
    #                   number.''')

    #########################################
    # OPTIMIZATION
    parser.add_option('--optimize', dest='optimize', default=False,
                      action='store_true',
                      help='Do the optimization of the region.')
    optimo.add_option('--maxdist', action='store', metavar="LIST",
                      default='400', dest='maxdist',
                      help='''range of numbers to be test as optimal maxdist
                      value, i.e. 400:1000:100 -- Can also pass only one
                      number.''')
    optimo.add_option('--upfreq', dest='upfreq', metavar="LIST",
                      default='0',
                      help='''range of numbers to be test as optimal upfreq
                      value, i.e. 0:1.2:0.3 -- Can also pass only one
                      number.''')
    optimo.add_option('--lowfreq', dest='lowfreq', metavar="LIST",
                      default='0',
                      help='''range of numbers to be test as optimal lowfreq 
                      value, i.e. -1.2:0:0.3 -- Can also pass only one number.
                      ''')
    optimo.add_option('--nmodels_opt', dest='nmodels_opt', metavar="INT",
                      default='500',
                      help='''[%default] Number of models to generate for
                      optimization''')
    optimo.add_option('--nkeep_opt', dest='nkeep_opt', metavar="INT",
                      default='100',
                      help='''[%default] Number of models to keep for
                      optimization''')

    #########################################
    # DESCRIPTION
    descro.add_option('--species', dest='species', metavar="STRING",
                      help='Species name, with no spaces, i.e.: homo_sapiens')
    descro.add_option('--cell', dest='cell', metavar="STRING",
                      help='Cell type name')
    descro.add_option('--exp_type', dest='exp_type', metavar="STRING",
                      help='Experiment type name (i.e.: Hi-C)')
    descro.add_option('--assembly', dest='assembly', metavar="STRING",
                      default=None,
                      help='''NCBI ID of the original assembly
                      (i.e.: NCBIM37 for mouse or NCBI36 for human)''')
    descro.add_option('--enzyme', dest='enzyme', metavar="STRING",
                      default=None,
                      help='''Name of the enzyme used to digest
                      chromatin (i.e. HindIII)''')
    descro.add_option('--identifier', dest='identifier', metavar="STRING",
                      default=None,
                      help='''NCBI identifier of the experiment.''')
    descro.add_option('--project', dest='project', metavar="STRING",
                      default=None,
                      help='''Project name.''')

    parser.add_option_group(optimo)
    parser.add_option_group(modelo)
    parser.add_option_group(descro)
    opts = parser.parse_args()[0]

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
            
        for key in new_opts:
            if not key in opts.__dict__:
                warn('WARNING: parameter "%s" not recognized' % (key))
                continue
            # if opts.__dict__[key]:
            #     continue
            opts.__dict__[key] = new_opts[key]
            print '  setting %12s to %s' % (key, new_opts[key])
    
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
    if opts.view_optimization and not opts.load_optimization:
        print ('ERROR: view_optimization needs to be used together with ' +
               'load_optimization')
        exit(parser.print_help())
    if not opts.load_optimization and not opts.load_optimization :
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
