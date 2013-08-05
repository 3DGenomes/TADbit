"""
04 Jul 2013
"""
from pytadbit import Chromosome
from pytadbit.imp.imp_model import generate_3d_models
from optparse   import OptionParser
from pytadbit.imp.CONFIG import CONFIG


def parse_zscores(zscore_f):
    """
    """
    zscores = {}
    for line in open(zscore_f):
        x, y, z = line.split()
        z = float(z)
        zscores.setdefault(x, {})
        zscores[x][y] = z
        zscores.setdefault(y, {})
        zscores[y][x] = z
    return zscores
    

def main():

    opts, params = get_options()
    if opts.inabc:
        zscores = parse_zscores(opts.inabc)
        models = generate_3d_models(zscores, opts.resolution, start=1,
                                    n_models=opts.nmodels,
                                    n_keep=opts.nkeep, n_cpus=opts.ncpus,
                                    keep_all=False, verbose=False,
                                    outfile=None,
                                    config=params)
    
    else:
        crm  = 'crm'
        xnam = 'X'
        crmbit=Chromosome(crm)
        crmbit.add_experiment(xnam, resolution=opts.resolution, xp_handler=opts.incrm)
        exp = crmbit.experiments[xnam]
        models = exp.model_region(start=opts.start, end=opts.end,
                                  n_models=opts.nmodels,
                                  n_keep=opts.nkeep, n_cpus=opts.ncpus,
                                  keep_all=False, verbose=False,
                                  config=params)

    if opts.save:
        models.save_models('%s/models_%s_%s.pik' % (opts.out, opts.start,
                                                    opts.start + opts.nmodels))
    for i in xrange(int(opts.cmm)):
        models.write_cmm(i, opts.out)

    if opts.full_report:
        
        models.cluster_models(dcutoff=200)
        models.cluster_analysis_dendrogram(n_best_clusters=10)
        models.model_consistency()

   

def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        usage=("%prog [options] file [options] file [options] " +
               "file [options [file ...]]"))
    parser.add_option('--incrm', dest='incrm', metavar="PATH",
                      default=None,
                      help='''path to a hi-c matrix file.''')
    parser.add_option('--inabc', dest='inabc', metavar="PATH",
                      default=None,
                      help='''path to a zscore file.''')
    parser.add_option('--out', dest='out', metavar="PATH",
                      default=None,
                      help='''Required: path to a directory store results''')
    parser.add_option('--save_models', dest='save', metavar="PATH",
                      default=False, action='store_true',
                      help='''Writes models to file.''')
    parser.add_option('--full_report', dest='full_report',default=False,
                      action='store_true',
                      help='''Do a fully automatic analysis on the models
                      generated in order to evaluate the quality of the
                      data, plots stored in out directory.''')
    parser.add_option('--cmm', dest='cmm', metavar="PATH",
                      default=0,
                      help='''[%default] number of CMM files to generate,
                      one per each of the best models.
                      (CMM files are input files for Chimera program).''')
    parser.add_option('--start', dest='start',
                      default=1,
                      help='''Start position for modelling (bin number
                      of the Hi-C matrix).''')
    parser.add_option('--end', dest='end',
                      default=None,
                      help='''End  position for modelling (bin number
                      of the Hi-C matrix).''')
    parser.add_option('--nmodels', dest='nmodels',
                      default=1,
                      help='''Number of models to generate''')
    parser.add_option('--nkeep', dest='nkeep',
                      default=1,
                      help='''Number of best models to keep from
                      the ones generated (usually 20%)''')
    parser.add_option('--resolution', dest='resolution',
                      default=1,
                      help='''resolution of the Hi-C experiment
                      (nucleotides/bin)''')
    parser.add_option('--ncpus', dest='ncpus',
                      default=1,
                      help='''number of cpus to use.''')
    parser.add_option('--param', dest='param',
                      default='dmel_01', type='choice',
                      choices=['dmel_01'],
                      help='''Precalculated set of parameters to use
                      (available: dmel_01)''')
    parser.add_option('--describe_param', dest='describe_param',
                      default=None, type='choice',
                      choices=['dmel_01'],
                      help='''Show a given configuration and exit (available:
                      dmel_01).''')
    parser.add_option('--consdist', dest='consdist',
                      default=None,
                      help="""Maximum experimental contact distance""")
    parser.add_option('--kforce', dest='kforce',
                      default=None,
                      help="""Force applied to the restraints inferred to
                      neighbor particles""")
    parser.add_option('--lowrdist', dest='lowrdist',
                      default=None,
                      help= """Minimum distance between two non-bonded
                      particles""")
    parser.add_option('--upfreq', dest='upfreq',
                      default=None,
                      help='''Maximum thresholds used to decide which
                      experimental values have to be included in the
                      computation of restraints. Z-score values bigger
                      than upfreq and less that lowfreq will be include,
                      whereas all the others will be rejected''')
    parser.add_option('--lowfreq', dest='lowfreq',
                      default=None,
                      help="""Minimum thresholds used to decide which
                      experimental values have to be included in the
                      computation of restraints. Z-score values bigger
                      than upfreq and less that lowfreq will be include,
                      whereas all the others will be rejected""")


                      
    opts = parser.parse_args()[0]
    if opts.describe_param:
        print 'parameters for %s:\n' % (opts.describe_param)
        for key, val in CONFIG[opts.describe_param].iteritems():
            print '%10s    %-10s' % (key, val)
        exit()
    if not opts.incrm or not opts.inabc:
        exit(parser.print_help())

    if not opts.param:
        perso = {}
        perso['consdist'] = opts.consdist
        perso['kforce']   = opts.kforce
        perso['lowrdist'] = opts.lowrdist
        perso['upfreq']   = opts.upfreq
        perso['lowfreq']  = opts.lowfreq
    else:
        perso = CONFIG[opts.describe_param]

    if not opts.out:
        print 'Should provide output path'
        exit(parser.print_help())
        
    return opts, perso


if __name__ == "__main__":
    exit(main())
