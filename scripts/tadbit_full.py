"""
21 Nov 2013

This script contains the main analysis that can be done using TADbit:
  * visualization of Hi-C data
  * detection of TADs
  * alignment of TADs
  * optimization of IMP parameters
  * building of models of chromatin structure
  
"""

from optparse                  import OptionParser, Option
from pytadbit import Chromosome
from warnings import warn


def load_experiments(opts):
    crm = Chromosome(opts.crm)
    for i, xpr in enumerate(opts.hic_files):
        if opts.exp_names:
            name = opts.exp_names[i]
        else:
            name = ''.join(xpr.split('/')[-1].split('.')[:-1])
        if opts.verbose:
            print ' Reading Hi-C datafile #%s (%s)' % (i+1, name)
        crm.add_experiment(name, hic_data=xpr,
                           resolution=int(opts.resolution))
        if opts.verbose:
            print '     loaded as: %s\n' % (crm.experiments[name])
    return crm


def search_tads(crm, opts):
    for xpr in crm.experiments:
        crm.find_tad(xpr, n_cpus=int(opts.ncpus), verbose=True)
        




def main():
    """
    main function
    """
    opts = get_options()

    if opts.verbose:
        print '\n\nStarting Analysis\n' + '=' * 17 + '\n'

    crm = load_experiments(opts)
    
    if opts.tads:
        search_tads(crm, opts)

    if opts.align:
        ali = crm.align_experiments()
    


def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        usage=("%prog [options] file [options] file [options] " +
               "file [options [file ...]]"), option_class=MultipleOption)
    parser.add_option('-c', '--chr_name', action='store', default='C',
                      dest='crm', help='[%default] chromosome name',
                      metavar='str')
    parser.add_option('-f', '--hic_files', dest='hic_files', metavar="PATH",
                      action='extend', default=None,
                      help='''path to a file containing Hi-C data in matrix
                      format. option can be called multiple times, for multiple
                      Hi-C data files.''')
    parser.add_option('-e', '--exp_names', dest='exp_names', metavar="str",
                      action='extend', default=None,
                      help='''name of Hi-C experiments passed using the
                      hic_files option. Should be in the same order. If not
                      used, names will be automatically generated using file
                      names''')
    parser.add_option('-r', '--resolution', dest='resolution', metavar="int",
                      action='store', default=None,
                      help='''resolution of Hi-C experiments passed using the
                      hic_files option. Should be the same for all experiments.''')
    parser.add_option('-a', '--align', action='store_true', dest='align',
                      default=False,
                      help='''[%default] If called, align all Hi-C experiments
                      passed.''')
    parser.add_option('--outdir', dest='outd', action="store", metavar="PATH", 
                      default='', help=
                      '''path to out-directory where to store generated outfiles
                      , images, text files, tables...''')
    parser.add_option('--tads', dest='tads', action='store_true',
                      help='Detect TADs in each experiments separately')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true',
                      default=False, help=
                      '''display progress messages.''')
    parser.add_option('--ncpus', dest='ncpus', metavar="int",
                      action='store', default=1,
                      help='''Number of CPUs to be used for the analysis''')
    opts = parser.parse_args()[0]
    if not opts.hic_files:
        exit(parser.print_help())
    if not opts.resolution:
        warn('ERROR: should provide resolution')
        exit(parser.print_help())
    return opts


class MultipleOption(Option):
    ACTIONS = Option.ACTIONS + ("extend",)
    STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
    TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
    ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

    def take_action(self, action, dest, opt, value, values, parser):
        if action == "extend":
            values.ensure_value(dest, []).append(value)
        else:
            Option.take_action(self, action, dest, opt, value, values, parser)


if __name__ == "__main__":
    exit(main())
