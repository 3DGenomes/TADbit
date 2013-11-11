"""
script to calculate the accessibilty, per particle, of a given model.
"""

from pytadbit.imp.impmodel import load_impmodel_from_cmm
from pytadbit.imp.impmodel import load_impmodel_from_xyz
from optparse                  import OptionParser



def main():

    opts = get_options()
    if opts.model_path.endswith('.cmm'):
        model = load_impmodel_from_cmm(opts.model_path)
    elif opts.model_path.endswith('.xyz'):
        model = load_impmodel_from_xyz(opts.model_path)

    nump   = int(opts.nump)
    radius = int(opts.radius)
    (accesible_parts,
     number_parts,
     accessible_area,
     total_area, by_part) = model.accessible_surface(
        radius, nump=nump, verbose=False, write_cmm_file=opts.outfcmm,
        include_edges=opts.slow)

    print '\n GLOBAL STATS:\n -------------\n'
    print '  - accesible dots (of the mesh)     : %s' % accesible_parts
    print '  - total number dots calculated     : %s' % number_parts
    print '  - proportion of accesible dots     : %.2f%%' % (100*float(accesible_parts)/number_parts)
    print '  - accessible area                  : %.4f square micrometers' % accessible_area
    print '  - total area (straight chromatin)  : %.4f square micrometers' % total_area     
    print '  - proportion of accesible surface  : %.2f%%' % (100*float(accessible_area)/total_area)
    print ''

    out = open(opts.outf, 'w')
    for part, acc, ina in by_part:
        if opts.burry:
            burried = 1 if (100*float(ina) / (acc + ina)) > float(opts.burry) else 0
            out.write('%s\t%.2f\n' % (part, burried))
        else:
            out.write('%s\t%.2f\n' % (part, 100*float(ina) / (acc + ina)))
    out.close()


def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        usage=("%prog [options] file [options] file [options] " +
               "file [options [file ...]]"))
    parser.add_option('--infile', dest='model_path', metavar="PATH",
                      action='store', default=None,
                      help='''path to a cmm or xyz file''')
    parser.add_option('--radius', dest='radius', default=300,
                      help='''[%default] radius of the object to fit''')
    parser.add_option('--dens', dest='nump', default=100, 
                      help='''[%default] density of the mesh, 100 mean 100
                      dots checked around each particle''')
    parser.add_option('--burry', dest='burry', default=None, 
                      help='''[%default] density of the mesh, 100 mean 100 dots
                      checked around each particle. 70, means 70% of the
                      particle is not accessible''')
    parser.add_option('--fast', dest='slow', default=True, action='store_false',
                      help='''skip computation of the accessibility of edges''')
    parser.add_option('--outfile', dest='outf', action="store", metavar="PATH", 
                      help=
                      '''path to outfile, 2 columns, particle number, and
                      percentage of burriement (if burriement value is passed,
                      than only 0 or 1)''')
    parser.add_option('--outfcmm', dest='outfcmm', metavar="PATH", default=None,
                      help='''path to second outfile, to visualize with chimera
                      the result. green=accessible, red=inaccessible''')
    opts = parser.parse_args()[0]
    if not opts.model_path:
        exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
