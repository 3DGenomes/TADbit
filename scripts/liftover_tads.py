"""
17 May 2013


Liftover (1) wrapper applied to the comparison of topologically associated
domains.


This script allows to compare Hi-C experiments (mainly align TAD boundaries)
done with different assemblies (e.g.: NCBI36 and GRCh37 for human genome), or
in different species.


INSTALL:

 - liftover tool needs to be downloaded from
   (http://hgdownload.cse.ucsc.edu/admin/exe/), installed, and appended to the
   path.
  - depending on the data a 'chain' file may also be downloaded. For example
    from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/


(1) Fujita, P. A., Rhead, B., Zweig, A. S., Hinrichs, A. S., Karolchik, D.,
    Cline, M. S., Goldman, M., et al. (2011).
    The UCSC Genome Browser database: update 2011.
    Nucleic Acids Research, 39(Database issue), D876-82. doi:10.1093/nar/gkq963
    
"""

from subprocess import Popen, PIPE
from os import system
from random import random
from optparse import OptionParser
from pytadbit import load_chromosome, Experiment


def liftover(coords, tmp_path, lft_path, chn_path):
    tmp = open(tmp_path + '/tmp', 'w')
    for coord in coords:
        tmp.write('chr{}\t{}\t{}\n'.format(coord[1], coord[2], coord[2]+1))
        tmp.write('chr{}\t{}\t{}\n'.format(coord[1], coord[3], coord[3]+1))
    tmp.close()
    _, err = Popen((lft_path + 'liftOver ' +
                    tmp_path + '/tmp ' +
                    chn_path + ' ' +
                    tmp_path + '/out ' +
                    tmp_path + '/out.log'
                    ), shell =True, stdout=PIPE, stderr=PIPE).communicate()
    if (not 'Reading' in err) or (not 'Mapping' in err):
        raise Exception(err)
    founds = [int(l.split()[0:1]) for l in open(tmp_path + '/out').readlines()]
    missed = [int(l.split()[1]) for l in open(tmp_path + '/out.log').readlines()\
              if not l.startswith('#')]
    system('rm -f {}/tmp'.format(tmp_path))
    relaunch = {}
    mapped = {}
    j = 0
    for k, (i, crm, beg, end) in enumerate(coords):
        mapped[i] = None
        if beg in missed:
            relaunch[i] = (i, crm, beg, end, 0)
            j += 2
            continue
        if end in missed:
            relaunch[i] = (i, crm, beg, end, 1)
            j += 2
            continue
        mapped[i] = (founds[0][3:], founds[1][k * 2 - j], founds[1][k * 2 + 1 - j])
    return mapped, relaunch
    
    
def remap_chr(crm_obj, crm, tmp, lft_path, chain_path, genome=None):
    missed = 0
    found  = 0
    max_dist = 20
    genome = {} or genome
    for exp in crm_obj.experiments:
        coords = []
        res = exp.resolution
        for t in exp.tads:
            coords.append((t, crm,
                           int(exp.tads[t]['start'] * res),
                           int(exp.tads[t]['end']   * res)))
        mapped, relaunch = liftover(coords, tmp, lft_path, chain_path)
        if not relaunch:
            continue
        for i in xrange(2, max_dist):
            new_coords = []
            for key, _, _, _, m in relaunch.values():
                _, crm, beg, end = coords[key]
                sign = -1 if i % 2 else 1
                if m == 1:
                    coord = (key, crm, beg, end + res/10 * (i/2) * sign)
                else:
                    coord = (key, crm, beg + res/10 * (i/2) * sign, end)
                new_coords.append(coord)
            new_mapped, relaunch = liftover(new_coords, tmp, lft_path, chain_path)
            mapped.update(new_mapped)
            if not relaunch:
                break
        for t in mapped:
            try:
                crm, beg, end = mapped[t]
                if crm not in genome:
                    genome[crm] = {}
                if random() < .5:
                    print exp.tads[t]['start'], float(int(beg / res + 0.5))
                genome[crm].setdefault('start', []).append(float(int(beg / res + 0.5)))
                genome[crm].setdefault('end', []).append(float(int(end / res + 0.5)))
                if exp.tads[t]['brk'] >= 0:
                    genome[crm].setdefault('brk', []).append(genome['end'][-1])
                else:
                    genome[crm].setdefault('brk', []).append(None)
                found += 1
            except TypeError:
                missed += 1
                genome[crm].setdefault('brk', []).append(None)
    print 'missed: {}, found: {}'.format(missed, found)
    return genome



def main():
    """
    main function
    """
    opts = get_options()
    crm_a = load_chromosome(opts.crm1)
    crm_b = load_chromosome(opts.crm2)
    system('mkdir -p ' + opts.tmp_path)
    genome = {}
    genome = remap_chr(crm_b, opts.crm_name if opts.crm_name else  crm_b.name,
                       opts.tmp_path, opts.lft_path, opts.chain_path, genome)
    for crm in genome:
        if crm == opts.crm_name:
            tads = []
            for t in genome[crm]:
                crm_a.add_experiment(exp)
    crm_a.save_chromosome(opts.out_path)


def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        usage=("%prog [options] file [options] file [options] " +
               "file [options [file ...]]"))
    parser.add_option('--crm1', dest='crm1', metavar="PATH", 
                      help='''path to input file, a chromosome saved through
                      tadbit (required)''')
    parser.add_option('--crm2', dest='crm2', metavar="PATH", 
                      help='''path to second input file, a chromosome saved
                      through tadbit (required)''')
    parser.add_option('--chain', dest='chain_path', action="store", \
                      help=
                      '''path to UCSC chain file (required)''')
    parser.add_option('-o', dest='out_path', metavar="PATH",
                      default='./',
                      help='''path to out file where merged tadbit chromosome
                      will be stored''')
    parser.add_option('--crm_name', dest='crm_name',
                      default=None,
                      help='''Chromosome name for crm1 (e.g. 21).''')
    parser.add_option('--tmp', dest='tmp_path', metavar="PATH",
                      default='./',
                      help='''path to temporary directory to store liftover
                      outfiles''')
    parser.add_option('--liftover_path', 
                      dest='lft_path', default='/usr/local/bin',\
                      help=
                      '''path to liftover binary
                      ''')
    opts = parser.parse_args()[0]
    if not opts.crm1 or not opts.crm2 or not opts.chain_path:
        exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
