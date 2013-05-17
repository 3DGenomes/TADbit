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
    Nucleic Acids Research, 39(Database issue), D876–82. doi:10.1093/nar/gkq963
    
"""

from subprocess import Popen, PIPE
from os import system
from random import random
from optparse import OptionParser

TMP_PATH = '/tmp/'
LFT_PATH = '/home/fransua/Tools/liftOver/'
CHN_PATH = '/home/fransua/Tools/liftOver/'


def liftover(coords):
    tmp = open(TMP_PATH + 'tmp', 'w')
    for coord in coords:
        tmp.write('chr{}\t{}\t{}\n'.format(coord[1], coord[2], coord[2]+1))
        tmp.write('chr{}\t{}\t{}\n'.format(coord[1], coord[3], coord[3]+1))
    tmp.close()
    _, err = Popen((LFT_PATH + 'liftOver ' +
                      TMP_PATH + 'tmp ' +
                      CHN_PATH + 'hg18ToHg19.over.chain ' +
                      TMP_PATH + '/out ' +
                      TMP_PATH + '/out.log'
                      ), shell =True, stdout=PIPE, stderr=PIPE).communicate()
    if (not 'Reading' in err) or (not 'Mapping' in err):
        raise Exception(err)
    system('rm -f {}tmp'.format(TMP_PATH))
    founds = [int(l.split()[1]) for l in open(TMP_PATH + 'out').readlines()]
    missed = [int(l.split()[1]) for l in open(TMP_PATH + 'out.log').readlines()\
              if not l.startswith('#')]
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
        mapped[i] = (crm, founds[k * 2 - j], founds[k * 2 + 1 - j])
    return mapped, relaunch
    
    
def remap_chr(crm_obj, crm):
    missed = 0
    found  = 0
    max_dist = 20
    for exp in crm_obj.experiments:
        coords = []
        res = exp.resolution
        for t in exp.tads:
            coords.append((t, crm,
                           int(exp.tads[t]['start'] * res),
                           int(exp.tads[t]['end']   * res)))
        mapped, relaunch = liftover(coords)
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
            new_mapped, relaunch = liftover(new_coords)
            mapped.update(new_mapped)
            if not relaunch:
                break
        for t in mapped:
            try:
                _, beg, end = mapped[t]
                if random() < .5:
                    print exp.tads[t]['start'], float(int(beg/res + 0.5))
                exp.tads[t]['start'] = float(int(beg/res + 0.5))
                exp.tads[t]['end'] = float(int(end/res + 0.5))
                if exp.tads[t]['brk'] >= 0:
                    exp.tads[t]['brk'] = exp.tads[t]['end']
                found += 1
            except TypeError:
                missed += 1
                exp.tads[t]['brk'] = None
        exp.brks = [t['brk'] for t in exp.tads.values() if t['brk']]
        yield exp
    print 'missed: {}, found: {}'.format(missed, found)



def main():
    """
    main function
    """
    pass


def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        usage="%prog [options] file [options [file ...]]",
        )
    parser.add_option('--crm1', dest='crm1', metavar="PATH", 
                      help='path to input file, a chromosome saved through tadbit')
    parser.add_option('--crm2', dest='crm2', metavar="PATH", 
                      help='path to second input file, a chromosome saved through tadbit')
    parser.add_option('--tmp', dest='tmp_path', metavar="PATH",
                      default='./',
                      help='path to temporary directory to store liftover outfiles')
    parser.add_option('--liftover_path', 
                      dest='lft_path', default='/usr/local/bin',\
                      help=
                      '''path to liftover binary
                      ''')
    parser.add_option('--chain', dest='chain_path', action="store", \
                      help=
                      '''path to UCSC chain file''')
    opts = parser.parse_args()[0]
    if not opts.crm1 or not opts.crm2 or not opts.chain_path:
        exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
