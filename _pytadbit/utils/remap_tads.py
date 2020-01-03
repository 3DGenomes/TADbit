"""
19 Sep 2013


LiftOver (1) wrapper applied to the comparison of topologically associated
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
from __future__ import print_function



from subprocess import Popen, PIPE
from os         import system


def liftover(coords, tmp_path, lft_path, chn_path):
    """
    LiftOver wrapper

    :param coords: a list of TADs each tad being another list with its name,
       its chromosome, and its position (tad['end'])
    :param tmp_path: path where to create temporary files for the LiftOver
       program
    :param chn_path: path to a 'chain' file needed by LiftOver to convert
       coordinates.

    :returns: a dictionary which keys are TAD names, and values are the new
       coordinates (chromosome name, genomic position).
    """
    tmp = open(tmp_path + '/tmp', 'w')
    for coord in coords:
        tmp.write('chr{}\t{}\t{}\n'.format(coord[1], coord[2], coord[2]+1))
    tmp.close()
    _, err = Popen((lft_path + 'liftOver ' +
                    tmp_path + '/tmp ' +
                    chn_path + ' ' +
                    tmp_path + '/out ' +
                    tmp_path + '/out.log'
                    ), shell =True, stdout=PIPE, stderr=PIPE).communicate()
    if (not 'Reading' in err) or (not 'Mapping' in err):
        raise Exception(err)
    founds = [[l.split()[0],
               int(l.split()[1])] for l in open(tmp_path + '/out').readlines()]
    missed = [int(l.split()[1]) for l in open(tmp_path + '/out.log').readlines()\
              if not l.startswith('#')]
    system('rm -f {}/tmp'.format(tmp_path))
    mapped = {}
    j = 0
    for k, (i, _, end) in enumerate(coords):
        mapped[i] = None
        if end in missed:
            j += 1
            continue
        mapped[i] = (founds[k - j][0][3:], founds[k - j][1])
    return mapped
    
    
def remap_chr(crm_obj, crm, tmp, lft_path, chain_path,
              wnt_exp=None, genome=None):
    """
    """
    missed = 0
    found  = 0
    genome = {} or genome
    for exp in crm_obj.experiments:
        if wnt_exp:
            if not exp.name in wnt_exp:
                continue
        print(exp)
        coords = []
        res = exp.resolution
        for t in exp.tads:
            coords.append((t, crm,
                           int(exp.tads[t]['end']   * res)))
        mapped = liftover(coords, tmp, lft_path, chain_path)
        for t in mapped:
            try:
                crm2, end = mapped[t]
                genome.setdefault(exp.name, {})
                if crm2 not in genome[exp.name]:
                    genome[exp.name][crm2] = {'end': [],
                                             'score': []}
                print(' -> crm{:<2}:{:>10} | crm{:<2}:{:>10}'.format(
                    crm, int(exp.tads[t]['end']) * 100000,
                    crm2, int(end / res + 0.5) * 100000))
                genome[exp.name][crm2]['end'  ].append(float(int(end / res + 0.5)))
                genome[exp.name][crm2]['score'].append(exp.tads[t]['score'])
                found += 1
            except TypeError:
                missed += 1
    print('missed: {}, found: {}'.format(missed, found))
    return genome


def reorder(genome):
    for exp in genome:
        for crm in genome[exp]:
            srt = sorted(genome[exp][crm]['end'])
            order = [genome[exp][crm]['end'].index(elt) for elt in srt]
            genome[exp][crm]['end']   = [genome[exp][crm]['end'  ][i] for i in order]
            genome[exp][crm]['score'] = [genome[exp][crm]['score'][i] for i in order]
            todel = []
            for brk in range(1, len(genome[exp][crm]['end'])):
                if genome[exp][crm]['end'][brk] == genome[exp][crm]['end'][brk-1]:
                    todel.append(brk-1)
                    genome[exp][crm]['score'][brk] = max(
                        genome[exp][crm]['score'][brk], genome[exp][crm]['score'][brk-1])
            for brk in todel[::-1]:
                genome[exp][crm]['end'].pop(brk)
                genome[exp][crm]['score'].pop(brk)
            genome[exp][crm]['start'] = []
            for brk in range(len(genome[exp][crm]['end']) + 1):
                genome[exp][crm]['start'].append(genome[exp][crm]['end'][brk - 1] + 1 \
                                                 if brk > 0 else 0)
                    
