#! /usr/bin/python

"""
reads a BAM file extracts valid pairs, and from them, computes an array of
biases (1 round ICE) and array of normalized expected counts according to
genomic distance, and an array of columns with poor_bins signal.
"""

import sys
import os
import multiprocessing  as mu
from collections                     import OrderedDict
from cPickle                         import dump, load
from argparse                        import ArgumentParser, SUPPRESS

from pysam                           import AlignmentFile
import numpy as np

from pytadbit.utils.file_handling    import mkdir
from pytadbit.utils.extraviews       import nicer
from pytadbit.mapping.filter         import MASKED
from pytadbit.parsers.hic_bam_parser import printime, print_progress
from pytadbit.utils.hic_filtering    import filter_by_cis_percentage


def read_bam_frag(inbam, filter_exclude, all_bins, sections,
                  resolution, outdir, region, start, end):
    bamfile = AlignmentFile(inbam, 'rb')
    refs = bamfile.references
    try:
        dico = {}
        for r in bamfile.fetch(region=region,
                               start=start - (1 if start else 0), end=end,  # coords starts at 0
                               multiple_iterators=True):
            if r.flag & filter_exclude:
                continue
            crm1 = r.reference_name
            pos1 = r.reference_start + 1
            crm2 = refs[r.mrnm]
            pos2 = r.mpos + 1
            try:
                pos1 = sections[(crm1, pos1 / resolution)]
                pos2 = sections[(crm2, pos2 / resolution)]
            except KeyError:
                continue  # not in the subset matrix we want
            try:
                dico[(pos1, pos2)] += 1
            except KeyError:
                dico[(pos1, pos2)] = 1
        cisprc = {}
        for (i, j), v in dico.iteritems():
            # out.write('%d\t%d\t%d\n' % (i, j, v))
            try:
                if all_bins[i][0] == all_bins[j][0]:
                    cisprc[i][0] += v
                cisprc[i][1] += v
            except KeyError:
                if all_bins[i][0] == all_bins[j][0]:
                    cisprc[i] = [v, v]
                else:
                    cisprc[i] = [0, v]
        out = open(os.path.join(outdir,
                                'tmp_%s:%d-%d.pickle' % (region, start, end)), 'w')
        dump(dico, out)
        out.close()
        out = open(os.path.join(outdir,
                                'tmp_bins_%s:%d-%d.pickle' % (region, start, end)), 'w')
        dump(cisprc, out)
        out.close()
    except Exception, e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print e
        print(exc_type, fname, exc_tb.tb_lineno)


def read_bam(inbam, filter_exclude, resolution, min_count=2500,
             sigma=2, ncpus=8, factor=1, outdir='.', check_sum=False):
    bamfile = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references,
                               [x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm] + 1
    bins = []
    for crm in sections:
        len_crm = sections[crm]
        bins.extend([(crm, i) for i in xrange(len_crm + 1)])

    start_bin = 0
    end_bin   = len(bins) + 1
    total = len(bins)

    total = end_bin - start_bin + 1
    regs = []
    begs = []
    ends = []
    njobs = min(total, 100) + 1
    nbins = total / njobs + 1
    for i in range(start_bin, end_bin, nbins):
        if i + nbins > end_bin:  # make sure that we stop at the right place
            nbins = end_bin - i
        try:
            (crm1, beg1), (crm2, end2) = bins[i], bins[i + nbins - 1]
        except IndexError:
            (crm1, beg1), (crm2, end2) = bins[i], bins[-1]
        if crm1 != crm2:
            end1 = sections[crm1]
            beg2 = 0
            regs.append(crm1)
            regs.append(crm2)
            begs.append(beg1 * resolution)
            begs.append(beg2 * resolution)
            ends.append(end1 * resolution + resolution)  # last nt included
            ends.append(end2 * resolution + resolution - 1)  # last nt not included (overlap with next window)
        else:
            regs.append(crm1)
            begs.append(beg1 * resolution)
            ends.append(end2 * resolution + resolution - 1)
    ends[-1] += 1  # last nucleotide included

    # print '\n'.join(['%s %d %d' % (a, b, c) for a, b, c in zip(regs, begs, ends)])
    printime('\n  - Parsing BAM (%d chunks)' % (len(regs)))
    bins_dict = dict([(j, i) for i, j in enumerate(bins)])
    pool = mu.Pool(ncpus)
    procs = []
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        procs.append(pool.apply_async(
            read_bam_frag, args=(inbam, filter_exclude, bins, bins_dict,
                                 resolution, outdir, region, start, end,)))
    pool.close()
    print_progress(procs)
    pool.join()

    ## COLLECT RESULTS
    verbose = True
    cisprc = {}
    for countbin, (region, start, end) in enumerate(zip(regs, begs, ends)):
        if verbose:
            if not countbin % 10 and countbin:
                sys.stdout.write(' ')
            if not countbin % 50 and countbin:
                sys.stdout.write(' %9s\n     ' % ('%s/%s' % (countbin , len(regs))))
            sys.stdout.write('.')
            sys.stdout.flush()

        fname = os.path.join(outdir,
                             'tmp_bins_%s:%d-%d.pickle' % (region, start, end))
        tmp_cisprc = load(open(fname))
        cisprc.update(tmp_cisprc)
    if verbose:
        print '%s %9s\n' % (' ' * (54 - (countbin % 50) - (countbin % 50) / 10),
                            '%s/%s' % (len(regs),len(regs)))

    # out = open(os.path.join(outdir, 'dicos_%s.pickle' % (
    #     nicer(resolution).replace(' ', ''))), 'w')
    # dump(cisprc, out)
    # out.close()
    # bad columns
    def func_gen(x, *args):
        cmd = "zzz = " + func_restring % (args)
        exec(cmd) in globals(), locals()
        #print cmd
        try:
            return np.lib.asarray_chkfinite(zzz)
        except:
            # avoid the creation of NaNs when invalid values for power or log
            return x
    print '  - Removing columns with too few or too much interactions'
    if not min_count:

        badcol = filter_by_cis_percentage(
            cisprc, sigma=sigma, verbose=True,
            savefig=os.path.join(outdir + 'filtered_bins_%s.png' % (
                nicer(resolution).replace(' ', ''))))
    else:
        print '      -> too few  interactions defined as less than %9d interactions' % (
            min_count)
        for k in cisprc:
            cisprc[k] = cisprc[k][1]
        badcol = {}
        countL = 0
        countZ = 0
        for c in xrange(total):
            if cisprc.get(c, 0) < min_count:
                badcol[c] = cisprc.get(c, 0)
                countL += 1
                if not c in cisprc:
                    countZ += 1
        print '      -> removed %d columns (%d/%d null/high counts) of %d (%.1f%%)' % (
            len(badcol), countZ, countL, total, float(len(badcol)) / total * 100)

    printime('  - Rescaling biases')
    size = len(bins)
    biases = [cisprc.get(k, 1.) for k in range(size)]
    mean_col = float(sum(biases)) / len(biases)
    biases = dict([(k, b / mean_col * mean_col**0.5)
                   for k, b in enumerate(biases)])

    # collect subset-matrices and write genomic one
    # out = open(os.path.join(outdir,
    #                         'hicdata_%s.abc' % (nicer(resolution).replace(' ', ''))), 'w')
    pool = mu.Pool(ncpus)
    procs = []
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        fname = os.path.join(outdir, 'tmp_%s:%d-%d.pickle' % (region, start, end))
        procs.append(pool.apply_async(sum_nrm_matrix, args=(fname, biases, )))
    pool.close()
    print_progress(procs)
    pool.join()

    # to correct biases
    sumnrm = sum(p.get() for p in procs)

    target = (sumnrm / float(size * size * factor))**0.5
    biases = dict([(b, biases[b] * target) for b in biases])

    # check the sum
    if check_sum:
        pool = mu.Pool(ncpus)
        procs = []
        for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
            fname = os.path.join(outdir, 'tmp_%s:%d-%d.pickle' % (region, start, end))
            procs.append(pool.apply_async(sum_nrm_matrix, args=(fname, biases, )))
        pool.close()
        print_progress(procs)
        pool.join()

        # to correct biases
        sumnrm = sum(p.get() for p in procs)
        print 'SUM:', sumnrm

    printime('  - Rescaling decay')
    # normalize decay by size of the diagonal, and by Vanilla correction
    # (all cells must still be equals to 1 in average)

    pool = mu.Pool(ncpus)
    procs = []
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        fname = os.path.join(outdir,
                             'tmp_%s:%d-%d.pickle' % (region, start, end))
        procs.append(pool.apply_async(sum_dec_matrix,
                                      args=(fname, biases, badcol, bins)))
    pool.close()
    print_progress(procs)
    pool.join()

    # collect results
    sumdec = {}
    for proc in procs:
        for k, v in proc.get().iteritems():
            try:
                sumdec[k] += v
            except KeyError:
                sumdec[k]  = v

    # count the number of cells per diagonal
    # TODO: parallelize
    # find larget chromsome
    len_big = max(section_pos[c][1] - section_pos[c][0] for c in section_pos)
    # initialize dictionary
    ndiags = dict((k, 0) for k in xrange(len_big))
    for crm in section_pos:
        beg_chr, end_chr = section_pos[crm][0], section_pos[crm][1]
        chr_size = end_chr - beg_chr
        thesebads = [b for b in badcol if beg_chr <= b <= end_chr]
        for dist in xrange(1, chr_size):
            ndiags[dist] += chr_size - dist
            # from this we remove bad columns
            # bad columns will only affect if they are at least as distant from
            # a border as the distance between the longest diagonal and the
            # current diagonal.
            bad_diag = set()  # 2 bad rows can point to the same bad cell in diagonal
            maxp = end_chr - dist
            minp = beg_chr + dist
            for b in thesebads:
                if b <= maxp:
                    bad_diag.add(b)
                if b >= minp:
                    bad_diag.add(b - dist)
            ndiags[dist] -= len(bad_diag)
        # chr_sizeerent behavior for longest diagonal:
        ndiags[0] += chr_size - len(thesebads)

    # normalize sum per diagonal by total number of cells in diagonal
    for k in sumdec:
        try:
            sumdec[k] /= ndiags[k]
        except ZeroDivisionError:  # all columns at this distance are "bad"
            pass

    return biases, sumdec, badcol


def sum_dec_matrix(fname, biases, badcol, bins):
    dico = load(open(fname))
    sumdec = {}
    for (i, j), v in dico.iteritems():
        # different chromosome
        if bins[i][0] != bins[j][0]:
            continue
        if i < j:
            continue
        if i in badcol or j in badcol:
            continue
        k = i - j
        val = v / biases[i] / biases[j]
        try:
            sumdec[k] += val
        except KeyError:
            sumdec[k]  = val
    os.system('rm -f %s' % (fname))
    return sumdec


def sum_nrm_matrix(fname, biases):
    dico = load(open(fname))
    sumnrm = sum(v / biases[i] / biases[j]
                 for (i, j), v in dico.iteritems())
    return sumnrm


def main():
    opts          = get_options()

    inbam          = opts.inbam
    resolution     = opts.reso
    filter_exclude = opts.filter
    min_count      = opts.min_count
    ncpus          = opts.cpus
    factor         = 1
    outdir         = opts.outdir
    sigma          = 2

    mkdir(outdir)

    sys.stdout.write('\nNormalization of full genome\n')

    biases, decay, badcol = read_bam(inbam, filter_exclude, resolution,
                                     min_count=min_count, ncpus=ncpus, sigma=sigma,
                                     factor=factor, outdir=outdir, check_sum=opts.check_sum)

    printime('  - Saving biases and badcol columns')
    # biases
    out = open(os.path.join(outdir, 'biases_%s.pickle' % (
        nicer(resolution).replace(' ', ''))), 'w')

    dump({'biases'    : biases,
          'decay'     : decay,
          'badcol'    : badcol,
          'resolution': resolution}, out)
    out.close()

    # hic_data.write_matrix('chr_names%s_%d-%d.mat' % (region, start, end), focus=())
    printime('\nDone.')


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbam', metavar='',
                        required=True, default=False, help='input HiC-BAM file.')
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='',
                        required=True, default=False, help='output directory.')
    parser.add_argument('-r', '--resolution', dest='reso', type=int, metavar='',
                        required=True, help='''wanted resolution form the
                        generated matrix''')
    parser.add_argument('--check_sum', dest='check_sum',
                        action='store_true', default=False, help=SUPPRESS
                        # '''print the sum_dec_matrix of the normalized matrix and exit'''
    )
    parser.add_argument('--min_count', dest='min_count', type=int, metavar='',
                        default=None,
                        help='''[%(default)s] minimum number of interactions
                        perc_zero bin. By default this value is optimized.''')
    parser.add_argument('-C', '--cpus', dest='cpus', metavar='', type=int,
                        default=8, help='''[%(default)s] number of cpus to be
                        used for parsing the HiC-BAM file''')
    parser.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 8, 9, 10],
                        choices = range(1, 11),
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))

    opts = parser.parse_args()
    # convert filters to binary for samtools
    opts.filter = filters_to_bin(opts.filter)

    return opts


def filters_to_bin(filters):
    return sum((k in filters) * 2**(k-1) for k in MASKED)


if __name__=='__main__':
    exit(main())
