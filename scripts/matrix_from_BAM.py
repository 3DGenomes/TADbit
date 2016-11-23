#! /usr/bin/python

"""
extract subset-matrix from a BAM file, and evantually normalizes it using
 precomputed biases
"""

from argparse                     import ArgumentParser
from pytadbit.utils.file_handling import mkdir
from cPickle                      import load, dump
from time                         import sleep, time
from collections                  import OrderedDict
from pytadbit.utils.extraviews    import nicer
import pysam
import datetime
import sys, os
import multiprocessing as mu

def printime(msg):
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')


def print_progress(procs):
    sys.stdout.write('     ')
    prev_done = done = i = 0
    while done < len(procs):
        sleep(2)
        done = sum(p.ready() for p in procs)
        for i in range(prev_done, done):
            if not i % 10 and i:
                sys.stdout.write(' ')
            if not i % 50 and i:
                sys.stdout.write(' %9s\n     ' % ('%s/%s' % (i , len(procs))))
            sys.stdout.write('.')
            sys.stdout.flush()
        prev_done = done    
    print '%s %9s\n' % (' ' * (54 - (i % 50) - (i % 50) / 10),
                        '%s/%s' % (len(procs),len(procs)))


def read_bam_frag(inbam, filter_exclude, sections, 
                  resolution, outdir, region, start, end, half=False):
    bamfile = pysam.AlignmentFile(inbam, 'rb')
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
        if half:
            for i, j in dico.keys():
                if i < j:
                    del(dico[(i,j)])
        out = open(os.path.join(outdir, 'tmp_%s:%d-%d.pickle' % (
            region, start, end)), 'w')
        dump(dico, out)
        out.close()
    except Exception, e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print e
        print(exc_type, fname, exc_tb.tb_lineno)


def read_bam(inbam, filter_exclude, resolution, biases, ncpus=8,
             region=None, start=None, end=None, outdir='.'):
    
    bamfile = pysam.AlignmentFile(inbam, 'rb')
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
    if region:
        regions = [region]
        start_bin = [i for i, b in enumerate(bins) if b[0] == region][0]
        end_bin   = [i for i, b in enumerate(bins[start_bin:], start_bin) if b[0] == region][-1]
    else:
        regions = bamfile.references
        total = len(bins)
        if start or end:
            raise Exception('ERROR: Cannot use start/end1 without region')

    if start:
        start_bin = section_pos[region][0] + start / resolution
    else:
        start = 0
    if end:
        end_bin = section_pos[region][0] + end / resolution
    else:
        end = len(bins)

    total = end_bin - start_bin + 1
    regs = []
    begs = []
    ends = []
    njobs = min(total, 100) + 1
    nbins = total / njobs + 1
    for i in xrange(start_bin, end_bin, nbins):
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
    
    # reduce dictionaries
    bins = []
    for crm in regions:
        beg_crm = section_pos[crm][0]
        if len(regions) == 1:
            start = start_bin - beg_crm
            end   = end_bin   - beg_crm
        else:
            start = 0
            end   = section_pos[crm][1] - section_pos[crm][0] + 1
        bins.extend([(crm, i) for i in xrange(start, end)])
    pool = mu.Pool(ncpus)
    bins_dict = dict([(j, i) for i, j in enumerate(bins)])
    ## RUN!
    printime('\n  - Parsing BAM (%d chunks)' % (len(regs)))
    procs = []
    for i, (region, b, e) in enumerate(zip(regs, begs, ends)):
        if ncpus == 1:
            read_bam_frag(inbam, filter_exclude, bins_dict, 
                          resolution, outdir, region, b, e,)
        else:
            procs.append(pool.apply_async(
                read_bam_frag, args=(inbam, filter_exclude, bins_dict, 
                                     resolution, outdir, region, b, e,)))
    pool.close()
    print_progress(procs)
    pool.join()

    printime('  - Writing matrices')
    bias  = dict((k - start_bin, v)
                 for k, v in biases.get('biases', {}).iteritems()
                 if start_bin <= k <= end_bin)
    decay = biases.get('decay' , {})
    bads  = dict((k - start_bin, v)
                 for k, v in biases.get('badcol', {}).iteritems()
                 if start_bin <= k <= end_bin)
    # hic_data = HiC_data((), len(bins_dict), sections,
    #                     bins_dict, resolution=resolution)
    if len(regions) == 1:
        name = '%s:%d-%d' % (region, start, end)
    else:
        name = 'full'
    out_raw = open(os.path.join(outdir, 'matrix_raw_%s_%s.abc' % (
        name, nicer(resolution).replace(' ', ''))), 'w')
    out_raw.write('# %s resolution:%d\n' % (name, resolution))
    out_raw.write('# BADS %s\n' % (','.join([str(b) for b in bads])))
    if biases:
        out_nrm = open(os.path.join(outdir, 'matrix_nrm_%s_%s.abc' % (
            name, nicer(resolution).replace(' ', ''))), 'w')
        out_nrm.write('# %s resolution:%d\n' % (name, resolution))
        out_nrm.write('# BADS %s\n' % (','.join([str(b) for b in bads])))
        out_dec = open(os.path.join(outdir, 'matrix_dec_%s_%s.abc' % (
            name, nicer(resolution).replace(' ', ''))), 'w')
        out_dec.write('# %s resolution:%d\n' % (
            name, resolution))
        out_dec.write('# BADS %s\n' % (','.join([str(b) for b in bads])))

    def write2matrix(a, b, c):
        out_raw.write('%d\t%d\t%d\n' % (a, b, c))
    def write2matrices(a, b, c):
        out_raw.write('%d\t%d\t%d\n' % (a, b, c))
        out_nrm.write('%d\t%d\t%f\n' % (a, b, c / (bias[a] * bias[b])))
        out_dec.write('%d\t%d\t%f\n' % (a, b, c / (bias[a] * bias[b] * decay[abs(a-b)])))
    def write2matrices_err(a, b, c):
        out_raw.write('%d\t%d\t%d\n' % (a, b, c))
        out_nrm.write('%d\t%d\t%f\n' % (a, b, c / (bias[a] * bias[b])))
        try:
            out_dec.write('%d\t%d\t%f\n' % (a, b, c / (bias[a] * bias[b] * decay[abs(a-b)])))
        except KeyError:  # different chromsomes
            out_dec.write('%d\t%d\t%s\n' % (a, b, 'nan'))

    if biases:
        if len(regions) == 1:
            write = write2matrices
        else:
            write = write2matrices_err
    else:
        write = write2matrix

    sys.stdout.write('     ')
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        if not i % 10 and i:
            sys.stdout.write(' ')
        if not i % 50 and i:
            sys.stdout.write(' %9s\n     ' % ('%s/%s' % (i , len(regs))))
        sys.stdout.write('.')
        sys.stdout.flush()

        fname = os.path.join(outdir, 'tmp_%s:%d-%d.pickle' % (region, start, end))
        dico = load(open(fname))
        for (j, k), v in dico.iteritems():
            if j in bads or k in bads:
                continue
            write(j, k, v)
        os.system('rm -f %s' % (fname))
    out_raw.close()
    if biases:
        out_nrm.close()
        out_dec.close()
    print '%s %9s\n' % (' ' * (54 - (i % 50) - (i % 50) / 10),
                        '%s/%s' % (len(regs),len(regs)))


def main():
    opts          = get_options()
    inbam          = opts.inbam
    resolution     = opts.reso
    filter_exclude = opts.filter
    ncpus          = opts.cpus
    if opts.biases:
        biases     = load(open(opts.biases))
    else:
        biases     = {}
    outdir         = opts.outdir
    coord          = opts.coord

    if biases['resolution'] != resolution:
        raise Exception('ERROR: different resolution in bias file (you want %d,'
                        ' there is %d).\n' % (resolution, biases['resolution']))
    
    if not coord:
        region = None
        start = None
        end   = None
    else:
        try:
            crm, pos   = coord.split(':')
            region     = crm
            start, end = pos.split('-')
            start      = int(start)
            end        = int(end)
        except ValueError:
            region = coord
            start  = None
            end    = None

    
    mkdir(outdir)
    if region:
        sys.stdout.write('\nExtraction of %s' % (region))
        if start:
            sys.stdout.write(':%s-%s\n' % (start, end))
        else:
            sys.stdout.write(' (full chromosome)\n')
    else:
        sys.stdout.write('\nExtraction of full genome\n')
    
    read_bam(inbam, filter_exclude, resolution, biases,
             region=region, start=start, end=end,
             ncpus=ncpus, outdir=outdir)
    
    printime('\nDone.')


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbam', metavar='',
                        required=True, default=False, help='input HiC-BAM file.')
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='',
                        required=True, default=False, help='output directory.')
    parser.add_argument('-r', '--resolution', dest='reso', type=int, metavar='',
                        required=True, help='''wanted resolution form th generated matrix''')
    parser.add_argument('-c', '--coord', dest='coord',  metavar='', 
                        default=None, help='''Coordinate of the region to 
                        retrieve. By default all genome, arguments can be 
                        either one chromosome name, or the_shape coordinate in 
                        the form "chr3:110000000-120000000"''')
    parser.add_argument('-b', '--biases', dest='biases', metavar='',
                        help='''path to pickle file with array of biases''')
    parser.add_argument('-C', '--cpus', dest='cpus', metavar='', type=int,
                        default=8, help='''[%(default)s] number of cpus to be 
                        used for parsing the HiC-BAM file''')
    parser.add_argument('-F', '--filter', dest='filter', metavar='', type=int,
                        default=391, help='''[%(default)s] binary code to 
                        exclude filtered reads [OPTION TO BE IMPROVED]''')
    
    opts = parser.parse_args()
    
    return opts


if __name__=='__main__':
    exit(main())
