"""
parser provided by Enrique Vidal <enrique.vidal@crg.eu> to read write 2D beds
into compressed BAM format.

"""

from pytadbit.utils.file_handling import mkdir
from cPickle                      import load, dump
from time                         import sleep, time
from collections                  import OrderedDict
from pytadbit.utils.extraviews    import nicer
from warnings                     import warn
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


def read_bam_frag(inbam, filter_exclude, sections1, sections2,
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
                pos1 = sections1[(crm1, pos1 / resolution)]
                pos2 = sections2[(crm2, pos2 / resolution)]
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
             region1=None, start1=None, end1=None, verbose=False,
             region2=None, start2=None, end2=None, outdir=None,
             tmpdir='/tmp/', normalized=False, by_decay=False,
             get_all_data=False, use_bads=False):
    """
    Extracts a (normalized) submatrix at wanted resolution from pseudo-BAM file

    :param inbam: path to pseudoBAM file
    :param filter_exclude:
    :param resolution:
    :param biases: path to pickle file with biases and low-coverage columns
    :param 8 ncpus:
    :param None region1: chromosome name of region 1
    :param None start1: start genomic coordinate of region 1
    :param None end1: end genomic coordinate of region 1
    :param None region1: chromosome name of region 2 (if not given use region1)
    :param None start1: start genomic coordinate of region 2 (if not given use region1)
    :param None end1: end genomic coordinate of region 2 (if not given use region1)
    :param False normalized: returns the dictionary of Vanilla normalized matrix
    :param False decay: returns the dictionary of Decay normalized matrix (decay
       option can not be used at the same time as normalized option)
    :param False get_all_data:

    returns: dictionary of interactions. If get_all_data is set to True, returns
       a dictionary with all biases used and bads1 columns (keys of the 
       dicitionary are: matrix, bias1, bias2, bads1, bads1, decay).
    """
    if outdir:
        mkdir(outdir)
    mkdir(tmpdir)
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
    if region1:
        regions = [region1]
        start_bin = [i for i, b in enumerate(bins) if b[0] == region1][0]
        end_bin   = [i for i, b in enumerate(bins[start_bin:], start_bin)
                     if b[0] == region1][-1]
    else:
        regions = bamfile.references
        total = len(bins)
        if start1 or end1:
            raise Exception('ERROR: Cannot use start/end1 without region')

    if start1:
        start_bin = section_pos[region1][0] + start1 / resolution
    else:
        start1 = 0
    if end1:
        end_bin = section_pos[region1][0] + end1 / resolution
    else:
        end = len(bins)
        end1 = (section_pos[region1][1] - section_pos[region1][0]) * resolution

    total = end_bin - start_bin + 1
    regs  = []
    begs  = []
    ends  = []
    njobs = min(total, 100) + 1
    nbins = total / njobs + 1
    for i in xrange(start_bin, end_bin, nbins):
        if i + nbins > end_bin:  # make sure that we stop at the right place
            nbins = end_bin - i
        try:
            (crm1, beg1), (crm2, fin2) = bins[i], bins[i + nbins - 1]
        except IndexError:
            (crm1, beg1), (crm2, fin2) = bins[i], bins[-1]
        if crm1 != crm2:
            fin1 = sections[crm1]
            beg2 = 0
            regs.append(crm1)
            regs.append(crm2)
            begs.append(beg1 * resolution)
            begs.append(beg2 * resolution)
            ends.append(fin1 * resolution + resolution)  # last nt included
            ends.append(fin2 * resolution + resolution - 1)  # last nt not included (overlap with next window)
        else:
            regs.append(crm1)
            begs.append(beg1 * resolution)
            ends.append(fin2 * resolution + resolution - 1)            
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
    bins_dict1 = dict([(j, i) for i, j in enumerate(bins)])
    if region2:
        bins = []
        beg_crm = section_pos[region2][0]
        if start2 is not None:
            start_bin2 = section_pos[region2][0] + start2 / resolution
            end_bin2   = section_pos[region2][0] + end2   / resolution
        else:
            start2     = 0
            start_bin2 = 0
            end_bin2   = section_pos[region2][1]
            end2       = sections[region2] * resolution
        start = start_bin2 - beg_crm
        end   = end_bin2   - beg_crm
        bins = [(region2, i) for i in xrange(start, end)]
        bins_dict2 = dict([(j, i) for i, j in enumerate(bins)])
    else:
        bins_dict2 = bins_dict1
    pool = mu.Pool(ncpus)
    ## RUN!
    if verbose:
        printime('\n  - Parsing BAM (%d chunks)' % (len(regs)))
    procs = []
    for i, (region, b, e) in enumerate(zip(regs, begs, ends)):
        if ncpus == 1:
            read_bam_frag(inbam, filter_exclude,
                          bins_dict1, bins_dict2,
                          resolution, tmpdir, region, b, e,)
        else:
            procs.append(pool.apply_async(
                read_bam_frag, args=(inbam, filter_exclude,
                                     bins_dict1, bins_dict2,
                                     resolution, tmpdir, region, b, e,)))
    pool.close()
    if verbose:
        print_progress(procs)
    pool.join()

    if verbose:
        printime('  - Writing matrices')
    bias1  = dict((k - start_bin, v)
                  for k, v in biases.get('biases', {}).iteritems()
                  if start_bin <= k < end_bin)
    if region2:
        bias2  = dict((k - start_bin2, v)
                      for k, v in biases.get('biases', {}).iteritems()
                      if start_bin2 <= k < end_bin2)
    else:
        bias2 = bias1
    decay = biases.get('decay' , {})
    bads1  = dict((k - start_bin, v)
                  for k, v in biases.get('badcol', {}).iteritems()
                  if start_bin <= k < end_bin)
    if region2:
        bads2  = dict((k - start_bin2, v)
                      for k, v in biases.get('badcol', {}).iteritems()
                      if start_bin2 <= k < end_bin2)
    else:
        bads2 = bads1
    if use_bads:
        bads2 = bads1 = {}
    # hic_data = HiC_data((), len(bins_dict), sections,
    #                     bins_dict, resolution=resolution)
    if len(regions) == 1:
        if region2:
            name = '%s:%d-%d_%s:%d-%d' % (region1, start1 / resolution, end1 / resolution,
                                          region2, start2 / resolution, end2 / resolution)
        else:
            name = '%s:%d-%d' % (region1, start1 / resolution, end1 / resolution)
    else:
        name = 'full'
    if outdir:
        out_raw = open(os.path.join(outdir, 'matrix_raw_%s_%s.abc' % (
            name, nicer(resolution).replace(' ', ''))), 'w')
        out_raw.write('# %s resolution:%d\n' % (name, resolution))
        if region2:
            out_raw.write('# BADROWS %s\n' % (','.join([str(b) for b in bads1])))
            out_raw.write('# BADCOLS %s\n' % (','.join([str(b) for b in bads2])))
        else:
            out_raw.write('# BADS %s\n' % (','.join([str(b) for b in bads1])))
        if biases:
            out_nrm = open(os.path.join(outdir, 'matrix_nrm_%s_%s.abc' % (
                name, nicer(resolution).replace(' ', ''))), 'w')
            out_nrm.write('# %s resolution:%d\n' % (name, resolution))
            if region2:
                out_nrm.write('# BADROWS %s\n' % (','.join([str(b) for b in bads1])))
                out_nrm.write('# BADCOLS %s\n' % (','.join([str(b) for b in bads2])))
            else:
                out_nrm.write('# BADS %s\n' % (','.join([str(b) for b in bads1])))
            out_dec = open(os.path.join(outdir, 'matrix_dec_%s_%s.abc' % (
                name, nicer(resolution).replace(' ', ''))), 'w')
            out_dec.write('# %s resolution:%d\n' % (
                name, resolution))
            if region2:
                out_dec.write('# BADROWS %s\n' % (','.join([str(b) for b in bads1])))
                out_dec.write('# BADCOLS %s\n' % (','.join([str(b) for b in bads2])))
            else:
                out_dec.write('# BADS %s\n' % (','.join([str(b) for b in bads1])))

        def write2matrix(a, b, c):
            out_raw.write('%d\t%d\t%d\n' % (a, b, c))
        def write2matrices(a, b, c):
            out_raw.write('%d\t%d\t%d\n' % (a, b, c))
            out_nrm.write('%d\t%d\t%f\n' % (a, b, c / (bias1[a] * bias2[b])))
            out_dec.write('%d\t%d\t%f\n' % (a, b, c / (bias1[a] * bias2[b] *
                                                       decay[abs(a-b)])))
        def write2matrices_2reg(a, b, c):
            out_raw.write('%d\t%d\t%d\n' % (a, b, c))
            out_nrm.write('%d\t%d\t%f\n' % (a, b, c / (bias1[a] * bias2[b])))
            out_dec.write('%d\t%d\t%f\n' % (a, b, c / (bias1[a] * bias2[b] *
                                                       decay[abs((a + start_bin) -
                                                                 (b + start_bin2))])))
        def write2matrices_err(a, b, c):
            out_raw.write('%d\t%d\t%d\n' % (a, b, c))
            out_nrm.write('%d\t%d\t%f\n' % (a, b, c / (bias1[a] * bias2[b])))
            try:
                out_dec.write('%d\t%d\t%f\n' % (a, b, c / (bias1[a] * bias2[b] *
                                                           decay[abs(a-b)])))
            except KeyError:  # different chromsomes
                out_dec.write('%d\t%d\t%s\n' % (a, b, 'nan'))

        if biases:
            if len(regions) == 1:
                if region2:
                    write = write2matrices_2reg
                else:
                    write = write2matrices
            else:
                write = write2matrices_err
        else:
            write = write2matrix
    
    if verbose:
        sys.stdout.write('     ')
    dico = {}
    for i, (region, start, end) in enumerate(zip(regs, begs, ends)):
        if not i % 10 and i:
            if verbose:
                sys.stdout.write(' ')
        if not i % 50 and i:
            if verbose:
                sys.stdout.write(' %9s\n     ' % ('%s/%s' % (i , len(regs))))
        if verbose:
            sys.stdout.write('.')
            sys.stdout.flush()
        fname = os.path.join(tmpdir, 'tmp_%s:%d-%d.pickle' % (region, start, end))
        if outdir:
            dico = load(open(fname))
            for (j, k), v in dico.iteritems():
                if j in bads1 or k in bads2:
                    continue
                write(j, k, v)
        else:
            dico.update(load(open(fname)))
        os.system('rm -f %s' % (fname))
    if outdir:
        out_raw.close()
        if biases:
            out_nrm.close()
            out_dec.close()
    if verbose:
        print '%s %9s\n' % (' ' * (54 - (i % 50) - (i % 50) / 10),
                            '%s/%s' % (len(regs),len(regs)))
    if normalized and by_decay:
        warn('WARNING: choose either normalized or by_decay. Using decay normalization')
    if not outdir:
        if by_decay:
            if region2:
                for i, j in dico:
                    if i in bads1 or j in bads2:
                        continue
                    try:
                        dico[(i, j)] /= bias1[i] * bias2[j] * decay[abs((i + start_bin) -
                                                                        (j + start_bin2))]
                    except KeyError:
                        dico[(i, j)] = float('nan')  # no value in decay
            else:
                for i, j in dico:
                    dico[(i, j)] /= bias1[i] * bias2[j] * decay[abs(i - j)]
        elif normalized:
            for i, j in dico:
                dico[(i, j)] /= bias1[i] * bias2[j]
        if get_all_data:
            return {'matrix': dico,
                    'bias1' : bias1,
                    'bias2' : bias2,
                    'bads1' : bads1,
                    'bads2' : bads2,
                    'decay' : decay}
        return dico
        

def bam_to_hic_data(inbam, resolution_list, filter_exclude, filter_include):
    """
    Load hacked BAM file, into list of hic_data objects (one perc_zero resolution)

    :param inbam: path to a BAM file
    :param resolution_list: list of reolutions
    :param filter_exclude: filters to exclude expects a number, which in binary would
       correspond to the presence/absence of the filters in the corresponding 
       order:
         - self_circle
         - dangling end
         - error
         - extra dangling_end
         - too close from RES
         - too short
         - too large
         - over represented
         - duplicated
         - random breaks
         - trans
    :param filter_include: filters to be included (see doc of filter_exclude param)
 
    :returns: a list of HiC_data objects
    """
    # open bam file
    bamfile = pysam.AlignmentFile(inbam, "rb")
    # get sections
    sections = OrderedDict(zip(bamfile.references,
                               bamfile.lengths))
    # init HiC_data objects
    dat_list = [init_hic_data(bamfile, resolution) for resolution in resolution_list]
    bin_list = [get_bins(bamfile, resolution) for resolution in resolution_list]
    # close bam file    
    bamfile.close()
    # access bam file per chromosome
    for chrom in sections.keys():
        j = 0
        print "Processing chromosome " + str(chrom)
        # read line by line
        for line in pysam.view("-F", str(filter_exclude),
                               "-f", str(filter_include),
                               inbam,
                               chrom):
            # get info
            info = line.strip().split("\t")
            chrom1 = info[2]
            pos1 = int(info[3])
            chrom2 = info[6]
            pos2 = int(info[7])
            if chrom2 == "=":
                chrom2 = chrom1
            # get bins and add counts
            try:
                for i in xrange(len(resolution_list)):
                    b1 = bin_list[i][(chrom1, int(pos1) / resolution_list[i])]
                    b2 = bin_list[i][(chrom2, int(pos2) / resolution_list[i])]
                    dat_list[i][b1, b2] += 1
            except KeyError:
                pass
            j += 1
        print str(j) + " lines processed"
    return dat_list, bin_list


def bam_to_2Dbed(inbam, outbed, resolution, filter_exclude=None, filter_include=None):
    """
    Load hacked BAM file, into list of hic_data objects (one perc_zero resolution)

    :param inbam: path to a BAM file
    :param outbed: path to output 2Dbed file
    :param resolution_list: list of reolutions
    :param None filter_exclude: filters to exclude expects a number, which in binary would
       correspond to the presence/absence of the filters in the corresponding 
       order:
         - self_circle
         - dangling end
         - error
         - extra dangling_end
         - too close from RES
         - too short
         - too large
         - over represented
         - duplicated
         - random breaks
         - trans
    :param None filter_include: filters to be included (see doc of filter_exclude param)
 
    :returns: a list of HiC_data objects
    """
    # open bam file
    bamfile = pysam.AlignmentFile(inbam, "rb")
    # get sections
    sections = OrderedDict(zip(bamfile.references,
                               bamfile.lengths))
    outbed = open(outbed, 'w')
    for sec in sections:
        outbed.write('# CRM %s %d\n' % (sec, sections[sec]))
    # close bam file    
    bamfile.close()
    # access bam file per chromosome
    j = 0
    print "Processing chromosome " + str(chrom)
    # read line by line
    for line in pysam.view("-F", str(filter_exclude),
                           "-f", str(filter_include),
                           inbam):
        # get info
        info = line.strip().split("\t")
        chrom1 = info[2]
        pos1 = int(info[3])
        chrom2 = info[6]
        pos2 = int(info[7])
        if chrom2 == "=":
            chrom2 = chrom1
        # get bins and add counts
        try:
            for i in xrange(len(resolution_list)):
                b1 = bin_list[i][(chrom1, int(pos1) / resolution_list[i])]
                b2 = bin_list[i][(chrom2, int(pos2) / resolution_list[i])]
                dat_list[i][b1, b2] += 1
        except KeyError:
            pass
        j += 1
    print str(j) + " lines processed"
    return dat_list, bin_list




