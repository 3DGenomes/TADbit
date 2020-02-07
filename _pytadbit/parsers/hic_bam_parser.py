"""
parser provided by Enrique Vidal <enrique.vidal@crg.eu> to read write 2D beds
into compressed BAM format.

"""
from __future__ import print_function
from builtins   import next

from future import standard_library
standard_library.install_aliases()
from pickle                      import load, dump
from time                         import sleep, time
from collections                  import OrderedDict
from subprocess                   import Popen, PIPE
from random                       import getrandbits
from tarfile                      import open as taropen
from io                     import StringIO
from shutil                       import copyfile
import datetime
from sys                          import stdout, stderr, exc_info, modules
from distutils.version            import LooseVersion
import os
import multiprocessing as mu

try:
    from lockfile                 import LockFile
except ImportError:
    pass  # silently pass, very specific need

from pysam                        import AlignmentFile

from pytadbit.utils.file_handling   import mkdir, which
from pytadbit.utils.extraviews      import nicer
from pytadbit.mapping.filter        import MASKED
try:
    from pytadbit.parsers.cooler_parser import cooler_file
except ImportError:
    pass

try:
    basestring
except NameError:
    basestring = str

def filters_to_bin(filters):
    return sum((k in filters) * 2**(k-1) for k in MASKED)


def printime(msg):
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')


def print_progress(procs):
    stdout.write('     ')
    prev_done = done = i = 0
    while done < len(procs):
        sleep(0.1)
        done = sum(p.ready() for p in procs)
        for i in range(prev_done, done):
            if not i % 10 and i:
                stdout.write(' ')
            if not i % 50 and i:
                stdout.write(' %9s\n     ' % ('%s/%s' % (i , len(procs))))
            stdout.write('.')
            stdout.flush()
        prev_done = done
    print('%s %9s\n' % (' ' * (54 - (i % 50) - (i % 50) // 10),
                        '%s/%s' % (len(procs),len(procs))))


def _map2sam_short(line, flag):
    """
    translate map + flag into hic-sam (two lines per contact)
    lose information
    58% of the size using RE sites, and 68% of the generation time
    """
    (qname,
     rname, pos, _, _, _, _,
     rnext, pnext, _) = line.strip().split('\t', 9)

    # multicontact?
    try:
        tc = qname.split('#')[1].split('/')[1]
    except IndexError:
        tc = '1'
    # trans contact?
    if rname != rnext:
        flag += 1024 # filter_keys['trans-chromosomic'] = 2**10
    r1r2 = ('{0}\t{1}\t{2}\t{3}\t0\t1P\t{4}\t{5}\t*\t*\t*\t'
            'TC:i:{6}\n'
            '{0}\t{1}\t{4}\t{5}\t0\t1S\t{2}\t{3}\t*\t*\t*\t'
            'TC:i:{6}\n'
           ).format(
               qname,               # 0
               flag,                # 1
               rname,               # 2
               pos,                 # 3
               rnext,               # 4
               pnext,               # 5
               tc)                  # 6
    return r1r2


def _map2sam_mid(line, flag):
    """
    translate map + flag into hic-sam (two lines per contact)
    only loses RE sites, that can be added back later
    63% of the size using RE sites, and 72% of the generation time
    """
    (qname,
     rname, pos, s1, l1, _, _,
     rnext, pnext, s2, l2, _) = line.strip().split('\t', 11)

    # multicontact?
    try:
        tc = qname.split('#')[1].split('/')[1]
    except IndexError:
        tc = '1'
    # trans contact?
    if rname != rnext:
        flag += 1024 # filter_keys['trans-chromosomic'] = 2**10
    r1r2 = ('{0}\t{1}\t{2}\t{3}\t0\t{4}P\t{6}\t{7}\t{5}\t*\t*\t'
            'TC:i:{8}\t'
            'S1:i:{9}\tS2:i:{10}\n'
            '{0}\t{1}\t{6}\t{7}\t0\t{5}S\t{2}\t{3}\t{4}\t*\t*\t'
            'TC:i:{8}\t'
            'S1:i:{9}\tS2:i:{10}\n').format(
                qname,               # 0
                flag,                # 1
                rname,               # 2
                pos,                 # 3
                l1,                  # 4
                l2,                  # 5
                rnext,               # 6
                pnext,               # 7
                tc,                  # 8
                s1,                  # 9
                s2)                  # 10
    return r1r2


def _map2sam_long(line, flag):
    """
    translate map + flag into hic-sam (two lines per contact)
    """
    (qname,
     rname, pos, s1, l1, e1, e2,
     rnext, pnext, s2, l2, e3, e4) = line.strip().split('\t')

    # multicontact?
    try:
        tc = qname.split('#')[1].split('/')[1]
    except IndexError:
        tc = '1'
    # trans contact?
    if rname != rnext:
        flag += 1024 # filter_keys['trans-chromosomic'] = 2**10
    r1r2 = ('{0}\t{1}\t{2}\t{3}\t0\t{4}P\t{6}\t{7}\t{5}\t*\t*\t'
            'TC:i:{8}\tS1:i:{13}\tS2:i:{14}\t'
            'E1:i:{9}\tE2:i:{10}\tE3:i:{11}\tE4:i:{12}\n'
            '{0}\t{1}\t{6}\t{7}\t0\t{5}S\t{2}\t{3}\t{4}\t*\t*\t'
            'TC:i:{8}\tS1:i:{14}\tS2:i:{13}\t'
            'E3:i:{11}\tE4:i:{12}\tE1:i:{9}\tE2:i:{10}\n').format(
                qname,               # 0
                flag,                # 1
                rname,               # 2
                pos,                 # 3
                l1,                  # 4
                l2,                  # 5
                rnext,               # 6
                pnext,               # 7
                tc,                  # 8
                e1,                  # 9
                e2,                  # 10
                e3,                  # 11
                e4,                  # 12
                s1,                  # 13
                s2)                  # 14
    return r1r2


def bed2D_to_BAMhic(infile, valid, ncpus, outbam, frmt, masked=None, samtools='samtools'):
    """
    function adapted from Enrique Vidal <enrique.vidal@crg.eu> scipt to convert
    2D beds into compressed BAM format.

    Gets the *_both_filled_map.tsv contacts from TADbit (and the corresponding
    filter files) and outputs a modified indexed BAM with the following fields:

       - read ID
       - filtering flag (see codes in header)
       - chromosome ID of the first pair of the contact
       - genomic position of the first pair of the contact
       - MAPQ set to 0
       - pseudo CIGAR with sequence length and info about current copy (P: first copy, S: second copy)
       - chromosome ID of the second pair of the contact
       - genomic position of the second pair of the contact
       - mapped length of the second pair of the contact
       - sequence is missing (*)
       - quality is missing (*)
       - TC tag indicating single (1) or multi contact (3 6 ... number being the number of times a given sequenced fragment is involved in a pairwise contact)
       - S1 and S2 tags are the strand orientation of the left and right read-end

    Each pair of contacts produces two lines in the output BAM
    """
    samtools = which(samtools)
    if not samtools:
        raise Exception('ERROR: samtools is needed to save a compressed '
                        'version of the results. Check '
                        'http://samtools.sourceforge.net/ \n')

    # define filter codes
    filter_keys = OrderedDict()
    for k in MASKED:
        filter_keys[MASKED[k]['name'].replace(' ', '-')] = 2 ** (k - 1)

    output = ''

    # write header
    output += ("\t".join(("@HD" ,"VN:1.5", "SO:queryname")) + '\n')
    fhandler = open(infile)
    line = next(fhandler)
    # chromosome lengths
    pos_fh = 0

    while line.startswith('#'):
        (_, _, cr, ln) = line.replace("\t", " ").strip().split(" ")
        output += ("\t".join(("@SQ", "SN:" + cr, "LN:" + ln)) + '\n')
        pos_fh += len(line)
        line = next(fhandler)

    # filter codes
    for i in filter_keys:
        output += ("\t".join(("@CO", "filter:" + i, "flag:" + str(filter_keys[i]))) + '\n')

    # tags
    output += ("\t".join(("@CO" ,"TC:i", "Number of time a sequenced fragment is involved in a pairwise contact\n")))
    output += ("\t".join(("@CO" ,("Each read is duplicated: once starting with the "
                                  "left read-end, once with the right read-end\n"))))
    output += ("\t".join(("@CO" , (" the order of RE sites and strands changes consequently "
                                   "depending on which read-end comes first ("
                                   "when right end is first: E3 E4 E1 E2)\n"))))
    output += ("\t".join(("@CO" ,(" CIGAR code contains the length of the "
                                  "1st read-end mapped and 'P' or 'S' "
                                  "if the copy is the first or the second\n"))))
    output += ("\t".join(("@CO" ,"E1:i", "Position of the left RE site of 1st read-end\n")))
    output += ("\t".join(("@CO" ,"E2:i", "Position of the right RE site of 1st read-end\n")))
    output += ("\t".join(("@CO" ,"E3:i", "Position of the left RE site of 2nd read-end\n")))
    output += ("\t".join(("@CO" ,"E4:i", "Position of the right RE site of 2nd read-end\n")))
    output += ("\t".join(("@CO" ,"S1:i", "Strand of the 1st read-end (1: positive, 0: negative)\n")))
    output += ("\t".join(("@CO" ,"S2:i", "Strand of the 2nd read-end  (1: positive, 0: negative)\n")))

    # open and init filter files
    if not valid:
        filter_line, filter_handler = get_filters(infile, masked)
    fhandler.seek(pos_fh)
    # check samtools version number and modify command line
    version = LooseVersion([l.split()[1]
                            for l in Popen(samtools, stderr=PIPE,
                                           universal_newlines=True).communicate()[1].split('\n')
                            if 'Version' in l][0])
    pre = '-o' if version >= LooseVersion('1.3') else ''

    proc = Popen(samtools + ' view -Shb -@ %d - | samtools sort -@ %d - %s %s' % (
        ncpus, ncpus, pre,
        outbam + '.bam' if  version >= LooseVersion('1.3') else ''),  # in new version '.bam' is no longer added
                 shell=True, stdin=PIPE, universal_newlines=True)
    proc.stdin.write(output)
    if frmt == 'mid':
        map2sam = _map2sam_mid
    elif frmt == 'long':
        map2sam = _map2sam_long
    else:
        map2sam = _map2sam_short

    if valid:
        for line in fhandler:
            flag = 0
            # get output in sam format
            proc.stdin.write(map2sam(line, flag))
    else:
        for line in fhandler:
            flag = 0
            # check if read matches any filter
            rid = line.split("\t")[0]
            for i in filter_line:
                if filter_line[i] == rid:
                    flag += filter_keys[i]
                    try:
                        filter_line[i] = next(filter_handler[i]).strip()
                    except StopIteration:
                        pass
            # get output in sam format
            proc.stdin.write(map2sam(line, flag))
    proc.stdin.close()
    proc.wait()

    # Index BAM
    _ = Popen(samtools + ' index %s.bam' % (outbam), shell=True, universal_newlines=True).communicate()

    # close file handlers
    fhandler.close()
    if not valid:
        for i in filter_handler:
            filter_handler[i].close()


def get_filters(infile, masked):
    """
    get all filters
    """
    basename = os.path.basename(infile)
    dirname = os.path.dirname(infile)
    if masked:  # if providedusemasked dictionary generate by filtering function
        filter_files = {}
        for i in masked:
            try:
                filter_files[masked[i]['name'].replace(' ', '_')] = masked[i]['fnam']
            except KeyError:
                if i != 11:  # error is expected for filter 11,which are trans reads.. a bit uggly... I know
                    raise Exception()
                continue
    else:  # otherwise search in the directory for the files
        filter_files = {}
        stderr.write('Using filter files:\n')
        for fname in os.listdir(dirname):
            if fname.startswith(basename + "_"):
                key = fname.replace(basename + "_", "").replace(".tsv", "")
                filter_files[key] = dirname + "/" + fname
            stderr.write('   - %-20s %s\n' %(key, fname))
    filter_handler = {}
    filter_line = {}
    for i in filter_files:
        filter_handler[i.replace('_', '-')] = open(filter_files[i])
        try:
            filter_line[i.replace('_', '-')] = next(filter_handler[i.replace('_', '-')]).strip()
        except StopIteration:
            filter_line[i.replace('_', '-')] = ""
    return filter_line, filter_handler


def _read_bam_frag(inbam, filter_exclude, all_bins, sections1, sections2,
                   rand_hash, resolution, tmpdir, region, start, end,
                   half=False, sum_columns=False):
    bamfile = AlignmentFile(inbam, 'rb')
    refs = bamfile.references
    bam_start = start - 2
    bam_start = max(0, bam_start)
    try:
        dico = {}
        for r in bamfile.fetch(region=region,
                               start=bam_start, end=end,  # coords starts at 0
                               multiple_iterators=True):
            if r.flag & filter_exclude:
                continue
            crm1 = r.reference_name
            pos1 = r.reference_start + 1
            crm2 = refs[r.mrnm]
            pos2 = r.mpos + 1
            try:
                pos1 = sections1[(crm1, pos1 // resolution)]
                pos2 = sections2[(crm2, pos2 // resolution)]
            except KeyError:
                continue  # not in the subset matrix we want
            crm = crm1 * (crm1 == crm2)
            try:
                dico[(crm, pos1, pos2)] += 1
            except KeyError:
                dico[(crm, pos1, pos2)] = 1
            # print '%-50s %5s %9s %5s %9s' % (r.query_name,
            #                                  crm1, r.reference_start + 1,
            #                                  crm2, r.mpos + 1)
        if half:
            for c, i, j in dico:
                if i < j:
                    del dico[(c, i, j)]
        out = open(os.path.join(tmpdir, '_tmp_%s' % (rand_hash),
                                '%s:%d-%d.tsv' % (region, start, end)), 'w')
        out.write(''.join('%s\t%d\t%d\t%d\n' % (c, a, b, v)
                          for (c, a, b), v in dico.items()))
        out.close()
        if sum_columns:
            sumcol = {}
            cisprc = {}
            for (c, i, j), v in dico.items():
                # out.write('%d\t%d\t%d\n' % (i, j, v))
                try:
                    sumcol[i] += v
                    cisprc[i][all_bins[i][0] == all_bins[j][0]] += v
                except KeyError:
                    sumcol[i]  = v
                    cisprc[i]  = [0, 0]
                    cisprc[i][all_bins[i][0] == all_bins[j][0]] += v
            return sumcol, cisprc
    except Exception as e:
        exc_type, exc_obj, exc_tb = exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(e)
        print(exc_type, fname, exc_tb.tb_lineno)


def read_bam(inbam, filter_exclude, resolution, ncpus=8,
             region1=None, start1=None, end1=None,
             region2=None, start2=None, end2=None, nchunks=100,
             tmpdir='.', verbose=True, normalize=False, max_size=None,
             chr_order=None):

    bamfile = AlignmentFile(inbam, 'rb')
    bam_refs = bamfile.references
    bam_lengths = bamfile.lengths
    if chr_order:
        bam_refs_idx = [bam_refs.index(chr_ord)
                        for chr_ord in chr_order if chr_ord in bam_refs]
        if not bam_refs_idx :
            raise Exception('''ERROR: Wrong number of chromosomes in chr_order.
                Found %s in bam file \n''' % (' '.join(bam_refs)))
        bam_refs = [bam_ref for bam_ref in [bam_refs[bam_ref_idx]
                                              for bam_ref_idx in bam_refs_idx]]
        bam_lengths = [bam_len for bam_len in [bam_lengths[bam_ref_idx]
                                                 for bam_ref_idx in bam_refs_idx]]
    sections = OrderedDict(list(zip(bam_refs,
                               [x // resolution + 1 for x in bam_lengths])))
    # get chromosomes and genome sizes
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    # define genomic bins
    bins = []
    for crm in sections:
        len_crm = sections[crm]
        bins.extend([(crm, i) for i in range(len_crm)])
    if not bins:
        raise Exception('ERROR: Chromosome %s smaller than bin size\n' % (crm))

    # define start, end position of region to grab
    start_bin1 = 0
    end_bin1   = len(bins) + 1
    regions = bam_refs
    if region1:
        regions = [region1]
        if region2:
            regions.append(region2)
    else:
        total = len(bins)
        if start1 is not None or end1:
            raise Exception('ERROR: Cannot use start/end1 without region')

    if start1 is not None:
        start_bin1 = section_pos[region1][0] + start1 // resolution
    else:
        if region1:
            start_bin1 = section_pos[region1][0]
        else:
            start_bin1 = 0
        start1 = 0
    if end1 is not None:
        end_bin1 = section_pos[region1][0] + end1 // resolution
    else:
        if region1:
            end_bin1 = section_pos[region1][1]
            end1 = sections[region1] * resolution
        else:
            end_bin1 = total
            end1 = total * resolution

    # define chunks, using at most 100 sub-divisions of region1
    total = end_bin1 - start_bin1
    regs  = []
    begs  = []
    ends  = []

    njobs = min(total, nchunks) + 1

    nbins = total // njobs + 1
    for i in range(start_bin1, end_bin1, nbins):
        if i + nbins > end_bin1:  # make sure that we stop at the right place
            nbins = end_bin1 - i
        try:
            (crm1, beg1), (crm2, fin2) = bins[i], bins[i + nbins - 1]
        except IndexError:
            (crm1, beg1), (crm2, fin2) = bins[i], bins[-1]
        if crm1 != crm2:
            fin1 = sections[crm1]
            regs.append(crm1)
            begs.append(beg1 * resolution)
            ends.append(fin1 * resolution + resolution)  # last nt included
            # be sure we don't miss regions in between the start and end bins
            start_chunk = i+fin1-beg1
            end_chunk = i + nbins - 1 if i + nbins - 1 < len(bins) else len(bins) - 1
            (crm1, beg1) = bins[start_chunk]
            fin1 = beg1
            for j in range(start_chunk, end_chunk+1):
                (crm2, beg2) = bins[j]
                if crm1 == crm2:
                    fin1 = beg2
                    continue
                regs.append(crm1)
                begs.append(beg1 * resolution)
                fin1 = sections[crm1]
                ends.append(fin1 * resolution + resolution - 1)  # last nt not included (overlap with next window)
                (crm1, beg1) = (crm2, beg2)
            regs.append(crm1)
            begs.append(beg1 * resolution)
            ends.append(fin1 * resolution + resolution - 1)  # last nt not included (overlap with next window)
        else:
            regs.append(crm1)
            begs.append(beg1 * resolution)
            ends.append(fin2 * resolution + resolution - 1)
    ends[-1] += 1  # last nucleotide included
    # reduce dictionaries
    all_bins = []
    seenbins = set()
    for crm in regions:
        beg_crm = section_pos[crm][0]
        if region1:
            start = start_bin1 - beg_crm
            end   = end_bin1   - beg_crm
        else:
            start = 0
            end   = section_pos[crm][1] - section_pos[crm][0]
        all_bins.extend([(crm, i) for i in range(start, end)
                          if not (crm, i) in seenbins])
        seenbins = set(all_bins)
    del(seenbins)

    bins_dict1 = dict((j, i) for i, j in enumerate(all_bins))
    if region2:
        if not region2 in section_pos:
            raise Exception('ERROR: chromosome %s not found' % region2)
        bins = []
        beg_crm = section_pos[region2][0]
        if start2 is not None:
            start_bin2 = section_pos[region2][0] + start2 // resolution
        else:
            start_bin2 = section_pos[region2][0]
            start2 = 0
        if end2 is not None:
            end_bin2   = section_pos[region2][0] + end2   // resolution
        else:
            end_bin2   = section_pos[region2][1]
            end2       = sections[region2] * resolution
        start = start_bin2 - beg_crm
        end   = end_bin2   - beg_crm
        bins = [(region2, i) for i in range(start, end)]
        bins_dict2 = dict([(j, i) for i, j in enumerate(bins)])
    else:
        start_bin2 = start_bin1
        end_bin2 = end_bin1
        bins_dict2 = bins_dict1

    size1 = end_bin1 - start_bin1
    size2 = end_bin2 - start_bin2
    if verbose:
        printime('\n  (Matrix size %dx%d)' % (size1, size2))
    if max_size and max_size < size1 * size2:
        raise Exception(('ERROR: matrix too large ({0}x{1}) should be at most '
                         '{2}x{2}').format(size1, size2, int(max_size**0.5)))

    pool = mu.Pool(ncpus)
    # create random hash associated to the run:
    rand_hash = "%016x" % getrandbits(64)

    ## RUN!
    if verbose:
        printime('\n  - Parsing BAM (%d chunks)' % (len(regs)))
    mkdir(os.path.join(tmpdir, '_tmp_%s' % (rand_hash)))
    # empty all_bins array if we are not going to normalize
    if not normalize:
        all_bins = []
    procs = []
    for i, (region, b, e) in enumerate(zip(regs, begs, ends)):
        if ncpus == 1:
            _read_bam_frag(inbam, filter_exclude, all_bins,
                           bins_dict1, bins_dict2, rand_hash,
                           resolution, tmpdir, region, b, e,)
        else:
            procs.append(pool.apply_async(
                _read_bam_frag, args=(inbam, filter_exclude, all_bins,
                                      bins_dict1, bins_dict2, rand_hash,
                                      resolution, tmpdir, region, b, e,)))
    pool.close()
    if verbose:
        print_progress(procs)
    pool.join()
    bin_coords = start_bin1, end_bin1, start_bin2, end_bin2
    chunks = regs, begs, ends
    return regions, rand_hash, bin_coords, chunks


def _iter_matrix_frags(chunks, tmpdir, rand_hash, clean=False, verbose=True,
                       include_chunk_count=False):
    if verbose:
        stdout.write('     ')
    countbin = 0
    for countbin, (region, start, end) in enumerate(zip(*chunks)):
        if verbose:
            if not countbin % 10 and countbin:
                stdout.write(' ')
            if not countbin % 50 and countbin:
                stdout.write(' %9s\n     ' % ('%s/%s' % (countbin , len(chunks[0]))))
            stdout.write('.')
            stdout.flush()

        fname = os.path.join(tmpdir, '_tmp_%s' % (rand_hash),
                             '%s:%d-%d.tsv' % (region, start, end))
        for l in open(fname):
            c, a, b, v = l.split('\t')
            if include_chunk_count:
                yield countbin, c, int(a), int(b), int(v)
            else:
                yield c, int(a), int(b), int(v)
        if clean:
            os.system('rm -f %s' % fname)
    if verbose:
        print('%s %9s\n' % (' ' * (54 - (countbin % 50) - (countbin % 50) // 10),
                            '%s/%s' % (len(chunks[0]),len(chunks[0]))))


def get_biases_region(biases, bin_coords, check_resolution=None):
    """
    Retrieve biases, decay, and bad bins from a dictionary, and re-index it
    according to a region of interest.
    """
    start_bin1, end_bin1, start_bin2, end_bin2 = bin_coords
    # load decay
    if isinstance(biases, basestring):
        biases = load(open(biases, 'rb'))
    resolution = biases.get('resolution', float('NaN'))
    if check_resolution is not None:
        if check_resolution != resolution:
            raise ValueError('ERROR: resolution not matching with wanted '
                             'resolution (wanted: {} vs found: {})'.format(
                                 check_resolution, resolution))
    decay = biases.get('decay' , {})
    # load biases and bad columns
    bias1  = dict((k - start_bin1, v)
                  for k, v in biases.get('biases', {}).items()
                  if start_bin1 <= k < end_bin1)
    bads1  = dict((k - start_bin1, v)
                  for k, v in biases.get('badcol', {}).items()
                  if start_bin1 <= k < end_bin1)
    if start_bin1 != start_bin2:
        bias2  = dict((k - start_bin2, v)
                      for k, v in biases.get('biases', {}).items()
                      if start_bin2 <= k < end_bin2)
        bads2  = dict((k - start_bin2, v)
                      for k, v in biases.get('badcol', {}).items()
                      if start_bin2 <= k < end_bin2)
    else:
        bias2 = bias1
        bads2 = bads1
    return bias1, bias2, decay, bads1, bads2


def get_matrix(inbam, resolution, biases=None,
               filter_exclude=(1, 2, 3, 4, 6, 7, 8, 9, 10),
               region1=None, start1=None, end1=None,
               region2=None, start2=None, end2=None, dico=None, clean=False,
               return_headers=False, tmpdir='.', normalization='raw', ncpus=8,
               nchunks=100, verbose=False, max_size=None, chr_order=None):
    """
    Get matrix from a BAM file containing interacting reads. The matrix
    will be extracted from the genomic BAM, the genomic coordinates of this
    matrix will be at the intersection of two regions defined byt the parameters
    region1, start1, end1 and region2, start2, end2. If the wanted matrix is
    than the second coodinate can be skipped.

    :param inbam: path to BAM file (generated byt TADbit)
    :param resolution: resolution at which we want to write the matrix
    :param biases: path to a file with biases
    :param (1, 2, 3, 4, 6, 7, 8, 9, 10) filter exclude: filters to define the
       set of valid pair of reads.
    :param 'raw' normalization: normalizations to use, can be 'decay',
       'norm' or 'raw'.
    :param None region1: chromosome name of the first region from which to
       extract the matrix
    :param None region1: chromosome name of the first region from which to
       extract the matrix
    :param None start1: start coordinate of the first region from which to
       extract the matrix
    :param None end1: end coordinate of the first region from which to
       extract the matrix
    :param None region2: chromosome name of the second region from which to
       extract the matrix
    :param None start2: start coordinate of the second region from which to
       extract the matrix
    :param None end2: end coordinate of the second region from which to
       extract the matrix
    :param '.' tmpdir: where to write temporary files
    :param 8 ncpus: number of cpus to use to read the BAM file
    :param True verbose: speak
    :param 100 nchunks: maximum number of chunks into which to cut the BAM
    :param None max_size: maximum size of matrix to read
    :param None chr_order: chromosome order

    :returns: dictionary with keys being tuples of the indexes of interacting
       bins: dico[(bin1, bin2)] = interactions
    """

    if not isinstance(filter_exclude, int):
        filter_exclude = filters_to_bin(filter_exclude)

    regions, rand_hash, bin_coords, chunks = read_bam(
        inbam, filter_exclude, resolution, ncpus=ncpus,
        region1=region1, start1=start1, end1=end1,
        region2=region2, start2=start2, end2=end2,
        tmpdir=tmpdir, nchunks=nchunks, verbose=verbose,
        max_size=max_size, chr_order=chr_order)

    if region1:
        regions = [region1]
        if region2:
            regions.append(region2)

    if biases:
        bias1, bias2, decay, bads1, bads2 = get_biases_region(biases, bin_coords)
    elif normalization != 'raw':
        raise Exception('ERROR: should provide path to file with biases (pickle).')
    else:
        bads1 = bads2 = {}
    start_bin1, start_bin2 = bin_coords[::2]

    if verbose:
        printime('  - Getting matrices')

    def transform_value_raw(_, __ , ___, v):
        return v
    def transform_value_norm(_, a, b, v):
        return v / bias1[a] / bias2[b]
    def transform_value_decay(c, a, b, v):
        return v / bias1[a] / bias2[b] / decay[c][abs(a-b)]
    def transform_value_decay_2reg(c, a, b, v):
        return v / bias1[a] / bias2[b] / decay[c][abs((a + start_bin1) - (b + start_bin2))]

    if normalization == 'raw':
        transform_value = transform_value_raw
    elif normalization == 'norm':
        transform_value = transform_value_norm
    elif normalization == 'decay':
        if start_bin1 == start_bin2:
            transform_value = transform_value_decay
        else:
            transform_value = transform_value_decay_2reg
    else:
        raise NotImplementedError(('ERROR: %s normalization not implemented '
                                   'here') % normalization)

    return_something = False
    if dico is None:
        return_something = True
        dico = dict(((i, j), transform_value(c, i, j, v))
                    for c, i, j, v in _iter_matrix_frags(
                        chunks, tmpdir, rand_hash, clean=clean, verbose=verbose)
                    if i not in bads1 and j not in bads2)
        # pull all sub-matrices and write full matrix
    else: # dico probably an HiC data object
        for _, i, j, v in _iter_matrix_frags(
                chunks, tmpdir, rand_hash,
                clean=clean, verbose=verbose):
            if i not in bads1 and j not in bads2:
                dico[i, j] = v

    if clean:
        os.system('rm -rf %s' % (os.path.join(tmpdir, '_tmp_%s' % (rand_hash))))
    if return_something:
        if return_headers:
            # define output file name
            name = _generate_name(regions, (start1, start2), (end1, end2), resolution,
                                  chr_order)
            return dico, bads1, bads2, regions, name, bin_coords
        return dico


def _generate_name(regions, starts, ends, resolution, chr_order=None):
    """
    Generate file name for write_matrix and get_matrix functions
    """
    name = []
    if len(regions) in [1, 2]:
        for i, region in enumerate(regions):
            try:
                name.append('%s:%d-%d' % (region, starts[i] // resolution,
                                        ends[i] // resolution))
            except TypeError: # all chromosomes
                name.append('%s' % (region))
        name = '_'.join(name)
    else:
        name = 'full'
    return name


def write_matrix(inbam, resolution, biases, outdir,
                 filter_exclude=(1, 2, 3, 4, 6, 7, 8, 9, 10),
                 normalizations=('decay',),
                 region1=None, start1=None, end1=None, clean=True,
                 region2=None, start2=None, end2=None, extra='',
                 half_matrix=True, nchunks=100, tmpdir='.', append_to_tar=None,
                 ncpus=8, cooler=False, row_names=False, chr_order=None, verbose=True):
    """
    Writes matrix file from a BAM file containing interacting reads. The matrix
    will be extracted from the genomic BAM, the genomic coordinates of this
    matrix will be at the intersection of two regions defined byt the parameters
    region1, start1, end1 and region2, start2, end2. If the wanted matrix is
    than the second coodinate can be skipped.

    :param inbam: path to BAM file (generated byt TADbit)
    :param resolution: resolution at which we want to write the matrix
    :param biases: path to a file with biases
    :param outdir: path to a folder where to write output files
    :param (1, 2, 3, 4, 6, 7, 8, 9, 10) filter exclude: filters to define the
       set of valid pair of reads.
    :param ('decay',) normalization: tuple with normalizations to use, can be 'decay',
       'norm' or/and 'raw'. One file per normalization will be created.
    :param None region1: chromosome name of the first region from which to
       extract the matrix
    :param None region1: chromosome name of the first region from which to
       extract the matrix
    :param None start1: start coordinate of the first region from which to
       extract the matrix
    :param None end1: end coordinate of the first region from which to
       extract the matrix
    :param None region2: chromosome name of the second region from which to
       extract the matrix
    :param None start2: start coordinate of the second region from which to
       extract the matrix
    :param None end2: end coordinate of the second region from which to
       extract the matrix
    :param True half_matrix: writes only half of the matrix (and the diagonal)
    :param '.' tmpdir: where to write temporary files
    :param None append_to_tar: path to a TAR file were generated matrices will
       be written directly
    :param 8 ncpus: number of cpus to use to read the BAM file
    :param True verbose: speak
    :param False row_names: Writes geneomic coocrdinates instead of bins.
       WARNING: results in two extra columns
    :param None chr_order: chromosome order
    :param 100 nchunks: maximum number of chunks into which to cut the BAM

    :returns: path to output files
    """
    if start1 is not None and end1:
        if end1 - start1 < resolution:
            raise Exception('ERROR: region1 should be at least as big as resolution')
    if start2 is not None and end2:
        if end2 - start2 < resolution:
            raise Exception('ERROR: region2 should be at least as big as resolution')

    if isinstance(normalizations, list):
        normalizations = tuple(normalizations)
    elif isinstance(normalizations, basestring):
        normalizations = tuple([normalizations])

    if not isinstance(filter_exclude, int):
        filter_exclude = filters_to_bin(filter_exclude)

    regions, rand_hash, bin_coords, chunks = read_bam(
        inbam, filter_exclude, resolution, ncpus=ncpus,
        region1=region1, start1=start1, end1=end1,
        region2=region2, start2=start2, end2=end2,
        tmpdir=tmpdir, nchunks=nchunks, chr_order=chr_order,
        verbose=verbose)

    if region1:
        regions = [region1]
        if region2:
            regions.append(region2)

    bamfile = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(list(zip(bamfile.references,
                               [x for x in bamfile.lengths])))

    start_bin1, end_bin1, start_bin2, end_bin2 = bin_coords

    section_pos1 = {}
    section_pos2 = {}
    totals = OrderedDict()
    total_num = 0
    for c in sections:
        totals[c] = total_num
        total_num += sections[c] / resolution + 1

    if row_names:
        if len(regions) in [1, 2]:
            offset = start_bin1 - totals[regions[0]]
            section_pos1 = dict((i, (region1, resolution * (i + offset)))
                               for i in range(end_bin1 - start_bin1))
            if region2:
                offset = start_bin2 - totals[regions[1]]
                section_pos2 = dict((i, (region2, resolution * (i + offset)))
                                    for i in range(end_bin2 - start_bin2))
        else:
            section_pos1 = dict((v + i, (c, i))
                                for c, v in totals.iteritems()
                                for i in range(sections[c] // resolution + 1))
            section_pos2 = section_pos1

    if biases:
        bias1, bias2, decay, bads1, bads2 = get_biases_region(biases, bin_coords)
    elif normalizations != ('raw', ):
        raise Exception('ERROR: should provide path to file with biases (pickle).')
    else:
        bads1 = bads2 = {}


    if verbose:
        printime('  - Writing matrices')
    # define output file name
    name = _generate_name(regions, (start1, start2), (end1, end2), resolution,
                          chr_order)

    # prepare file header
    outfiles = []
    if cooler:
        if 'h5py' not in modules:
            raise Exception('ERROR: cooler output is not available. Probably ' +
                            'you need to install h5py\n')
        if 'decay' in normalizations or 'raw&decay' in normalizations:
            raise Exception('ERROR: decay and raw&decay matrices cannot be exported '
                'to cooler format. Cooler only accepts weights per column/row')
        fnam = 'raw_%s_%s%s.mcool' % (name,
                                    nicer(resolution).replace(' ', ''),
                                    ('_' + extra) if extra else '')
        if os.path.exists(os.path.join(outdir, fnam)):
            os.remove(os.path.join(outdir, fnam))
        out_raw = cooler_file(os.path.join(outdir, fnam), resolution, sections, regions)
        out_raw.create_bins()
        out_raw.prepare_matrix(start_bin1, start_bin2)
        outfiles.append((os.path.join(outdir, fnam), fnam))
    else:
        if 'raw' in normalizations:
            fnam = 'raw_%s_%s%s.abc' % (name,
                                        nicer(resolution).replace(' ', ''),
                                        ('_' + extra) if extra else '')
            if append_to_tar:
                out_raw = StringIO()
                outfiles.append((out_raw, fnam))
            else:
                out_raw = open(os.path.join(outdir, fnam), 'w')
                outfiles.append((os.path.join(outdir, fnam), fnam))
            for reg in regions:
                out_raw.write('# CRM %s\t%d\n' % (reg, sections[reg]))

            out_raw.write('# %s resolution:%d\n' % (name, resolution))
            if region2:
                out_raw.write('# BADROWS %s\n' % (','.join([str(b) for b in bads1])))
                out_raw.write('# BADCOLS %s\n' % (','.join([str(b) for b in bads2])))
            else:
                out_raw.write('# MASKED %s\n' % (','.join([str(b) for b in bads1])))

        # write file header
        if 'norm' in normalizations:
            fnam = 'nrm_%s_%s%s.abc' % (name,
                                        nicer(resolution).replace(' ', ''),
                                        ('_' + extra) if extra else '')
            if append_to_tar:
                out_nrm = StringIO()
                outfiles.append((out_nrm, fnam))
            else:
                out_nrm = open(os.path.join(outdir, fnam), 'w')
                outfiles.append((os.path.join(outdir, fnam), fnam))
            for reg in regions:
                out_nrm.write('# CRM %s\t%d\n' % (reg, sections[reg]))

            out_nrm.write('# %s resolution:%d\n' % (name, resolution))
            if region2:
                out_nrm.write('# BADROWS %s\n' % (','.join([str(b) for b in bads1])))
                out_nrm.write('# BADCOLS %s\n' % (','.join([str(b) for b in bads2])))
            else:
                out_nrm.write('# MASKED %s\n' % (','.join([str(b) for b in bads1])))
        if 'decay' in normalizations or 'raw&decay' in normalizations:
            fnam = 'dec_%s_%s%s.abc' % (name,
                                        nicer(resolution).replace(' ', ''),
                                        ('_' + extra) if extra else '')
            if append_to_tar:
                out_dec = StringIO()
                outfiles.append((out_dec, fnam))
            else:
                out_dec = open(os.path.join(outdir, fnam), 'w')
                outfiles.append((os.path.join(outdir, fnam), fnam))
            for reg in regions:
                out_dec.write('# CRM %s\t%d\n' % (reg, sections[reg]))

            out_dec.write('# %s resolution:%d\n' % (
                name, resolution))
            if region2:
                out_dec.write('# BADROWS %s\n' % (','.join([str(b) for b in bads1])))
                out_dec.write('# BADCOLS %s\n' % (','.join([str(b) for b in bads2])))
            else:
                out_dec.write('# MASKED %s\n' % (','.join([str(b) for b in bads1])))

    # functions to write lines of pairwise interactions
    def write_raw(func=None):
        def writer2(c, a, b, v):
            func(c, a, b, v)
            out_raw.write('{}\t{}\n'.format(get_name(a, b), v))
        def writer(_, a, b, v):
            out_raw.write('{}\t{}\n'.format(get_name(a, b), v))
        return writer2 if func else writer

    def write_bias(func=None):
        def writer2(c, a, b, v):
            func(c, a, b, v)
            out_nrm.write('{}\t{}\n'.format(get_name(a, b), v / bias1[a] / bias2[b]))
        def writer(_, a, b, v):
            out_nrm.write('{}\t{}\n'.format(get_name(a, b), v / bias1[a] / bias2[b]))
        return writer2 if func else writer

    def write_expc(func=None):
        def writer2(c, a, b, v):
            func(c, a, b, v)
            out_dec.write('{}\t{}\n'.format(
                get_name(a, b), v / bias1[a] / bias2[b] / decay[c][abs(a-b)]))
        def writer(c, a, b, v):
            out_dec.write('{}\t{}\n'.format(
                get_name(a, b), v / bias1[a] / bias2[b] / decay[c][abs(a-b)]))
        return writer2 if func else writer

    def write_expc_2reg(func=None):
        def writer2(c, a, b, v):
            func(c, a, b, v)
            out_dec.write('{}\t{}\n'.format(
                get_name(a, b), v / bias1[a] / bias2[b] / decay[c][abs((a + start_bin1) - (b + start_bin2))]))
        def writer(c, a, b, v):
            out_dec.write('{}\t{}\n'.format(
                get_name(a, b), v / bias1[a] / bias2[b] / decay[c][abs((a + start_bin1) - (b + start_bin2))]))
        return writer2 if func else writer

    def write_expc_err(func=None):
        def writer2(c, a, b, v):
            func(c, a, b, v)
            try:
                out_dec.write('{}\t{}\n'.format(
                    get_name(a, b), v / bias1[a] / bias2[b] / decay[c][abs(a-b)]))
            except KeyError:  # different chromosomes
                out_dec.write('{}\t{}\n'.format(get_name(a, b), 'nan'))
        def writer(c, a, b, v):
            try:
                out_dec.write('{}\t{}\n'.format(
                    get_name(a, b), v / bias1[a] / bias2[b] / decay[c][abs(a-b)]))
            except KeyError:  # different chromosomes
                out_dec.write('{}\t{}\n'.format(get_name(a, b), 'nan'))
        return writer2 if func else writer

    def write_raw_and_expc(func=None):
        def writer2(c, a, b, v):
            func(c, a, b, v)
            try:
                out_dec.write('{}\t{}\n'.format(
                    get_name(a, b), v, v / bias1[a] / bias2[b] / decay[c][abs(a-b)]))
            except KeyError:  # different chromosomes
                out_dec.write('{}\t{}\n'.format(
                    get_name(a, b), v, v / bias1[a] / bias2[b]))
        def writer(c, a, b, v):
            try:
                out_dec.write('{}\t{}\n'.format(
                    get_name(a, b), v, v / bias1[a] / bias2[b] / decay[c][abs(a-b)]))
            except KeyError:  # different chromosomes
                out_dec.write('{}\t{}\n'.format(
                    get_name(a, b), v, v / bias1[a] / bias2[b]))
        return writer2 if func else writer

    def get_row_name(a, b):
        return '{}\t{}\t{}\t{}\t'.format(*(section_pos1[a] + section_pos2[b]))

    def get_bin_name(a, b):
        return '{}\t{}\t'.format(a, b)

    get_name = get_row_name if row_names else get_bin_name

    write = None
    if 'raw'   in normalizations:
        write = write_raw(write)
    if 'norm'  in normalizations and not cooler:
        write = write_bias(write)
    if 'decay' in normalizations and not cooler:
        if len(regions) in [1, 2]:
            if region2:
                write = write_expc_2reg(write)
            else:
                write = write_expc(write)
        else:
            write = write_expc_err(write)
    if 'raw&decay' in normalizations and not cooler:
        write = write_raw_and_expc(write)

    # pull all sub-matrices and write full matrix
    if region2 is not None:  # already half-matrix in this case
        half_matrix = False

    if cooler:
        for ichunk, c, j, k, v in _iter_matrix_frags(chunks, tmpdir, rand_hash,
                                                     verbose=verbose, clean=clean,
                                                     include_chunk_count=True):
            if j > k:
                continue
            if j not in bads1 and k not in bads2:
                out_raw.write_iter(ichunk, j, k, v)
        out_raw.close()
    else:
        if half_matrix:
            for c, j, k, v in _iter_matrix_frags(chunks, tmpdir, rand_hash,
                                                 verbose=verbose, clean=clean):
                if k > j:
                    continue
                if j not in bads1 and k not in bads2:
                    write(c, j, k, v)
        else:
            for c, j, k, v in _iter_matrix_frags(chunks, tmpdir, rand_hash,
                                                 verbose=verbose, clean=clean):
                if j not in bads1 and k not in bads2:
                    write(c, j, k, v)

    fnames = {}
    if append_to_tar:
        lock = LockFile(append_to_tar)
        with lock:
            archive = taropen(append_to_tar, "a:")
            for fobj, fnam in outfiles:
                fobj.seek(0)
                info = archive.tarinfo(name=fnam)
                info.size=len(fobj.buf)
                archive.addfile(tarinfo=info, fileobj=fobj)
            archive.close()
    else:
        if cooler:
            fnames['RAW'] = out_raw.name
            if 'norm' in normalizations:
                fnam = 'nrm_%s_%s%s.mcool' % (name,
                                            nicer(resolution).replace(' ', ''),
                                            ('_' + extra) if extra else '')
                copyfile(outfiles[0][0], os.path.join(outdir, fnam))
                out_nrm = cooler_file(os.path.join(outdir, fnam), resolution, sections, regions)
                bias_data_row = [1./b if b > 0 else 0 for b in bias1]
                bias_data_col = [1./b if b > 0 else 0 for b in bias2]
                out_nrm.write_weights(bias_data_row, bias_data_col, *bin_coords)
                outfiles.append((os.path.join(outdir, fnam), fnam))
                fnames['NRM'] = os.path.join(outdir, fnam)
        else:
            if 'raw' in normalizations:
                out_raw.close()
                fnames['RAW'] = out_raw.name
            if 'norm' in normalizations:
                out_nrm.close()
                fnames['NRM'] = out_nrm.name
            if 'decay' in normalizations:
                out_dec.close()
                fnames['DEC'] = out_dec.name
            if 'raw&decay' in normalizations:
                out_dec.close()
                fnames['RAW&DEC'] = out_dec.name

    # this is the last thing we do in case something goes wrong
    if clean:
        os.system('rm -rf %s' % (os.path.join(tmpdir, '_tmp_%s' % (rand_hash))))

    return fnames
