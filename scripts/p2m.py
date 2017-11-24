from argparse                     import ArgumentParser
import pysam
from collections import OrderedDict
from datetime import datetime
import subprocess
import glob
import multiprocessing as mu
import functools
import os
from collections import defaultdict
from cPickle import dump,load
import itertools


def extract_coordinates(peak_list,resolution,section_pos,tmpdir,chromosome,badcols):
    '''Chunk file into multiple, and write them in parallel per file write coord of 10,000 peak pairs
    Write a dictionary depending of pairs of peak per target'''

    unique_filename, peak_list = peak_list
    w = open(tmpdir+chromosome+'_{}.tmp'.format(unique_filename),'wa')
    for line in peak_list:
        chr1, beg1, end1, beg2, end2 = line.split()
        pos = section_pos[chr1][0]
        if beg1 < beg2:
            start_bin1 = pos + (int(beg1) / resolution)
            end_bin1 = pos + (int(end1) / resolution)

            start_bin2 = pos + (int(beg2) / resolution)
            end_bin2 = pos + (int(end2) / resolution)
        else:
            start_bin1 = pos + (int(beg2) / resolution)
            end_bin1 = pos + (int(end2) / resolution)

            start_bin2 = pos + (int(beg1) / resolution)
            end_bin2 = pos + (int(end1) / resolution)
        for x, p1 in enumerate(xrange(start_bin1, end_bin1+1)):
            for y, p2 in enumerate(xrange(start_bin2, end_bin2+1)):
                if p1 in badcols or p2 in badcols:
                    continue
                w.write('{}\t{}\t{}\t{}\n'.format(p1,p2,x,y))
    w.close()


def eq_pos(pos1, pos2):
    return pos1 == pos2


def greater_pos(pos1, pos2):
    return pos1 > pos2


def readfiles(file1,file2,chromosome):
    def split_line1(l):
        a, b, c, d = l.split()
        return (int(a), int(b)), int(c), float(d)
    def split_line2(l):
        a, b, c, d = map(int, l.split())
        return (a, b), c, d
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Reading BAM and peaks ...'
    fh1 = open(file1)
    fh2 = open(file2)
    pos1, raw, nrm = split_line1(fh1.next())
    pos2, x, y = split_line2(fh2.next())
    try:
        while True:
            if eq_pos(pos1, pos2):
                avg_raw[(x,y)] += raw
                avg_nrm[(x,y)] += nrm
                avg_pass[(x,y)] += 1
                pos2_ = pos2
                pos2, x, y = split_line2(fh2.next())
                if pos2_ != pos2:  # some cells in the peak file are repeated
                    pos1, raw, nrm = split_line1(fh1.next())
            elif greater_pos(pos1, pos2):
                avg_pass[(x,y)] += 1
                pos2, x, y = split_line2(fh2.next())
            else:
                pos1, raw, nrm = split_line1(fh1.next())

    except StopIteration:
        fh1.close()
        fh2.close()
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Finished'


def read_line(line):
   c, p1, p2 = line.split()[:3]
   return c, (int(p1) + int(p2)) / 2


def binning_bed(peak_file, resolution, windows_span, max_dist, outdir, name,
                chrom_sizes):
    '''Input bed file of ChIP peaks and bin into desired resolution of Hi-C'''

    peaks = open(peak_file,"r")

    peak_coord = set((c, p / resolution)
                     for c, p in map(read_line, peaks)
                     if p > windows_span)  # take into account to add windows span both sides

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Total of different bin coordinates: ', len(bin_coordinate)

    wsp = ((windows_span * 2) / resolution) + 1
    mdr = max_dist / resolution

    # get all combination of peaks
    # - same chromosomes
    # - inside window span
    # - below max_dist
    pairs = ((a[0], a[1], b[1]) for a, b in itertools.combinations(peak_coord, 2)
             if a[0] == b[0] and wsp <= abs(b[1]-a[1]) <= mdr)

    windows = ((255000  , 1000000),
               (1000000 , 2500000),
               (2500000 , 5000000),
               (5000000 , 10000000),
               (10000000, 15000000),
               (15000000, 20000000),
               (20000000, 30000000))

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '- Writing list of coordinates...'
    intervals = dict(((b, e), []) for b, e in (windows))

    for chromosome, bs1, bs2 in pairs:
       distance = abs(bs1 - bs2)
       for (lower,upper) in windows:
           if lower / resolution < distance <= upper / resolution:
               intervals[(lower, upper)].append((chromosome, bs1, bs2))

    adding = windows_span / resolution

    for beg, end in intervals:
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Writing interval: ', beg, end
        w = open('%s%s_%d_%d.tsv'%(outdir, name, beg, end),'wa')
        for c, s1, s2 in intervals[(beg,end)]:
            start1, end1 = s1 - adding, s1 + adding
            start2, end2 = s2 - adding, s2 + adding
            #  check chromosome length
            new_start1, new_end1 = start1 * resolution, end1 * resolution
            new_start2, new_end2 = start2 * resolution, end2 * resolution
            if new_end1 > chrom_sizes[c] or new_end2 > chrom_sizes[c]:
                continue
            else:
                w.write('%s\t%d\t%d\t%d\t%d\n' % (
                    c, new_start1, new_end1, new_start2, new_end2))
        w.close()


def main():
    opts         = get_options()

    inbam        = opts.inbam
    resolution   = opts.resolution
    peak_file    = opts.peak_file
    tmpdir       = opts.tmpdir
    outdir       = opts.outdir
    ncpus        = opts.ncpus
    name         = opts.name
    biases       = opts.biases
    mats         = opts.mats
    windows_span = opts.windows_span
    max_dist     = opts.max_dist

    ## peaks file sorted per chromosome
    bamfile  = pysam.AlignmentFile(inbam,'rb')
    sections = OrderedDict(zip(bamfile.references,[x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    chrom_sizes = OrderedDict(zip(bamfile.references, [x for x in bamfile.lengths]))

    binning_bed(peak_file, resolution, windows_span, max_dist, outdir, name,
                chrom_sizes)

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),'Sublists written!'

    #get sublists in directory
    sublists = glob.glob(outdir + '*.tsv')
    for l in sublists:
        label = l.split('/')[-1].split('.')[0]
        print  datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Getting finger: ', l
        #split file  peaks per chromosome
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Splitting peak pairs per chromosome...'
        fh = subprocess.Popen("awk -F: '{print >> " + '"' +  tmpdir + '"' +  " $1; close($1)}' %s " % (l),
                              shell=True)
        fh.communicate() # wait until the process finishes
        badcols = load(open(biases))['badcol']
        chromosomes_file = glob.glob(tmpdir + '/*')
        global avg_raw
        avg_raw = defaultdict(int)
        global avg_nrm
        avg_nrm = defaultdict(float)
        global avg_pass
        avg_pass = defaultdict(int)

        for peak_list in chromosomes_file:
            chromosome = peak_list.split('/')[-1]
            peakfile = open(peak_list, 'r')
            total_lines = peakfile.readlines()
            total_lines.sort()
            starts = [x for x in range(0, len(total_lines), 10000)]
            lines = {}
            print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), ' - Generate: ', len(starts), ' tmp files for ', chromosome
            for n, i in enumerate(starts):
                if i < starts[-1]:
                    lines[n] = total_lines[i:starts[n + 1]]
                else:
                    lines[n] = total_lines[i:len(total_lines)]

            # parallel job to write coordinates, splitting peak file, write tmp files
            print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Writing coordinates files...'
            pool = mu.Pool(ncpus)
            pool.map(functools.partial(
                extract_coordinates, chromosome=chromosome,
                resolution=resolution, section_pos=section_pos, tmpdir=tmpdir,
                badcols=badcols), lines.iteritems())
            pool.close()
            pool.join()
            print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Written tmps', chromosome
            tmp_chr = [tmpdir + '%s_%d.tmp' % (chromosome, n) for n in lines])

            out = open(tmp + "%s%s_sorted" % (tmpdir, chromosome), 'w')
            for tmpf in tmp_chr:
                out.write(''.join(l for l in tmpf))
                os.system("rm -f " + tmpf)
            out.close()

            # read bam chromosome and peak file same time
            file1 = mats+'%s_bam_5kb.tsv'%(chromosome)
            file2 = '%s%s_sorted'%(tmpdir,chromosome)
            readfiles(file1, file2, chromosome)
            os.system("rm %s%s_sorted"%(tmpdir,chromosome))
            os.system("rm %s%s"%(tmpdir,chromosome))

        out_raw=open(outdir+'raw_%s.pickle'%(label),'wb')
        out_nrm = open(outdir+'nrm_%s.pickle'%(label),'wb')
        out_pas = open(outdir+'pass_%s.pickle'%(label),'wb')
        dump(avg_raw,out_raw)
        dump(avg_nrm,out_nrm)
        dump(avg_pass,out_pas)
        out_raw.close()
        out_nrm.close()
        out_pas.close()

def get_options():
    parser = ArgumentParser(usage="-i Peaks -r INT [options]")

    parser.add_argument('-i','--peak', dest='peak_file',required=True, default=False,
                        help='''Pairwise peaks to compute average submatrix (norm and raw)''')
    parser.add_argument('-bam','--bam',dest='inbam',required=True, default=False,
                        help= 'Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-t','--tmp',dest='tmpdir',required=True, default=False,
                        help='Tmpdir to store coordinates files')
    parser.add_argument('-o','--outdir',dest='outdir',default=True,help='output directory')
    parser.add_argument('-C','--cpus',dest='ncpus',type=int,default=8,
                        help='''[%(default)s] number of cpus to be used for parsing the HiC-BAM file''')
    parser.add_argument('-n','--name',dest='name',default=True, help = 'Output name')
    parser.add_argument('-b','--biases',dest='biases',default=True, help = 'Biases')
    parser.add_argument('-mats','--mats',dest='mats',default=True, help = 'Folder where matrices are located')
    parser.add_argument('-w','-w', dest='windows_span',required=True, default=False,type=int,
                        help='''Windows span around center of the peak''')
    parser.add_argument('-m','-max', dest='max_dist',required=True, default=False,type=int,
                        help='''Max dist between center peaks''')
    opts = parser.parse_args()

    return opts

if __name__=='__main__':
    exit(main())
