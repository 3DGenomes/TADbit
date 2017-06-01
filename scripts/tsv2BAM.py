#!/usr/bin/env python

"""
script provided by Enrique Vidal <enrique.vidal@crg.eu> to convert 2D beds
into compressed BAM format.

Gets the *_both_filled_map.tsv contacts from TADbit (and the corresponding filter files)
 and outputs a modified SAM with the following fields:

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

from argparse    import ArgumentParser
from collections import OrderedDict
from subprocess  import Popen, PIPE
import sys
import os

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
        flag += 1024 # filter_keys['trans'] = 2**10
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
        flag += 1024 # filter_keys['trans'] = 2**10
    r1r2 = ('{0}\t{1}\t{2}\t{3}\t0\t{4}P\t{6}\t{7}\t{5}\t*\t*\t'
            'TC:i:{8}\t'
            'S1:i:{9}\tS2:i:{10}\n'
            '{0}\t{1}\t{6}\t{7}\t0\t{5}S\t{2}\t{3}\t{4}\t*\t*\t'
            'TC:i:{8}\t'
            'S2:i:{9}\tS1:i:{10}\n').format(
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
        flag += 1024 # filter_keys['trans'] = 2**10
    r1r2 = ('{0}\t{1}\t{2}\t{3}\t0\t{4}P\t{6}\t{7}\t{5}\t*\t*\t'
            'TC:i:{8}\tS1:i:{13}\tS2:i:{14}\t'
            'E1:i:{9}\tE2:i:{10}\tE3:i:{11}\tE4:i:{12}\n'
            '{0}\t{1}\t{6}\t{7}\t0\t{5}S\t{2}\t{3}\t{4}\t*\t*\t'
            'TC:i:{8}\tS2:i:{14}\tS1:i:{13}\t'
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


def generate_BAM(infile, valid, ncpus, outbam, frmt):
    # define filter codes
    filter_keys = OrderedDict()
    filter_keys['self-circle']        = 2 ** 0
    filter_keys['dangling-end']       = 2 ** 1
    filter_keys['error']              = 2 ** 2
    filter_keys['extra-dangling-end'] = 2 ** 3
    filter_keys['too-close-from-RES'] = 2 ** 4
    filter_keys['too-short']          = 2 ** 5
    filter_keys['too-large']          = 2 ** 6
    filter_keys['over-represented']   = 2 ** 7
    filter_keys['duplicated']         = 2 ** 8
    filter_keys['random-breaks']      = 2 ** 9
    filter_keys['trans-chromosomic']  = 2 ** 10

    output = ''

    # write header
    output += ("\t".join(("@HD" ,"VN:1.5", "SO:queryname")) + '\n')

    fhandler = open(infile)
    line = fhandler.next()
    # chromosome lengths
    pos_fh = 0

    while line.startswith('#'):
        (_, _, cr, ln) = line.replace("\t", " ").strip().split(" ")
        output += ("\t".join(("@SQ", "SN:" + cr, "LN:" + ln)) + '\n')
        pos_fh += len(line)
        line = fhandler.next()

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
    output += ("\t".join(("@CO" ,"S1:i", "Strand of the 1st read-end (1: positive, 0: negative)")))
    output += ("\t".join(("@CO" ,"S2:i", "Strand of the 2nd read-end  (1: positive, 0: negative)")))

    # open and init filter files
    if not valid:
        filter_line, filter_handler = get_filters(infile)
    fhandler.seek(pos_fh)
    proc = Popen('samtools view -Shb -@ %d - | samtools sort -@ %d - %s' % (
        ncpus, ncpus, outbam),
                 shell=True, stdin=PIPE)
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
                        filter_line[i] = filter_handler[i].next().strip()
                    except StopIteration:
                        pass
            # get output in sam format
            proc.stdin.write(map2sam(line, flag))
    proc.stdin.close()
    proc.wait()

    # Index BAM
    _ = Popen('samtools index %s.bam' % (outbam), shell=True).communicate()

    # close file handlers
    fhandler.close()
    if not valid:
        for i in filter_handler:
            filter_handler[i].close()


def get_filters(infile):
    """
    get all filters
    """
    basename = os.path.basename(infile)
    dirname = os.path.dirname(infile)

    filter_files = {}
    sys.stderr.write('Using filter files:\n')
    for fname in os.listdir(dirname):
        if fname.startswith(basename + "_"):
            key = fname.replace(basename + "_", "").replace(".tsv", "")
            filter_files[key] = dirname + "/" + fname
            sys.stderr.write('   - %-20s %s\n' %(key, fname))
    filter_handler = {}
    filter_line = {}
    for i in filter_files:
        filter_handler[i.replace('_', '-')] = open(filter_files[i])
        try:
            filter_line[i.replace('_', '-')] = filter_handler[i.replace('_', '-')].next().strip()
        except StopIteration:
            filter_line[i.replace('_', '-')] = ""
    return filter_line, filter_handler


def main():
    opts = get_options()
    infile = os.path.realpath(opts.inbed)
    generate_BAM(infile, opts.valid, opts.ncpus, opts.outbam, opts.format)


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbed', metavar='',
                        required=True, default=False,
                        help="input TADbit's 2D bed.")
    parser.add_argument('-o', '--outfile', dest='outbam', metavar='',
                        required=True, default=False,
                        help="output TADbit's BAM")
    parser.add_argument('--format', dest='format', metavar='',
                        default='mid', choices=['short', 'mid', 'long'],
                        help=("[%(default)s] output format, in terms of number "
                              "of extra tags (can be any of: %(choices)s)"))
    parser.add_argument('--cpus', dest='ncpus', metavar='',
                        default=8,
                        help="Number of threads for compressing/sorting BAM")
    parser.add_argument('--valid', dest='valid', action='store_true',
                        default=False, help='input already filtered')
    opts = parser.parse_args()

    return opts

if __name__ == '__main__':
    exit(main())
