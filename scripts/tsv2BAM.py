#!/usr/bin/env python

"""
script provided by Enrique Vidal <enrique.vidal@crg.eu> to convert 2D beds
into compressed BAM format.

"""
#
# Gets the *_both_filled_map.tsv contacts from TADBIT (and the corresponding filter files)
#  and outputs a modified SAM with the following fields:
#
#    - read ID
#    - filtering flag (see codes in header)
#    - chromosome ID of the first pair of the contact
#    - genomic position of the first pair of the contact
#    - mapped length of the first pair of the contact
#    - pseudo CIGAR with strand orientation of the first pair (0M = -, 1M = +)
#    - chromosome ID of the second pair of the contact
#    - genomic position of the second pair of the contact
#    - SIGNED mapped length of the first pair of the contact (sign ~ orientation)
#    - sequence is missing (*)
#    - quality is missing (*)
#    - TC tag indicating single (1) or multi (2) contact
#
# Each pair of contacts porduces two lines in the output SAM
#


# dependencies

from argparse                     import ArgumentParser
import sys
import os
import collections

# get arguments

def map2sam (line, flag):
    """
    translate map + flag into hic-sam (two lines per contact)
    """
    (qname,
     rname, pos, s1, l1, e1, e2,
     rnext, pnext, s2, l2, e3, e4) = line.strip().split('\t')

    # store mapped length and strand in a single value
    mapq = ('-' * (s1 == '0')) + l1
    tlen = ('-' * (s2 == '0')) + l2

    # multicontact?
    tc = str(int("~" in qname) + 1)
    # trans contact?
    if rname != rnext:
        flag += 1024 # filter_keys['trans'] = 2**10
    r1r2 = ('{0}\t{1}\t{2}\t{3}\t{4}\t1M\t{6}\t{7}\t{8}\t*\t*\t'
            'TC:i:{10}\tE1:i:{11}\tE2:i:{12}\tE3:i:{13}\tE4:i:{14}\n'
            '{0}\t{1}\t{6}\t{7}\t{5}\t1M\t{2}\t{3}\t{9}\t*\t*\t'
            'TC:i:{10}\tE3:i:{13}\tE4:i:{14}\tE1:i:{11}\tE2:i:{12}\n').format(
                qname,               # 0
                flag,                # 1
                rname,               # 2
                pos,                 # 3
                l1,                  # 4
                l2,                  # 5 
                rnext,               # 6 
                pnext,               # 7 
                tlen,                # 8
                mapq,                # 9
                tc,                  # 10 
                e1,                  # 11
                e2,                  # 12
                e3,                  # 13
                e4)                  # 14
    return r1r2

def main():
    opts = get_options()
    infile = os.path.realpath(opts.inbed)
    
    # define filter codes
    filter_keys = collections.OrderedDict()
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
    filter_keys['trans']              = 2 ** 10
        
    # write header
    print "\t".join(("@HD" ,"VN:1.5", "SO:queryname"))
    
    fhandler = open(infile)
    line = fhandler.next()
    # chromosome lengths
    pos_fh = 0
    
    while line.startswith('#'):
        (_, _, cr, ln) = line.replace("\t", " ").strip().split(" ")
        print "\t".join(("@SQ", "SN:" + cr, "LN:" + ln))
        pos_fh += len(line)
        line = fhandler.next()
    
    # filter codes
    for i in filter_keys:
        print "\t".join(("@CO", "filter:" + i, "flag:" + str(filter_keys[i])))
    
    # tags
    print "\t".join(("@CO" ,"TC:i", "Multicontact? 0 = no 1 = yes"))
    print "\t".join(("@CO" ,"E1:i", "Position of the left RE site of first read"))
    print "\t".join(("@CO" ,"E2:i", "Position of the right RE site of first read"))
    print "\t".join(("@CO" ,"E3:i", "Position of the left RE site of second read"))
    print "\t".join(("@CO" ,"E4:i", "Position of the right RE site of second read"))    
    
    # open and init filter files
    if not opts.valid:
        filter_line, filter_handler = get_filters(infile)
    
    fhandler.seek(pos_fh)
    if opts.valid:
        for line in fhandler:
            flag = 0
            # get output in sam format
            sys.stdout.write(map2sam(line, flag))
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
            sys.stdout.write(map2sam(line, flag))
    
    # close file handlers
    fhandler.close()
    if not opts.valid:
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


def get_options():
    parser = ArgumentParser(usage="%(prog)s -i PATH -r INT [options]")

    parser.add_argument('-i', '--infile', dest='inbed', metavar='',
                        required=True, default=False,
                        help="input TADbit's 2D bed.")
    parser.add_argument('--valid', dest='valid', action='store_true',
                        default=False, help='input already filtered')    
    opts = parser.parse_args()
    
    return opts

if __name__ == '__main__':
    exit(main())
