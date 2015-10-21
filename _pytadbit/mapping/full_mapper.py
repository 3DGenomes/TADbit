"""
09 Jun 2015
"""

import os
from pytadbit.utils.file_handling import mkdir
from warnings import warn
from pytadbit.utils.file_handling import magic_open, get_free_space_mb
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, religated
from tempfile import gettempdir, mkstemp
try:
    import gem
except ImportError:
    warn('WARNING: GEMTOOLS not found')

def transform_fastq(fastq_path, out_fastq, trim=None, r_enz=None, add_site=True,
                    min_seq_len=15, fastq=True, verbose=True):
    """
    Given a FASTQ file it can split it into chunks of a given number of reads,
    trim each read according to a start/end positions or split them into
    restriction enzyme fragments

    :param True add_site: when splitting the sequence by ligated sites found,
       removes the ligation site, and put back the original RE site.

    """
    ## define local funcitons to process reads and sequences
    def _get_fastq_read(rlines):
        """
        returns header and sequence of 1 FASTQ entry
        Note: header also contains the sequence
        """
        rlines = rlines.rstrip('\n').split()[0][1:]
        seq = fhandler.next()
        _   = fhandler.next()  # lose qualities header but not needed
        qal = fhandler.next()  # lose qualities but not needed
        # header now also contains original read
        return (rlines + ' ' + seq.strip() + ' ' + qal.strip(),
                seq.strip(), qal.strip())

    def _get_map_read(line):
        header = line.split('\t', 1)[0]
        seq, qal    = header.rsplit(' ', 2)[-2:]
        return header, seq, qal
        
    def _split_read_re(seq, qal, pattern, max_seq_len=None, site='', cnt=0):
        """
        Recursive generator that splits reads according to the
        predefined restriction enzyme.
        RE fragments yielded are followed by the RE site if a ligation
        site was found after the fragment.
        The RE site before the fragment is added outside this function
        """
        try:
            cnt += 1
            pos = seq.index(pattern)
            if pos < min_seq_len:
                split_read(seq[pos + len_relg:], qal[pos + len_relg:],
                           pattern, max_seq_len, cnt=cnt)
            else:
                yield seq[:pos] + site, qal[:pos] + ('H' * len(site)), cnt
            for subseq, subqal, cnt in split_read(seq[pos + len_relg:],
                                             qal[pos + len_relg:],
                                             pattern,
                                             max_seq_len, cnt=cnt):
                yield subseq, subqal, cnt
        except ValueError:
            if len(seq) == max_seq_len:
                raise ValueError
            if len(seq) > min_seq_len:
                yield seq, qal, cnt

    # Define function for stripping lines according to ficus
    if isinstance(trim, tuple):
        beg, end = trim
        beg -= 1
        strip_line = lambda x: x[beg:end]
    else:
        strip_line = lambda x: x

    # define function to split reads according to restriction enzyme sites
    if isinstance(r_enz, str):
        enzyme = RESTRICTION_ENZYMES[r_enz].replace('|', '')
        enz_pattern = religated(r_enz)
        sub_enz_pattern = enz_pattern[:len(enz_pattern) / 2]
        len_relg = len(enz_pattern)
        print '  - splitting into restriction enzyme (RE) fragments using ligation sites'
        print '  - ligation sites are replaced by RE sites to match the reference genome'
        print '    * enzyme: %s, ligation site: %s, RE site: %s' % (r_enz, enz_pattern, enzyme)
        split_read = _split_read_re
    else:
        enz_pattern = ''
        split_read = lambda x, y, z, after_z, after_after_z: (yield x, y , 1)

    # function to yield reads from input file
    get_seq = _get_fastq_read if fastq else _get_map_read

    ## Start processing the input file
    if verbose:
        print 'Preparing %s file' % ('FASTQ' if fastq else 'MAP')
        if fastq:
            print '  - conversion to MAP format'
        if trim:
            print '  - trimming reads %d-%d' % tuple(trim)
            
    # open input file
    fhandler = magic_open(fastq_path)
    # create output file
    out_name = out_fastq
    out = open(out_fastq, 'w')
    # iterate over reads and strip them
    site = '' if add_site else enzyme
    for header in fhandler:
        header, seq, qal = get_seq(header)
        # trim on wanted region of the read
        seq = strip_line(seq)
        qal = strip_line(qal)
        # get the generator of restriction enzyme fragments
        iter_frags = split_read(seq, qal, enz_pattern, len(seq), site)
        # the first fragment should not be preceded by the RE site
        try:
            seq, qal, cnt = iter_frags.next()
        except StopIteration:
            # read full of ligation events, fragments not reaching minimum
            continue
        except ValueError:
            # or not ligation site found, in which case we try with half
            # ligation site in case there was a sequencing error (half ligation
            # site is a RE site or nearly, and thus should not be found anyway)
            iter_frags = split_read(seq, qal, sub_enz_pattern, len(seq), '')
            try:
                seq, qal, cnt = iter_frags.next()
            except ValueError:
                continue
            except StopIteration:
                continue
        out.write(_map2fastq('\t'.join((insert_mark(header, cnt),
                                        seq, qal, '0', '-\n'))))
        # the next fragments should be preceded by the RE site
        # continue
        for seq, qal, cnt in  iter_frags:
            out.write(_map2fastq('\t'.join((insert_mark(header, cnt),
                                            seq + site, qal + 'H' * (len(site)),
                                            '0', '-\n'))))
    out.close()
    return out_name

def insert_mark(header, num):
    if num == 1 :
        return header
    h, s, q = header.rsplit(' ', 2)
    return '%s~%d~ %s %s' % (h, num, s, q)

def _map2fastq(read):
    return '@{0}\n{1}\n+\n{2}\n'.format(*read.split('\t', 3)[:-1])


def _gem_filter(fnam, unmap_out, map_out):
    """
    Divides reads in a map file in two categories: uniquely mapped, and not.
    Writes them in two files
    
    Notes:
       - GEM unique-maps can not be used as it gets rid of reads like 1:0:0:5
       - not feasible with gt.filter
    """
    fhandler = magic_open(fnam) if isinstance(fnam, str) else fnam
    unmap_out = open(unmap_out, 'w')
    map_out   = open(map_out  , 'w')
    def _strip_read_name(line):
        """
        remove original sequence from read name when read is mapped uniquely
        """
        header, line = line.split('\t', 1)
        return '\t'.join((header.rsplit(' ', 2)[0], line))
    for line in fhandler:
        matches = line.rsplit('\t', 2)[1]
        bad = False
        if matches != '1':
            for m in matches.replace('+', ':').split(':'):
                if m == '0':
                    continue
                if  m != '1':
                    bad = True
                    unmap_out.write(line)
                    break
                break
            else:
                bad = True
                unmap_out.write(line)
        if not bad:
            map_out.write(_strip_read_name(line))
    unmap_out.close()


def gem_mapping(gem_index_path, fastq_path, out_map_path, **kwargs):
    """
    :param None focus: trims the sequence in the input FASTQ file according to a
       (start, end) position, or the name of a restriction enzyme. By default it
       uses the full sequence.
    :param 33 quality: set it to 'ignore' in order to speed-up the mapping
    """
    gem_index_path    = os.path.abspath(os.path.expanduser(gem_index_path))
    fastq_path        = os.path.abspath(os.path.expanduser(fastq_path))
    out_map_path      = os.path.abspath(os.path.expanduser(out_map_path))
    nthreads          = kwargs.get('nthreads'            , 8)
    max_edit_distance = kwargs.get('max_edit_distance'   , 0.04)
    mismatches        = kwargs.get('mismatches'          , 0.04)
    quality           = kwargs.get('quality'             , 33)

    # check kwargs
    for kw in kwargs:
        if not kw in ['nthreads', 'max_edit_distance',
                      'mismatches', 'max_reads_per_chunk',
                      'out_files', 'temp_dir']:
            warn('WARNING: %s not is usual keywords, misspelled?' % kw)

    # input
    inputf = gem.files.open(fastq_path)

    # mapping
    print 'TO GEM', fastq_path
    return gem.mapper(inputf, gem_index_path, min_decoded_strata=0,
                      max_decoded_matches=1, unique_mapping=False,
                      max_edit_distance=max_edit_distance,
                      mismatches=mismatches, quality=quality,
                      output=out_map_path,
                      threads=nthreads)

def full_mapping(gem_index_path, fastq_path, out_map_dir, r_enz=None, frag_map=True,
                 min_seq_len=15, windows=((None, None),), add_site=True, clean=False,
                 **kwargs):
    """
    Do the mapping

    :param gem_index_path: path to index file created from a reference genome
       using gem-index tool
    :param fastq_path: PATH to fastq file, either compressed or not.
    :param out_map_dir: path to a directory where to store mapped reads in MAP
       format .
    :param None r_enz: name of the restriction enzyme used in the experiment e.g.
       HindIII. This is optional if frag_map option is False
    :param True frag_map: two step mapper, first full length is mapped, then
       remaining, unmapped reads, are divided into restriction-enzyme fragments
       andeach is mapped.
    :param True add_site: when splitting the sequence by ligated sites found,
       removes the ligation site, and put back the original RE site.
    :param 15 min_seq_len: minimum size of a fragment to map
    :param ((None, None),) windows: tuple of ranges for begining and end of the
       mapping. This parameter allows to do classical iterative mapping.
    :param False clean: remove intermedite files created in temp_dir
    :param 4 nthreads: number of threads to use for mapping (number of CPUs)
    :param 0.04 max_edit_distance: The maximum number of edit operations allowed
       while verifying candidate matches by dynamic programming.
    :param 0.04 mismatches: The maximum number of nucleotide substitutions
       allowed while mapping each k-mer. It is always guaranteed that, however
       other options are chosen, all the matches up to the specified number of
       substitutions will be found by the program.
    :param /tmp temp_dir: important to change. Intermediate FASTQ files will be
       written there.

    :returns: a list of paths to generated outfiles. To be passed to 
       :func:`pytadbit.parsers.map_parser.parse_map`
    """
    outfiles = []
    temp_dir = os.path.abspath(os.path.expanduser(
        kwargs.get('temp_dir', gettempdir())))
    # create directories
    for rep in [temp_dir, out_map_dir]:
        mkdir(rep)
    # check space
    if get_free_space_mb(temp_dir, div=3) < 50:
        warn('WARNING: less than 50 Gb left on tmp_dir: %s\n' % temp_dir)

    # iterative mapping
    base_name = os.path.split(fastq_path)[-1].replace('.gz', '')
    base_name = base_name.replace('.fastq', '')
    input_reads = fastq_path
    for beg, end in windows:
        # Prepare the FASTQ file and iterate over them
        curr_map = transform_fastq(input_reads, 
                                   mkstemp(prefix=base_name + '_',
                                           dir=temp_dir)[1],
                                   fastq=(input_reads.endswith('.fastq')
                                          or input_reads.endswith('.fastq.gz')),
                                   min_seq_len=min_seq_len, trim=(beg, end))
        # clean
        if input_reads != fastq_path and clean:
            print '   x removing original input %s' % input_reads
            os.system('rm -f %s' % (input_reads))
        # First mapping, full length
        out_map_path = curr_map + '_full_%s-%s.map' % (beg, end)
        if end:
            print 'Mapping reads in window %s-%s...' % (beg, end)
        else:
            print 'Mapping full reads...', curr_map
        map_file = gem_mapping(gem_index_path, curr_map, out_map_path, **kwargs)
        map_file.close()

        # parse map file to extract not uniquely mapped reads
        print 'Parsing result...'
        _gem_filter(out_map_path, curr_map + '_filt_%s-%s.map' % (beg, end),
                    os.path.join(out_map_dir,
                                 base_name + '_full_%s-%s.map' % (beg, end)))
        # clean
        if clean:
            print '   x removing GEM input %s' % curr_map
            os.system('rm -f %s' % (curr_map))
            print '   x removing map %s' % out_map_path
            os.system('rm -f %s' % (out_map_path))
        # for next round, we will use remaining unmapped reads
        input_reads = curr_map + '_filt_%s-%s.map' % (beg, end)
        outfiles.append(os.path.join(out_map_dir,
                                     base_name + '_full_%s-%s.map' % (beg, end)))

    # map again splitting unmapped reads into RE fragments
    # (no need to trim this time)
    if frag_map:
        if not r_enz:
            raise Exception('ERROR: need enzyme name to fragment.')
        frag_map = transform_fastq(input_reads,
                                   mkstemp(prefix=base_name + '_',
                                           dir=temp_dir)[1],
                                   min_seq_len=min_seq_len, trim=(beg, end),
                                   fastq=False, r_enz=r_enz, add_site=add_site)
        out_map_path = frag_map + '_frag.map'
        print 'Mapping fragments of remaining reads...'
        map_file = gem_mapping(gem_index_path, frag_map, out_map_path,
                               **kwargs)
        map_file.close()
        print 'Parsing result...'
        _gem_filter(out_map_path, curr_map + '_fail.map',
                    os.path.join(out_map_dir, base_name + '_frag.map'))
        outfiles.append(os.path.join(out_map_dir, base_name + '_frag.map'))
    return outfiles

def main():

    fastq          = '/scratch/db/FASTQs/hsap/dixon_2012/dixon-2012_200bp.fastq'
    fastq          = 'short_dixon-2012_200bp.fastq'
    # fastq        = '/scratch/test/sample_dataset/FASTQs/sample_hsap_HindIII.fastq'
    gem_index_path = '/scratch/db/index_files/Homo_sapiens-79/Homo_sapiens.gem'
    out_map_dir1   = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/'
    out_map_dir2   = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/'
    temp_dir1      = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/tmp1/'
    temp_dir2      = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/tmp2/'

    print 'read 1'
    outfiles1 = full_mapping(gem_index_path, fastq, out_map_dir1, 'HindIII',
                             temp_dir=temp_dir1, windows=((1,100),), add_site=True)
    print 'read 2'
    outfiles2 = full_mapping(gem_index_path, fastq, out_map_dir2, 'HindIII',
                             temp_dir=temp_dir2, windows=((101, 200),), add_site=True)
    # print 'read 1'
    # outfiles1 = mapping(gem_index_path, fastq, out_map_dir1, 'HindIII',
    #                     temp_dir=temp_dir1,
    #                     windows=(zip(*([0] * len(range(25, 105, 5)),
    #                                    range(25,105,5)))))
    # print 'read 2'
    # outfiles2 = mapping(gem_index_path, fastq, out_map_dir2, 'HindIII',
    #                     temp_dir=temp_dir2,
    #                     windows=(zip(*([100] * len(range(125, 205, 5)),
    #                                            range(125,205,5)))))
    
    print outfiles1
    print 'xcmvnkljnv'
    print outfiles2
    
    from pytadbit.parsers.map_parser import parse_map
    from pytadbit.parsers.genome_parser import parse_fasta
    from pytadbit.mapping.mapper import get_intersection
    from pytadbit.mapping.filter import filter_reads, apply_filter
    
    read1, read2 = 'read1.tsv', 'read2.tsv',
    parse_map(outfiles1, outfiles2, out_file1=read1, out_file2=read2,
              genome_seq=parse_fasta('/scratch/db/index_files/Homo_sapiens-79/Homo_sapiens.fa'),
              re_name='HindIII', verbose=True)

    reads = 'both_reads.tsv'
    get_intersection(read1, read2, reads)

    masked = filter_reads(reads)
    freads = 'filtered_reads.tsv'
    apply_filter(reads, freads, masked)

if __name__ == "__main__":
    exit(main())


"""
outfiles1 = ['/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-25.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-30.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-35.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-40.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-45.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-50.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-55.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-60.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-65.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-70.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-75.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-80.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-85.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-90.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-95.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_full_0-100.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read1/short_dixon-2012_200bp_frag.map']
outfiles2 = ['/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-125.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-130.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-135.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-140.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-145.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-150.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-155.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-160.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-165.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-170.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-175.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-180.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-185.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-190.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-195.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_full_100-200.map',
             '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/read2/short_dixon-2012_200bp_frag.map']
"""


