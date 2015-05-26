"""
18 may 2015
"""
import os
from warnings import warn
from pytadbit.utils.file_handling import magic_open, get_free_space_mb
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, religated
from tempfile import gettempdir, mkstemp
try:
    import gem
except ImportError:
    warn('WARNING: GEMTOOLS not found')

def transform_fastq(fastq_path, out_fastq, trim=None, r_enz=None,
                    min_seq_len=20, fastq=True, verbose=True):
    """
    Given a FASTQ file it can split it into chunks of a given number of reads,
    trim each read according to a start/end positions or split them into
    restriction enzyme fragments
    """
    ## define local funcitons to process reads and sequences
    def _get_fastq_read(rlines):
        """
        returns header and sequence of 1 FASTQ entry
        """
        rlines = rlines.rstrip('\n')
        line = fhandler.next()
        _ = fhandler.next()  # lose qualities but not needed
        _ = fhandler.next()  # lose qualities but not needed
        return rlines, line.strip()

    def _split_read_re(x, max_seq_len=None):
        """
        Recursive generator that splits reads according to the
        predefined restriction enzyme.
        RE fragments yielded are followed by the RE site if a ligation
        site was found after the fragment.
        The RE site before the fragment is added outside this function
        """
        try:
            pos = x.index(enz_pattern)
            if pos < min_seq_len:
                split_read(x[pos + len_relg:], max_seq_len)
            else:
                yield x[:pos] + enzyme
            for x in split_read(x[pos + len_relg:], max_seq_len):
                yield x
        except ValueError:
            if len(x) > min_seq_len:
                if len(x) == max_seq_len:
                    raise StopIteration
                yield x

    # Define function for stripping lines according to ficus
    if isinstance(trim, tuple):
        beg, end = trim
        strip_line = lambda x: x[beg:end]
    else:
        strip_line = lambda x: x

    # define function to split reads according to restriction enzyme sites
    if isinstance(r_enz, str):
        enzyme = RESTRICTION_ENZYMES[r_enz].replace('|', '')
        enz_pattern = religated(r_enz)
        len_relg = len(enz_pattern)
        print '  - splitting into restriction enzyme (RE) fragments using ligation sites'
        print '  - ligation sites are replaced by RE sites to match the reference genome'
        print '    * enzyme: %s, ligation site: %s, RE site: %s' % (r_enz, enz_pattern, enzyme)
        split_read = _split_read_re
    else:
        split_read = lambda x, y: (yield x)

    # function to yield reads from input file
    get_seq = _get_fastq_read if fastq else lambda x: x.split('\t', 2)[:2]

    ## Start processing the input file
    if verbose:
        print 'Preparing %s file' % ('FASTQ' if fastq else 'MAP')
        if fastq:
            print '  - conversion to MAP format'
        if trim:
            print '  - triming reads %d-%d' % tuple(trim)


    # open input file
    fhandler = magic_open(fastq_path)
    # create output file
    out_name = out_fastq
    out = open(out_fastq, 'w')
    # iterate over reads and strip them
    for header in fhandler:
        header, line = get_seq(header)
        # trim on wanted region of the read
        line = strip_line(line)
        # get the generator of restriction enzyme fragments
        iter_frags = split_read(line, len(line))
        # the first fragment should not be preceded by the RE site
        try:
            frag = iter_frags.next()
        except StopIteration:
            # read full of ligation events, fragments not reaching minimum
            continue
        out.write('\t'.join((header, frag, 'H' * len(frag), '0', '-\n')))
        # the next fragments should be preceded by the RE site
        for frag in  iter_frags:
            out.write('\t'.join((header, frag + enzyme,
                                 'H' * (len(frag) + len(enzyme)), '0', '-\n')))
    out.close()
    return out_name
        

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
    map_out = open(map_out, 'w')
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
            map_out.write(line)
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
    return gem.mapper(inputf, gem_index_path, min_decoded_strata=0,
                      max_decoded_matches=1, unique_mapping=False,
                      max_edit_distance=max_edit_distance,
                      mismatches=mismatches, quality=quality,
                      output=out_map_path,
                      threads=nthreads)

def fragment_based_mapping(gem_index_path, fastq_path, out_map_dir, r_enz,
                           trim=None, min_seq_len=20, **kwargs):
    """
    Do the mapping
    """
    temp_dir = os.path.abspath(os.path.expanduser(
        kwargs.get('temp_dir', gettempdir())))
    # create directories
    for rep in [temp_dir, out_map_dir]:
        try:
            os.mkdir(rep)
        except OSError, error:
            if error.strerror != 'File exists':
                warn('ERROR: problem loading file, probable problem with the ' +
                     'use of relative path')
                raise error
    # check space
    if get_free_space_mb(temp_dir, div=3) < 50:
        warn('WARNING: less than 50 Gb left on tmp_dir: %s\n' % temp_dir)
    base_name = os.path.split(fastq_path)[-1].replace('.gz', '')
    base_name = base_name.replace('.fastq', '')
    # Prepare the FASTQ file and iterate over them
    full_map = transform_fastq(fastq_path,
                               mkstemp(prefix='fastq_', dir=temp_dir)[1],
                               min_seq_len=min_seq_len, trim=trim)
    # First mapping, full length
    out_map_path = full_map + '_full.map'
    print 'Mapping full reads...'
    map_file = gem_mapping(gem_index_path, full_map, out_map_path, **kwargs)
    map_file.close()

    # parse map file to extract not uniquely mapped reads
    print 'Parsing result...'
    _gem_filter(out_map_path, full_map + '_filt.map',
                out_map_dir + base_name + '_full.map')

    # map again splitting unmapped reads into RE fragments
    # (no need to trim this time)
    frag_map = transform_fastq(out_map_path,
                               mkstemp(prefix='fastq_', dir=temp_dir)[1],
                               min_seq_len=min_seq_len,
                               fastq=False, r_enz=r_enz)
    out_map_path = frag_map + '_frag.map'
    print 'Mapping fragments of remaining reads...'
    map_file = gem_mapping(gem_index_path, frag_map, out_map_path,
                           **kwargs)
    map_file.close()
    print 'Parsing result...'
    _gem_filter(out_map_path, full_map + '_fail.map',
                out_map_dir + base_name + '_frag.map')

    
def main():

    fastq = '/scratch/test/sample_dataset/FASTQs/sample_hsap_HindIII.fastq'
    gem_index_path = '/scratch/db/index_files/Homo_sapiens-79/Homo_sapiens.gem'
    out_map_dir = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/tmp/'
    temp_dir = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/tmp/'
    fragment_based_mapping(gem_index_path, fastq, out_map_dir, 'HindIII',
                           temp_dir=temp_dir, trim=(0,100))


if __name__ == "__main__":
    exit(main())

# from pytadbit.parsers.genome_parser import parse_fasta

# genome_seq = parse_fasta('/scratch/db/index_files/Homo_sapiens-79/Homo_sapiens.fa')

# frag_seq = {}

# for crm in genome_seq:
#     for i in xrange(0, len(genome_seq[crm])):
#         if genome_seq[crm][i:i+4] == 'GATC':
#             frag_seq[crm+'_'+str(i)] = genome_seq[crm][max(0, i-500):i+504]

# out = open('frag_based_genome.fa', 'w')
# for frag in frag_seq:
#     out.write('>%s\n%s\n' % (frag, frag_seq[frag]))
# out.close()
