"""
18 may 2015
"""
import os
from warnings import warn
from pytadbit.utils.file_handling import magic_open, get_free_space_mb
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, religated
from itertools import chain, islice
from tempfile import gettempdir, mkstemp
try:
    import gem
except ImportError:
    warn('WARNING: GEMTOOLS not found')

def chunks(iterable, n):
    """
    chunks(ABCDE,2) => AB CD E
    From stackexchange answer:
    http://codereview.stackexchange.com/questions/57395/split-large-file-into-smaller-files
    """
    while True:
        yield chain([next(iterable)], islice(iterable, n - 1))

def transform_fastq(fastq_path, out_fastq, trim=None, r_enz=None,
                    min_seq_len=20, max_reads_per_chunk=None, fastq=True):
    """
    Given a FASTQ file it can split it into chunks of a given number of reads,
    trim each read according to a start/end positions or split them into
    restriction enzyme fragments
    """
    fhandler = magic_open(fastq_path)
    for chunk, lines in enumerate(chunks(fhandler, max_reads_per_chunk) if
                                  max_reads_per_chunk else [fhandler]):
        if max_reads_per_chunk:
            print 'Preparing FASTQ file chunk', chunk + 1
        else:
            print 'Preparing FASTQ file'
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

            def split_read(x, max_seq_len=None):
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
        else:
            split_read = lambda x, y: (yield x)
        # iterate over reads in current chunk and strip them
        if max_reads_per_chunk:
            out_name = out_fastq + 'chunk_%d' % chunk
            out = open(out_fastq + 'chunk_%d' % chunk, 'w')
        else:
            out_name = out_fastq
            out = open(out_fastq, 'w')
        if fastq:
            def get_seq(rlines):
                rlines = rlines.rstrip('\n')
                line = fhandler.next()
                _ = fhandler.next()  # lose qualities but not needed
                _ = fhandler.next()  # lose qualities but not needed
                return rlines, line.strip()
        else:
            get_seq = lambda x: x.split('\t', 2)[:2]
        for header in lines:
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
            out.write('\t'.join((header, frag, '!' * len(frag), '0', '-\n')))
            # the next fragments should be preceded by the RE site
            for frag in  iter_frags:
                out.write('\t'.join((header, frag + enzyme,
                                     '!' * (len(frag) + len(enzyme)), '0', '-\n')))
        out.close()
        yield out_name
        

def _gem_filter(fnam, unmap_out):
    """
    not feasible with gt.filter
    """
    fhandler = magic_open(fnam) if isinstance(fnam, str) else fnam
    unmap_out = open(unmap_out, 'w')
    # ok = open('tmp/of.lala', 'w')
    for line in fhandler:
        # bad = False
        matches = line.rsplit('\t', 2)[1]
        if matches != '1':
            for m in matches.replace('+', ':').split(':'):
                if m == '0':
                    continue
                if  m != '1':
                    # bad = True
                    unmap_out.write(line)
                    break
                break
            else:
                # bad = True
                unmap_out.write(line)
        # if not bad:
        #     ok.write(line)                
    unmap_out.close()


def gem_mapping(gem_index_path, fastq_path, out_map_path, **kwargs):
    """
    :param None focus: trims the sequence in the input FASTQ file according to a
       (start, end) position, or the name of a restriction enzyme. By default it
       uses the full sequence.
    """
    gem_index_path    = os.path.abspath(os.path.expanduser(gem_index_path))
    fastq_path        = os.path.abspath(os.path.expanduser(fastq_path))
    out_map_path      = os.path.abspath(os.path.expanduser(out_map_path))
    nthreads          = kwargs.get('nthreads'            , 8)
    max_edit_distance = kwargs.get('max_edit_distance'   , 0.04)
    mismatches        = kwargs.get('mismatches'          , 0.04)

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
                      mismatches=mismatches, quality='ignore',
                      output=out_map_path,
                      threads=nthreads)


def adaptive_mapping(gem_index_path, fastq_path, out_map_dir, r_enz, trim=None,
                     max_reads_per_chunk=None, min_seq_len=20, **kwargs):
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
    # Prepare the FASTQ file and iterate over them
    for full_map in transform_fastq(fastq_path,
                                    mkstemp(prefix='fastq_', dir=temp_dir)[1],
                                    min_seq_len=min_seq_len, trim=trim,
                                    max_reads_per_chunk=max_reads_per_chunk):
        # First mapping, full length
        out_map_path = full_map + '_full.map'
        map_file = gem_mapping(gem_index_path, full_map, out_map_path, **kwargs)
        map_file.close()
        # parse map file to extract not uniquely mapped reads
        _gem_filter(out_map_path, full_map + '_filt.map')
        # map again splitting unmapped reads into RE fragments
        for frag_map in transform_fastq(out_map_path,
                                        mkstemp(prefix='fastq_', dir=temp_dir)[1],
                                        min_seq_len=min_seq_len, trim=trim,
                                        max_reads_per_chunk=max_reads_per_chunk,
                                        fastq=False, r_enz=r_enz):
            out_map_path = frag_map + '_frag.map'
            map_file = gem_mapping(gem_index_path, frag_map, out_map_path,
                                   **kwargs)
            map_file.close()
    
def main():

    fastq = '/scratch/test/sample_dataset/FASTQs/sample_hsap_HindIII.fastq'
    gem_index_path = '/scratch/db/index_files/Homo_sapiens-79/Homo_sapiens.gem'
    out_map_dir = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/tmp/'
    temp_dir = '/home/fransua/Box/tadbits/tadbit/_pytadbit/mapping/tmp/'
    adaptive_mapping(gem_index_path, fastq, out_map_dir, 'HindIII',
                     temp_dir=temp_dir, trim=(0,100))
    
if __name__ == "__main__":
    exit(main())
