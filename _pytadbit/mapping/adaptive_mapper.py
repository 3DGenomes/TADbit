"""
18 may 2015
"""
import os
from warnings import warn
from pytadbit.utils.file_handling import magic_open, get_free_space_mb
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, religated
from itertools import chain, islice
from tempfile import gettempdir
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
                    min_seq_len=20, max_reads_per_chunk=None):
    fhandler = magic_open(fastq_path)
    """
    Given a FASTQ file it can split it into chunks of a given number of reads,
    trim each read according to a start/end positions or split them into
    restriction enzyme fragments
    """
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

            def split_read(x):
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
                        split_read(x[pos + len_relg:])
                    else:
                        yield x[:pos] + enzyme
                    for x in split_read(x[pos + len_relg:]):
                        yield x
                except ValueError:
                    if len(x) > min_seq_len:
                        yield x
        else:
            split_read = lambda x: [x]

        # iterate over reads in current chunk and strip them
        if max_reads_per_chunk:
            out_name = out_fastq + 'chunk_%d' % chunk
            out = open(out_fastq + 'chunk_%d' % chunk, 'w')
        else:
            out_name = out_fastq
            out = open(out_fastq, 'w')
        lenz = len(enzyme) if isinstance(r_enz, str) else None
        for header in lines:
            line = fhandler.next().rstrip('\n')
            _ = fhandler.next()  # lose qualities but not needed
            _ = fhandler.next()  # lose qualities but not needed
            # trim on wanted region of the read
            line = strip_line(line)
            # get the generator of restriction enzyme fragments
            iter_frags = split_read(line)
            # the first fragment should not be preceded by the RE site
            try:
                frag = iter_frags.next()
            except StopIteration:
                # read full of ligation events, fragments not reaching minimum
                continue
            out.write(''.join((header, frag, '\n+\n', '!' * len(frag),'\n')))
            # the next fragments should be preceded by the RE site
            for frag in  iter_frags:
                out.write(''.join((header, enzyme, frag, '\n+\n',
                                   '!' * (lenz + len(frag)),'\n')))
        out.close()
        yield out_name


def gem_mapping(gem_index_path, fastq_path, out_sam_path, **kwargs):
    """
    :param None focus: trims the sequence in the input FASTQ file according to a
       (start, end) position, or the name of a restriction enzyme. By default it
       uses the full sequence.
    """
    gem_index_path      = os.path.abspath(os.path.expanduser(gem_index_path))
    fastq_path          = os.path.abspath(os.path.expanduser(fastq_path))
    out_sam_path        = os.path.abspath(os.path.expanduser(out_sam_path))
    nthreads            = kwargs.get('nthreads'            , 4)
    max_edit_distance   = kwargs.get('max_edit_distance'   , 0.04)
    mismatches          = kwargs.get('mismatches'          , 0.04)
    max_reads_per_chunk = kwargs.get('max_reads_per_chunk' , -1)
    temp_dir = os.path.abspath(os.path.expanduser(
        kwargs.get('temp_dir', tempfile.gettempdir())))

    # check kwargs
    for kw in kwargs:
        if not kw in ['nthreads', 'max_edit_distance',
                      'mismatches', 'max_reads_per_chunk',
                      'out_files', 'temp_dir']:
            warn('WARNING: %s not is usual keywords, misspelled?' % kw)
    
    # create directories
    for rep in [temp_dir, os.path.split(out_sam_path)[0]]:
        try:
            os.mkdir(rep)
        except OSError, error:
            if error.strerror != 'File exists':
                raise error


    # prepare the FASTQ
    

def adaptive_mapping(gem_index_path, fastq_path, out_sam_path, r_enz, trim=None,
                     max_reads_per_chunk=None, min_seq_len=20, **kwargs):
    """
    Do the mapping
    """
    temp_dir = os.path.abspath(os.path.expanduser(
        kwargs.get('temp_dir', gettempdir())))
    # check space
    if get_free_space_mb(temp_dir, div=3) < 50:
        warn('WARNING: less than 50 Gb left on tmp_dir: %s\n' % temp_dir)
    # Prepare the file
    for fastq in transform_fastq(fastq_path, temp_dir, trim=trim,
                                 max_reads_per_chunk=max_reads_per_chunk):
    
        gem_mapping(gem_index_path, fastq, out_sam_path, **kwargs)
    
    # First mapping, full length
    transform_fastq(fastq_path, temp_dir, trim=trim, r_enz=r_enz)
    
    


