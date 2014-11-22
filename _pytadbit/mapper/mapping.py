"""
12 Dec 2013


"""
from os.path import abspath, expanduser, join as joinpath, split as splitpath
from os      import remove
import os
import tempfile
import gem


def trimming(raw_seq_len, seq_start, min_seq_len):
    return seq_start, raw_seq_len - seq_start - min_seq_len


def iterative_mapping_gem(gem_index_path, fastq_path, out_sam_path,
                          min_seq_len, len_step, **kwargs):
    '''Map raw HiC reads iteratively with bowtie2.
    http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Iterative mapping accounts for the modification of fragments sequences
    due to ligation.

    The algorithm of iterative correction:

    1. Truncate the sequences to the first N = **min_seq_len** base pairs,
       starting at the **seq_start** position.
    2. Map the sequences using bowtie2.
    3. Store the uniquely mapped sequences in a SAM file at **out_sam_path**.
    4. Go to the step 1, increase the truncation length N by **len_step** base
       pairs, and map the non-mapped and non-uniquely mapped sequences,
       ...
       Stop when the 3 end of the truncated sequence reaches the **seq_end**
       position.

    Parameters
    ----------

    bowtie_path : str
        The path to the bowtie2 executable.

    bowtie_index_path : str
        The path to the bowtie2 genome index. Since the index consists of
        several files with the different suffices (e.g., hg18.1.bt2,
        hg18.2.bt.2), provide only the common part (hg18).

    fastq_path : str
        The path to the input FASTQ or gzipped FASTQ file.

    out_sam_path : str
        The path to the output SAM file. If ends with .bam, then the output
        is converted to BAM.

    min_seq_len : int
        The truncation length at the first iteration of mapping.

    len_step : int
        The increase in truncation length at each iteration.

    seq_start, seq_end : int, optional
        Slice the FASTQ sequences at [seq_start:seq_end]. Default is [O:None].

    nthreads : int, optional
        The number of Bowtie2 threads. Default is 8

    bowtie_flags : str, optional
        Extra command-line flags for Bowtie2. Default is ''.

    max_reads_per_chunk : int, optional
        If positive then split input into several chunks with
        `max_reads_per_chunk` each and map them separately. Use for large
        datasets and low-memory machines.

    temp_dir : str, optional
        The path to the temporary folder. If not specified, this path is
        supplied by the OS.

    bash_reader : str, optional
        A bash application to convert the input to the FASTQ format. The 
        application have to save the output into stdout.
        The default value is None, that is the app is autodetected by the 
        extension (i.e. cat for .fastq, gunzip for .gz).

    '''
    gem_index_path = os.path.abspath(os.path.expanduser(gem_index_path))
    fastq_path = os.path.abspath(os.path.expanduser(fastq_path))
    out_sam_path = os.path.abspath(os.path.expanduser(out_sam_path))

    single_end = kwargs.get('single_end', False)
    seq_start = kwargs.get('seq_start', 0)
    seq_end = kwargs.get('seq_end', None)
    nthreads = kwargs.get('nthreads', 4)
    max_reads_per_chunk = kwargs.get('max_reads_per_chunk', -1)
    temp_dir = os.path.abspath(os.path.expanduser(
        kwargs.get('temp_dir', tempfile.gettempdir())))

    bash_reader = kwargs.get('bash_reader', None)
    if bash_reader is None:
        extension = fastq_path.split('.')[-1].lower()
        if extension == 'gz':
            bash_reader = 'gunzip -c'
        else:
            bash_reader = 'cat'

    reading_command = bash_reader.split() + [fastq_path, ]
    # print 'READING COMD' , reading_command

    # If bash reader is not 'cat', convert file to FASTQ first and
    # run iterative_mapping recursively on the converted file.
    if bash_reader != 'cat':
        converted_fastq = os.path.join(temp_dir, os.path.split(fastq_path)[1] + '.fastq')
        log.info('Bash reader is not trivial, read input with %s, store in %s',
                 ' '.join(reading_command), converted_fastq)
        converting_process = subprocess.Popen(
            reading_command,
            stdout=open(converted_fastq, 'w'))
        converting_process.wait()

        kwargs['bash_reader'] = 'cat'
        log.info('Run iterative_mapping recursively on %s', converted_fastq)
        iterative_mapping_gem(
            gem_index_path, converted_fastq,
            out_sam_path, min_seq_len, len_step,
            **kwargs)

        os.remove(converted_fastq)
        return

    output_is_bam = (out_sam_path.split('.')[-1].lower() == 'bam') or (out_sam_path.split('.')[-2].lower() == 'bam')


    # Split input files if required and apply iterative mapping to each
    # segment separately.
    if max_reads_per_chunk > 0:
        kwargs['max_reads_per_chunk'] = -1
        log.info('Split input file %s into chunks', fastq_path)
        chunked_files = _chunk_file(
            fastq_path,
            os.path.join(temp_dir, os.path.split(fastq_path)[1]),
            max_reads_per_chunk * 4)
        log.debug('%d chunks obtained', len(chunked_files))
        for i, fastq_chunk_path in enumerate(chunked_files):
            log.info('Run iterative_mapping recursively on %s', fastq_chunk_path)
            iterative_mapping_gem(
                gem_index_path, fastq_chunk_path,
                out_sam_path + '.%d' % (i + 1), min_seq_len, len_step,
                **kwargs)

            # Delete chunks only if the file was really chunked.
            if len(chunked_files) > 1:
                log.info('Remove the chunks: %s', ' '.join(chunked_files))
                os.remove(fastq_chunk_path)
        return

    # Convert input relative arguments to the absolute length scale.
    reading_process = subprocess.Popen(reading_command,
                                       stdout=subprocess.PIPE)
    reading_process.stdout.readline()
    raw_seq_len = len(reading_process.stdout.readline().strip())
    reading_process.terminate()

    if (seq_start < 0
        or seq_start > raw_seq_len
        or (seq_end and seq_end > raw_seq_len)):
        raise Exception('An incorrect trimming region is supplied: [%d, %d], '
                        'the raw sequence length is %d' % (
                            seq_start, seq_end, raw_seq_len))
    local_seq_end = min(raw_seq_len, seq_end) if seq_end else raw_seq_len

    if min_seq_len <= local_seq_end - seq_start:
        trim_5 = seq_start
        trim_3 = raw_seq_len - seq_start - min_seq_len
        local_out_sam = out_sam_path + '.' + str(min_seq_len)
        # mapping_command = [
        #     bowtie_path, '-x', bowtie_index_path, '-q', '-',
        #     '-5', str(trim_5), '-3', str(trim_3), '-p', str(nthreads)
        #     ] + bowtie_flags.split()

        pipeline = []
        cutoff = 0.04
        try:
            log.info('Reading command: %s', ' '.join(reading_command))
            # pipeline.append(
            #     subprocess.Popen(reading_command, stdout=subprocess.PIPE))
            if not os.path.exists(local_out_sam):
                sys.stderr.write('Mapping in ' + local_out_sam + '\n')
                import gem

                # inputf = gem.files.open(pipeline[-1].stdout)
                inputf = gem.files.open(fastq_path)
                trimmed = gem.filter.run_filter(
                    inputf, ['--hard-trim', '%d,%d' % (trim_5, trim_3)],
                    threads=nthreads, paired=not single_end)
                mapped = gem.mapper(trimmed, gem_index_path, min_decoded_strata=0,
                                    max_decoded_matches=2, unique_mapping=False,
                                    max_edit_distance=cutoff, mismatches=cutoff,
                                    # max_edit_distance=0.5, mismatches=0.02,
                                    # output=local_out_sam + '.map',#'/tmp/test.map',
                                    output=temp_dir + '/test.map',
                                    threads=nthreads)
                if output_is_bam:
                    sam = gem.gem2sam(mapped, index=gem_index_path, threads=nthreads,
                                      output=local_out_sam.replace('bam', 'sam'),
                                      single_end=single_end)
                    bam = gem.sam2bam(sam, output=local_out_sam, threads=nthreads)
                else:
                    sam = gem.gem2sam(mapped, index=gem_index_path, output=local_out_sam,
                                      threads=nthreads, single_end=single_end)
            else:
                sys.stderr.write('WARNING: file ' + local_out_sam + ' found, mapping skipped\n')
        finally:
            for process in pipeline:
                if process.poll() is None:
                    process.terminate()

        # Check if the next iteration is required.
        if len_step <= 0:
            return
        if min_seq_len + len_step > local_seq_end - seq_start:
            return

        # Recursively go to the next iteration.
        log.info('Save unique aligments and send the '
                     'non-unique ones to the next iteration')

        unmapped_fastq_path = os.path.join(
            temp_dir, os.path.split(fastq_path)[1] + '.%d' % min_seq_len)
        _filter_unmapped_fastq(fastq_path, local_out_sam, unmapped_fastq_path)

        iterative_mapping_gem(gem_index_path, unmapped_fastq_path,
                              out_sam_path,
                              min_seq_len=min_seq_len + len_step,
                              len_step=len_step, **kwargs)

        os.remove(unmapped_fastq_path)

    

def _filter_unmapped_fastq(in_fastq, in_sam, nonunique_fastq):
    '''Read raw sequences from **in_fastq** and alignments from
    **in_sam** and save the non-uniquely aligned and unmapped sequences
    to **unique_sam**.
    '''
    samfile = pysam.Samfile(in_sam)

    nonunique_ids = set()
    for read in samfile:
        tags_dict = dict(read.tags)
        read_id = read.qname
        # If exists, the option 'XS' contains the score of the second
        # best alignment. Therefore, its presence means non-unique alignment.
        if 'XS' in tags_dict or read.is_unmapped or (
            'NM' in tags_dict and int(tags_dict['NM']) > 1):
            nonunique_ids.add(read_id)

    _filter_fastq(nonunique_ids, in_fastq, nonunique_fastq)


def _filter_fastq(ids, in_fastq, out_fastq):
    '''Filter FASTQ sequences by their IDs.

    Read entries from **in_fastq** and store in **out_fastq** only those
    the whose ID are in **ids**.
    '''
    out_file = open(out_fastq, 'w')
    in_file = _gzopen(in_fastq)
    while True:
        line = in_file.readline()
        if not line:
            break

        if not line.startswith('@'):
            raise Exception(
                '{0} does not comply with the FASTQ standards.'.format(in_fastq))

        fastq_entry = [line, in_file.readline(),
                       in_file.readline(), in_file.readline()]
        read_id = line.split()[0][1:]
        if read_id.endswith('/1') or read_id.endswith('/2'):
            read_id = read_id[:-2]
        if read_id in ids:
            out_file.writelines(fastq_entry)

def _chunk_file(in_path, out_basename, max_num_lines):
    '''Slice lines from a large file.
    The line numbering is as in Python slicing notation.
    '''
    num_lines = _line_count(in_path)
    if num_lines <= max_num_lines:
        return [in_path, ]

    out_paths = []

    for i, line in enumerate(open(in_path)):
        if i % max_num_lines == 0:
            out_path = out_basename + '.%d' % (i // max_num_lines + 1)
            out_paths.append(out_path)
            out_file = file(out_path, 'w')
        out_file.write(line)

    return out_paths

def _gzopen(path):
    if path.endswith('.gz'):
        return gzip.open(path)
    else:
        return open(path)

