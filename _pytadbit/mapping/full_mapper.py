"""
09 Jun 2015
"""

import os
from pytadbit.utils.file_handling import mkdir, which
from warnings import warn
from pytadbit.utils.file_handling import magic_open, get_free_space_mb
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, religated
from tempfile import gettempdir, mkstemp
from subprocess import CalledProcessError, PIPE, Popen

def transform_fastq(fastq_path, out_fastq, trim=None, r_enz=None, add_site=True,
                    min_seq_len=15, fastq=True, verbose=True,
                    light_storage=False, **kwargs):
    """
    Given a FASTQ file it can split it into chunks of a given number of reads,
    trim each read according to a start/end positions or split them into
    restriction enzyme fragments

    :param True add_site: when splitting the sequence by ligated sites found,
       removes the ligation site, and put back the original RE site.

    """
    skip = kwargs.get('skip', False)
    ## define local funcitons to process reads and sequences
    def _get_fastq_read_heavy(rlines):
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

    def _get_fastq_read_light(rlines):
        """
        returns header and sequence of 1 FASTQ entry
        Note: header also contains the sequence
        """
        rlines = rlines.rstrip('\n').split()[0][1:]
        seq = fhandler.next()
        _   = fhandler.next()  # lose qualities header but not needed
        qal = fhandler.next()  # lose qualities but not needed
        return (rlines, seq.strip(), qal.strip())

    def _get_map_read_heavy(line):
        header = line.split('\t', 1)[0]
        seq, qal    = header.rsplit(' ', 2)[-2:]
        return header, seq, qal

    def _get_map_read_light(line):
        header, seq, qal, _ = line.split('\t', 3)
        return header, seq, qal

    def _split_read_re(seq, qal, pattern, max_seq_len=None, site='', cnt=0):
        """
        Recursive generator that splits reads according to the
        predefined restriction enzyme.
        RE fragments yielded are followed and preceded by the RE site if a
        ligation site was found after the fragment.

        EXAMPLE:

           seq = '-------oGATCo========oGATCGATCo_____________oGATCGATCo~~~~~~~~~~~~'
           qal = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

        should yield these fragments:
        
            -------oGATCo========oGATC
            xxxxxxxxxxxxxxxxxxxxxxHHHH

            GATCo_____________oGATC
            HHHHxxxxxxxxxxxxxxxHHHH

            GATCo~~~~~~~~~~~~
            HHHHxxxxxxxxxxxxx
        
        """
        cnt += 1
        try:
            pos = seq.index(pattern)
        except ValueError:
            if len(seq) == max_seq_len:
                raise ValueError
            if len(seq) > min_seq_len:
                yield seq, qal, cnt
            return
        xqal = ('H' * len(site))
        if pos < min_seq_len:
            split_read(site + seq[pos + len_relg:],
                       xqal + qal[pos + len_relg:],
                       pattern, max_seq_len, cnt=cnt)
        else:
            yield seq[:pos] + site, qal[:pos] + xqal, cnt
        new_pos = pos + len_relg
        for sseq, sqal, cnt in split_read(site + seq[new_pos:],
                                          xqal + qal[new_pos:], pattern,
                                          max_seq_len, site=site, cnt=cnt):
            yield sseq, sqal, cnt

    # Define function for stripping lines according to focus
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
        enzyme = ''
        enz_pattern = ''
        split_read = lambda x, y, z, after_z, after_after_z: (yield x, y , 1)

    # function to yield reads from input file
    if light_storage:
        get_seq = _get_fastq_read_light if fastq else _get_map_read_light
        insert_mark = insert_mark_light
    else:
        get_seq = _get_fastq_read_heavy if fastq else _get_map_read_heavy
        insert_mark = insert_mark_heavy
        
    ## Start processing the input file
    if verbose:
        print 'Preparing %s file' % ('FASTQ' if fastq else 'MAP')
        if fastq:
            print '  - conversion to MAP format'
        if trim:
            print '  - trimming reads %d-%d' % tuple(trim)
    counter = 0
    if skip:
        if fastq:
            print '    ... skipping, only counting lines'
            counter = sum(1 for _ in magic_open(fastq_path,
                                                cpus=kwargs.get('nthreads')))
            counter /= 4 if fastq else 1
            print '            ' + fastq_path, counter, fastq
        return out_fastq, counter
    # open input file
    fhandler = magic_open(fastq_path, cpus=kwargs.get('nthreads'))
    # create output file
    out_name = out_fastq
    out = open(out_fastq, 'w')
    # iterate over reads and strip them
    site = enzyme if add_site else ''
    for header in fhandler:
        header, seq, qal = get_seq(header)
        counter += 1
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
                                            seq, qal, '0', '-\n'))))
    out.close()
    return out_name, counter

def insert_mark_heavy(header, num):
    if num == 1 :
        return header
    h, s, q = header.rsplit(' ', 2)
    return '%s~%d~ %s %s' % (h, num, s, q)

def insert_mark_light(header, num):
    if num == 1 :
        return header
    return '%s~%d~' % (header, num)

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

def gem_mapping(gem_index_path, fastq_path, out_map_path,
                gem_binary='gem-mapper', **kwargs):
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

    # check that we have the GEM binary:
    gem_binary = which(gem_binary)
    if not gem_binary:
        raise Exception('\n\nERROR: GEM binary not found, install it from:'
                        '\nhttps://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%202/'
                        '\n - Download the GEM-binaries-Linux-x86_64-core_i3 if'
                        'have a recent computer, the '
                        'GEM-binaries-Linux-x86_64-core_2 otherwise\n - '
                        'Uncompress with "tar xjvf GEM-binaries-xxx.tbz2"\n - '
                        'Copy the binary gem-mapper to /usr/local/bin/ for '
                        'example (somewhere in your PATH).\n\nNOTE: GEM does '
                        'not provide any binary for MAC-OS.')

    # mapping
    print 'TO GEM', fastq_path
    kgt = kwargs.get
    gem_cmd = [
        gem_binary, '-I', gem_index_path,
        '-q'                        , kgt('q', 'offset-33'                    ),
        '-m'                        , kgt('m', str(max_edit_distance       )  ),
        '-s'                        , kgt('s', kgt('strata-after-best', '0')  ),
        '--allow-incomplete-strata' , kgt('allow-incomplete-strata', '0.00'   ),
        '--granularity'             , kgt('granularity', '10000'              ),
        '--max-decoded-matches'     , kgt('max-decoded-matches', kgt('d', '1')),
        '--min-decoded-strata'      , kgt('min-decoded-strata', kgt('D', '0') ),
        '--min-insert-size'         , kgt('min-insert-size', '0'              ),
        '--max-insert-size'         , kgt('max-insert-size', '0'              ),
        '--min-matched-bases'       , kgt('min-matched-bases', '0.8'          ),
        '--gem-quality-threshold'   , kgt('gem-quality-threshold', '26'       ),
        '--max-big-indel-length'    , kgt('max-big-indel-length', '15'        ),
        '--mismatch-alphabet'       , kgt('mismatch-alphabet', 'ACGT'         ),
        '-E'                        , kgt('E', '0.30'                         ),
        '--max-extendable-matches'  , kgt('max-extendable-matches', '20'      ),
        '--max-extensions-per-match', kgt('max-extensions-per-match', '1'     ),
        '-e'                        , kgt('e', str(mismatches)                ),
        '-T'                        , str(nthreads),
        '-i'                        , fastq_path,
        '-o', out_map_path.replace('.map', '')]

    if 'paired-end-alignment' in kwargs or 'p' in kwargs:
        gem_cmd.append('--paired-end-alignment')
    if 'map-both-ends' in kwargs or 'b' in kwargs:
        gem_cmd.append('--map-both-ends')
    if 'fast-mapping' in kwargs:
        gem_cmd.append('--fast-mapping')
    if 'unique-mapping' in kwargs:
        gem_cmd.append('--unique-mapping')
    if 'unique-pairing' in kwargs:
        gem_cmd.append('--unique-pairing')

    # check kwargs
    for kw in kwargs:
        if not kw in ['nthreads', 'max_edit_distance',
                      'mismatches', 'max_reads_per_chunk',
                      'out_files', 'temp_dir', 'skip', 'q', 'm', 's',
                      'strata-after-best', 'allow-incomplete-strata',
                      'granularity', 'max-decoded-matches',
                      'min-decoded-strata', 'min-insert-size',
                      'max-insert-size', 'min-matched-bases',
                      'gem-quality-threshold', 'max-big-indel-length',
                      'mismatch-alphabet', 'E', 'max-extendable-matches',
                      'max-extensions-per-match', 'e', 'paired-end-alignment',
                      'p', 'map-both-ends', 'fast-mapping', 'unique-mapping',
                      'unique-pairing', 'suffix']:
            warn('WARNING: %s not in usual keywords, misspelled?' % kw)

    print ' '.join(gem_cmd)
    try:
        # check_call(gem_cmd, stdout=PIPE, stderr=PIPE)
        out, err = Popen(gem_cmd, stdout=PIPE, stderr=PIPE).communicate()
    except CalledProcessError as e:
        print out
        print err
        raise Exception(e.output)

def full_mapping(gem_index_path, fastq_path, out_map_dir, r_enz=None, frag_map=True,
                 min_seq_len=15, windows=None, add_site=True, clean=False,
                 get_nread=False, **kwargs):
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
    :param None windows: tuple of ranges for begining and end of the
       mapping. This parameter allows to do classical iterative mapping, e.g.
         windows=((1,25),(1,30),(1,35),(1,40),(1,45),(1,50))
       A unique window can also be passed, for trimming, like this:
         windows=((1,101),)
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
    :param False get_nreads: returns a list of lists where each element contains
       a path and the number of reads processed

    :returns: a list of paths to generated outfiles. To be passed to 
       :func:`pytadbit.parsers.map_parser.parse_map`
    """

    skip = kwargs.get('skip', False)
    suffix = kwargs.get('suffix', '')
    suffix = ('_' * (suffix != '')) + suffix
    nthreads = kwargs.get('nthreads', 8)
    outfiles = []
    temp_dir = os.path.abspath(os.path.expanduser(
        kwargs.get('temp_dir', gettempdir())))
    # create directories
    for rep in [temp_dir, out_map_dir]:
        mkdir(rep)
    # check space
    fspace = int(get_free_space_mb(temp_dir, div=3))
    if fspace < 200:
        warn('WARNING: only %d Gb left on tmp_dir: %s\n' % (fspace, temp_dir))

    # iterative mapping
    base_name = os.path.split(fastq_path)[-1].replace('.gz', '')
    base_name = base_name.replace('.fastq', '')
    input_reads = fastq_path
    if windows is None:
        light_storage = True
        windows = (None, )
    elif isinstance(windows[0], int):
        # if windows starts at zero we do not need to store all the sequence
        # otherwise we need it because sequence can be trimmed two times
        # in fragment based mapping
        light_storage = True if not windows[0] else False
        windows = [tuple(windows)]
    else:
        # ensure that each element is a tuple, not a list
        windows = [tuple(win) for win in windows]
        # in this case we will need to keep the information about original
        # sequence at any point, light storage is thus not possible.
        light_storage = False
    for win in windows:
        # Prepare the FASTQ file and iterate over them
        curr_map, counter = transform_fastq(
            input_reads, mkstemp(prefix=base_name + '_', dir=temp_dir)[1],
            fastq=(   input_reads.endswith('.fastq'   )
                   or input_reads.endswith('.fastq.gz')
                   or input_reads.endswith('.fq.gz'   )
                   or input_reads.endswith('.dsrc'    )),
            min_seq_len=min_seq_len, trim=win, skip=skip, nthreads=nthreads,
            light_storage=light_storage)
        # clean
        if input_reads != fastq_path and clean:
            print '   x removing original input %s' % input_reads
            os.system('rm -f %s' % (input_reads))
        # First mapping, full length
        if not win:
            beg, end = 1, 'end'
        else:
            beg, end = win
        out_map_path = curr_map + '_full_%s-%s%s.map' % (beg, end, suffix)
        if end:
            print 'Mapping reads in window %s-%s%s...' % (beg, end, suffix)
        else:
            print 'Mapping full reads...', curr_map

        if not skip:
            gem_mapping(gem_index_path, curr_map, out_map_path, **kwargs)
            # parse map file to extract not uniquely mapped reads
            print 'Parsing result...'
            _gem_filter(out_map_path,
                        curr_map + '_filt_%s-%s%s.map' % (beg, end, suffix),
                        os.path.join(out_map_dir,
                                     base_name + '_full_%s-%s%s.map' % (
                                         beg, end, suffix)))
            # clean
            if clean:
                print '   x removing GEM input %s' % curr_map
                os.system('rm -f %s' % (curr_map))
                print '   x removing map %s' % out_map_path
                os.system('rm -f %s' % (out_map_path))
            # for next round, we will use remaining unmapped reads
            input_reads = curr_map + '_filt_%s-%s%s.map' % (beg, end, suffix)
        outfiles.append(
            (os.path.join(out_map_dir,
                          base_name + '_full_%s-%s%s.map' % (beg, end, suffix)),
             counter))

    # map again splitting unmapped reads into RE fragments
    # (no need to trim this time)
    if frag_map:
        if not r_enz:
            raise Exception('ERROR: need enzyme name to fragment.')
        frag_map, counter = transform_fastq(
            input_reads, mkstemp(prefix=base_name + '_', dir=temp_dir)[1],
            min_seq_len=min_seq_len, trim=win, fastq=False, r_enz=r_enz,
            add_site=add_site, skip=skip, nthreads=nthreads,
            light_storage=light_storage)
        # clean
        if clean:
            print '   x removing pre-GEM input %s' % input_reads
            os.system('rm -f %s' % (input_reads))
        if not win:
            beg, end = 1, 'end'
        else:
            beg, end = win
        out_map_path = frag_map + '_frag_%s-%s%s.map' % (beg, end, suffix)
        if not skip:
            print 'Mapping fragments of remaining reads...'
            gem_mapping(gem_index_path, frag_map, out_map_path, **kwargs)
            print 'Parsing result...'
            _gem_filter(out_map_path, curr_map + '_fail%s.map' % (suffix),
                        os.path.join(out_map_dir,
                                     base_name + '_frag_%s-%s%s.map' % (beg, end, suffix)))
        # clean
        if clean:
            print '   x removing GEM input %s' % frag_map
            os.system('rm -f %s' % (frag_map))
            print '   x removing failed to map ' + curr_map + '_fail%s.map' % (suffix)
            os.system('rm -f %s' % (curr_map + '_fail%s.map' % (suffix)))
            print '   x removing tmp mapped %s' % out_map_path
            os.system('rm -f %s' % (out_map_path))
        outfiles.append((os.path.join(out_map_dir,
                                      base_name + '_frag_%s-%s%s.map' % (beg, end, suffix)),
                         counter))
    if get_nread:
        return outfiles
    return [out for out, _ in outfiles]

