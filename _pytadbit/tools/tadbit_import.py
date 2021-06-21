"""

information needed

 - path working directory

"""
from future import standard_library
standard_library.install_aliases()
from argparse                        import HelpFormatter
from multiprocessing                 import cpu_count
from os                              import path, system, remove
from collections                     import OrderedDict
from subprocess                      import Popen, PIPE
import time
import logging
import sqlite3 as lite
from distutils.version               import LooseVersion

from pytadbit                        import HiC_data
from pytadbit.parsers.gzopen         import gzopen
from pytadbit.parsers.hic_parser     import autoreader
from pytadbit.parsers.cooler_parser  import parse_cooler, is_cooler, parse_header
from pytadbit.utils.file_handling    import mkdir, which
from pytadbit                        import get_dependencies_version
from pytadbit.utils.sqlite_utils     import add_path, get_path_id, print_db
from pytadbit.utils.sqlite_utils     import get_jobid, digest_parameters, retry
from pytadbit.parsers.hic_bam_parser import printime
from pytadbit.parsers.hic_bam_parser import _map2sam_mid
from pytadbit.mapping.filter         import MASKED

DESC = 'import Hi-C data to TADbit toy BAM'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()
    param_hash = digest_parameters(opts, extra=['quiet'])
    
    coord1         = opts.coord1

    if not coord1:
        region1 = None
        start1  = None
        end1    = None
    else:
        try:
            crm1, pos1   = coord1.split(':')
            start1, end1 = pos1.split('-')
            region1 = crm1
            start1  = int(start1)
            end1    = int(end1)
        except ValueError:
            region1 = coord1
            start1  = None
            end1    = None

    printime('Importing hic in %s format'%opts.format)
    if opts.format == 'matrix' or opts.format == 'text':
        with gzopen(opts.input) as f_thing:
            masked, chroms_gen, crm, beg, _, _ = read_file_header(f_thing)
        if not chroms_gen or (region1 and region1 not in chroms_gen):
            raise Exception('''ERROR: Chromosome size not included in import file.
                             Please include the chromosome sizes of the data that
                             you want to import in the header of the file. Example:
                             # CRM chr1    249250621''')
    elif opts.format == 'cooler':
        if is_cooler(opts.input, opts.reso if opts.reso > 1 else None):
            chroms_gen = parse_header(opts.input, opts.reso if opts.reso > 1 else None)
            if not chroms_gen or (region1 and region1 not in chroms_gen):
                raise Exception('''ERROR: Chromosome size not included in import file.
                                ''')
        else:
            raise Exception('''ERROR: The input file is not a cooler''')
        
    chroms = OrderedDict((crm,int(chroms_gen[crm]//opts.reso)+1) for crm in chroms_gen)
    sections = []
    if not region1:
        size = 0
        for crm in chroms:
            size += chroms[crm]
            sections.extend([(crm, i) for i in range(chroms[crm])])
    elif not start1:
        size = chroms[region1]
        sections.extend([(region1, i) for i in range(size)])
    else:
        #size = (end1 - start1)//opts.reso
        size = chroms[region1]
        sections.extend([(region1, i) for i in range(start1//opts.reso,
                                                     (end1//opts.reso))])
    dict_sec = dict([(j, i) for i, j in enumerate(sections)])
    
    if opts.format == 'text':
        with gzopen(opts.input) as f_thing:
            matrix = abc_reader(f_thing, size, start1//opts.reso if start1 else None)
        size_mat = size
    elif opts.format == 'matrix':
        with gzopen(opts.input) as in_f:
            matrix, size_mat, _, masked, _ = autoreader(in_f)
        if size != size_mat:
            raise Exception('''ERROR: The size of the specified region is
                            different from the data in the matrix''')
    elif opts.format == 'cooler':
        matrix, size, _, masked, sym = parse_cooler(opts.input,
                                                    opts.reso if opts.reso > 1 else None,
                                                    False)
        size_mat = size
    
    hic = HiC_data(matrix, size_mat, dict_sec=dict_sec,
                   chromosomes=chroms, masked=masked,
                   resolution=opts.reso)
    
    #from pytadbit.mapping.analyze import hic_map
    #hic_map(hic, normalized=False, focus='chr1', show=True, cmap='viridis')
    
    
    outbam = path.join(opts.workdir, '03_filtered_reads',
                       'intersection_%s' % param_hash)

    total_counts = create_BAMhic(hic, opts.cpus, outbam, chroms_gen,
                                 opts.reso, samtools=opts.samtools)

    finish_time = time.localtime()
    # save all job information to sqlite DB
    save_to_db(opts, total_counts, outbam + '.bam', 
               launch_time, finish_time)

def abc_reader(f, size, beg):
    """
    Read matrix stored in 3 column format (bin1, bin2, value)

    :param f: an iterable (typically an open file).

    :returns: An iterator to be converted in dictionary
    """
    fpos = f.tell()
    for line in f:
        if line[0] != '#':
            break
        fpos += len(line)
    f.seek(fpos)
    offset = (beg or 0) * (1 + size)
    def _disect(x):
        a, b, v = x.split()
        return (int(a) + int(b) * size + offset, int(v))
    items = tuple(_disect(line) for line in f)
    return items

def create_BAMhic(hic, ncpus, outbam, chromosomes, reso,
                  masked=None, samtools='samtools'):
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
    # chromosome lengths
    pos_fh = 0

    for chrom in chromosomes:
        output += ("\t".join(("@SQ", "SN:" + chrom, "LN:" + str(chromosomes[chrom]))) + '\n')

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
    map2sam = _map2sam_mid
    
    rownam = [(k[0],k[1] * reso + 1)
              for k in sorted(hic.sections,
                              key=lambda x: hic.sections[x])]
    total_counts=0
    iter_rows = hic.yield_matrix()
    for nrow, row in enumerate(rownam):
        line = next(iter_rows)
        iter_cols = iter(line[nrow:])
        for ncol in range(nrow,len(rownam)):
            col = rownam[ncol]
            val = int(next(iter_cols))
            total_counts += val
            if not val:
                continue
            readid = '%s.%d.%s.%d'% (row[0],nrow, col[0], ncol)
            for nval in range(val):
                line_out = '%s.%d\t%s\t%d\t.\t1\t.\t.\t%s\t%d\t.\t1\t.\t.' % (readid,nval,
                                                        row[0],
                                                        row[1],col[0],
                                                        col[1])
                flag = 0
                proc.stdin.write(map2sam(line_out, flag))    
    proc.stdin.close()
    proc.wait()

    # Index BAM
    _ = Popen(samtools + ' index %s.bam' % (outbam), shell=True, universal_newlines=True).communicate()
    
    return total_counts

@retry(lite.OperationalError, tries=20, delay=2)
def save_to_db(opts, count, outbam, launch_time, finish_time):
    if 'tmpdb' in opts and opts.tmpdb:
        # check lock
        while path.exists(path.join(opts.workdir, '__lock_db')):
            time.sleep(0.5)
        # close lock
        open(path.join(opts.workdir, '__lock_db'), 'a').close()
        # tmp file
        dbfile = opts.tmpdb
        try: # to copy in case read1 was already mapped for example
            copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
        except IOError:
            pass
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='JOBs'""")
        if not cur.fetchall():
            try:
                cur.execute("""
                create table PATHs
                (Id integer primary key,
                JOBid int, Path text, Type text,
                unique (Path))""")
            except lite.OperationalError:
                pass  # may append when mapped files cleaned
            cur.execute("""
            create table JOBs
               (Id integer primary key,
                Parameters text,
                Launch_time text,
                Finish_time text,
                Type text,
                Parameters_md5 text,
                unique (Parameters_md5))""")
        cur.execute("""SELECT name FROM sqlite_master WHERE
                       type='table' AND name='INTERSECTION_OUTPUTs'""")
        if not cur.fetchall():     
            cur.execute("""
                create table INTERSECTION_OUTPUTs
                   (Id integer primary key,
                    PATHid int,
                    Total_interactions int,
                    Multiple_interactions text,
                    Median_fragment_length,
                    MAD_fragment_length,
                    Max_fragment_length,
                    unique (PATHid))""")
            cur.execute("""
                create table FILTER_OUTPUTs
                   (Id integer primary key,
                    PATHid int,
                    Name text,
                    Count int,
                    Applied text,
                    JOBid int,
                    unique (PATHid))""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            param_hash = digest_parameters(opts, get_md5=True )
            cur.execute("""
    insert into JOBs
     (Id  , Parameters, Launch_time, Finish_time,    Type, Parameters_md5)
    values
     (NULL,       '%s',        '%s',        '%s', 'Import',           '%s')
     """ % (parameters,
            time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
            time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass

        jobid = get_jobid(cur)

        add_path(cur, outbam, 'HIC_BAM', jobid, opts.workdir)
        add_path(cur, outbam + '.bai', 'HIC_BAI', jobid, opts.workdir)
        try:
            cur.execute("""
            insert into INTERSECTION_OUTPUTs
            (Id  , PATHid, Total_interactions, Multiple_interactions, Median_fragment_length, MAD_fragment_length, Max_fragment_length)
            values
            (NULL,    %d,                  %d,                  '%s',                     %d,                  %d,                  %d)
            """ % (get_path_id(cur, outbam, opts.workdir),
                   count, '',
                   1, 1, 1))
        except lite.IntegrityError:
            print('WARNING: already filtered')
        try:
            cur.execute("""
        insert into FILTER_OUTPUTs
            (Id  , PATHid, Name, Count, Applied, JOBid)
        values
            (NULL,     %d, '%s',  '%s',    '%s',    %d)
            """ % (get_path_id(cur, outbam, opts.workdir),
                   'valid-pairs', count, '', jobid))
        except lite.IntegrityError:
            print('WARNING: already filtered')
        print_db(cur, 'PATHs')
        #print_db(cur, 'MAPPED_OUTPUTs')
        #print_db(cur, 'PARSED_OUTPUTs')
        print_db(cur, 'JOBs')
        print_db(cur, 'INTERSECTION_OUTPUTs')
        print_db(cur, 'FILTER_OUTPUTs')
    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
    # release lock
    try:
        remove(path.join(opts.workdir, '__lock_db'))
    except OSError:
        pass

def read_file_header(f):
    """
    Read file header, inside first commented lines of a file

    :returns masked dict, chromsomes orderedDict, crm, beg, end, resolution:
    """
    masked = {}
    chromosomes = OrderedDict()
    crm, beg, end, reso = None, None, None, None
    fpos = f.tell()
    for line in f:
        if line[0] != '#':
            break
        fpos += len(line)
        if line.startswith('# MASKED'):
            try:
                masked = dict([(int(n), True) for n in line.split()[-1].split(',')])
            except ValueError:  # nothing here
                pass
        elif line.startswith('# CRM'):
            _, _, crm, size = line.split()
            chromosomes[crm] = int(size)
        elif 'resolution:' in line:
            _, coords, reso = line.split()
            try:
                crm, pos = coords.split(':')
                beg, end = list(map(int, pos.split('-')))
            except ValueError:
                crm = coords
                beg, end = None, None
            reso = int(reso.split(':')[1])
    f.seek(fpos)
    if crm == 'full':
        crm = None
    return masked, chromosomes, crm, beg, end, reso

def check_options(opts):
    if not path.exists(opts.workdir):
        mkdir(opts.workdir)
        # write version log
        vlog_path = path.join(opts.workdir, 'TADbit_and_dependencies_versions.log')
        dependencies = get_dependencies_version()
        if not path.exists(vlog_path) or open(vlog_path).readlines() != dependencies:
            logging.info('Writing versions of TADbit and dependencies')
            vlog = open(vlog_path, 'w')
            vlog.write(dependencies)
            vlog.close()

    mkdir(path.join(opts.workdir, '03_filtered_reads'))
    
    # create empty DB if don't exists
    dbpath = path.join(opts.workdir, 'trace.db')
    open(dbpath, 'a').close()

    # for LUSTRE file system....
    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)
        try:
            copyfile(path.join(opts.workdir, 'trace.db'), opts.tmpdb)
        except IOError:
            pass

    # number of cpus
    if opts.cpus == 0:
        opts.cpus = cpu_count()
    else:
        opts.cpus = min(opts.cpus, cpu_count())

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    oblopt = parser.add_argument_group('Required options')
    glopts = parser.add_argument_group('General options')

    oblopt.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    oblopt.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int, required=True,
                        help='''resolution at which to output matrices''')
    
    oblopt.add_argument('--format', dest='format', action='store',
                        default='text', choices=['text', 'matrix', 'cooler'],
                        help='''[%(default)s] can be any of [%(choices)s]''')

    oblopt.add_argument('-i', '--input', dest='input', metavar="STR",
                        action='store', default=None, type=str, required=True,
                        help='''path to input file''')

    glopts.add_argument('-c', '--coord', dest='coord1',  metavar='',
                        default=None, help='''Coordinate of the region to
                        import. By default all genome, arguments can be
                        either one chromosome name, or the coordinate in
                        the form: "-c chr3:110000000-120000000"''')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument("-C", "--cpus", dest="cpus", type=int,
                        default=cpu_count(), help='''[%(default)s] Maximum
                        number of CPU cores  available in the execution host.
                        If higher than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    glopts.add_argument('--samtools', dest='samtools', metavar="PATH",
                        action='store', default='samtools', type=str,
                        help='''path samtools binary''')

