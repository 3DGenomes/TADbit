"""

information needed

 - path working directory with parsed reads

"""
from future import standard_library
standard_library.install_aliases()
from argparse                        import HelpFormatter
from os                              import path, remove, system, rename
from sys                             import stdout
from pickle                          import load
from multiprocessing                 import cpu_count
from collections                     import OrderedDict
from subprocess                      import Popen, PIPE

from pysam                           import AlignmentFile

from pytadbit.tools.tadbit_bin       import load_parameters_fromdb
from pytadbit.mapping.filter         import MASKED
from pytadbit.utils.file_handling    import mkdir
from pytadbit.parsers.hic_bam_parser import filters_to_bin, printime
from pytadbit.parsers.hic_bam_parser import write_matrix, get_matrix
from pytadbit.utils.sqlite_utils     import digest_parameters

DESC = 'export Hi-C data to other formats'

def run(opts):
    check_options(opts)
    param_hash = digest_parameters(opts, extra=['quiet'])
    opts.normalizations = ['norm' if opts.norm else 'raw']
    biases = None

    clean = True  # change for debug
    if opts.bam:
        mreads = path.realpath(opts.bam)
        if not opts.biases and opts.norm:
            raise Exception('ERROR: external BAM input, should provide path to'
                            ' biases file.')
    else:
        biases, mreads = load_parameters_fromdb(opts)
        mreads = path.join(opts.workdir, mreads)
        biases = path.join(opts.workdir, biases) if biases else None
    if opts.biases:
        biases = opts.biases

    coord1         = opts.coord1
    coord2         = opts.coord2

    if coord2 and not coord1:
        coord1, coord2 = coord2, coord1

    if not coord1:
        region1 = None
        start1  = None
        end1    = None
        region2 = None
        start2  = None
        end2    = None
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
        if coord2:
            try:
                crm2, pos2   = coord2.split(':')
                start2, end2 = pos2.split('-')
                region2 = crm2
                start2  = int(start2)
                end2    = int(end2)
            except ValueError:
                region2 = coord2
                start2  = None
                end2    = None
        else:
            region2 = None
            start2  = None
            end2    = None

    outdir = path.join(opts.workdir, '05_sub-matrices')
    mkdir(outdir)
    tmpdir = path.join(opts.workdir, '05_sub-matrices',
                       '_tmp_sub-matrices_%s' % param_hash)
    mkdir(tmpdir)

    if region1:
        if region1:
            if not opts.quiet:
                stdout.write('\nExtraction of %s' % (region1))
            if start1:
                if not opts.quiet:
                    stdout.write(':%s-%s' % (start1, end1))
            else:
                if not opts.quiet:
                    stdout.write(' (full chromosome)')
            if region2:
                if not opts.quiet:
                    stdout.write(' intersection with %s' % (region2))
                if start2:
                    if not opts.quiet:
                        stdout.write(':%s-%s\n' % (start2, end2))
                else:
                    if not opts.quiet:
                        stdout.write(' (full chromosome)\n')
            else:
                if not opts.quiet:
                    stdout.write('\n')
    else:
        if not opts.quiet:
            stdout.write('\nExtraction of %s genome\n'
                         %('partial' if opts.chr_name else 'full'))

    norm = 'norm' if opts.norm else 'raw'
    
    if opts.format == 'matrix' or opts.format == 'hic':
        bamfile = AlignmentFile(mreads, 'rb')
        bam_refs = bamfile.references
        bam_lengths = bamfile.lengths
        if opts.chr_name:
            bam_refs_idx = [bam_refs.index(chr_ord)
                            for chr_ord in opts.chr_name if chr_ord in bam_refs]
            if not bam_refs_idx :
                raise Exception('''ERROR: Wrong number of chromosomes in chr_order.
                    Found %s in bam file \n''' % (' '.join(bam_refs)))
            bam_refs = [bam_ref for bam_ref in [bam_refs[bam_ref_idx]
                                                  for bam_ref_idx in bam_refs_idx]]
            bam_lengths = [bam_len for bam_len in [bam_lengths[bam_ref_idx]
                                                     for bam_ref_idx in bam_refs_idx]]
        sections = OrderedDict(list(zip(bam_refs,
                                   [x for x in bam_lengths])))
        printime('Getting %s matrices' % norm)
        matrix, bads1, bads2, regions, name, bin_coords = get_matrix(
            mreads, opts.reso,
            load(open(biases, 'rb')) if biases and norm != 'raw' else None,
            normalization=norm, filter_exclude=opts.filter,
            region1=region1, start1=start1, end1=end1,
            region2=region2, start2=start2, end2=end2,
            tmpdir=tmpdir, ncpus=opts.cpus,
            return_headers=True,
            nchunks=opts.nchunks, verbose=not opts.quiet,
            clean=clean, chr_order=opts.chr_name)

        b1, e1, b2, e2 = bin_coords
        b1, e1 = 0, e1 - b1
        b2, e2 = 0, e2 - b2

        if opts.format == 'matrix':
            if opts.row_names:
                starts = [start1, start2]
                ends = [end1, end2]
                row_names = (
                    (reg, p + 1 , p + opts.reso)
                    for r, reg in enumerate(regions)
                    for p in range(starts[r] if r < len(starts)
                                   and starts[r] else 0,
                                   ends[r] if r < len(ends)
                                   and ends[r] else sections[reg],
                                   opts.reso))
    
            printime(' - Writing: %s' % norm)
            out = open(opts.out, 'w')
            for reg in regions:
                out.write('# CRM %s\t%d\n' % (reg, sections[reg]))
            if region2:
                out.write('# BADROWS %s\n' % (','.join([str(b) for b in bads1])))
                out.write('# BADCOLS %s\n' % (','.join([str(b) for b in bads2])))
            else:
                out.write('# MASKED %s\n' % (','.join([str(b) for b in bads1])))
            if opts.row_names:
                out.write('\n'.join('%s\t%d\t%d\t' % (next(row_names)) +
                                    '\t'.join(str(matrix.get((i, j), 0))
                                              for i in range(b1, e1))
                                    for j in range(b2, e2)) + '\n')
            else:
                out.write('\n'.join('\t'.join(str(matrix.get((i, j), 0))
                                              for i in range(b1, e1))
                                    for j in range(b2, e2)) + '\n')
            out.close()
        else:
            printime(' - Writing: %s' % norm)
            tmp_chromsize = path.join(tmpdir,'hic_%s.chrom.sizes'% param_hash)
            out = open(tmp_chromsize, 'w')
            for reg in regions:
                out.write('%s\t%d\n' % (reg, sections[reg]))
            out.close()
            tmpfl = path.join(tmpdir,'hic_export_%s.tsv'% param_hash)
            out = open(tmpfl, 'w')
            out_ln = '0\t%s\t%d\t0\t1\t%s\t%d\t1\t1%f'if opts.norm else '0\t%s\t%d\t0\t1\t%s\t%d\t1\t1%d'                
            if region1:
                starts = [start1, start2]
                ends = [end1, end2]
                row_names = [(reg, pos+1) for r, reg in enumerate(regions)
                             for pos in range(starts[r] if r < len(starts)
                                            and starts[r] else 0,
                                            ends[r] if r < len(ends)
                                            and ends[r] else sections[reg],
                                            opts.reso)]
                out.write('\n'.join(out_ln % (row_names[i][0],
                                          row_names[i][1],
                                          row_names[j][0],
                                          row_names[j][1],
                                          matrix.get((i,j),0))
                                    for i in range(b1, e1)
                                    for j in range(i, e2)))
            else:
                totals = OrderedDict()
                total_num = 0
                for c in sections:
                    totals[c] = (total_num, total_num + sections[c] // opts.reso + 1)
                    total_num += sections[c] // opts.reso + 1
                
                for crm1_id, crm1 in enumerate(sections):
                    b1, e1 = totals[crm1]
                    row_names1 = dict((b1+ipos,pos+1) for ipos,pos in enumerate(range(0,sections[crm1],opts.reso)))
                    for crm2 in list(sections.keys())[crm1_id:]:
                        b2, e2 = totals[crm2]
                        row_names2 = dict((b2+ipos,pos+1) for ipos,pos in enumerate(range(0,sections[crm2],opts.reso)))
                        out.write('\n'.join(out_ln % (crm1,row_names1[i],
                                                      crm2,row_names2[j],
                                                      matrix.get((i,j),0))
                                    for i in range(b1, e1)
                                    for j in range(max(b2,i), e2)))
            
            out.close()
            do_norm = '-n' if opts.norm else ''
            _ = Popen('java -Xmx32g -jar %s pre -j %d %s %s %s %s'%(opts.juicerjar,
                                                         opts.cpus,
                                                         do_norm,
                                                         tmpfl,
                                                         opts.out,
                                                         tmp_chromsize), shell=True, universal_newlines=True).communicate()
    elif opts.format == 'text':
        printime('Getting and writing matrix to text format')
        fnames = write_matrix(
            mreads, opts.reso,
            load(open(biases, 'rb')) if biases else None,
            outdir, filter_exclude=opts.filter,
            normalizations=[norm],
            region1=region1, start1=start1, end1=end1,
            region2=region2, start2=start2, end2=end2,
            tmpdir=tmpdir, append_to_tar=None, ncpus=opts.cpus,
            nchunks=opts.nchunks, verbose=not opts.quiet,
            extra=param_hash, cooler=False, clean=clean,
            chr_order=opts.chr_name)
        rename(list(fnames.values())[0],opts.out)
    elif opts.format == 'cooler':
        printime('Getting and writing matrix to cooler format')
        fnames = write_matrix(
            mreads, opts.reso,
            load(open(biases, 'rb')) if biases else None,
            outdir, filter_exclude=opts.filter,
            normalizations=[norm],
            region1=region1, start1=start1, end1=end1,
            region2=region2, start2=start2, end2=end2,
            tmpdir=tmpdir, append_to_tar=None, ncpus=opts.cpus,
            nchunks=opts.nchunks, verbose=not opts.quiet,
            extra=param_hash, cooler=True, clean=clean,
            chr_order=opts.chr_name)
        rename(fnames['NRM' if opts.norm else 'RAW'],opts.out)
        if 'NRM' in fnames and not opts.norm:
            remove(fnames['NRM'])
        if 'RAW' in fnames and opts.norm:
            remove(fnames['RAW'])

    if clean:
        printime('Cleaning')
        system('rm -rf %s '% tmpdir)


def check_options(opts):
    mkdir(opts.workdir)

    # transform filtering reads option
    opts.filter = filters_to_bin(opts.filter)

    if not path.exists(opts.workdir):
        raise IOError('ERROR: workdir not found.')

    if opts.format=='hic':
        if not opts.juicerjar:
            raise IOError('ERROR: juicer jar file needed for "hic" export.')
            
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
    rfiltr = parser.add_argument_group('Read filtering options')
    normpt = parser.add_argument_group('Normalization options')
    outopt = parser.add_argument_group('Output options')

    oblopt.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    oblopt.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int, required=True,
                        help='''resolution at which to output matrices''')
    
    oblopt.add_argument('--format', dest='format', action='store',
                        default='text', choices=['text', 'matrix', 'cooler', 'hic'],
                        help='''[%(default)s] can be any of [%(choices)s]''')

    oblopt.add_argument('-o', '--output', dest='out', metavar="STR",
                        action='store', default=None, type=str, required=True,
                        help='''path to output file''')

    glopts.add_argument('--bam', dest='bam', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to a TADbit-generated BAM file with
                        all reads (other wise the tool will guess from the
                        working directory database)''')

    glopts.add_argument('-j', '--jobid', dest='jobid', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')

    glopts.add_argument('--force', dest='force', action='store_true',
                        default=False,
                        help='overwrite previously run job')

    glopts.add_argument('-q', '--quiet', dest='quiet', action='store_true',
                        default=False,
                        help='remove all messages')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument('--nchunks', dest='nchunks', action='store', default=100,
                        type=int,
                        help='''maximum number of chunks into which to
                        cut the BAM''')

    glopts.add_argument("-C", "--cpus", dest="cpus", type=int,
                        default=cpu_count(), help='''[%(default)s] Maximum
                        number of CPU cores  available in the execution host.
                        If higher than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

    glopts.add_argument('--chr_name', dest='chr_name', metavar="STR", nargs='+',
                        default=None, type=str,
                        help='''[fasta header] chromosome name(s). Order of chromosomes
                        in the output matrices.''')

    glopts.add_argument('--juicerjar', dest='juicerjar', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to the juicer tools jar file needed to export
                        matrices to hic format (check https://github.com/aidenlab/juicer/wiki/Download).
                        Note that you also need java available in the path.''')

    outopt.add_argument('--rownames', dest='row_names', action='store_true',
                        default=False,
                        help='''To store row names in the output text matrix.
                        WARNING: when non-matrix, results in two extra columns''')

    outopt.add_argument('-c', '--coord', dest='coord1',  metavar='',
                        default=None, help='''Coordinate of the region to
                        retrieve. By default all genome, arguments can be
                        either one chromosome name, or the coordinate in
                        the form: "-c chr3:110000000-120000000"''')

    outopt.add_argument('-c2', '--coord2', dest='coord2',  metavar='',
                        default=None, help='''Coordinate of a second region to
                        retrieve the matrix in the intersection with the first
                        region.''')

    normpt.add_argument('--biases',   dest='biases', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to file with pre-calculated biases by
                        columns''')

    normpt.add_argument('--norm', dest='norm', action='store_true',
                        default=False,
                        help='export normalized matrix')

    rfiltr.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                        choices = list(range(0, 11)),
                        help=("""[%(default)s] Use filters to define a set of
                        valid pair of reads e.g.:
                        '--filter 1 2 3 4 8 9 10'. Where these numbers """ +
                              "correspond to: 0: nothing, %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))


