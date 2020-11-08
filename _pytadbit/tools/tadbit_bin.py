"""

information needed

 - path working directory with parsed reads

"""
from future import standard_library
standard_library.install_aliases()
from argparse                        import HelpFormatter
from os                              import path, remove, system
from sys                             import stdout
from shutil                          import copyfile
from string                          import ascii_letters
from random                          import random
from pickle                          import load
from warnings                        import warn
from multiprocessing                 import cpu_count
from collections                     import OrderedDict
import time

import sqlite3 as lite
from numpy                           import zeros_like
from numpy                           import array
from numpy                           import ma, log, log2
from matplotlib                      import pyplot as plt
from matplotlib.ticker               import FuncFormatter
from pysam                           import AlignmentFile

from pytadbit.mapping.filter         import MASKED
from pytadbit.utils.file_handling    import mkdir
from pytadbit.parsers.hic_bam_parser import filters_to_bin, printime
from pytadbit.parsers.hic_bam_parser import write_matrix, get_matrix
from pytadbit.parsers.tad_parser     import parse_tads
from pytadbit.utils.sqlite_utils     import already_run, digest_parameters
from pytadbit.utils.sqlite_utils     import add_path, get_jobid, print_db, retry
from pytadbit.utils.extraviews       import tadbit_savefig, nicer
from pytadbit.utils.extraviews       import plot_HiC_matrix


DESC = 'bin Hi-C data into matrices'

def run(opts):
    check_options(opts)
    launch_time = time.localtime()
    param_hash = digest_parameters(opts, extra=['quiet'])
    biases = None

    if opts.zrange:
        vmin = float(opts.zrange.split(',')[0])
        vmax = float(opts.zrange.split(',')[1])
    else:
        vmin = vmax = None

    if opts.figsize:
        opts.figsize = list(map(float, opts.figsize.split(',')))

    clean = True  # change for debug
    if opts.bam:
        mreads = path.realpath(opts.bam)
        if not opts.biases and all(v != 'raw' for v in opts.normalizations):
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

    if opts.plot and not opts.force_plot:
        if opts.interactive:
            max_size = 3500**2
        else:
            max_size = 5000**2
    else:
        max_size = None

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

    out_files = {}
    out_plots = {}

    if opts.matrix or opts.plot:
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
        total = 0
        section_pos = OrderedDict()
        for crm in sections:
            section_pos[crm] = (total, total + sections[crm])
            total += sections[crm]
        for norm in opts.normalizations:
            norm_string = ('RAW' if norm == 'raw' else 'NRM'
                           if norm == 'norm' else 'DEC')
            printime('Getting %s matrices' % norm)
            try:
                matrix, bads1, bads2, regions, name, bin_coords = get_matrix(
                    mreads, opts.reso,
                    load(open(biases, 'rb')) if biases and norm != 'raw' else None,
                    normalization=norm, filter_exclude=opts.filter,
                    region1=region1, start1=start1, end1=end1,
                    region2=region2, start2=start2, end2=end2,
                    tmpdir=tmpdir, ncpus=opts.cpus,
                    return_headers=True,
                    nchunks=opts.nchunks, verbose=not opts.quiet,
                    clean=clean, max_size=max_size,
                    chr_order=opts.chr_name)
            except NotImplementedError:
                if norm == "raw&decay":
                    warn('WARNING: raw&decay normalization not implemented '
                         'for matrices\n... skipping\n')
                    continue
                raise

            b1, e1, b2, e2 = bin_coords
            b1, e1 = 0, e1 - b1
            b2, e2 = 0, e2 - b2

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

            if opts.matrix:
                printime(' - Writing: %s' % norm)
                fnam = '%s_%s_%s%s.mat' % (norm, name,
                                           nicer(opts.reso, sep=''),
                                           ('_' + param_hash))
                out_files[norm_string] = path.join(outdir, fnam)
                out = open(path.join(outdir, fnam), 'w')
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
            if opts.plot:
                # transform matrix
                matrix = array([array([matrix.get((i, j), 0)
                                       for i in range(b1, e1)])
                                for j in range(b2, e2)])
                m = zeros_like(matrix)
                for bad1 in bads1:
                    m[:,bad1] = 1
                    for bad2 in bads2:
                        m[bad2,:] = 1
                matrix = ma.masked_array(matrix, m)
                printime(' - Plotting: %s' % norm)
                fnam = '%s_%s_%s%s%s.%s' % (
                    'nrm' if norm == 'norm' else norm[:3], name,
                    nicer(opts.reso, sep=''),
                    ('_' + param_hash), '_tri' if opts.triangular else '',
                    opts.format)
                out_plots[norm_string] = path.join(outdir, fnam)
                pltbeg1 = 0 if start1 is None else start1
                pltend1 = sections[regions[0]] if end1 is None else end1
                pltbeg2 = 0 if start2 is None else start2
                pltend2 = sections[regions[-1]] if end2 is None else end2
                xlabel = '{}:{:,}-{:,}'.format(
                    regions[0], pltbeg1 if pltbeg1 else 1, pltend1)
                ylabel = '{}:{:,}-{:,}'.format(
                    regions[-1], pltbeg2 if pltbeg2 else 1, pltend2)
                section_pos = OrderedDict((k, section_pos[k]) for k in section_pos
                                   if k in regions)
                transform = (log2 if opts.transform == 'log2' else
                             log if opts.transform == 'log' else lambda x: x)
                tads=None
                if opts.tad_def and not region2:
                    tads, _ = parse_tads(opts.tad_def)
                    if start1:
                        tads = dict([(t, tads[t]) for t in tads
                             if  (int(tads[t]['start']) >= start1 // opts.reso
                                   and int(tads[t]['end']) <= end1 // opts.reso)])
                        for tad in tads:
                            tads[tad]['start'] -= start1 // opts.reso
                            tads[tad]['end'] -= start1 // opts.reso
                ax1, _ = plot_HiC_matrix(
                    matrix, triangular=opts.triangular,
                    vmin=vmin, vmax=vmax, cmap=opts.cmap,
                    figsize=opts.figsize, transform=transform,
                    bad_color=opts.bad_color if norm != 'raw' else None,
                    tad_def=tads)
                ax1.set_title('Region: %s, normalization: %s, resolution: %s' % (
                    name, norm, nicer(opts.reso)), y=1.05)
                _format_axes(ax1, start1, end1, start2, end2, opts.reso,
                             regions, section_pos, sections,
                             opts.xtick_rotation, triangular=False)
                if opts.interactive:
                    plt.show()
                    plt.close('all')
                else:
                    tadbit_savefig(path.join(outdir, fnam))
    if not opts.matrix and not opts.only_plot:
        printime('Getting and writing matrices')
        out_files.update(write_matrix(
            mreads, opts.reso,
            load(open(biases, 'rb')) if biases else None,
            outdir, filter_exclude=opts.filter,
            normalizations=opts.normalizations,
            region1=region1, start1=start1, end1=end1,
            region2=region2, start2=start2, end2=end2,
            tmpdir=tmpdir, append_to_tar=None, ncpus=opts.cpus,
            nchunks=opts.nchunks, verbose=not opts.quiet,
            extra=param_hash, cooler=opts.cooler, clean=clean,
            chr_order=opts.chr_name))

    if clean:
        printime('Cleaning')
        system('rm -rf %s '% tmpdir)

    if not opts.interactive:
        printime('Saving to DB')
        finish_time = time.localtime()
        save_to_db(opts, launch_time, finish_time, out_files, out_plots)


def _format_axes(axe1, start1, end1, start2, end2, reso, regions,
                 section_pos, sections, xtick_rotation, triangular=False):
    if len(regions) <= 2:
        pltbeg1 = 0 if start1 is None else start1
        pltend1 = sections[regions[0]] if end1 is None else end1
        pltbeg2 = (pltbeg1 if len(regions) == 1 else
                   0 if start2 is None else start2)
        pltend2 = (pltend1 if len(regions) == 1 else
                   sections[regions[-1]] if end2 is None else end2)
        axe1.set_xlabel('{}:{:,}-{:,}'.format(
            regions[0] , pltbeg1 if pltbeg1 else 1, pltend1))
        if not triangular:
            axe1.set_ylabel('{}:{:,}-{:,}'.format(
                regions[-1], pltbeg2 if pltbeg2 else 1, pltend2))

        def format_xticks(tickstring, _=None):
            tickstring = int(tickstring * reso + pltbeg1)
            return nicer(tickstring if tickstring else 1,
                         comma=',', allowed_decimals=1)

        def format_yticks(tickstring, _=None):
            tickstring = int(tickstring * reso + pltbeg2)
            return nicer(tickstring if tickstring else 1,
                         comma=',', allowed_decimals=1)

        axe1.xaxis.set_major_formatter(FuncFormatter(format_xticks))
        axe1.yaxis.set_major_formatter(FuncFormatter(format_yticks))
        if triangular:
            axe1.set_yticks([])

        labels = axe1.get_xticklabels()
        plt.setp(labels, rotation=xtick_rotation,
                 ha='left' if xtick_rotation else 'center')
    else:
        vals = [0]
        keys = []
        total = 0
        for crm in section_pos:
            total += (section_pos[crm][1]-section_pos[crm][0]) // reso + 1
            vals.append(total)
            keys.append(crm)
        axe1.set_yticks(vals)
        axe1.set_yticklabels('')
        axe1.set_yticks([float(vals[i]+vals[i + 1]) / 2
                         for i in range(len(vals) - 1)],
                         minor=True)
        axe1.set_yticklabels(keys, minor=True)
        for t in axe1.yaxis.get_minor_ticks():
            t.tick1line.set_visible(False)
            t.tick2line.set_visible(False)

        axe1.set_xticks(vals)
        axe1.set_xticklabels('')
        axe1.set_xticks([float(vals[i]+vals[i+1])/2
                         for i in range(len(vals) - 1)],
                        minor=True)
        axe1.set_xticklabels(keys, minor=True)
        for t in axe1.xaxis.get_minor_ticks():
            t.tick1line.set_visible(False)
            t.tick2line.set_visible(False)
        axe1.set_xlabel('Chromosomes')
        if not triangular:
            axe1.set_ylabel('Chromosomes')


@retry(lite.OperationalError, tries=20, delay=2)
def save_to_db(opts, launch_time, finish_time, out_files, out_plots):
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
        try:
            parameters = digest_parameters(opts, get_md5=False, extra=['quiet'])
            param_hash = digest_parameters(opts, get_md5=True , extra=['quiet'])
            cur.execute("""
            insert into JOBs
            (Id  , Parameters, Launch_time, Finish_time, Type , Parameters_md5)
            values
            (NULL,       '%s',        '%s',        '%s', 'Bin',           '%s')
            """ % (parameters,
                   time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
                   time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass
        except lite.OperationalError:
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
            cur.execute("""
            insert into JOBs
            (Id  , Parameters, Launch_time, Finish_time, Type , Parameters_md5)
            values
            (NULL,       '%s',        '%s',        '%s', 'Bin',           '%s')
            """ % (parameters,
                   time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
                   time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        jobid = get_jobid(cur)
        for fnam in out_files:
            add_path(cur, out_files[fnam], fnam + '_MATRIX', jobid, opts.workdir)
        for fnam in out_plots:
            add_path(cur, out_plots[fnam], fnam + '_FIGURE', jobid, opts.workdir)
        if not opts.quiet:
            print_db(cur, 'JOBs')
            print_db(cur, 'PATHs')
    if 'tmpdb' in opts and opts.tmpdb:
        # copy back file
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)
    # release lock
    try:
        remove(path.join(opts.workdir, '__lock_db'))
    except OSError:
        pass


def check_options(opts):
    mkdir(opts.workdir)

    # transform filtering reads option
    opts.filter = filters_to_bin(opts.filter)

    # enlighten plotting parameter writing
    if opts.only_plot:
        opts.plot = True
    if opts.interactive:
        if opts.nox:
            raise Exception('ERROR: no screen no fun.\n'
                            'Interactive plot incompatible with noX option.')
        opts.plot = True
        opts.only_plot = True

    # check resume
    if not path.exists(opts.workdir):
        raise IOError('ERROR: workdir not found.')

    # check resume
    if opts.triangular and opts.coord2:
        raise NotImplementedError('ERROR: triangular is only available for '
                                  'symmetric matrices.')

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

    # check if job already run using md5 digestion of parameters
    try:
        if already_run(opts):
            if not opts.force:
                if 'tmpdb' in opts and opts.tmpdb:
                    remove(path.join(dbdir, dbfile))
                    exit('WARNING: exact same job already computed, see JOBs table above')
            else:
                warn('WARNING: exact same job already computed, overwriting...')
    except IOError:
        warn((""
              "\nWARNING:\n  new working directory created. It's ok... "
              "but next time use TADbit since the beginning!! :)"))


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
    pltopt = parser.add_argument_group('Plotting options')

    oblopt.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    glopts.add_argument('--noX', help='no display server (X screen)',
                        dest='nox', action='store_true')

    oblopt.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int, required=True,
                        help='''resolution at which to output matrices''')

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

    outopt.add_argument('--matrix', dest='matrix', action='store_true',
                        default=False,
                        help='''Write text matrix in multiple columns (square).
                        By defaults matrices are written in BED-like format (also
                        only way to get a raw matrix with all values including
                        the ones in masked columns).''')

    outopt.add_argument('--cooler', dest='cooler', action='store_true',
                        default=False,
                        help='''Write i,j,v matrix in cooler format instead of text.
                        ''')

    outopt.add_argument('--rownames', dest='row_names', action='store_true',
                        default=False,
                        help='''To store row names in the output text matrix.
                        WARNING: when non-matrix, results in two extra columns''')

    pltopt.add_argument('--plot', dest='plot', action='store_true',
                        default=False,
                        help='Plot matrix in desired format.')

    pltopt.add_argument('--force_plot', dest='force_plot', action='store_true',
                        default=False,
                        help=('Force plotting even with demoniacally big '
                              'matrices (more than 5000x5000, or 1500x1500'
                              'with interactive option).'))

    outopt.add_argument('--only_plot', dest='only_plot', action='store_true',
                        default=False,
                        help='[%(default)s] Skip writing matrix in text format.')

    outopt.add_argument('-i', '--interactive', dest='interactive', action='store_true',
                        default=False,
                        help='''[%(default)s] Open matplotlib interactive plot
                        (nothing will be saved).''')

    pltopt.add_argument('--triangular', dest='triangular', action='store_true',
                        default=False,
                        help='''[%(default)s] represents only half matrix. Note
                        that this also results in truly vectorial images of
                        matrix.''')

    pltopt.add_argument('--xtick_rotation', dest='xtick_rotation',
                        default=-25, type=int,
                        help='''[%(default)s] x-tick rotation''')

    pltopt.add_argument('--cmap', dest='cmap', action='store',
                        default='viridis',
                        help='[%(default)s] Matplotlib color map to use.')

    pltopt.add_argument('--bad_color', dest='bad_color', action='store',
                        default='white',
                        help='''[%(default)s] Matplotlib color to use on bins
                        filtered out (only used with normalized matrices, not
                        raw).''')

    pltopt.add_argument('--format', dest='format', action='store',
                        default='png',
                        help='[%(default)s] plot file format.')

    pltopt.add_argument('--zrange', dest='zrange', action='store',
                        default=None,
                        help='''Range, in log2 scale of the color scale.
                        i.e.: --zrange=-2,2''')

    pltopt.add_argument('--transform', dest='transform', action='store',
                        default='log2', choices=['log2', 'log', 'none'],
                        help='''[%(default)s] can be any of [%(choices)s]''')

    pltopt.add_argument('--figsize', dest='figsize', action='store',
                        default=None,
                        help='''Range, in log2 scale of the color scale.
                        default for triangular matrices: --figsize=16,10
                        and for square matrices:  --figsize=16,14''')

    pltopt.add_argument('--tad_def', dest='tad_def', action='store',
                        default=None,
                        help='''tsv file with tad definition, columns:
                        #    start    end    score    density''')

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

    normpt.add_argument('--norm', dest='normalizations', metavar="STR",
                        action='store', default=['raw'], type=str, nargs='+',
                        choices=['norm', 'decay', 'raw', 'raw&decay'],
                        help='''[%(default)s] normalization(s) to apply.
                        Choices are: [%(choices)s]''')

    rfiltr.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                        choices = list(range(0, 11)),
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: 0: nothing, %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))

    outopt.add_argument('--only_txt', dest='only_txt', action='store_true',
                        default=False,
                        help='Save only text file for matrices, not images')

    # rfiltr.add_argument('--valid', dest='only_valid', action='store_true',
    #                     default=False,
    #                     help='input BAM file contains only valid pairs (already filtered).')


def load_parameters_fromdb(opts):
    if opts.tmpdb:
        dbfile = opts.tmpdb
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        if not opts.jobid:
            # get the JOBid of the parsing job
            try:
                cur.execute("""
                select distinct Id from JOBs
                where Type = '%s'
                """ % ('Normalize' if opts.normalizations != ['raw'] else 'Filter'))
                jobids = cur.fetchall()
                parse_jobid = jobids[0][0]
            except IndexError:
                cur.execute("""
                select distinct Id from JOBs
                where Type = '%s'
                """ % ('Filter'))
                jobids = cur.fetchall()
                try:
                    parse_jobid = jobids[0][0]
                except IndexError:
                    parse_jobid = 1
            if len(jobids) > 1:
                found = False
                if opts.normalizations != ('raw', ):
                    cur.execute("""
                    select distinct JOBid from NORMALIZE_OUTPUTs
                    where Resolution = %d
                    """ % (opts.reso))
                    jobs = cur.fetchall()
                    try:
                        parse_jobid = jobs[0][0]
                        found = True
                    except IndexError:
                        found = False
                    if len(jobs ) > 1:
                        found = False
                if not found:
                    raise Exception('ERROR: more than one possible input found, use'
                                    '"tadbit describe" and select corresponding '
                                    'jobid with --jobid')
        else:
            parse_jobid = opts.jobid
        # fetch path to parsed BED files
        # try:
        biases = mreads = reso = None
        if len(opts.normalizations) > 1 or opts.normalizations[0] != 'raw':
            try:
                cur.execute("""
                select distinct Path from PATHs
                where paths.jobid = %s and paths.Type = 'BIASES'
                """ % parse_jobid)
                biases = cur.fetchall()[0][0]

                cur.execute("""
                select distinct Path from PATHs
                inner join NORMALIZE_OUTPUTs on PATHs.Id = NORMALIZE_OUTPUTs.Input
                where NORMALIZE_OUTPUTs.JOBid = %d;
                """ % parse_jobid)
                mreads = cur.fetchall()[0][0]

                cur.execute("""
                select distinct Resolution from NORMALIZE_OUTPUTs
                where NORMALIZE_OUTPUTs.JOBid = %d;
                """ % parse_jobid)
                reso = int(cur.fetchall()[0][0])
                if reso != opts.reso:
                    warn('WARNING: input resolution does not match '
                         'the one of the precomputed normalization')
            except IndexError:
                warn('WARNING: normalization not found')
                cur.execute("""
                select distinct path from paths
                inner join filter_outputs on filter_outputs.pathid = paths.id
                where filter_outputs.name = 'valid-pairs' and paths.jobid = %s
                """ % parse_jobid)
                try:
                    mreads = cur.fetchall()[0][0]
                except IndexError:
                    raise Exception('ERROR: missing normalization input file.')
        else:
            cur.execute("""
            select distinct path from paths
            inner join filter_outputs on paths.type = 'HIC_BAM'
            where filter_outputs.name = 'valid-pairs' and paths.jobid = %s
            """ % parse_jobid)
            fetched = cur.fetchall()
            if len(fetched) > 1:
                raise Exception('ERROR: more than one item in the database')
            mreads = fetched[0][0]
        return biases, mreads
