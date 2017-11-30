"""

information needed

 - path working directory with parsed reads

"""
from argparse                        import HelpFormatter
from os                              import path, rmdir, remove
from sys                             import stdout
from shutil                          import copyfile
from string                          import ascii_letters
from random                          import random
from cPickle                         import load
from warnings                        import warn
from multiprocessing                 import cpu_count
import time
import sqlite3 as lite

from pytadbit.mapping.filter         import MASKED
from pytadbit.utils.file_handling    import mkdir
from pytadbit.parsers.hic_bam_parser import filters_to_bin, write_matrix
from pytadbit.utils.sqlite_utils     import already_run, digest_parameters
from pytadbit.utils.sqlite_utils     import add_path, get_jobid, print_db


DESC = 'bin Hi-C data into matrices'


def run(opts):
    check_options(opts)
    launch_time = time.localtime()
    param_hash = digest_parameters(opts)

    if opts.bam:
        mreads = path.realpath(opts.bam)
        if not opts.biases and all(v !='raw' for v in opts.normalizations):
            raise Exception('ERROR: external BAM input, should provide path to'
                            ' biases file.')
        biases = opts.biases
    else:
        biases, mreads = load_parameters_fromdb(opts)
        mreads = path.join(opts.workdir, mreads)
        biases = path.join(opts.workdir, biases)
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
            stdout.write('\nExtraction of full genome\n')

    out_files = write_matrix(mreads, opts.reso,
                             load(open(biases)) if biases else None,
                             outdir, filter_exclude=opts.filter,
                             normalizations=opts.normalizations,
                             region1=region1, start1=start1, end1=end1,
                             region2=region2, start2=start2, end2=end2,
                             tmpdir='.', append_to_tar=None, ncpus=opts.cpus,
                             nchunks=opts.nchunks, verbose=not opts.quiet,
                             extra=param_hash, clean=False)

    rmdir(tmpdir)
    finish_time = time.localtime()

    save_to_db(opts, launch_time, finish_time, out_files)


def save_to_db(opts, launch_time, finish_time, out_files):
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
        jobid = get_jobid(cur)
        for fnam in out_files:
            add_path(cur, out_files[fnam], fnam + '_MATRIX', jobid, opts.workdir)
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

    # check resume
    if not path.exists(opts.workdir):
        raise IOError('ERROR: workdir not found.')

    # for lustre file system....
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
    if already_run(opts):
        if 'tmpdb' in opts and opts.tmpdb:
            remove(path.join(dbdir, dbfile))
        exit('WARNING: exact same job already computed, see JOBs table above')


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

    glopts.add_argument('--nchunks', dest='nchunks', action='store', default=None,
                        type=int,
                        help='''maximum number of chunks into which to
                        cut the BAM''')

    glopts.add_argument("-C", "--cpus", dest="cpus", type=int,
                        default=0, help='''[%(default)s] Maximum number of CPU
                        cores  available in the execution host. If higher
                        than 1, tasks with multi-threading
                        capabilities will enabled (if 0 all available)
                        cores will be used''')

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
                        help='''path to file with precalculated biases by
                        columns''')

    normpt.add_argument('--norm', dest='normalizations', metavar="STR",
                        action='store', default=['raw'], type=str, nargs='+',
                        choices=['norm', 'decay', 'raw'],
                        help='''[%(default)s] normalization(s) to apply.
                        Order matters. Choices: [%(choices)s]''')

    outopt.add_argument('--keep', dest='keep', action='store',
                        default=['intra', 'genome'], nargs='+',
                        choices = ['intra', 'inter', 'genome', 'none'],
                        help='''%(default)s Matrices to save, choices are
                        "intra" to keep intra-chromosomal matrices, "inter" to
                        keep inter-chromosomal matrices and "genome", to keep
                        genomic matrices.''')

    rfiltr.add_argument('-F', '--filter', dest='filter', nargs='+',
                        type=int, metavar='INT', default=[1, 2, 3, 4, 6, 7, 9, 10],
                        choices = range(1, 11),
                        help=("""[%(default)s] Use filters to define a set os
                        valid pair of reads e.g.:
                        '--apply 1 2 3 4 8 9 10'. Where these numbers""" +
                              "correspond to: %s" % (', '.join(
                                  ['%2d: %15s' % (k, MASKED[k]['name'])
                                   for k in MASKED]))))

    outopt.add_argument('--only_txt', dest='only_txt', action='store_true',
                        default=False,
                        help='Save only text file for matrices, not images')

    rfiltr.add_argument('--valid', dest='only_valid', action='store_true',
                        default=False,
                        help='input BAM file contains only valid pairs (already filtered).')


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
                """ % ('Normalize' if opts.normalizations != ('raw', ) else 'Filter'))
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
        if opts.normalizations != ('raw', ):
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
                mreads = cur.fetchall()[0][0]
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
