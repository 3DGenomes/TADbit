"""

information needed

 - path working directory with parsed reads

"""
from argparse                     import HelpFormatter
from pytadbit                     import load_hic_data_from_reads
from pytadbit.mapping.analyze     import correlate_matrices
from pytadbit.mapping             import merge_2d_beds
from pytadbit.mapping.analyze     import eig_correlate_matrices
from pytadbit.utils.sqlite_utils  import already_run, digest_parameters
from pytadbit.utils.sqlite_utils  import add_path, get_jobid, print_db
from pytadbit.utils.sqlite_utils  import get_path_id
from pytadbit.utils.file_handling import mkdir
from os                           import path, remove
from string                       import ascii_letters
from random                       import random
from shutil                       import copyfile
import sqlite3 as lite
import time

DESC = ('load two working directories with different Hi-C data samples and ' +
        'merges them into a new working directory generating some statistics')

def run(opts):
    check_options(opts)
    launch_time = time.localtime()

    param_hash = digest_parameters(opts)

    if opts.bed1:
        mreads1 = path.realpath(opts.bed1)
    else:
        mreads1 = path.join(opts.workdir1, load_parameters_fromdb(
            opts.workdir1, opts.jobid1, opts))

    if opts.bed2:
        mreads2 = path.realpath(opts.bed2)
    else:
        mreads2 = path.join(opts.workdir2, load_parameters_fromdb(
            opts.workdir2, opts.jobid2, opts))

    print 'loading first sample', mreads1
    hic_data1 = load_hic_data_from_reads(mreads1, opts.reso)

    print 'loading second sample', mreads2
    hic_data2 = load_hic_data_from_reads(mreads2, opts.reso)

    mkdir(path.join(opts.workdir, '00_merge'))

    decay_corr_dat = path.join(opts.workdir, '00_merge', 'decay_corr_dat_%s_%s.txt' % (opts.reso, param_hash))
    decay_corr_fig = path.join(opts.workdir, '00_merge', 'decay_corr_dat_%s_%s.png' % (opts.reso, param_hash))
    eigen_corr_dat = path.join(opts.workdir, '00_merge', 'eigen_corr_dat_%s_%s.txt' % (opts.reso, param_hash))
    eigen_corr_fig = path.join(opts.workdir, '00_merge', 'eigen_corr_dat_%s_%s.png' % (opts.reso, param_hash))

    print 'correlation between equidistant loci'
    corr, dist, bads = correlate_matrices(hic_data1, hic_data2, normalized=False,
                                          remove_bad_columns=True,
                                          savefig=decay_corr_fig,
                                          savedata=decay_corr_dat, get_bads=True)

    print 'correlation between eigenvectors'
    eig_corr = eig_correlate_matrices(hic_data1, hic_data2, normalized=False,
                                      remove_bad_columns=True,
                                      savefig=eigen_corr_fig,
                                      savedata=eigen_corr_dat)

    # merge inputs
    mkdir(path.join(opts.workdir, '03_filtered_reads'))
    outbed = path.join(opts.workdir, '03_filtered_reads', 'valid_r1-r2_intersection_%s.tsv' % (
        param_hash))

    nreads = merge_2d_beds(mreads1, mreads2, outbed)

    finish_time = time.localtime()
    save_to_db (opts, mreads1, mreads2, decay_corr_dat, decay_corr_fig,
                len(bads.keys()), len(hic_data1), nreads,
                eigen_corr_dat, eigen_corr_fig, outbed, 
                launch_time, finish_time)

def save_to_db(opts, mreads1, mreads2, decay_corr_dat, decay_corr_fig,
               nbad_columns, ncolumns, nreads,
               eigen_corr_dat, eigen_corr_fig, outbed,
               launch_time, finish_time):
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
                       type='table' AND name='NORMALIZE_OUTPUTs'""")
        if not cur.fetchall():
            cur.execute("""
            create table PATHs
               (Id integer primary key,
                JOBid int, Path text, Type text,
                unique (Path))""")
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
            create table FILTER_OUTPUTs
               (Id integer primary key,
                PATHid int,
                Name text,
                Count int,
                JOBid int,
                unique (PATHid))""")
            cur.execute("""
            create table MERGE_OUTPUTs
               (Id integer primary key,
                JOBid int,
                Wrkd1_PathId int,
                Wrkd2_PathId int,
                Bed1_PathId int,
                Bed2_PathId int,
                Merge_PathId int,
                decay_corr text,
                eigen_corr text,
                N_columns int,
                N_filtered int,
                Resolution int,
                unique (JOBid))""")
        try:
            parameters = digest_parameters(opts, get_md5=False)
            param_hash = digest_parameters(opts, get_md5=True )
            cur.execute("""
            insert into JOBs
            (Id  , Parameters, Launch_time, Finish_time, Type   , Parameters_md5)
            values
            (NULL,       '%s',        '%s',        '%s', 'Merge',           '%s')
            """ % (parameters,
                   time.strftime("%d/%m/%Y %H:%M:%S", launch_time),
                   time.strftime("%d/%m/%Y %H:%M:%S", finish_time), param_hash))
        except lite.IntegrityError:
            pass

        jobid = get_jobid(cur)
        add_path(cur, decay_corr_dat, 'CORR'   , jobid, opts.workdir)
        add_path(cur, decay_corr_fig, 'FIGURE' , jobid, opts.workdir)
        add_path(cur, eigen_corr_dat, 'CORR'   , jobid, opts.workdir)
        add_path(cur, eigen_corr_fig, 'FIGURE' , jobid, opts.workdir)

        add_path(cur, opts.workdir , 'WORKDIR' , jobid)
        add_path(cur, opts.workdir1, 'WORKDIR1', jobid, opts.workdir)
        add_path(cur, opts.workdir2, 'WORKDIR2', jobid, opts.workdir)
        add_path(cur, mreads1      , '2D_BED'  , jobid, opts.workdir)
        add_path(cur, mreads2      , '2D_BED'  , jobid, opts.workdir)
        add_path(cur, outbed       , '2D_BED'  , jobid, opts.workdir)
        
        cur.execute("select id from paths where path = '%s'" % (
            path.relpath(mreads1, opts.workdir)))
        bed1 = cur.fetchall()[0][0]
        cur.execute("select id from paths where path = '%s'" % (
            path.relpath(opts.workdir1, opts.workdir)))
        w1path = cur.fetchall()[0][0]
        cur.execute("select id from paths where path = '%s'" % (
            path.relpath(mreads2, opts.workdir)))
        bed2 = cur.fetchall()[0][0]
        cur.execute("select id from paths where path = '%s'" % (
            path.relpath(opts.workdir2, opts.workdir)))
        w2path = cur.fetchall()[0][0]
        cur.execute("select id from paths where path = '%s'" % (
            path.relpath(outbed, opts.workdir)))
        outbedid = cur.fetchall()[0][0]
        cur.execute("""
        insert into MERGE_OUTPUTs
        (Id  , JOBid, Wrkd1_PathId, Wrkd2_PathId, Bed1_PathId, Bed2_PathId, Merge_PathId, N_columns,   N_filtered, Resolution)
        values
        (NULL,    %d,           %d,           %d,          %d,          %d,           %d,        %d,           %d,         %d)
        """ % (jobid,          w1path,          w2path,        bed1,        bed2,      outbedid,  ncolumns, nbad_columns, opts.reso ))
            
        masked1 = {}
        if opts.workdir1:
            tmpcon = lite.connect(path.join(opts.workdir1, 'trace.db'))
            with tmpcon:
                tmpcur = tmpcon.cursor()
                tmpcur.execute("select Name, PATHid, Count from filter_outputs")
                for name, pathid, count in tmpcur.fetchall():
                    res = tmpcur.execute("select Path from PATHs where Id = %d" % (pathid))
                    tmppath = res.fetchall()[0][0]
                    masked1[name] = {'path': tmppath, 'count': count}

        masked2 = {}
        if opts.workdir2:
            tmpcon = lite.connect(path.join(opts.workdir2, 'trace.db'))
            with tmpcon:
                tmpcur = tmpcon.cursor()
                tmpcur.execute("select Name, PATHid, Count from filter_outputs")
                for name, pathid, count in tmpcur.fetchall():
                    res = tmpcur.execute("select Path from PATHs where Id = %d" % (pathid))
                    tmppath = res.fetchall()[0][0]
                    masked2[name] = {'path': tmppath, 'count': count}

        for f in masked1:
            if f  != 'valid-pairs':
                outmask = path.join(opts.workdir, '03_filtered_reads',
                                    'all_r1-r2_intersection_%s.tsv_%s.tsv' % (
                    param_hash, f))
                out = open(outmask, 'w')
                for line in open(path.join(opts.workdir1, masked1[f]['path'])):
                    out.write(line)
                for line in open(path.join(opts.workdir2, masked2[f]['path'])):
                    out.write(line)
                add_path(cur, outmask, 'FILTER', jobid, opts.workdir)
            else:
                outmask = outbed

            cur.execute("""
            insert into FILTER_OUTPUTs
            (Id  , PATHid, Name, Count, JOBid)
            values
            (NULL,     %d, '%s',  '%s',    %d)
            """ % (get_path_id(cur, outmask, opts.workdir),
                   f, masked1[f]['count'] + masked2[f]['count'], jobid))

        print_db(cur, 'PATHs')
        print_db(cur, 'JOBs')
        print_db(cur, 'MERGE_OUTPUTs')        
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

def load_parameters_fromdb(workdir, jobid, opts):
    if 'tmpdb' in opts and opts.tmpdb:
        dbfile = opts.tmpdb
    else:
        dbfile = path.join(workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        if not jobid:
            # get the JOBid of the parsing job
            cur.execute("""
            select distinct Id from JOBs
            where Type = 'Filter'
            """)
            jobids = cur.fetchall()
            if len(jobids) > 1:
                raise Exception('ERROR: more than one possible input found, use'
                                '"tadbit describe" and select corresponding '
                                'jobid with --jobid')
            parse_jobid = jobids[0][0]
        else:
            parse_jobid = opts.jobid
        # fetch path to parsed BED files
        cur.execute("""
        select distinct path from paths
        inner join filter_outputs on filter_outputs.pathid = paths.id
        where filter_outputs.name = 'valid-pairs' and paths.jobid = %s
        """ % parse_jobid)
        return cur.fetchall()[0][0]

def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to a new output folder''')

    glopts.add_argument('-w1', '--workdir1', dest='workdir1', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory of the first HiC data
                        sample to merge''')


    glopts.add_argument('-w2', '--workdir2', dest='workdir2', metavar="PATH",
                        action='store', default=None, type=str, required=True,
                        help='''path to working directory of the second HiC data
                        sample to merge''')

    glopts.add_argument('--bed1', dest='bed1', metavar="PATH",
                        action='store', default=None, type=str, 
                        help='''path to the first TADbit-generated BED file with
                        filtered reads (other wise the tool will guess from the
                        working directory database)''')

    glopts.add_argument('--bed2', dest='bed2', metavar="PATH",
                        action='store', default=None, type=str, 
                        help='''path to the second TADbit-generated BED file with
                        filtered reads (other wise the tool will guess from the
                        working directory database)''')

    glopts.add_argument('-r', '--resolution', dest='reso', metavar="INT",
                        action='store', default=None, type=int, required=True,
                        help='''resolution at which to do the comparison,
                        and generate the matrices.''')

    glopts.add_argument('--perc_zeros', dest='perc_zeros', metavar="FLOAT",
                        action='store', default=95, type=float, 
                        help=('[%(default)s%%] maximum percentage of zeroes '
                              'allowed per column.'))

    glopts.add_argument('--normalization', dest='resolution', metavar="STR",
                        action='store', default='ICE', nargs='+', type=str,
                        choices=['ICE', 'EXP'],
                        help='''[%(default)s] normalization(s) to apply.
                        Order matters.''')

    glopts.add_argument('--save', dest='save', metavar="STR",
                        action='store', default='genome', nargs='+', type=str,
                        choices=['genome', 'chromosomes'],
                        help='''[%(default)s] save genomic or chromosomic matrix.''')

    glopts.add_argument('--jobid1', dest='jobid1', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')    

    glopts.add_argument('--jobid2', dest='jobid2', metavar="INT",
                        action='store', default=None, type=int,
                        help='''Use as input data generated by a job with a given
                        jobid. Use tadbit describe to find out which.''')    

    glopts.add_argument('--force', dest='force', action='store_true',
                      default=False,
                      help='overwrite previously run job')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    parser.add_argument_group(glopts)


def check_options(opts):
    mkdir(opts.workdir)

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

    # check if job already run using md5 digestion of parameters
    if already_run(opts):
        if 'tmpdb' in opts and opts.tmpdb:
            remove(path.join(dbdir, dbfile))
        exit('WARNING: exact same job already computed, see JOBs table above')

def nice(reso):
    if reso >= 1000000:
        return '%dMb' % (reso / 1000000)
    return '%dkb' % (reso / 1000)


