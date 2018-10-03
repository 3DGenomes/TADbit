"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""

from os                          import path, system, remove
from string                      import ascii_letters
from random                      import random
from shutil                      import copyfile
from argparse                    import HelpFormatter
import sqlite3 as lite

from pytadbit.utils.sqlite_utils import delete_entries
from pytadbit.utils.sqlite_utils import update_wordir_path

DESC = "Delete jobs and results of a given list of jobids in a given directories"

def run(opts):
    # print summary of what will be removed
    # prmpt if sure?
    check_options(opts)
    if 'tmpdb' in opts and opts.tmpdb:
        dbfile = opts.tmpdb
        copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()

        # change working directory
        if opts.new_workdir:
            new_workdir = path.abspath(opts.new_workdir)
            update_wordir_path(cur, new_workdir)
            print 'Updated working directory to %s' % (new_workdir)
            return

        # get PATHids corresponding to JOBid:
        paths = []
        protected_types = ['INDEX', 'FASTA', 'MAPPED_FASTQ', 'WORKDIR']
        for jobid in opts.jobids:
            cur.execute("SELECT Id, Path, Type FROM PATHs WHERE JOBid=%s" % jobid)
            paths.extend([p for p in cur.fetchall()])

        # delete files and directories
        if opts.delete:
            print 'deleting %d files' % len(paths)
            for _, lpath, typ in paths:
                if typ in protected_types or typ.startswith('EXT_'):
                    continue
                # we do not remove anything outside the work directory
                if not path.realpath(path.join(opts.workdir, lpath)).startswith(
                        path.realpath(opts.workdir)):
                    continue
                print '  x ' + path.join(opts.workdir, lpath)
                system('rm -rf ' + path.join(opts.workdir, lpath))

        # remove entry
        cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        for table in [t[0] for t in cur.fetchall()]:
            print 'cleaning table: %s' % table
            for jobid in opts.jobids:
                if table.lower() == 'jobs':
                    col = 'Id'
                else:
                    col = 'JOBid'
                delete_entries(cur, table, col, jobid)

            for cpath, _, _ in paths:
                if table.lower() == 'mapped_outputs':
                    elt = 'BEDid'
                else:
                    elt = 'PATHid'
                for jobid in opts.jobids:
                    delete_entries(cur, table, elt, cpath)
            # remove empty tables:
            cur.execute("select count(*) from %s" % table)
            count = cur.fetchall()[0][0]
            if count == 0:
                cur.execute("drop table %s" % table)
                print ' * dropped table %s' % table
    if 'tmpdb' in opts and opts.tmpdb:
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)


def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class=lambda prog: HelpFormatter(prog, width=95,
                                                      max_help_position=27)

    glopts = parser.add_argument_group('General options')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to working directory (generated with the
                        tool tadbit mapper)''')

    glopts.add_argument('-j', '--jobids', dest='jobids', metavar="INT",
                        action='store', default=None, nargs='+', type=int,
                        help='jobids of the files and entries to be removed')

    glopts.add_argument('--delete', dest='delete', action="store_true",
                        default=False,
                        help='delete files, otherwise only DB entries.')

    glopts.add_argument('--compress', dest='compress', action="store_true",
                        default=False,
                        help='compress files and update paths accordingly')

    glopts.add_argument('--noX', action='store_true', help='no display server (X screen)')

    parser.add_argument('--change_workdir', dest='new_workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''In case folder was moved, input the new path''')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    parser.add_argument_group(glopts)

def check_options(opts):
    if not opts.workdir: raise Exception('ERROR: output option required.')

    dbpath = path.join(opts.workdir, 'trace.db')
    if not path.exists(dbpath):
        raise Exception('ERROR: DB file: %s not found.' % dbpath)

    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)
