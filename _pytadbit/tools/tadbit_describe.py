"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""

from argparse                    import HelpFormatter
from pytadbit.utils.sqlite_utils import print_db
import sqlite3 as lite
from os                          import path, remove
from string                      import ascii_letters
from random                      import random
from shutil                      import copyfile

DESC = "Describe jobs and results in a given working directory"

def run(opts):
    check_options(opts)
    if 'tmp' in opts and opts.tmp:
        dbfile = opts.tmp
        copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    with con:
        cur = con.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        for table in cur.fetchall():
            if table[0].lower() in opts.tables:
                print_db(cur, table[0])
    if 'tmp' in opts and opts.tmp:
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

    glopts.add_argument('-t', '--table', dest='tables', metavar='',
                        action='store', nargs='+', type=str,
                        default=[str(t) for t in range(1, 10)],
                        help='''[%(default)s] what tables to show, wrte either the sequence of
                        names or indexes, according to this list:
                        1: paths, 2: jobs, 3: mapped_outputs,
                        4: mapped_inputs, 5: parsed_outputs,
                        6: intersection_outputs, 7: filter_outputs,
                        8: normalize_outputs, 9: segment_outputs''')

    glopts.add_argument('--tmp', dest='tmp', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    parser.add_argument_group(glopts)

def check_options(opts):
    if not opts.workdir: raise Exception('ERROR: output option required.')

    choices=['1', 'paths', '2', 'jobs',
             '3', 'mapped_outputs',
             '4', 'mapped_inputs', '5', 'parsed_outputs',
             '6', 'intersection_outputs',
             '7', 'filter_outputs', '8', 'normalize_outputs',
             '9', 'segment_outputs']
    table_idx = {
        '1': 'paths',
        '2': 'jobs',
        '3': 'mapped_outputs',
        '4': 'mapped_inputs',
        '5': 'parsed_outputs',
        '6': 'intersection_outputs',
        '7': 'filter_outputs',
        '8': 'normalize_outputs',
        '9': 'segment_outputs'}
    recovered = []
    bads = []
    for t in range(len(opts.tables)):
        opts.tables[t] = opts.tables[t].lower()
        if not opts.tables[t] in choices:
            # check if the begining of the input string matches any of
            # the possible choices
            found = False
            for choice in table_idx.values():
                if choice.startswith(opts.tables[t]):
                    recovered.append(choice)
                    found = True
            if not found:
                print(('error: argument -t/--table: invalid choice: %s'
                       '(choose from %s )') % (opts.tables[t], str(choices)))
                exit()
            bads.append(t)
        opts.tables[t] = table_idx.get(opts.tables[t], opts.tables[t])
    for bad in bads[::-1]:
        del(opts.tables[bad])
    for rec in recovered:
        opts.tables.append(rec)

    if 'tmp' in opts and opts.tmp:
        dbdir = opts.tmp
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmp = path.join(dbdir, dbfile)
