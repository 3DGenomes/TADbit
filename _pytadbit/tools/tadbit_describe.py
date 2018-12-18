"""

information needed

 - path working directory with mapped reads or list of SAM/BAM/MAP files

"""

from argparse                    import HelpFormatter
import sqlite3 as lite
from os                          import path, remove
from string                      import ascii_letters
from random                      import random
from shutil                      import copyfile

from pytadbit.utils.sqlite_utils import print_db

DESC = "Describe jobs and results in a given working directory"

TABLE_IDX = {
    '1' : 'paths',
    '2' : 'jobs',
    '3' : 'mapped_outputs',
    '4' : 'mapped_inputs',
    '5' : 'parsed_outputs',
    '6' : 'intersection_outputs',
    '7' : 'filter_outputs',
    '8' : 'normalize_outputs',
    '9' : 'merge_stats',
    '10': 'merge_outputs',
    '11': 'segment_outputs',
    '12': 'models',
    '13': 'modeled_regions'}


def run(opts):
    check_options(opts)
    if 'tmpdb' in opts and opts.tmpdb:
        dbfile = opts.tmpdb
        copyfile(path.join(opts.workdir, 'trace.db'), dbfile)
    else:
        dbfile = path.join(opts.workdir, 'trace.db')
    con = lite.connect(dbfile)
    if opts.output and path.exists(opts.output):
        remove(opts.output)
    with con:
        cur = con.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        for table in cur.fetchall():
            if table[0].lower() in ['jobs', 'paths'] and opts.tsv:
                continue
            if table[0].lower() in opts.tables:
                print_db(cur, table[0],
                         savedata=opts.output,
                         append=True, tsv=opts.tsv,
                         no_print=['JOBid', 'Id', 'Input',
                                   '' if table[0] == 'MAPPED_OUTPUTs'
                                   else 'PATHid'] if opts.tsv else '',
                         jobids=opts.jobids, columns=opts.select,
                         **dict(f.split(',') for f in opts.where))
    if 'tmpdb' in opts and opts.tmpdb:
        copyfile(dbfile, path.join(opts.workdir, 'trace.db'))
        remove(dbfile)


def populate_args(parser):
    """
    parse option from call
    """
    parser.formatter_class = lambda prog: HelpFormatter(prog, width=95,
                                                        max_help_position=27)

    glopts = parser.add_argument_group('General options')

    glopts.add_argument('-w', '--workdir', dest='workdir', metavar="PATH",
                        action='store', default=None, type=str,
                        help='''path to working directory (generated with the
                        tool tadbit map)''')

    glopts.add_argument('--noX', action='store_true', help='no display server (X screen)')

    glopts.add_argument('-t', '--tables', dest='tables', metavar='',
                        action='store', nargs='+', type=str,
                        default=[str(t) for t in range(1, len(TABLE_IDX) + 1)],
                        help=('[%(default)s] what tables to show, write either '
                              'the sequence of names or indexes, according to '
                              'this list: {}').format(', '.join(
                                  ['%s: %s' % (k, v)
                                   for k, v in TABLE_IDX.iteritems()])))

    glopts.add_argument('-T', '--skip_tables', dest='skip_tables', metavar='',
                        action='store', nargs='+', type=str,
                        default=[],
                        help=('[%(default)s] what tables NOT to show, write either '
                              'the sequence of names or indexes, according to '
                              'this list: {}').format(', '.join(
                                  ['%s: %s' % (k, v)
                                   for k, v in TABLE_IDX.iteritems()])))

    glopts.add_argument('-j', '--jobids', dest='jobids', metavar="INT",
                        nargs='+', action='store', default=None, type=int,
                        help='''Display only items matching these jobids.''')

    glopts.add_argument('-W', '--where', dest='where', metavar="STR",
                        action='store', default={}, type=str, nargs='+',
                        help='''Select rows. List pairs of keywords (column header)
                        and values to filter results. For example to get only results
                        where "18" appears in the column "Chromosome", the option
                        should be set as: `-W Chromosome,18`''')

    glopts.add_argument('-s', '--select', dest='select', metavar="STR",
                        action='store', default={}, type=str, nargs='+',
                        help='''Select columns. List the keyword (column
                        headers) to be displayed. E.g.
                        to show only the colmns JobIds: `-s Jobids`''')

    glopts.add_argument('--tmpdb', dest='tmpdb', action='store', default=None,
                        metavar='PATH', type=str,
                        help='''if provided uses this directory to manipulate the
                        database''')

    glopts.add_argument('--tsv', dest='tsv', default=False,
                        action='store_true',
                        help='''Print output in tab separated format''')

    glopts.add_argument('-o', '--output', dest='output', default=None,
                        action='store',
                        help='''Writes output in specified file.''')
    parser.add_argument_group(glopts)


def check_options(opts):
    if not opts.workdir:
        raise Exception('ERROR: output option required.')

    choices = reduce(lambda x, y: x + y,
                     [kv for kv in sorted(TABLE_IDX.iteritems(),
                                          key=lambda x: int(x[0]))])

    recovered = []
    bads = []
    for t in range(len(opts.tables)):
        opts.tables[t] = opts.tables[t].lower()
        if not opts.tables[t] in choices:
            # check if the beginning of the input string matches any of
            # the possible choices
            found = False
            for choice in TABLE_IDX.values():
                if choice.startswith(opts.tables[t]):
                    recovered.append(choice)
                    found = True
            if not found:
                print(('ERROR: argument -t/--tables: invalid choice: %s. '
                       'Should be one of:\n  - %s') % (opts.tables[t].strip(),
                                                       '\n  - '.join(choices)))
                exit()
            bads.append(t)
        opts.tables[t] = TABLE_IDX.get(opts.tables[t], opts.tables[t])

    for t in range(len(opts.skip_tables)):
        opts.skip_tables[t] = opts.skip_tables[t].lower()
        if not opts.skip_tables[t] in choices:
            # check if the beginning of the input string matches any of
            # the possible choices
            for choice in TABLE_IDX.values():
                if choice.startswith(opts.skip_tables[t]):
                    break
            else:
                print(('ERROR: argument -T/--skip_tables: invalid choice: %s. '
                       'Should be one of:\n  - %s') % (opts.tables[t].strip(),
                                                       '\n  - '.join(choices)))
                exit()
            opts.skip_tables[t] = choice
        else:
            opts.skip_tables[t] = TABLE_IDX.get(opts.skip_tables[t], opts.skip_tables[t])
    for bad in bads[::-1]:
        del(opts.tables[bad])
    for rec in recovered:
        opts.tables.append(rec)
    for t in opts.tables[:]:
        if t in opts.skip_tables:
            opts.tables.pop(opts.tables.index(t))

    if 'tmpdb' in opts and opts.tmpdb:
        dbdir = opts.tmpdb
        # tmp file
        dbfile = 'trace_%s' % (''.join([ascii_letters[int(random() * 52)]
                                        for _ in range(10)]))
        opts.tmpdb = path.join(dbdir, dbfile)
