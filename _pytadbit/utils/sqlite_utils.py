"""
some utils relative to sqlite
"""
import sqlite3 as lite
from os.path import abspath, relpath, join
from hashlib import md5


def digest_parameters(opts, get_md5=True, extra=None):
    """
    digestion is truncated at 10 characters. More than a 1000 runs are not
    expected on a same sample (usually ~10).

    :param None extra: extra parameter to remove from digestion
    """
    extra = extra or []
    if get_md5:
        # print 'MD5', ' '.join(
        #     ['%s:%s' % (k, int(v) if isinstance(v, bool) else v)
        #      for k, v in sorted(opts.__dict__.iteritems())
        #      if k not in ['force', 'workdir', 'func', 'tmp',
        #                   'skip', 'keep_tmp', 'tmpdb'] + extra])
        param_hash = md5(' '.join(
            ['%s:%s' % (k, int(v) if isinstance(v, bool) else v)
             for k, v in sorted(opts.__dict__.iteritems())
             if k not in ['force', 'workdir', 'func', 'tmp',
                          'skip', 'keep_tmp', 'tmpdb'] + extra])).hexdigest()[:10]
        return param_hash
    parameters = ' '.join(
        ['%s:%s' % (k, int(v) if isinstance(v, bool) else v)
         for k, v in opts.__dict__.iteritems()
         if k not in ['fastq', 'index', 'renz', 'iterative', 'workdir',
                      'skip', 'func', 'tmp', 'keep_tmp'] + extra and v is not None])
    parameters = parameters.replace("'", "")
    return parameters


def update_wordir_path(cur, new_path):
    cur.execute("update paths set path='%s' where type='WORKDIR'" % (
        new_path))

def delete_entries(cur, table, col, val):
    try:
        cur.execute("select %s from %s where %s.%s = %d" % (
            col, table, table, col, val))
        elts = [e[0] for e in cur.fetchall()]
        if len(elts) == 0:
            return
    except lite.OperationalError:
        return
    try:
        cur.execute("delete from %s where %s.%s = %d" % (
            table, table, col, val))
        print ' - deleted %d elements with %s = %s' % (len(elts), col, val)
    except lite.OperationalError:
        pass


def already_run(opts):
    if 'tmpdb' in opts and 'tmp' in opts and opts.tmp and opts.tmpdb:
        dbpath = opts.tmpdb
    else:
        dbpath = join(opts.workdir, 'trace.db')
    con = lite.connect(dbpath)
    try:
        with con:
            # check if table exists
            cur = con.cursor()
            param_hash = digest_parameters(opts, get_md5=True)
            cur.execute("select * from JOBs where Parameters_md5 = '%s'" % param_hash)
            found = len(cur.fetchall()) == 1
            if found:
                print_db(cur, 'JOBs')
    except lite.OperationalError:
        return False
    return found


def get_jobid(cur=None, workdir=None):
    try:
        if cur:
            cur.execute("select Id from JOBs where Id = (select max(id)  from JOBs)")
            return cur.fetchall()[0][0]
        con = lite.connect(join(abspath(workdir), 'trace.db'))
        with con:
            # check if table exists
            cur = con.cursor()
            cur.execute("select Id from JOBs where Id = (select max(id)  from JOBs)")
        return cur.fetchall()[0][0]
    except lite.OperationalError:
        return 0


def get_path_id(cur, path, workdir=None):
    path = abspath(path)
    if workdir:
        workdir = abspath(workdir)
        path    = relpath(path, workdir)
    cur.execute('SELECT Id from PATHs where Path="%s"' % path)
    return cur.fetchall()[0][0]


def add_path(cur, path, typ, jobid, workdir=None):
    if not path:  # case where path is None
        return
    path    = abspath(path)
    if workdir:
        workdir = abspath(workdir)
        path    = relpath(path, workdir)
    try:
        cur.execute("""
        insert into PATHs (Id  , Path, Type, JOBid)
        values (NULL, '%s', '%s', '%s')""" % (path, typ, jobid))
    except lite.IntegrityError:
        pass


def print_db(cur, name, no_print='', jobids=None, savedata=None, append=False,
             **kwargs):
    """
    print the content of a table to stdout in ascii/human-friendly format,
    or write it to a file in tab separated format (suitable for excel).

    :param cur: sqlite cursor
    :param name: DB name
    :param '' no_print: columns to skip
    :param None jobids: limit to a given list of jobids
    :param None savedata: path to file. If provided than output will be saved in
       tsv format,otherwise, it will be written to stdout in ascii format
    :param False append: whether to append to file,or to overwrite it.
    :param kwargs: dictionary with column names and values to use as filter
    """
    if jobids:
        cur.execute('select * from %s where %s in(%s)' % (
            name, 'JOBID' if name != 'JOBs' else 'ID',
            ','.join(map(str, jobids))))
    elif kwargs:
        filterstr = ' AND '.join("%s='%s'" % (k, 'NULL' if v == 'None' else v)
                                 for k, v in kwargs.iteritems())
        try:
            cur.execute("select * from %s where %s" % (name, filterstr))
        except lite.OperationalError:
            return
    else:
        cur.execute('select * from %s' % name)
    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    if not rows:
        return
    if isinstance(no_print, str):
        no_print = [no_print]
    for nop in no_print:
        if nop in names:
            bad = names.index(nop)
            names.pop(bad)
            rows = [[v for i, v in enumerate(row) if i != bad] for row in rows]
    if savedata is None:
        cols = [max(vals) for vals in zip(*[[len(str(v)) for v in row]
                                            for row in rows + [names]])]
        _ascii_print_db(name, names, cols, rows)
    else:
        if name == 'FILTER_OUTPUTs':
            _rev_tsv_print_db(names, rows, savedata, append)
        else:
            _tsv_print_db(name, names, rows, savedata, append)


def _ascii_print_db(name, names, cols, rows):
    print ',-' + '-' * len(name) + '-.'
    print '| ' + name + ' |'
    print ',-' + '-.-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-.'
    print '| ' + ' | '.join([('%{}s'.format(cols[i])) % str(v)
                             for i, v in enumerate(names)]) + ' |'
    print '|-' + '-+-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-|'
    print '| ' + '\n| '.join([' | '.join([('%{}s'.format(cols[i])) % str(v)
                                          for i, v in enumerate(row)]) + ' |'
                              for row in rows])
    print "'-" + '-^-'.join(['-' * cols[i] for i, v in enumerate(names)]) + "-'"


def _tsv_print_db(name, names, rows, savedata, append):
    out = open(savedata, 'a' if append else 'w')
    out.write('\t' + '\t'.join(names) + '\n')
    out.write(name)
    for row in rows:
        out.write('\t' + '\t'.join([str(v) for v in row]) + '\n')
    out.write('\n')
    out.close()


def _rev_tsv_print_db(names, rows, savedata, append):
    out = open(savedata, 'a' if append else 'w')
    pos = names.index('Name')
    out.write('\t'.join([row[pos] for row in rows]) + '\n')
    for col in xrange(len(rows[pos])):
        if col == pos:
            continue
        out.write('\t'.join([str(row[col]) for row in rows]) + '\n')
    out.write('\n')
    out.close()
