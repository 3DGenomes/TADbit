"""
some utils relative to sqlite
"""
import sqlite3 as lite
from os.path import abspath, relpath, join
from hashlib import md5

def already_run(opts):
    con = lite.connect(join(opts.workdir, 'trace.db'))
    with con:
        # check if table exists
        cur = con.cursor()
        param_hash = md5(' '.join(
            ['%s:%s' % (k, int(v) if isinstance(v, bool) else v)
             for k, v in sorted(opts.__dict__.iteritems())
             if not k in ['workdir', 'func', 'tmp', 'keep_tmp']])).hexdigest()
        cur.execute("select * from JOBs where Parameters_md5 = '%s'" % param_hash)
        found = len(cur.fetchall()) == 1
        if found:
            print_db(cur, 'JOBs')
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

def print_db(cur, name, no_print='Parameters_md5'):
    cur.execute('select * from %s' % name)
    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    if no_print in names:
        bad = names.index(no_print)
        names.pop(bad)
        rows = [[v for i, v in enumerate(row) if i != bad] for row in rows]
    cols = [max(vals) for vals in zip(*[[len(str(v)) for v in row]
                                        for row in rows + [names]])]
    print ',-' + '-' * len(name) + '-.'
    print '| ' + name + ' |'
    print ',-' + '-.-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-.'
    print '| ' + ' | '.join([('%{}s'.format(cols[i])) % str(v) for i, v in enumerate(names)]) + ' |'
    print '|-' + '-+-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-|'
    print '| ' + '\n| '.join([' | '.join([('%{}s'.format(cols[i])) % str(v)
                                        for i, v in enumerate(row)]) + ' |'
                              for row in rows])
    print "'-" + '-^-'.join(['-' * cols[i] for i, v in enumerate(names)]) + "-'"
