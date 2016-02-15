"""
some utils relative to sqlite
"""
import sqlite3 as lite
from os.path import abspath, relpath, join

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

def print_db(cur, name):
    cur.execute('select * from %s' % name)
    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    cols = [max(vals) for vals in zip(*[[len(str(v)) for v in row]
                                        for row in rows + [names]])]
    print ',-' + '-' * len(name) + '-.'
    print '| ' + name + ' |'
    print ',-' + '-.-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-.'
    print '| ' + ' | '.join([('%{}s'.format(cols[i])) % str(v) for i, v in enumerate(names)]) + ' |'
    print '|-' + '-+-'.join(['-' * cols[i] for i, v in enumerate(names)]) + '-|'
    print '| ' + '\n| '.join([' | '.join([('%{}s'.format(cols[i])) % str(v)
                                        for i, v in enumerate(row)]) + ' |'  for row in rows])
    print "'-" + '-^-'.join(['-' * cols[i] for i, v in enumerate(names)]) + "-'"
