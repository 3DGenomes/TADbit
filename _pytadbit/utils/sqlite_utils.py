"""
some utils relative to sqlite
"""
import sqlite3 as lite

def get_path_id(cur, name):
    cur.execute('SELECT Id from PATHs where Path="%s"' % name)
    return cur.fetchall()[0][0]

def add_path(cur, path, typ):
    try:
        cur.execute("insert into PATHs (Id  , Path, Type) values (NULL, '%s', '%s')" % (
            path, typ))
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
