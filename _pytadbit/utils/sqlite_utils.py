"""
some utils relative to sqlite
"""
from __future__ import print_function
from sys        import stdout
from os.path    import abspath, relpath, join, exists
from hashlib    import md5
from functools  import wraps
from time       import sleep
import sqlite3 as lite


try:
    basestring
except NameError:
    basestring = str

def digest_parameters(opts, get_md5=True, extra=None):
    """
    digestion is truncated at 10 characters. More than a 1000 runs are not
    expected on a same sample (usually ~10).

    :param None extra: extra parameter to remove from digestion
    """
    extra = extra or []
    if get_md5:
        param_hash = md5(' '.join(
            ['%s:%s' % (k, int(v) if isinstance(v, bool) else v)
             for k, v in sorted(opts.__dict__.items())
             if k not in ['force', 'workdir', 'func', 'tmp',
                          'skip', 'keep_tmp', 'tmpdb'] + extra]).encode('utf-8')).hexdigest()[:10]
        return param_hash
    parameters = ' '.join(
        ['%s:%s' % (k, int(v) if isinstance(v, bool) else v)
         for k, v in opts.__dict__.items()
         if k not in ['fastq', 'index', 'renz', 'iterative', 'workdir',
                      'skip', 'func', 'tmp', 'keep_tmp'] + extra and v is not None])
    parameters = parameters.replace("'", "")
    parameters = parameters.replace("|", "_")
    parameters = parameters.replace(";", "_")
    parameters = parameters.replace("#", "_")
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
        print(' - deleted %d elements with %s = %s' % (len(elts), col, val))
    except lite.OperationalError:
        pass


def already_run(opts, extra=None):
    if 'tmpdb' in opts and 'tmp' in opts and opts.tmp and opts.tmpdb:
        dbpath = opts.tmpdb
    else:
        dbpath = join(opts.workdir, 'trace.db')
    if not exists(dbpath):
        raise IOError('ERROR: DB file: %s not found.' % dbpath)
    con = lite.connect(dbpath)
    try:
        with con:
            # check if table exists
            cur = con.cursor()
            param_hash = digest_parameters(opts, get_md5=True, extra=extra)
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
             columns=None, tsv=False, **kwargs):
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
    if not savedata:
        savedata = stdout
    if not columns:
        columns = ['*']
    columns = [c.lower() for c in columns]
    if jobids:
        cur.execute('select %s from %s where %s in(%s)' % (
            ','.join(columns), name, 'JOBID' if name != 'JOBs' else 'ID',
            ','.join(map(str, jobids))))
    elif kwargs:
        filterstr = ' AND '.join("%s='%s'" % (k, 'NULL' if v == 'None' else v)
                                 for k, v in kwargs.items())
        try:
            cur.execute("select %s from %s where %s" % (
                ','.join(columns), name, filterstr))
        except lite.OperationalError:
            return
    else:
        cur.execute('select %s from %s' % (','.join(columns), name))
    names = [x[0] for x in cur.description]
    rows = cur.fetchall()
    rows = [['{:,}'.format(v) if isinstance(v, int) else v for v in vals]
            for vals in rows]
    if not rows:
        return
    if isinstance(no_print, basestring):
        no_print = [no_print]
    for nop in no_print:
        if nop.lower() in columns:
            continue
        if nop in names:
            bad = names.index(nop)
            names.pop(bad)
            rows = [[v for i, v in enumerate(row) if i != bad] for row in rows]
    if not tsv:
        cols = [max(vals) for vals in zip(*[[len(str(v)) for v in row]
                                            for row in rows + [names]])]
        _ascii_print_db(name, names, cols, rows, savedata)
    else:
        if name == 'FILTER_OUTPUTs':
            _rev_tsv_print_db(names, rows, savedata, append)
        else:
            _tsv_print_db(name, names, rows, savedata, append)


def _ascii_print_db(name, names, cols, rows, savedata):
    if isinstance(savedata, basestring):
        out = open(savedata, 'w')
        to_close = True
    else:
        out = savedata
        to_close = False
    chars = ''
    chars += ',-' + '-' * len(name) + '-.\n'
    chars += '| ' + name + ' |\n'
    chars += ',-' + '-.-'.join('-' * cols[i] for i, v in enumerate(names)) + '-.\n'
    chars += '| ' + ' | '.join(('%{}s'.format(cols[i])) % str(v)
                             for i, v in enumerate(names)) + ' |\n'
    chars += '|-' + '-+-'.join('-' * cols[i] for i, v in enumerate(names)) + '-|\n'
    chars += '| ' + '\n| '.join(' | '.join(('%{}s'.format(cols[i])) % str(v)
                                          for i, v in enumerate(row)) + ' |'
                              for row in rows) + '\n'
    chars += "'-" + '-^-'.join('-' * cols[i] for i, v in enumerate(names)) + "-'\n"
    out.write(chars)
    if to_close:
        out.close()


def _tsv_print_db(name, names, rows, savedata, append):
    if isinstance(savedata, basestring):
        out = open(savedata, 'a' if append else 'w')
        to_close = True
    else:
        out = savedata
        to_close = False
    out.write('\t' + '\t'.join(names) + '\n')
    out.write(name)
    for row in rows:
        out.write('\t' + '\t'.join([str(v) for v in row]) + '\n')
    out.write('\n')
    if to_close:
        out.close()


def _rev_tsv_print_db(names, rows, savedata, append):
    if isinstance(savedata, basestring):
        out = open(savedata, 'a' if append else 'w')
        to_close = True
    else:
        out = savedata
        to_close = False
    pos = names.index('Name')
    out.write('\t'.join([row[pos] for row in rows]) + '\n')
    for col in range(len(rows[pos])):
        if col == pos:
            continue
        out.write('\t'.join([str(row[col]) for row in rows]) + '\n')
    out.write('\n')
    if to_close:
        out.close()


def retry(exceptions, tries=4, delay=3, backoff=2, silent=False, logger=None):
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    :param exceptions: the exception(s) to check. may be a tuple of
        exceptions to check.
    :type exceptions: Exception type, exception instance, or tuple containing
        any number of both (eg. IOError, IOError(errno.ECOMM), (IOError,), or
        (ValueError, IOError(errno.ECOMM))
    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param silent: If set then no logging will be attempted
    :type silent: bool
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    """
    try:
        len(exceptions)
    except TypeError:
        exceptions = (exceptions,)
    all_exception_types = tuple(set(x if type(x) == type else x.__class__ for x in exceptions))
    exception_types = tuple(x for x in exceptions if type(x) == type)
    exception_instances = tuple(x for x in exceptions if type(x) != type)

    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except all_exception_types as e:
                    if (not any(x for x in exception_types if isinstance(e, x))
                        and not any(x for x in exception_instances if type(x) == type(e) and x.args == e.args)):
                        raise
                    msg = "%s, Retrying in %d seconds..." % (str(e) if str(e) != "" else repr(e), mdelay)
                    if not silent:
                        if logger:
                            logger.warning(msg)
                        else:
                            print(msg)
                    sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry
