"""
02 Dec 2012


"""

import bz2, gzip, zipfile, tarfile

def magic_open(filename, verbose=False):
    """
    To read uncompressed zip gzip bzip2 or tar.xx files

    :param filename: either a path to a file, or a file handler

    :returns: opened file ready to be iterated
    """
    if isinstance(filename, str):
        fhandler = file(filename, 'rb')
        inputpath = True
        if tarfile.is_tarfile(filename):
            print 'tar'
            thandler = tarfile.open(filename)
            if len(thandler.members) != 1:
                raise NotImplementedError(
                    'Not exactly one file in this tar archieve.')
            return magic_open(thandler.extractfile(thandler.getnames()[0]))
    else:
        fhandler = filename
        filename = fhandler.name
        inputpath = False
        start_of_file = ''
    if inputpath:
        start_of_file = fhandler.read(1024)
        fhandler.seek(0)
    if start_of_file.startswith('\x50\x4b\x03\x04'):
        if verbose:
            print 'zip'
        zhandler = zipfile.ZipFile(fhandler)
        if len(zhandler.NameToInfo) != 1:
            raise NotImplementedError(
                'Not exactly one file in this zip archieve.')
        return zhandler.open(zhandler.NameToInfo.keys()[0])
    if start_of_file.startswith('\x42\x5a\x68'):
        if verbose:
            print 'bz2'
        fhandler.close()
        return bz2.BZ2File(filename)
    if start_of_file.startswith('\x1f\x8b\x08'):
        if verbose:
            print 'gz'
        return gzip.GzipFile(fileobj=fhandler)
    if verbose:
        print 'text'
    return fhandler
