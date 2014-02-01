"""
31 Jan 2014


"""

import inspect
import pytadbit
from pytadbit.tad_clustering import tad_cmo
from string import uppercase
import re
from sys import stderr


def condition(member):
    return type(member.__doc__) != type(None)

def get_all(module, all_members=None):
    if not all_members:
        all_members = {}
    members = inspect.getmembers(module, condition)
    if not members:
        return
    for name, member in members:
        if name.startswith('_'):
            continue
        if 'OLD' in name:
            continue
        if name.startswith('im_'):
            continue
        if name in all_members:
            continue
        try:
            if not 'pytadbit' in inspect.getsourcefile(member):
                continue
        except TypeError:
            continue
        if not hasattr(member, '__module__'):
            get_all(member, all_members)
            continue
        if inspect.ismodule(member):
            get_all(member, all_members)
            continue
        all_members[name] = {'dady': module.__name__, 'son': member}
        get_all(member, all_members)
    return all_members


def print_doc(member, header=1, indent=3, offset=45):
    doc = member.__doc__
    dochead = doc.split(':param')[0].split(':returns')[0]
    dochead = dochead.split('::')[0].split('..')[0].strip()
    dochead = re.sub(r':[a-z:]+:`(?:[a-zA-Z\.-_0-9]+\.)?([a-zA-Z0-9_-]+)`', '\\1', dochead)
    dochead = dochead.replace(']_', ']')
    dochead = re.sub('\. [^\.]+:.*', '.', dochead)
    if header == 1:
        title = ' '.join(member.split('.')[1:]).capitalize()
        # out =  '-' * (7 + len(title)) + '\n'
        out = title + ' module\n'
        out += '-' * (7 + len(title)) + '\n'
        return out
    elif header == 2:
        title = member.__name__
        indent=0
        # out =  (' ' * indent) + ('+' * (6 + len(title))) + '\n'
        out = (' ' * indent) + title + ' class\n'
        out += (' ' * indent) + ('+' * (6 + len(title))) + '\n'
        out += ' ' * (offset / 2)
        out += ('\n' + ' ' * (offset / 2)).join(
            [line.strip() for line in dochead.split('\n')])
        return out + '\n'
    elif header == 3:
        title = '**' + member.__name__ + '**'
        args = inspect.getargspec(member).args
        extra = []
        if 'savefig' in args or 'axe' in args:
            extra.append('1')
        if 'savedata' in args or 'outfile' in args or 'directory' in args:
            extra.append('2')
        title += ' (%s): ' % (','.join(extra)) if extra else ': '
        out  =  (' ' * indent) + '- ' + title
        out +=  ' ' * (offset - len(title) - 2 - indent)
        reout = ''
        for line in dochead.split('\n'):
            reout += out + line.strip() # + '\n'
            out = ' ' * offset
        return reout + '\n'


def main():
    """
    main function
    """

    all_members = get_all(pytadbit)
    get_all(pytadbit.utils, all_members)
    get_all(tad_cmo, all_members)
    modules = set([all_members[m]['son'].__module__ for m in all_members])

    numclasses = 0
    nummodules = len(modules)
    nfunctions = 0
    
    print '======================================='
    print 'Summary of TADbit classes and functions'
    print '=======================================\n'

    print '**Note**:'
    print '  - *1: for functions generating plots*'
    print '  - *2: for functions generating out files*'
    print ''
    for module in modules:
        print print_doc(module, header=1)
        
        submodules = [m for m in all_members
                      if all_members[m]['son'].__module__ == module]
        dadies = set([all_members[m]['dady'] for m in submodules
                      if all_members[m]['dady'][0] in uppercase])
        for member in sorted(submodules,
                             key=lambda x:all_members[x]['son'].__name__[0] in uppercase):
            if all_members[member]['dady'] in dadies or member in dadies:
                continue
            if all_members[member]['son'].__name__[0] in uppercase:
                print print_doc(all_members[member]['son'], header=2, indent=3)
                numclasses += 1
            else:
                nfunctions += 1
                print print_doc(all_members[member]['son'], header=3, indent=3)
        for dady in dadies:
            numclasses += 1
            print print_doc(all_members[dady]['son'], offset=9, header=2)
            for member in sorted(submodules):
                if all_members[member]['dady'] != dady:
                    continue
                nfunctions += 1
                print print_doc(all_members[member]['son'], header=3,
                                indent=6)
    stderr.write('Reporting %s modules, %s classes and %s functions\n' %(
        nummodules, numclasses, nfunctions))
    
if __name__ == "__main__":
    exit(main())
