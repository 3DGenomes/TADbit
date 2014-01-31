"""
31 Jan 2014


"""

import inspect
import pytadbit
from string import uppercase

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
        # if not inspect.isclass(member):
        #     print '   ' + name + ' '*(27-len(name)),
        # else:
        #     print name + ' '*(30-len(name)),
        # print ('\n'+' '*23).join(member.__doc__.split(':param')[0].split(':returns')[0].split('::')[0].split('..')[0].strip().split('\n'))
        # print '-'*105
        get_all(member, all_members)
    return all_members


def print_doc(member, header=1, offset=34):
    headers = {1: 5, 2: 7}
    print ' '*(offset - len(member.__name__) - headers[header]),
    print ('\n'+' '*offset).join([line.strip() for line in member.__doc__.split(':param')[0].split(':returns')[0].split('::')[0].split('..')[0].strip().split('\n')])


def main():
    """
    main function
    """
    all_members = get_all(pytadbit)

    modules = set([all_members[m]['son'].__module__ for m in all_members])

    for module in modules:
        print '\n' + '*' * 100
        print '\n>'+module
        submodules = [m for m in all_members
                      if all_members[m]['son'].__module__ == module]
        dadies = set([all_members[m]['dady'] for m in submodules
                      if all_members[m]['dady'][0] in uppercase])
        for dady in dadies:
            print '\n * ' + dady,
            print_doc(all_members[dady]['son'], offset=25)
            print '-' * (34)
            for member in submodules:
                if all_members[member]['dady'] != dady:
                    continue
                print '   - ' + member,
                print_doc(all_members[member]['son'], header=2)
        if dadies:
            print '=' * 100
        for member in submodules:
            if all_members[member]['dady'] in dadies or member in dadies:
                continue
            print ' * ' + member,
            print_doc(all_members[member]['son'], header=1, offset=32)
        


if __name__ == "__main__":
    exit(main())
