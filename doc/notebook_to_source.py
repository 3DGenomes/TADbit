"""
24 Oct 2013


"""

import os, re

PATH = 'notebooks'

def main():
    """
    main function
    """
    
    os.chdir(PATH)
    for fname in os.listdir('.'):
        if not fname.endswith('.ipynb'):
            continue
        os.system('ipython nbconvert --format=rst ' + fname)
        out = open('../source/' + fname[:-6] + '.rst', 'w')
        for line in open(fname[:-6] + '.rst'):
            line = re.sub('In\[[0-9 ]+\]:\n', '', line)
            line = re.sub(fname[:-6] + '_files/', 'pictures/', line)
            # line = re.sub('~', '^', line)
            # print line,
            out.write(line)
        out.close()
        for iname in os.listdir(fname[:-6] + '_files'):
            if not iname.endswith('.png'):
                continue
            os.system('cp '+ fname[:-6] + '_files/' + iname + ' ../source/pictures/')
            
    os.system('cd ..; make html')
    
if __name__ == "__main__":
    exit(main())
