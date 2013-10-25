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
    os.system('mkdir -p ../source/nbpictures')
    for fname in os.listdir('.'):
        if not fname.endswith('.ipynb'):
            continue
        # os.system('ipython nbconvert --format=rst ' + fname)
        os.system('ipython nbconvert --to rst ' + fname)
        extra = (fname.split('_')[0] + '/') if '_' in fname else ''
        out = open('../source/' + extra + fname[:-6] + '.rst', 'w')
        passing = False
        for line in open(fname[:-6] + '.rst'):
            if '## REMOVE' in line:
                passing=True
            elif '## STOP REMOVE' in line:
                passing=False
                continue
            if passing:
                continue
            line = re.sub('In\[[0-9 ]+\]:\n', '', line)
            if extra:
                line = re.sub(fname[:-6] + '_files/', '../nbpictures/', line)
            else:
                line = re.sub(fname[:-6] + '_files/', 'nbpictures/', line)
            out.write(line)
        out.close()
        os.system('rm -f ' + fname[:-6] + '.rst')
        try:
            for iname in os.listdir(fname[:-6] + '_files'):
                if not iname.endswith('.png'):
                    continue
                os.system('cp '+ fname[:-6] + '_files/' + iname + ' ../source/nbpictures/')
            os.system('rm -rf '+ fname[:-6] + '_files')
        except OSError:
            pass
    os.system('cd ..; make html')
    
if __name__ == "__main__":
    exit(main())
