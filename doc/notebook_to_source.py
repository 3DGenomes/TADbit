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
        passing = False
        lines = []
        for line in open(fname[:-6] + '.rst'):
            if '## REMOVE' in line:
                lines = lines[:-2]
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
            lines.append(line)
        out = open('../source/' + extra + fname[:-6] + '.rst', 'w')
        out.write(''.join(lines))
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
