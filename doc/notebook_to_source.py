"""
24 Oct 2013
"""

import os, re

CURPATH = os.path.abspath(os.path.split(os.path.realpath(__file__))[0])
PATH = os.path.realpath(os.path.join(CURPATH, 'notebooks'))

def main():
    """
    main function
    """    
    global CURPATH, PATH
    os.system('mkdir -p ' + os.path.join(CURPATH, 'source', 'nbpictures'))
    for fname in os.listdir(PATH):
        if not fname.endswith('.ipynb'):
            continue
        os.system('jupyter nbconvert --to rst %s --output %s' % (
            os.path.join(PATH, fname), os.path.join(PATH, fname[:-6] + '.rst')))
        extra = (fname.split('_')[0] + '/') if '_' in fname else ''
        passing = False
        lines = []
        for line in open(os.path.join(PATH, fname[:-6] + '.rst')):
            if '## REMOVE' in line:
                lines = lines[:-2]
                passing=True
            elif '## STOP REMOVE' in line:
                passing=False
                continue
            if passing:
                continue
            line = re.sub('In\[[0-9 ]+\]:\n', '', line)
            # these two next are just to display alignment in colors
            line = re.sub('\[m', '[0m', line)
            line = re.sub('parsed-literal', 'ansi-block', line)
            if extra:
                line = re.sub(fname[:-6] + '_files/', '../nbpictures/', line)
                line = re.sub(os.path.join(PATH), '../nbpictures/', line)
            else:
                line = re.sub(fname[:-6] + '_files/', 'nbpictures/', line)
                line = re.sub(os.path.join(PATH), 'nbpictures/', line)
            lines.append(line)
        out = open(os.path.join(CURPATH, 'source', extra + fname[:-6] + '.rst'), 'w')
        out.write(''.join(lines))
        out.close()
        os.system('rm -f ' + os.path.join(PATH, fname[:-6] + '.rst'))
        try:
            for iname in os.listdir(os.path.join(PATH, fname[:-6] + '_files')):
                if not iname.endswith('.png'):
                    continue
                os.system('cp %s %s' % (
                    os.path.join(
                        PATH, fname[:-6] + '_files', iname),
                    os.path.join(
                        CURPATH, 'source', 'nbpictures')))
            os.system('rm -rf ' + os.path.join(PATH, fname[:-6] + '_files'))
        except OSError:
            pass
        
    PATH = os.path.realpath('notebooks')
    os.system('cd %s; make html' % CURPATH)

if __name__ == "__main__":
    exit(main())
