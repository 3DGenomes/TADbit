"""
This is to quickly test eqv_rms_drms.rmsdRMSD_wrapper function.
just compile with:
g++-4.4 -shared eqv_rms_drms_py.cpp -I/usr/include/python2.7 -lm -lpthread  -fPIC -g -O3 -Wall -o eqv_rms_drms.so
and run this test as:
python test.py
or through gdb:
gdb python -ex 'run test.py'
"""

import eqv_rms_drms

def lala(lo):
    print 'avant'
    la = eqv_rms_drms.rmsdRMSD_wrapper(list1, list2, len(list1), 200, lo)
    print la
    for i, j in enumerate(la):
        print i+1, j
    print 'apres'

list1 = []
list2 = []
for line in open ('/tmp/tmp_cons/model_1_rnd57.xyz'):
    _, _, x, y, z = line.split()
    list1.append((float(x), float(y), float(z)))

    
for line in open ('/tmp/tmp_cons/model_9_rnd4.xyz'):
    _, _, x, y, z = line.split()
    list2.append((float(x), float(y), float(z)))


lala(1)
lala(1)
lala(1)
lala(1)


print 'finished test'

