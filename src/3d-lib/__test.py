"""
This is to quickly test eqv_rms_drms.rmsdRMSD_wrapper function.
just compile with:
g++-4.4 -shared eqv_rms_drms_py.cpp -I/usr/include/python2.7 -lm -lpthread  -fPIC -g -O3 -Wall -o eqv_rms_drms.so
g++-4.4 -shared centroid_py.cpp -I/usr/include/python2.7 -lm -lpthread  -fPIC -g -O3 -Wall -o centroid.so
and run this test as:
python test.py
or through gdb:
gdb python -ex 'run test.py'
"""

import eqv_rms_drms


models = ['2', '3', '4', '5', '8']

lists = {}
for model in models:
    lists[model] = []
    for line in open ('/scratch/shared/TADs/T0/models/T0_10_1065_1545/models/model.{}.xyz'.format(model)):
        _, _, x, y, z = line.split()
        lists[model].append((float(x), float(y), float(z)))
    lists[model] = [list(i) for i in zip(*lists[model])]
    
for model1 in ['2', '3', '4', '5', '8']:
    for model2 in ['2', '3', '4', '5', '8']:
        print model1, model2, eqv_rms_drms.rmsdRMSD_wrapper(lists[model1][0], lists[model1][1], lists[model1][2], lists[model2][0], lists[model2][1], lists[model2][2], len(lists[model1][0]),
                                                            200., 0)

print 'finished test'


import centroid

from itertools import permutations

for m1, m2, m3, m4, m5 in permutations(models, 5):
    
    print m1, m2, m3, m4, m5
    lala = centroid.centroid_wrapper([lists[m1][0], lists[m2][0], lists[m3][0], lists[m4][0], lists[m5][0]],
                                     [lists[m1][1], lists[m2][1], lists[m3][1], lists[m4][1], lists[m5][1]],
                                     [lists[m1][2], lists[m2][2], lists[m3][2], lists[m4][2], lists[m5][2]],
                                     len(lists['2'][0]), 5, 1)
    print '   ====>>>> ', models[lala]


lala = centroid.centroid_wrapper([lists['2'][0], lists['3'][0], lists['4'][0], lists['5'][0], lists['8'][0]],
                                 [lists['2'][1], lists['3'][1], lists['4'][1], lists['5'][1], lists['8'][1]],
                                 [lists['2'][2], lists['3'][2], lists['4'][2], lists['5'][2], lists['8'][2]],
                                 len(lists['2'][0]), 5, 1)
print models[lala]


lala = centroid.centroid_wrapper([lists['4'][0], lists['3'][0], lists['5'][0], lists['8'][0], lists['2'][0]],
                                 [lists['4'][1], lists['3'][1], lists['5'][1], lists['8'][1], lists['2'][1]],
                                 [lists['4'][2], lists['3'][2], lists['5'][2], lists['8'][2], lists['2'][2]],
                                 len(lists['2'][0]), 5, 1)
print models[lala]


lala = centroid.centroid_wrapper([lists['2'][0], lists['8'][0], lists['5'][0], lists['4'][0], lists['3'][0]],
                                 [lists['2'][1], lists['8'][1], lists['5'][1], lists['4'][1], lists['3'][1]],
                                 [lists['2'][2], lists['8'][2], lists['5'][2], lists['4'][2], lists['3'][2]],
                                 len(lists['2'][0]), 5, 1)
print models[lala]




