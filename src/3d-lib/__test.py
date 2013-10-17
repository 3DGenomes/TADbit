"""
This is to quickly test eqv_rms_drms.rmsdRMSD_wrapper function.
just compile with:
g++-4.4 -shared eqv_rms_drms_py.cpp -I/usr/include/python2.7 -lm -lpthread  -fPIC -g -O3 -Wall -o eqv_rms_drms.so
g++-4.4 -shared centroid_py.cpp -I/usr/include/python2.7 -lm -lpthread  -fPIC -g -O3 -Wall -o centroid_py.so
and run this test as:
python test.py
or through gdb:
gdb python -ex 'run test.py'
"""

import eqv_rms_drms


lists = {}
for model in ['2', '3', '4', '5', '8']:
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


import centroid_py

lala = centroid_py.centroid_wrapper([lists['2'][0], lists['3'][0], lists['4'][0], lists['5'][0], lists['8'][0]],
                                    [lists['2'][1], lists['3'][1], lists['4'][1], lists['5'][1], lists['8'][1]],
                                    [lists['2'][2], lists['3'][2], lists['4'][2], lists['5'][2], lists['8'][2]],
                                    len(lists['2'][0]), 5)
print lala

