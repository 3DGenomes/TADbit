"""
This is to quickly test eqv_rms_drms.rmsdRMSD_wrapper function.
just compile with:
g++ -shared eqv_rms_drms_py.cpp -I/usr/include/python2.7 -lm -lpthread  -fPIC -g -O3 -Wall -o eqv_rms_drms.so
g++ -shared centroid_py.cpp -I/usr/include/python2.7 -lm -lpthread  -fPIC -g -O3 -Wall -o centroid.so
and run this test as:
python test.py
or through gdb:
gdb python -ex 'run test.py'
"""

import eqv_rms_drms
from pytadbit.imp.impmodel import IMPmodel


models = ['2', '3', '4', '5', '8']

lists = {}
for model in models:
    lists[model] = []
    for line in open ('/scratch/shared/TADs/T0/models/T0_10_1065_1545/models/model.{}.xyz'.format(model)):
    # for line in open ('test/model.{}.xyz'.format(model)):
        _, _, x, y, z = line.split()
        lists[model].append((float(x), float(y), float(z)))
    lists[model] = [list(i) for i in zip(*lists[model])]

lala =  eqv_rms_drms.rmsdRMSD_wrapper([lists[m][0] for m in models], [lists[m][1] for m in models], [lists[m][2] for m in models],
                                    len(lists[models[0]][0]),
                                    200., 0, models, len(models), 0)

for model1 in ['2', '3', '4', '5', '8']:
    for model2 in ['2', '3', '4', '5', '8']:
        print model1, model2, eqv_rms_drms.rmsdRMSD_wrapper(lists[model1][0], lists[model1][1], lists[model1][2], lists[model2][0], lists[model2][1], lists[model2][2],
                                                            len(lists[model1][0]),
                                                            200., 0)

print 'finished test'


import centroid

from itertools import permutations

for m1, m2, m3, m4, m5 in permutations(models, 5):
    
    # print m1, m2, m3, m4, m5
    models2 = [m1, m2, m3, m4, m5]
    lala = centroid.centroid_wrapper([lists[m1][0], lists[m2][0], lists[m3][0], lists[m4][0], lists[m5][0]],
                                     [lists[m1][1], lists[m2][1], lists[m3][1], lists[m4][1], lists[m5][1]],
                                     [lists[m1][2], lists[m2][2], lists[m3][2], lists[m4][2], lists[m5][2]],
                                     len(lists['2'][0]), 5, 0, 0)
    print '   ====>>>> ', models2[lala]


avg = centroid.centroid_wrapper([lists[m1][0], lists[m2][0], lists[m3][0], lists[m4][0], lists[m5][0]],
                                [lists[m1][1], lists[m2][1], lists[m3][1], lists[m4][1], lists[m5][1]],
                                [lists[m1][2], lists[m2][2], lists[m3][2], lists[m4][2], lists[m5][2]],
                                len(lists['2'][0]), 5, 0, 1)

avgmodel = IMPmodel((('x', avg[0]), ('y', avg[1]), ('z', avg[2]),
                     ('rand_init', 'avg'), ('objfun', None),('radius',0.1)))
avgmodel.write_cmm('test/')


head = '<marker_set name="{}">\n'
mark = '<marker id="{0}" x="{1}" y="{2}" z="{3}" r="0.1" g="0" b="1" radius="0.1" note="{0}"/>\n'
link = '<link id1="{}" id2="{}" r="1" g="1" b="1" radius="0.05"/>\n'
end = '</marker_set>\n'
for model in models:
    out = open('test/model.{}.cmm'.format(model), 'w')
    out.write(head.format(model))
    for i, (x, y, z) in enumerate(zip(*lists[model])):
        out.write(mark.format(i+1, x,y,z))
    for i in range(len(lists[model][0])-1):
        out.write(link.format(i+1, i+2))
    out.write(end)
    out.close()


exit()

