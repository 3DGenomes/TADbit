"""
02 Dec 2012


"""

from os.path import isfile

def parse_tads(f_name, max_size=3000000, bin_size=1):
    tads = {}
    forbidden = []
    if type(f_name) is dict:
        for pos in xrange(len(f_name['end'])):
            start = float(f_name['start'][pos])
            end   = float(f_name['end'][pos])
            try:
                score = float(f_name['score'][pos])
            except TypeError:
                score = None
            diff  = end - start
            tads[pos] = {'start': start,
                         'end'  : end,
                         'brk'  : end,
                         'score': score}
            if diff * bin_size > max_size:
                forbidden += range(int(start), int(end+1))
                tads[pos]['brk'] = None
    elif isfile(f_name):
        for line in open(f_name):
            if line.startswith('#'): continue
            pos, start, end, score = line.split()
            start = float(start)
            end   = float(end)
            pos   = int(pos)
            try:
                score = float(score)
            except ValueError:
                score = None
            diff  = end - start
            tads[pos] = {'start': start,
                         'end'  : end,
                         'brk'  : end,
                         'score': score}
            if diff * bin_size > max_size:
                forbidden += range(int(start), int(end+1))
                tads[pos]['brk'] = None
    else:
        raise Exception('File {} not found\n'.format(f_name))
    return tads, set(forbidden)
