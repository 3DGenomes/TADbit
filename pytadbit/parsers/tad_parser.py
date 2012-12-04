"""
02 Dec 2012


"""

def parse_tads(f_name, max_size=3000000, bin_size=1):
    tads = {}
    forbidden = []
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
    return tads, set(forbidden)
