"""
02 Dec 2012


"""

from os.path import isfile


def parse_tads(handler):
    """
    Parse a tab separated value file that contains the list of TADs of a given
    experiment. This file might have been generated whith the
    :func:`pytadbit.tadbit.print_result_R` or with the R binding for tadbit

    :param handler: path to file
    :param 1 bin_size: resolution of the experiment

    :returns: list of TADs and list of weights, each TAD being a dict of type:

    ::

      {TAD_num: {'start': start,
                 'end'  : end,
                 'brk'  : end,
                 'score': score}}
    """
    tads = {}
    weights = None
    if isinstance(handler, tuple):
        handler, weights = handler
    if isinstance(handler, dict):
        try:
            for pos in xrange(len(handler['end'])):
                start = float(handler['start'][pos])
                end   = float(handler['end'][pos])
                try:
                    score = float(handler['score'][pos])
                except KeyError: # missing
                    score = 10.0
                except TypeError: # last one
                    score = 10.0
                tads[pos + 1] = {'start': start,
                                 'end'  : end,
                                 'brk'  : end,
                                 'score': score}
        except KeyError: # the other dict format
            for pos in handler:
                if not ('start' in handler[pos]
                        or not 'end' in handler[pos]
                        or not 'brk' in handler[pos]
                        or not 'score' in handler[pos]):
                    raise Exception('ERROR: bad format\n')
            tads = handler
    elif isfile(handler):
        for line in open(handler):
            if line.startswith('#'): continue
            try:
                pos, start, end, score = line.split()
                dens = float('nan')
            except ValueError:
                pos, start, end, score, dens = line.split()
                dens  = float(dens)
            start = float(start) - 1
            end   = float(end)   - 1
            pos   = int(pos)
            try:
                score = float(score)
            except ValueError: # last one
                score = 10.0
            tads[pos] = {'start' : start,
                         'end'   : end,
                         'brk'   : end,
                         'score' : score,
                         'height': dens}
    else:
        raise Exception('File %s not found\n' % (handler))
    return tads, weights
