"""
02 Dec 2012


"""

from os.path import isfile


def parse_tads(handler):
    """
    Parse a tsv file that contains the list of TADs of a given experiment.
    This file might have been generated whith the
    :func:`pytadbit.tadbit.print_result_R` or with the R binding for tadbit

    :param handler: path to file
    :param 1 bin_size: resolution of the experiment

    :returns: list of TADs, each TAD being a dict of type:

    ::
    
      {TAD_num: {'start': start,
                 'end'  : end,
                 'brk'  : end,
                 'score': score}}
    """
    tads = {}
    weights = None
    if type(handler) is tuple:
        handler, weights = handler
    if type(handler) is dict:
        for pos in xrange(len(handler['end'])):
            start = float(handler['start'][pos])
            end   = float(handler['end'][pos])
            try:
                score = float(handler['score'][pos])
            except TypeError: # last one
                score = 10.0
            tads[pos] = {'start': start,
                         'end'  : end,
                         'brk'  : end,
                         'score': score}
    elif isfile(handler):
        for line in open(handler):
            if line.startswith('#'): continue
            pos, start, end, score = line.split()
            start = float(start)
            end   = float(end)
            pos   = int(pos)
            try:
                score = float(score)
            except ValueError: # last one
                score = 10.0
            tads[pos] = {'start': start,
                         'end'  : end,
                         'brk'  : end,
                         'score': score}
    else:
        raise Exception('File %s not found\n' % (handler))
    return tads, weights
