"""
02 Dec 2012


"""

from os.path import isfile


def parse_tads(handler, max_size=3000000, bin_size=1):
    """
    Parse a tsv file that contains the list of TADs of a given experiment.
    This file might have been generated whith the
    :func:`pytadbit.tadbit.print_result_R` or with the R binding for tadbit

    :param handler: path to file
    :param 3000000 max_size: maximum size allowed for a TAD
    :param 1 bin_size: resolution of the experiment

    :returns: list of TADs, each TAD being a dict of type:

    ::
    
      {TAD_num: {'start': start,
                 'end'  : end,
                 'brk'  : end,
                 'score': score}}
    """
    tads = {}
    if type(handler) is dict:
        for pos in xrange(len(handler['end'])):
            start = float(handler['start'][pos])
            end   = float(handler['end'][pos])
            try:
                score = float(handler['score'][pos])
            except TypeError:
                score = None
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
            except ValueError:
                score = None
            tads[pos] = {'start': start,
                         'end'  : end,
                         'brk'  : end,
                         'score': score}
    else:
        raise Exception('File {} not found\n'.format(handler))
    return tads
