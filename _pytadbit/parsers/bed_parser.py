"""
simple BED and BEDgraph parser
"""

from pytadbit.utils.file_handling import magic_open

def _bed_float(line):
    crm, beg, end, _, val, _ =  line.split('\t', 5)
    return crm, int(beg), int(end), float(val)

def _bed_one(line):
    crm, beg, end, _ =  line.split('\t', 3)
    return crm, int(beg), int(end), 1

def _bedgraph_float(line):
    crm, beg, end, val =  line.split()
    return crm, int(beg), int(end), float(val)

def _3_col(line):
    crm, beg, end =  line.split()
    return crm, int(beg), int(end), 1

def _2_col(line):
    crm, beg =  line.split()
    beg = int(beg)
    return crm, beg, beg, 1


def parse_bed(fnam, resolution=1):
    """
    simple BED and BEDgraph parser that only checks for the fields 1, 2, 3 and 5
       (or 1, 2 and 3 if 5 not availbale).

    .. note::

        2 or 3 columns files can also be passed and will be interpreted,
        respectively, as chromosome/begin and chromosome/begin/end


    :param fnam: path to BED file
    :param 1 resolution: to bin the resulting dictionary

    :returns: a dictionnary with a count of number of entries found per bin. In
       case column 5 is present the values used tyo weight entries, otherwise
       each entry will weight 1.

    """

    fhandler = magic_open(fnam)
    line = fhandler.next()
    fpos = len(line)
    while (line.startswith('#')     or
           line.startswith('track') or
           line.startswith('browser')):
        fpos += len(line)
        line = fhandler.next()
    ##################
    # check file type
    try:
        # classic BED
        _, _, _, _, val, _ =  line.split('\t', 5)
        try:
            float(val)
            parse_line = _bed_float
        except ValueError:
            parse_line = _bed_one
    except ValueError:
        try:
            # BEDgraph
            _, _, _, val =  line.split('\t', 5)
            parse_line = _bedgraph_float
        except ValueError:
            try:
                # BEDgraph with no values
                _, _, _ =  line.split()
                parse_line = _3_col
            except ValueError:
                # only chromosome and begin position available
                parse_line = _2_col

    ####################################
    # go back to first informative line
    # parse
    dico = {}
    fhandler.seek(fpos)
    for line in fhandler:
        crm, beg, end, val = parse_line(line)
        pos = (beg + end - beg) / resolution
        dico.setdefault(crm, {})
        dico[crm].setdefault(pos, 0)
        dico[crm][pos] += val

    return dico
