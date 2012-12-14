"""
02 Dec 2012


"""
from pytadbit.tads_aligner.globally import needleman_wunsch


def consensusize(ali1, ali2):
    consensus = []
    for pos in xrange(len(ali1)):
        if ali1[pos] != ali2[pos]:
            try:
                bound = (ali1[pos]+ali2[pos])/2
            except TypeError:
                bound = ali1[pos] if type(ali1[pos]) is float else ali2[pos]
        else:
            bound = ali1[pos]
        consensus.append(bound)
    return consensus


def align(sequences, method='global', **kwargs):
    """
    Align Topologically associated domains. Supports multiple alignment by
    building a consensus TAD and aligning each TAD to it.
    Note: as long as we are using multiple alignments in an iterative way,
    the order of sequences we be relevant. Here TADs are sorted in order to try
    to reduce this problem.

    :argument global method: method used to align.
    """
    if method=='global':
        aligner = needleman_wunsch
    if len(sequences) == 2:
        return aligner(sequences[0], sequences[1], **kwargs)
    if len(sequences) > 2:
        dico = {}
        reference = None
        others = []
        for j, (i, seq) in enumerate(sorted(enumerate(sequences),
                                            key=lambda x: x[1])):
            if not reference:
                reference = seq
            else:
                others.append(seq)
            dico[j] = {'sort':i,
                       'seq' :seq}
        aligneds = []
        scores = 0
        for other in xrange(1, len(sequences)):
            [align1, align2], score = needleman_wunsch(reference, dico[other]['seq'],
                                                       **kwargs)
            scores += score
            if len(reference) != len(align1):
                for pos in xrange(len(align1)):
                    try:
                        if align1[pos] != '-':
                            continue
                        for ali in aligneds:
                            ali.insert(pos, '-')
                    except IndexError:
                        for ali in aligneds:
                            ali.append('-')
            if not aligneds:
                aligneds.append(align1)
            aligneds.append(align2)
            reference = consensusize(align1, align2)
        sort_alis = [[] for _ in xrange(len(dico))]
        for seq in xrange(len(dico)):
            sort_alis[dico[seq]['sort']] = aligneds[seq][:]
        return sort_alis, scores

