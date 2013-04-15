"""
02 Dec 2012


"""
from pytadbit.boundary_aligner.globally     import needleman_wunsch
from pytadbit.boundary_aligner.reciprocally import reciprocal


def consensusize(ali1, ali2, passed):
    """
    :param ali1: first aligned sequence
    :param ali1: second aligned sequence
    :param passed: in case first aligned sequence is already a consensus, it
       might be weighted
    :returns: consensus sequence corresponding to ali1 and ali2
    """
    consensus = []
    for pos in xrange(len(ali1)):
        if ali1[pos] != ali2[pos]:
            try:
                bound = (ali1[pos] * passed + ali2[pos]) / (1 + passed)
            except TypeError:
                bound = ali1[pos] if type(ali1[pos]) is not str else ali2[pos]
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

    :param global method: method used to align.
    """
    if method == 'global':
        aligner = needleman_wunsch
    elif method == 'reciprocal':
        aligner = reciprocal
    else:
        raise NotImplementedError(('Only "global" and "reciprocal" are ' +
                                   'implemented right now.\n'))
    if len(sequences) > 2:
        dico = {}
        reference = None
        for j, (i, seq) in enumerate(sorted(enumerate(sequences),
                                            key=lambda x: x[1])):
            reference = reference or seq
            dico[j] = {'sort':i,
                       'seq' :seq}
        aligneds = []
        scores = 0
        for other in xrange(1, len(sequences)):
            [align1, align2], score = aligner(reference, dico[other]['seq'],
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
            reference = consensusize(align1, align2, other)
            
        sort_alis = [[] for _ in xrange(len(dico))]
        for seq in xrange(len(dico)):
            sort_alis[dico[seq]['sort']] = aligneds[seq][:]
        return sort_alis, scores
    return aligner(sequences[0], sequences[1], **kwargs)

