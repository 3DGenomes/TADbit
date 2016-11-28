"""
02 Dec 2012


"""
from pytadbit.boundary_aligner.globally     import needleman_wunsch
from pytadbit.boundary_aligner.reciprocally import reciprocal


def consensusize(ali1, ali2, passed):
    """
    Given two alignments returns a consensus alignment. Used for the generation
    of multiple alignments
    
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
                bound = ali1[pos] if not isinstance(ali1[pos],
                                                    str) else ali2[pos]
        else:
            bound = ali1[pos]
        consensus.append(bound)
    return consensus


def align(sequences, method='reciprocal', **kwargs):
    """
    Align Topologically Associating Domain borders. Supports multiple alignment
    by building a consensus TAD sequence and aligning each experiment to it.

    .. note::

      as long as we are using multiple alignments in an iterative way,
      the order of sequences will be relevant. Here experiments are sorted
      according to the value of the first boundary found in order to try to
      reduce this problem.

    :param reciprocal method: method used to align
    :returns: the result of the aligner used
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
        consensus = None
        for j, (i, seq) in enumerate(sorted(enumerate(sequences),
                                            key=lambda x: x[1])):
            consensus = consensus or seq
            dico[j] = {'sort':i,
                       'seq' :seq}
        aligneds = []
        scores = 0
        perc1 = 0
        perc2 = 0
        for other in xrange(1, len(sequences)):
            try:
                ([align1, align2], score,
                 p1, p2) = aligner(consensus, dico[other]['seq'], **kwargs)
                perc1 += p1
                perc2 += p2
            except ValueError:
                [align1, align2], score = aligner(consensus, dico[other]['seq'],
                                                  **kwargs)
            scores += score
            if len(consensus) != len(align1):
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
            consensus = consensusize(align1, align2, other)
            
        sort_alis = [[] for _ in xrange(len(dico))]
        for seq in xrange(len(dico)):
            sort_alis[dico[seq]['sort']] = aligneds[seq][:]
        return (sort_alis, scores,
                perc1 / (len(sequences) - 1.),
                perc2 / (len(sequences) - 1.)), consensus
    ([align1, align2], score, p1, p2) = aligner(sequences[0], sequences[1], **kwargs) 
    consensus = consensusize(align1, align2, sequences[1])
    return ([align1, align2], score, p1, p2), consensus

