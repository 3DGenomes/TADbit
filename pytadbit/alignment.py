"""
14 May 2013


"""

from copy import deepcopy

class Alignment(object):
    def __init__(self, name, alignment, experiments, score=None):
        self.name = name
        self.__experiments = experiments
        self.__values = []        
        self.__keys = []
        self.__len = None
        for exp in alignment:
            self.add_aligned_exp(exp, alignment[exp])
        print self.__values
        self.score = score


    def __len__(self):
        return self.__len


    def __str__(self):
        """
        TODO: the crm.print_alignment
        """
        pass


    def __repr__(self):
        """
        TODO: description
        """
        pass


    def __getitem__(self):
        pass


    def __setitem__(self):
        pass


    def __delitem__(self):
        pass


    def __iter__(self):
        for key in self.__keys:
            yield key


    def get_column(self, cond1, cond2=None, min_num=None):
        """
        TODO: with score?
        """
        if type(cond1) is int:
            val = cond1
            cond1 = lambda x: x['pos'] == val
        if not cond2:
            cond2 = lambda x: True
        cond = lambda x: cond1(x) and cond2(x)
        min_num = min_num or len(self.__keys)
        column = []
        for pos in xrange(len(self)):
            col = []
            for exp in xrange(len(self.__keys)):
                if cond(self.__values[exp][pos]):
                    col.append(self.__values[exp][pos])
            if len(col) >= min_num:
                column.append(col)
        return column


    def add_aligned_exp(self, name, seq):
        p = 0
        scores = []
        exp = self.__experiments[name]
        for i, pos in enumerate(seq):
            if pos == '-':
                scores.append(TAD((('brk', None),
                                   ('end', 0.0),
                                   ('score', 0.0),
                                   ('start', 0.0)), i, name))
                continue
            try:
                while exp.tads[p]['brk'] < 0:
                    p += 1
            except KeyError:
                continue
            if exp.tads[p]['score'] < 0:
                exp.tads[p]['score'] = 10
            scores.append(TAD(exp.tads[p], i, name))
            p += 1
        if not self.__len:
            self.__len = len(scores)
        elif self.__len != len(scores):
            raise AssertionError('ERROR: alignments of different lengths\n')
        self.__values.append(scores)
        self.__keys.append(name)


class TAD(dict):
    def __init__(self, thing, i, exp):
        super(TAD, self).__init__(thing)
        self.update(dict((('pos', i),('exp', exp))))
        
    def __repr__(self):
        return '>' + str(int(self['end']))
