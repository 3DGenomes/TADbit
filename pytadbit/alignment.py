"""
14 May 2013


"""

from sys import stdout
from pytadbit.utils import colorize


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
        calls Alignment.print_alignment
        """
        self.print_alignment(string=False)


    def __repr__(self):
        """
        Alignment description
        """
        print ('Alignment of boundaries (length: %s, ' +
               'number of experiments: %s)') % (len(self), len(self.__keys))


    def __getitem__(self, i):
        try:
            return self.__values[i]
        except TypeError:
            try:
                i = self.__keys.index(i)
                return self[i]
            except ValueError:
                raise ValueError('ERROR: %s not in alignment' % (i))


    def __setitem__(self, i, j):
        if i in self.__keys:
            self.__values[self.__keys[i]] = j
        else:
            self.__values.append(j)
            self.__keys.append(i)


    def __delitem__(self, i):
        try:
            self.__values.pop(i)
            self.__keys.pop(i)
        except TypeError:
            try:
                i = self.__keys.index(i)
                self.__values.pop(i)
                self.__keys.pop(i)
            except ValueError:
                raise ValueError('ERROR: %s not in alignment' % (i))


    def __iter__(self):
        for key in self.__keys:
            yield key


    def iteritems(self):
        """
        Iterate over experiment names and aligned boundaries
        """
        for i in xrange(len(self)):
            yield (self.__keys[i], self.__values[i])


    def itervalues(self):
        """
        Iterate over experiment names and aligned boundaries
        """
        for i in xrange(len(self)):
            yield self.__values[i]


    def itercolumns(self):
        """
        Iterate over columns in the alignment
        """
        for pos in xrange(len(self)):
            col = []
            for exp in xrange(len(self.__keys)):
                col.append(self[exp][pos])
            yield col


    def print_alignment(self, name=None, xpers=None, string=False,
                        title='', ftype='ansi'):
        """
        Print alignment of TAD boundaries between different experiments.
           Alignment are displayed with colors according to the tadbit
           confidence score for each boundary.
        
        :param None names: if None print all experiments
        :param None xpers: if None print all experiments
        :param False string: return string instead of printing
        :param ansi ftype: display colors in 'ansi' or 'html' format
        """
        if xpers:
            xpers = [self.__experiments(n.name) for n in xpers]
        else:
            xpers = self.__experiments
        if not name:
            name = self.__keys[0]
        length = max([len(n.name) for n in xpers])
        if ftype == 'html':
            out = '''<?xml version="1.0" encoding="UTF-8" ?>
            <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
            <!-- This file was created with the aha Ansi HTML Adapter. http://ziz.delphigl.com/tool_aha.php -->
            <html xmlns="http://www.w3.org/1999/xhtml">
            <head>
            <meta http-equiv="Content-Type" content="application/xml+xhtml; charset=UTF-8" />
            <title>stdin</title>
            </head>
            <h1>{}</h1>
            <body>
            <pre>'''.format(title)
        elif ftype == 'ansi':
            out = title + '\n' if title else ''
        else:
            raise NotImplementedError('Only ansi and html ftype implemented.\n')
        out += 'Alignment shown in Kb (%s experiments) (' % (len(xpers))
        out += 'scores: {})\n'.format(' '.join(
            [colorize(x, x, ftype) for x in range(11)]))
        for xpr in xpers:
            if not xpr.name in self.__keys:
                continue
            res = xpr.resolution / 1000
            out += '{1:{0}}:'.format(length, xpr.name)
            i = 0 # we are working with the end position
            for x in self[xpr.name]:
                if x['end'] == 0.0:
                    out += '| ' + '-' * 4 + ' '
                    continue
                try:
                    while xpr.tads[i]['brk'] < 0:
                        i += 1
                except KeyError:
                    continue
                cell = str(int(x['end'] * res))
                out += ('|' + ' ' * (6 - len(cell)) +
                        colorize(cell, xpr.tads[i]['score'], ftype))
                i += 1
            out += '\n'
        if ftype == 'html':
            out += '</pre></body></html>'
        if string:
            return out
        print out



    def get_column(self, cond1, cond2=None, min_num=None):
        """
        Get a list of column responding to a given characteristic.

        :param cond1: can be either a column number or a function to search for
           a given condition that may fulfill boundaries.
        :param None cond2: a second function to search for a given condition
           that may fulfill boundaries.
        :param None min_num: Number of boundaries in a given column that have to
           pass the defined conditions (all boundaries may pass them by default).

        :returns: a list of column, each column represented by a tuple which
           first element is the column number in the alignment and second
           element is the list of boundaries in this column.

        Example
        -------

        ::

          ali = Alignment('ALL', crm.get_alignment('ALL'),
                          crm.experiments, score)
          ali.get_column(3)
          # [(3, [>55<, ,>56<, >55<, >58<])]

          # now we want boundaries with high scores
          cond1 = lambda x: x['score'] > 5
          # and we want boundaries to be detected in Experiments exp1 and exp2
          cond2=lambda x: x['exp'] in ['exp1', 'exp2']
          
          ali.get_column(cond1=cond1, cond2=cond2, min_num=2)
          # [(33, [>268<, >192<, >-<, >-<]),	 
          #  (46, [>324<, >323<, >335<, >329<]), 
          #  (51, [>348<, >335<, >357<, >347<]), 
          #  (56, [>374<, >358<, >383<, >-<]),	 
          #  (64, [>397<, >396<, >407<, >-<]),	 
          #  (77, [>444<, >442<, >456<, >-<])]   
          
          
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
            cnt = 0
            for exp in xrange(len(self.__keys)):
                col.append(self.__values[exp][pos])
                if cond(self.__values[exp][pos]):
                    cnt += 1
            if cnt >= min_num:
                column.append((pos, col))
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
        return '>' + (str(int(self['end'])) if int(self['end']) else '-') + '<'
