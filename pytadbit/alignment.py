"""
14 May 2013


"""

from pytadbit.utils import colorize
from random import random, shuffle
from sys import stdout
from pytadbit.boundary_aligner.aligner import align
try:
    from scipy.interpolate import interp1d
except ImportError:
    from pytadbit.utils import Interpolate as interp1d


class Alignment(object):
    def __init__(self, name, alignment, experiments, score=None):
        self.name = name
        self.__experiments = experiments
        self.__values = []        
        self.__keys = []
        self.__len = None
        for seq, exp in zip(alignment, experiments):
            self.add_aligned_exp(exp.name, seq)
        self.score = score


    def __len__(self):
        return self.__len


    def __str__(self):
        """
        calls Alignment.print_alignment
        """
        return self.write_alignment(string=True)


    def __repr__(self):
        """
        Alignment description
        """
        return ('Alignment of boundaries (length: %s, ' +
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


    def write_alignment(self, name=None, xpers=None, string=False,
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
            for x in self[xpr.name]:
                if x['end'] == 0.0:
                    out += '| ' + '-' * 4 + ' '
                    continue
                cell = str(int(x['end'] * res))
                out += ('|' + ' ' * (6 - len(cell)) +
                        colorize(cell, x['score'], ftype))
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

        **Example**

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
                                   ('start', 0.0)), i,
                                  self.__experiments[name]))
                continue
            try:
                while exp.tads[p]['brk'] < 0:
                    p += 1
            except KeyError:
                continue
            scr = exp.tads[p]['score'] if exp.tads[p]['score'] >= 0 else 10
            scores.append(TAD(exp.tads[p], i, self.__experiments[name]))
            scores[-1]['score'] = scr
            p += 1
        # print name, len(scores)
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
        return '>' + (str(int(self['end']) * self['exp'].resolution / 1000) \
                      if int(self['end']) else '-') + '<'


def _interpolation(experiments):
    """
    Calculate the distribution of TAD lengths, and interpolate it using
    interp1d function from scipy.

    :param experiments: names of experiments included in the
        distribution
    :return: function to interpolate a given TAD length according to a
        probability value
    """
    # get all TAD lengths and multiply it by bin size of the experiment
    norm_tads = []
    for tad in experiments:
        for brk in tad.tads.values():
            if not brk['brk']:
                continue
            norm_tads.append((brk['end'] - brk['start']) * tad.resolution)
    win = [0.0]
    cnt = [0.0]
    max_d = max(norm_tads)
    # we ask here for a mean number of 20 values per bin 
    bin_s = int(max_d / (float(len(norm_tads)) / 20))
    for bee in range(0, int(max_d) + bin_s, bin_s):
        win.append(len([i for i in norm_tads
                        if bee < i <= bee + bin_s]) + win[-1])
        cnt.append(bee + float(bin_s))
    win = [float(v) / max(win) for v in win]
    ## to see the distribution and its interpolation
    #distr = interp1d(x, y)
    #from matplotlib import pyplot as plt
    #plt.plot([distr(float(i)/1000) for i in xrange(1000)],
    #         [float(i)/1000 for i in xrange(1000)])
    #plt.hist(norm_tads, normed=True, bins=20, cumulative=True)
    #plt.show()
    return interp1d(win, cnt)


def randomization_test(xpers, score=None, num=1000, verbose=False,
                       method='interpolate', r_size=None):
    """
    Return the probability that original alignment is better than an
    alignment of randomized boundaries.

    :param tads: original TADs of each experiment to align
    :param distr: the function to interpolate TAD lengths from probability
    :param None score: just to print it when verbose
    :param 1000 num: number of random alignment to generate for comparison
    :param False verbose: to print something nice
    :param interpolate method: how to generate random tads (alternative is
       'shuffle').
    """
    if not method in ['interpolate', 'shuffle']:
        raise Exception('method should be either "interpolate" or ' +
                        '"shuffle"\n')
    if method == 'interpolate' and not r_size:
        raise Exception('should provide Chromosome r_size if interpolate\n')
    tads = []
    for xpr in xpers:
        if not xpr.tads:
            raise Exception('No TADs defined, use find_tad function.\n')
        tads.append([(t['end'] - t['start']) * \
                     xpr.resolution for t in xpr.tads.values()])
    rnd_distr = []
    # rnd_len = []
    distr = _interpolation(xpers) if method is 'interpolate' else None
    rnd_exp = lambda : tads[int(random() * len(tads))]
    for val in xrange(num):
        if verbose:
            val = float(val)
            if not val / num * 100 % 5:
                stdout.write('\r' + ' ' * 10 + 
                             ' randomizing: '
                             '{:.2%} completed'.format(val/num))
                stdout.flush()
        if method is 'interpolate':
            rnd_tads = [generate_rnd_tads(r_size, distr)
                        for _ in xrange(len(tads))]
            # rnd_len.append(float(sum([len(r) for r in rnd_tads]))
            #                / len(rnd_tads))
        else:
            rnd_tads = [generate_shuffle_tads(rnd_exp())
                        for _ in xrange(len(tads))]
            # rnd_len.append(len(tads))
        rnd_distr.append(align(rnd_tads, verbose=False)[1])
        # aligns, sc = align(rnd_tads, verbose=False)
        # rnd_distr.append(sc)
        # for xpr in aligns:
        #     print sc, '|'.join(['%5s' % (str(x/1000)[:-2] \
        # if x!='-' else '-' * 4)\
        #                         for x in xpr])
        # print ''
    pval = float(len([n for n in rnd_distr if n > score])) / len(rnd_distr)
    if verbose:
        stdout.write('\n {} randomizations finished.'.format(num))
        stdout.flush()
        print '  Observed alignment score: {}'.format(score)
        # print '  Mean number of boundaries: {}; observed: {}'.format(
        #     sum(rnd_len)/len(rnd_len),
        #     str([len(self.experiments[e].brks)
        #          for e in self.experiments]))
        print 'Randomized scores between {} and {}; observed: {}'.format(
            min(rnd_distr), max(rnd_distr), score)
        print 'p-value: {}'.format(pval if pval else '<{}'.format(1./num))
    return pval


def generate_rnd_tads(chromosome_len, distr, start=0):
    """
    Generates random TADs over a chromosome of a given size according to a given
    distribution of lengths of TADs.
    
    :param chromosome_len: length of the chromosome
    :param distr: function that returns a TAD length depending on a p value
    :param bin_size: size of the bin of the Hi-C experiment
    :param 0 start: starting position in the chromosome
    
    :returns: list of TADs
    """
    pos = start
    tads = []
    while True:
        pos += distr(random())
        if pos > chromosome_len:
            break
        tads.append(float(pos))
    return tads


def generate_shuffle_tads(tads):
    """
    Returns a shuffle version of a given list of TADs

    :param tads: list of TADs

    :returns: list of shuffled TADs
    """
    rnd_tads = tads[:]
    shuffle(rnd_tads)
    tads = []
    for tad in rnd_tads:
        if tads:
            tad += tads[-1]
        tads.append(tad)
    return tads
