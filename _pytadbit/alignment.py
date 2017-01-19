"""
14 May 2013


"""

from pytadbit.utils.extraviews         import colorize, tadbit_savefig
from pytadbit.utils.extraviews         import _tad_density_plot
from random                            import random, shuffle
from sys                               import stdout
from pytadbit.boundary_aligner.aligner import align


try:
    from scipy.interpolate import interp1d
except ImportError:
    from pytadbit.utils.tadmaths import Interpolate as interp1d


class Alignment(object):
    """
    Alignment of TAD borders
    """
    def __init__(self, name, alignment, experiments, consensus, score=None):
        self.name = name
        self.__experiments = experiments
        self.__values = []        
        self.__keys = []
        self.__len = None
        self.consensus = consensus
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


    def write_alignment(self, names=None, xpers=None, string=False,
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
            xpers = [self.__experiments[n.name] for n in xpers]
        else:
            xpers = self.__experiments
        if not names:
            names = [n.name for n in xpers]
        length = max([len(n) for n in names])
        if ftype == 'html':
            out = '''<?xml version="1.0" encoding="UTF-8" ?>
            <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
            <!-- This file was created with the aha Ansi HTML Adapter. http://ziz.delphigl.com/tool_aha.php -->
            <html xmlns="http://www.w3.org/1999/xhtml">
            <head>
            <meta http-equiv="Content-Type" content="application/xml+xhtml; charset=UTF-8" />
            <title>stdin</title>
            </head>
            <h1>%s</h1>
            <body>
            <pre>''' % (title)
        elif ftype == 'ansi':
            out = title + '\n' if title else ''
        else:
            raise NotImplementedError('Only ansi and html ftype implemented.\n')
        out += 'Alignment shown in %s Kb (%s experiments) (' % (
            int(xpers[0].resolution / 1000), len(xpers))
        out += 'scores: %s)\n' % (' '.join(
            [colorize(x, x, ftype) for x in range(11)]))
        for i, xpr in enumerate(xpers):
            if not xpr.name in self.__keys:
                continue
            # res = xpr.resolution / 1000
            out += ('%' + str(length) + 's') % (names[i])
            out += ':'
            for x in self[xpr.name]:
                if x['end'] == 0.0:
                    out += '| ' + '-' * 4 + ' '
                    continue
                cell = str(int(x['end']) + 1) # * res # TODO: +1
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
        if isinstance(cond1, int):
            val = cond1 - 1
            cond1 = lambda x: x['pos'] == val
        elif isinstance(cond1, list):
            val = [v - 1 for v in cond1]
            cond1 = lambda x: x['pos'] in val
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
        p = 1
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
            scr = abs(exp.tads[p]['score'])
            scores.append(TAD(exp.tads[p], i, self.__experiments[name]))
            scores[-1]['score'] = scr
            p += 1
        # print name, len(scores)
        if not self.__len:
            self.__len = len(scores)
        elif self.__len != len(scores):
            raise AssertionError('ERROR: alignments of different lengths ' +
                                 '(%s != %s)\n' % (self.__len, len(scores)))
        self.__values.append(scores)
        self.__keys.append(name)


    def draw(self, focus=None, extras=None, ymax=None, ali_colors=('grey',),
             normalized=True, savefig=None, shape='ellipse'):
        """
        Draw alignments as a plot.
        
        :param None focus: can pass a tuple (bin_start, bin_stop) to display the
           alignment between these genomic bins
        :param None extras: list of coordinates (genomic bin) where to draw a
           red cross
        :param None ymax: limit the y axis up to a given value
        :param ('grey', ): successive colors for alignment
        :param True normalized: normalized Hi-C count are plotted instead of raw
           data.
        :param 'ellipse' shape: which kind of shape to use as schematic
           representation of TADs. Implemented: 'ellipse', 'rectangle',
           'triangle'
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        """
        from matplotlib.cm import jet
        from matplotlib import pyplot as plt
        experiments = self.__experiments

        maxres = max([e.resolution for e in experiments])
        facts = [maxres / e.resolution for e in experiments]
        

        siz = experiments[0].size
        if focus:
            figsiz = 4 + (focus[1] - focus[0]) / 30
        else:
            figsiz = 4 + siz / 30
        fig, axes = plt.subplots(nrows=len(experiments),
                                 sharex=True, sharey=True,
                                 figsize=(figsiz, 1 + len(experiments) * 1.8))
        fig.subplots_adjust(hspace=0)

        maxys = []
        for iex, xpr in enumerate(experiments):
            if not xpr.name in self:
                continue
            _tad_density_plot(xpr, maxys=maxys, normalized=normalized,
                              fact_res=facts[iex], axe=axes[iex],
                              extras=extras, shape=shape, focus=focus)
        # draw alignment columns
        start = focus[0] if focus else 1
        end = focus[1] if focus else xpr.tads[max(xpr.tads)]['end']
        maxy = (ymax or max(maxys)) + 0.4
        maxxs = []
        for iex in range(len(experiments)):
            starting = focus[0] if focus else 1
            ending = (focus[1] if focus
                      else experiments[iex].tads.values()[-1]['end'])
            axes[iex].hlines(1, 1, end, 'k', lw=1.5)
            axes[iex].set_ylim((0, maxy))
            maxxs.append(ending / facts[iex])
            axes[iex].text(starting + 1, float(maxy) / 20,
                           experiments[iex].name, {'ha': 'left',
                                                   'va': 'bottom'})
            axes[iex].set_yticks([float(i) / 2
                                  for i in range(1, int(maxy + .5) * 2)])
            if ymax:
                axes[iex].set_ylim((0, ymax))
            axes[iex].set_xlim(xmin=starting, xmax=max(maxxs))

        pos = {'ha': 'center', 'va': 'bottom'}
        for i, col in enumerate(self.itercolumns()):
            ends = sorted([(t['end'], j) for j, t in enumerate(col) if t['end']])
            beg = (ends[0][0] + 0.9) / facts[ends[0][1]]
            end = (ends[-1][0] + 1.1) / facts[ends[-1][1]]
            if focus:
                if beg < focus[0] or end > focus[1]:
                    continue
            axes[0].text(beg + float(end - beg) / 2, maxy + float(maxy) / 20,
                         str(i + 1), pos,
                         rotation=90, size='small')
            for iex, tad in enumerate(col):
                if not tad['end']:
                    continue
                axes[iex].axvspan(beg-.2, end+.2, alpha=0.2,
                                  color=ali_colors[i%(len(ali_colors))])
        axes[iex].set_xlabel('Genomic bin')
        tit1 = fig.suptitle("TAD borders' alignment", size='x-large')
        tit2 = axes[0].set_title("Alignment column number")
        tit2.set_y(1.3)
        plt.subplots_adjust(top=0.76)
        # This was for color bar instead of legend
        # ax1 = fig.add_axes([0.9 + 0.3/figsiz, 0.05, 0.2/figsiz, 0.9])
        # cb1 = colorbar.ColorbarBase(ax1, cmap=jet,
        #                             norm=colors.Normalize(vmin=0., vmax=1.))
        # cb1.set_label('Border prediction score')
        # cb1.ax.set_yticklabels([str(i)for i in range(1, 11)])
        fig.set_facecolor('white')
        plots = []
        for scr in xrange(1, 11):
            plots += plt.plot((100,),(100,), marker=6, ms=9,
                              color=jet(float(scr) / 10), mec='none')
        try:
            axes[-1].legend(plots,
                            [str(scr) for scr in xrange(1, 11)],
                            numpoints=1, title='Border scores',
                            fontsize='small', loc='lower left',
                            bbox_to_anchor=(1, 0.5))
        except TypeError:
            axes[-1].legend(plots,
                            [str(scr) for scr in xrange(1, 11)],
                            numpoints=1, title='Border scores',
                            loc='lower left',
                            bbox_to_anchor=(1, 0.5))
        if savefig:
            tadbit_savefig(savefig)
        else:
            plt.show()


class TAD(dict):
    """
    Specific class of TADs, used only within Alignment objects.
    It is directly inheriting from python dict.
    a TAD these keys:

     - 'start': position of the TAD
     - 'end': position of the TAD
     - 'score': of the prediction of boundary
     - 'brk': same as 'end'
     - 'pos': in the alignment (column number)
     - 'exp': Experiment this TAD belongs to
     - 'index': of this TAD within all TADs in the Experiment
    
    """
    def __init__(self, thing, i, exp):
        super(TAD, self).__init__(thing)
        idx = [t for t in exp.tads if exp.tads[t]['start']==self['start']][0]
        self.update(dict((('pos', i),('exp', exp), ('index', idx))))
        
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


def randomization_test(xpers, score=None, num=1000, verbose=False, max_dist=100000,
                       rnd_method='interpolate', r_size=None, method='reciprocal'):
    """
    Return the probability that original alignment is better than an
    alignment of randomized boundaries.

    :param tads: original TADs of each experiment to align
    :param distr: the function to interpolate TAD lengths from probability
    :param None score: just to print it when verbose
    :param 1000 num: number of random alignment to generate for comparison
    :param False verbose: to print something nice
    :param interpolate method: how to generate random tads (alternative is
       'shuffle'). 'interpolate' will calculate the distribution of TAD lengths,
       and generate a random set of TADs according to this distribution (see
       :func:`pytadbit.alignment.generate_rnd_tads`). In contrast, the 'shuffle'
       method uses directly the set of observed TADs and shuffle them (see
       :func:`pytadbit.alignment.generate_shuffle_tads`).
    """
    if not rnd_method in ['interpolate', 'shuffle']:
        raise Exception('method should be either "interpolate" or ' +
                        '"shuffle"\n')
    if rnd_method == 'interpolate' and not r_size:
        raise Exception('should provide Chromosome r_size if interpolate\n')
    tads = []
    for xpr in xpers:
        if not xpr.tads:
            raise Exception('No TADs defined, use find_tad function.\n')
        tads.append([(t['end'] - t['start']) * \
                     xpr.resolution for t in xpr.tads.values()])
    rnd_distr = []
    # rnd_len = []
    distr = _interpolation(xpers) if rnd_method is 'interpolate' else None
    rnd_exp = lambda : tads[int(random() * len(tads))]
    for val in xrange(num):
        if verbose:
            val = float(val)
            if not val / num * 100 % 5:
                stdout.write('\r' + ' ' * 10 + 
                             ' randomizing: '
                             '%.2f completed' % (100 * val/num))
                stdout.flush()
        if rnd_method is 'interpolate':
            rnd_tads = [generate_rnd_tads(r_size, distr)
                        for _ in xrange(len(tads))]
            # rnd_len.append(float(sum([len(r) for r in rnd_tads]))
            #                / len(rnd_tads))
        else:
            rnd_tads = [generate_shuffle_tads(rnd_exp())
                        for _ in xrange(len(tads))]
            # rnd_len.append(len(tads))
        rnd_distr.append(align(rnd_tads, verbose=False, method=method,
                               max_dist=max_dist)[0][1])
        # aligns, sc = align(rnd_tads, verbose=False)
        # rnd_distr.append(sc)
        # for xpr in aligns:
        #     print sc, '|'.join(['%5s' % (str(x/1000)[:-2] \
        # if x!='-' else '-' * 4)\
        #                         for x in xpr])
        # print ''
    pval = float(len([n for n in rnd_distr if n > score])) / len(rnd_distr)
    if verbose:
        stdout.write('\n %s randomizations finished.' % (num))
        stdout.flush()
        print '  Observed alignment score: %s' % (score)
        # print '  Mean number of boundaries: {}; observed: {}'.format (
        #     sum(rnd_len)/len(rnd_len),
        #     str([len(self.experiments[e].brks)
        #          for e in self.experiments]))
        print 'Randomized scores between %s and %s; observed: %s' % (
            min(rnd_distr), max(rnd_distr), score)
        print 'p-value: %s' % (pval if pval else '<%s' % (1./num))
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
