"""
26 Nov 2012


"""

from math import sqrt, log
from sys import stdout
from pytadbit.parsers.tad_parser import parse_tads
from pytadbit.tads_aligner.aligner import align
from pytadbit.parsers.hic_parser import read_matrix
from pytadbit import tadbit


try:
    from scipy.interpolate import interp1d
    from numpy import log2
    from matplotlib import pyplot as plt
except ImportError:
    from pytadbit.utils import Interpolate as interp1d

from random import random


class Chromosome():
    """
    Chromosome object designed to deal with Topologically associated domains
    predictions from different experiments, in different cell types and compare
    them.
    """
    def __init__(self, name, resolution, experiments=None, 
                 experiment_names=None,
                 max_tad_size=3000000, chr_len=None):
        """
        :argument name: name of the chromosome
        :argument resolution: resolution of the experiments. All experiments may have\
        the same resolution
        :argument None experiments: list of paths to files containing the definition\
        of TADs corresponding to different experiments (or output of tadbit)
        :argument None experiment_names: list of names for each experiment
        :argument 3000000 max_tad_size: maximum size of TAD allowed. TADs longer than\
        this will not be considered, and relative chromosome size will be reduced\
        accordingly
        :argument None chr_len: size of the chromosome in bp. By default it will be\
        inferred from the distribution of TADs.

        :return: Chromosome object
        """
        self.name = name
        self.max_tad_size = max_tad_size
        self.size = self.r_size = chr_len
        self.resolution = resolution
        self.forbidden = {}
        self.experiments = {}
        experiments = experiments or []
        for i, f_name in enumerate(experiments):
            name = experiment_names[i] if experiment_names else None
            if type(f_name) is dict:
                self.add_tad_def(f_name, name)
            else:
                self.add_experiment(f_name, name)


    def align_experiments(self, names=None, verbose=False, randomize=False,
                          **kwargs):
        """
        Align prediction of boundaries of two different experiments

        :argument None names: list of names of experiments to align. If None\
        align all.
        :argument experiment1: name of the first experiment to align
        :argument experiment2: name of the second experiment to align
        :argument -0.1 penalty: penalty of inserting a gap in the alignment
        :argument 500000 max_dist: Maximum distance between 2 boundaries allowing match
        :argument False verbose: print somethings
        :argument False randomize: check alignment quality by comparing randomization\
        of boundaries over chromosomes of same size. This will return a extra value,\
        the p-value of accepting that observed alignment is not better than random\
        alignment
        """
        experiments = names or self.experiments.keys()
        tads = []
        for e in experiments:
            if not self.experiments[e]['tads']:
                raise Exception('No TADs defined, use find_tad function.\n')
            tads.append(self.experiments[e]['brks'])
        aligneds, score = align(tads, bin_size=self.resolution,
                                max_dist=self.max_tad_size, **kwargs)
        for e, ali in zip(experiments, aligneds):
            self.experiments[e]['align'] = ali
            self.experiments[e]['align'] = ali
        if verbose:
            self.print_alignment(experiments)
        if not randomize:
            return self.get_alignment(names), score
        #mean, std = self._get_tads_mean_std(experiments)
        #print 'mean', mean, 'std', std, self.r_size, self.r_size/mean
        #p_value = randomization_test(len(experiments), mean, std, score,
        #                             self.r_size, self.resolution,
        #                             verbose=verbose, **kwargs)
        distr = self.interpolation(experiments)
        p_value = self.randomization_test(len(experiments), distr, score,
                                          verbose=verbose, **kwargs)
        return score, p_value


    def print_alignment(self, names=None, string=False):
        """
        print alignment
        :argument None names: if None print all experiments
        :argument False string: return string instead of printing
        """
        names = names or self.experiments.keys()
        length = max([len(n) for n in names])
        out = 'Alignment (%s TADs)\n' % (len(names))
        for e in names:
            if not 'align' in self.experiments[e]:
                continue
            out += '{1:{0}}:'.format(length, e)
            out += '|'.join(['%5s' % (str(x)[:-2] if x!='-' else '-' * 4)\
                             for x in self.experiments[e]['align']]) + '\n'
        if string:
            return out
        print out


    def get_alignment(self, names=None):
        """
        Return dictionary corresponding to alignment of experiments
        :argument None names: if None print all experiments
        :return: alignment as dict
        """
        names = names or self.experiments.keys()
        return dict([(e, self.experiments[e]['align']) \
                     for e in names if 'align' in self.experiments[e]])
                    

    def add_experiment(self, f_name, name, force=False):
        """
        Add Hi-C experiment to Chromosome
        """
        nums, size = read_matrix(f_name)
        if name in self.experiments:
            if 'hi-c' in self.experiments[name] and not force:
                raise Exception(\
                    '''Hi-C data already loaded under the name: {}.
                    Force loading or use an other name.\n'''.format(name))
            self.experiments[name]['hi-c'] = nums
            self.experiments[name]['size'] = size
        else:
            self.experiments[name] = {'hi-c': nums, 'size': size,
                                      'tads': None, 'brks': None,
                                      'wght': None}


    def find_tad(self, experiments, n_cpus=None, verbose=True,
                 max_tad_size="auto", no_heuristic=False):
        """
        Call tadbit function to calculate the position of Topologically associated
        domains
        
        :argument experiment: A square matrix of interaction counts in hi-C data or a list of\
        such matrices for replicated experiments. The counts must be evenly sampled\
        and not normalized.\
        'experiment' might be either a list of list, a path to a file or a file handler
        :argument None n_cpus: The number of CPUs to allocate to tadbit. The value default\
        is the total number of CPUs minus 1.
        :argument auto max_tad_size: an integer defining maximum size of TAD.\
        Default (auto) defines it to the number of rows/columns.
        :argument False no_heuristic: whether to use or not some heuristics
        
        """
        for experiment in experiments:
            result, weights = tadbit(self.experiments[experiment]['hi-c'],
                                     n_cpus=n_cpus, verbose=verbose,
                                     max_tad_size=max_tad_size,
                                     no_heuristic=no_heuristic,
                                     get_weights=True)
            self.add_tad_def(result, name=experiment, weights=weights)
        

    def visualize(self, name, tad=None, paint_tads=False, ax=None):
        """
        Visualize the matrix of Hi-C interactions

        :argument name: name of the experiment to visualize
        
        """
        vmin = log2(min(self.experiments[name]['hi-c'][0]) or 1)
        vmax = log2(max(self.experiments[name]['hi-c'][0]))
        size = self.experiments[name]['size']
        if not ax:
            ax = plt.subplot(111)
        if tad:
            matrix = [[self.experiments[name]['hi-c'][0][i+size*j] \
                       for i in xrange(int(tad['start']), int(tad['end']))] \
                      for j in xrange(int(tad['start']), int(tad['end']))]
        else:
            matrix = [[self.experiments[name]['hi-c'][0][i+size*j] \
                       for i in xrange(size)] \
                      for j in xrange(size)]
        ax.imshow(log2(matrix), origin='lower', vmin=vmin, vmax=vmax,
                  interpolation="nearest")
        if not paint_tads:            
            return
        for i, tad in self.experiments[name]['tads'].iteritems():
            ax.hlines(tad['start'], tad['start'], tad['end'], colors='k')
            ax.hlines(tad['end'], tad['start'], tad['end'], colors='k')
            ax.vlines(tad['start'], tad['start'], tad['end'], colors='k')
            ax.vlines(tad['end'], tad['start'], tad['end'], colors='k')
            ax.text(tad['start'] + abs(tad['start']-tad['end'])/2 - 1,
                    tad['start'] + abs(tad['start']-tad['end'])/2 - 1, str(i))

    def add_tad_def(self, f_name, name=None, weights=None):
        """
        Add Topologically Associated Domains defintinion detection to chromosome

        :argument f_name: path to file
        :argument None name: name of the experiment, if None f_name will be used:

        """
        name = name or f_name
        tads, forbidden = parse_tads(f_name, max_size=self.max_tad_size,
                                     bin_size=self.resolution)
        brks = [t['brk'] for t in tads.values() if t['brk']]
        if not name in self.experiments:
            self.experiments[name] = {'hi-c': None, 'size': None,
                                      'tads': None, 'brks': None,
                                      'wght': None}
        self.experiments[name]['tads'] = tads
        self.experiments[name]['brks'] = brks
        if weights:
            self.experiments[name]['wght'] = weights
        if not self.forbidden:
            self.forbidden = dict([(f, None) for f in forbidden])
        else:
            self.forbidden = dict([(f, None) for f in \
                                   forbidden.intersection(self.forbidden)])
        if not self.size:
            self.size = tads[max(tads)]['end'] * self.resolution
        self.r_size = self.size - len(self.forbidden) * self.resolution


    def interpolation(self, experiments):
        """
        calculate the distribution of TAD lengths, and interpolate it using\
        interp1d function from scipy.
        :arguments experiments: names of experiments included in the\
        distribution
        :return: function to interpolate a given TAD length according to a\
        probability value
        """
        # get all TAD lengths and multiply it by bin size of the experiment
        norm_tads = []
        for tad in experiments:
            for brk in self.experiments[tad]['tads'].values():
                if not brk['brk']:
                    continue
                norm_tads.append((brk['end'] - brk['start']) * self.resolution)
        x = [0.0]
        y = [0.0]
        max_d = max(norm_tads)
        # we ask here for a mean number of 20 values per bin 
        bin_s = int(max_d / (float(len(norm_tads)) / 20))
        for b in range(0, int(max_d) + bin_s, bin_s):
            x.append(len([i for i in norm_tads if b < i <= b + bin_s]) + x[-1])
            y.append(b + float(bin_s))
        x = [float(v) / max(x) for v in x]
        ## to see the distribution and its interpolation
        #distr = interp1d(x, y)
        #from matplotlib import pyplot as plt
        #plt.plot([distr(float(i)/1000) for i in xrange(1000)],
        #         [float(i)/1000 for i in xrange(1000)])
        #plt.hist(norm_tads, normed=True, bins=20, cumulative=True)
        #plt.show()
        return interp1d(x, y)
        

    def _get_tads_mean_std(self, experiments):
        """
        returns mean and standard deviation of TAD lengths. Value is for both
        TADs.

        Note: no longer used in core functions
        """
        norm_tads = []
        for tad in experiments:
            for brk in self.experiments[tad]['tads'].values():
                if not brk['brk']:
                    continue
                norm_tads.append(log((brk['end'] - brk['start']) * self.resolution))
        length = len(norm_tads)
        mean   = sum(norm_tads)/length
        std    = sqrt(sum([(t-mean)**2 for t in norm_tads])/length)
        return mean, std


    def randomization_test(self, num_sequences, distr, score=None,
                           num=1000, verbose=False):
        """
        Return the probability that original alignment is better than an
        alignment of randomized boundaries.
        :argument num_sequences: number of sequences aligned
        :argument distr: the function to interpolate TAD lengths from\
        probability
        :argument None score: just to print it when verbose
        :argument 1000 num: number of random alignment to generate for\
        comparison
        :argument False verbose: to print something nice
        """
        rnd_distr = []
        rnd_len = []
        for val in xrange(num):
            if verbose:
                val = float(val)
                if not val / num * 100 % 5:
                    stdout.write('\r' + ' ' * 10 + \
                                 ' randomizing: '
                                 '{:.2%} completed'.format(val/num))
                    stdout.flush()
            rnd_tads = [generate_rnd_tads(self.r_size, distr, self.resolution) \
                           for _ in xrange(num_sequences)]
            rnd_len.append(float(sum([len(r) for r in rnd_tads])) \
                           / len(rnd_tads))
            rnd_distr.append(align(rnd_tads,
                                    bin_size=self.resolution,
                                    verbose=False)[1])
        pval = float(len([n for n in rnd_distr if n > score])) / len(rnd_distr)
        if verbose:
            stdout.write('\n {} randomizations finished.'.format(num))
            stdout.flush()
            print '  Observed alignment score: {}'.format(score)
            print '  Mean number of boundaries: {}; observed: {}'.format(\
                sum(rnd_len)/len(rnd_len),
                str([len(self.experiments[e]['brks']) \
                     for e in self.experiments]))
            print 'Randomized scores between {} and {}; observed: {}'.format(\
                min(rnd_distr), max(rnd_distr), score)
            print 'p-value: {}'.format(pval if pval else '<{}'.format(1./num))
        return pval


    def get_tad_hic(self, tad, x_name, normed=True):
        """
        Retrieve the Hi-C data matrix corresponding to ma given TAD.
        
        :argument tad: a given TAD -> dict
        :argument x_name: name og the experiment
        :argument True normed: if Hi-C data has to be normalized
        
        :returns: Hi-C data matrix for the given TAD
        """
        beg, end = int(tad['start']), int(tad['end'])
        xpr = self.experiments[x_name]
        size = xpr['size']
        matrix = [[[] for _ in xrange(beg, end)]\
                  for _ in xrange(beg, end)]
        for i, ii in enumerate(xrange(beg - 1, end - 1)):
            ii = ii * size
            for j, jj in enumerate(xrange(beg, end)):
                matrix[j][i] = float(xpr['hi-c'][0][ii + jj])
                if not normed:
                    continue
                try:
                    matrix[j][i] = matrix[j][i] / xpr['wght'][0][ii + jj]
                except ZeroDivisionError:
                    matrix[j][i] = 0.0
        return matrix


    def iter_tads(self, x_name, normed=True):
        """
        Iterate over TADs corresponding to a given experiment.
        
        :argument x_name: name of the experiment
        :argument True normed: normalize Hi-C data returned
        
        :yields: Hi-C data corresponding to each TAD
        """
        if not self.experiments[x_name]['hi-c']:
            raise Exception('No Hi-c data for {} experiment\n'.format(x_name))
        for tad in self.experiments[x_name]['tads']:
            yield self.get_tad_hic(self.experiments[x_name]['tads'][tad],
                                   x_name, normed=normed)


def generate_rnd_tads(chr_len, distr, bin_size, start=0):
    """
    Generates random TADs over a chromosome of a given size according to \
    a given distribution of lengths of TADs.
    
    :argument chr_len: length of the chromosome
    :argument distr: function that returns a TAD length depending on a p value
    :argument bin_size: size of the bin of the Hi-C experiment
    :argument 0 start: starting position in the chromosome
    
    :returns: list of TADs
    """
    pos = start
    tads = []
    while True:
        pos += distr(random())
        if pos > chr_len:
            break
        tads.append(float(int(pos / bin_size + .5)))
    return tads


def randomization_test_old(num_sequences, mean, std, score, chr_len, bin_size,
                       num=1000, verbose=False):
    """
    No longer used
    
    TODO: use Guillaume idea.
    """
    rand_distr = []
    rand_len = []
    best = (None, None)
    for n in xrange(num):
        if verbose:
            n = float(n)
            if not n / num * 100 % 5:
                stdout.write('\r' + ' ' * 10 + \
                             ' randomizing: {:.2%} completed'.format(n/num))
                stdout.flush()
        random_tads = [generate_rnd_tads(chr_len, mean,
                                            std, bin_size) \
                       for _ in xrange(num_sequences)]
        rand_len.append(float(sum([len(r) for r in random_tads]))/len(random_tads))
        rand_distr.append(align(random_tads,
                                bin_size=bin_size, chr_len=chr_len,
                                verbose=False)[1])
        if rand_distr[-1] > best[0]:
            best = rand_distr[-1], random_tads
    p_value = float(len([n for n in rand_distr if n > score]))/len(rand_distr)
    if verbose:
        stdout.write('\n {} randomizations finished.'.format(num))
        stdout.flush()
        align(best[-1], bin_size=bin_size, chr_len=chr_len, verbose=True)
        print 'Observed alignment score: {}'.format(score)
        print '  Randomized scores between {} and {}'.format(min(rand_distr),
                                                             max(rand_distr))
        print 'p-value: {}'.format(p_value if p_value else '<{}'.format(1./num))
        print sum(rand_len)/len(rand_len)
    return p_value    

