"""
26 Nov 2012


"""

from math import sqrt, log
from sys import stdout
from pytadbit.parsers.tad_parser import parse_tads
from pytadbit.tads_aligner.aligner import align
from pytadbit.parsers.hic_parser import read_matrix
from pytadbit import tadbit
from string import ascii_lowercase as letters
from warnings import warn
from copy import deepcopy as copy
from cPickle import load, dump


try:
    from scipy.interpolate import interp1d
    from numpy import log2
    from matplotlib import pyplot as plt
except ImportError:
    from pytadbit.utils import Interpolate as interp1d

from random import random, shuffle



def load_slice(in_f):
    """
    Load Slice from file. Slice might have been saved through the
    :func:`Slice.save_slice`.
    
    :para in_f: path to a saved Slice file
    
    :returns: Slice object
    """
    dico = load(open(in_f))
    slicy = Slice(dico['name'], resolution=dico['resolution'])
    slicy.experiments     = dico['experiments']
    slicy.resolution      = dico['resolution']
    slicy._ori_resolution = dico['_ori_resolution']
    slicy.size            = dico['size']
    slicy.r_size          = dico['r_size']
    slicy.max_tad_size    = dico['max_tad_size']
    slicy.forbidden       = dico['forbidden']
    slicy._centromere     = dico['_centromere']
    return slicy


class Slice():
    """
    Slice object designed to deal with Topologically Associating Domains
    predictions from different experiments, in different cell types for a given
    slice of DNA, and compare them.

    :param name: name of the slice (might be a chromosome name for example).
    :param None resolution: resolution of the experiments. All experiments may
        have the same resolution
    :param None experiments: :py:func:`list` of paths to files containing the
        definition of TADs corresponding to different experiments (or output
        of tadbit)
    :param None experiment_names: :py:func:`list` of names for each experiment
    :param None wanted_resoultion: change the resolution of the experiment
       loaded. This with :func:`set_resolution`
    :param 3000000 max_tad_size: maximum size of TAD allowed. TADs longer than
        this will not be considered, and relative slice size will be reduced
        accordingly
    :param None chr_len: size of the DNA slice in bp. By default it will be
        inferred from the distribution of TADs.

    :return: Slice object

    """
    def __init__(self, name, resolution=None, experiments=None, 
                 experiment_names=None, wanted_resoultion=None,
                 max_tad_size=3000000, chr_len=None, parser=None):
        self.name             = name
        self.max_tad_size     = max_tad_size
        self.size             = self.r_size = chr_len
        if not resolution:
            raise Exception ('ERROR: please provide a resolution value\n')
        self.resolution       = resolution
        self._ori_resolution  = resolution
        self.forbidden        = {}
        self.experiments      = {}
        self._ori_experiments = {}
        self._centromere      = None
        experiments           = experiments or []
        for i, f_name in enumerate(experiments):
            name = experiment_names[i] if experiment_names else None
            if type(f_name) is dict:
                self.add_tad_def(f_name, name)
            else:
                self.add_experiment(f_name, name, parser=parser)
        if wanted_resoultion:
            self.set_resolution(wanted_resoultion)


    def save_slice(self, out_f):
        """
        Save Slice object to file (it uses :py:func:`pickle.load` from the
        :py:mod:`cPickle`). Once saved, the object may be loaded through
        :func:`load_slice`.

        :param out_f: path to file to dump the :py:mod:`cPickle` object.
        """
        dico = {}
        dico['experiments']     = self.experiments
        dico['name']            = self.name
        dico['resolution']      = self.resolution
        dico['_ori_resolution'] = self._ori_resolution
        dico['size']            = self.size
        dico['r_size']          = self.r_size
        dico['max_tad_size']    = self.max_tad_size
        dico['forbidden']       = self.forbidden
        dico['_centromere']     = self._centromere
        out = open(out_f, 'w')
        dump(dico, out)
        out.close()


    def set_resolution(self, resolution, keep_original=True):
        """
        Set a new value for resolution. copy original data into
        Slice.ori_experiments and replaces the Slice.experiments
        with the data corresponding to new data.

        :param resolution: an integer, representing resolution. This number must
            be a multiple of the original resolution, and higher than it.
        :param True keep_original: either to keep or not the original data
        """
        if resolution < self._ori_resolution:
            raise Exception('New resolution might be higher than original.')
        if resolution % self._ori_resolution:
            raise Exception('New resolution might be a mulitple original...\n' +
                            '  otherwise it is too complicated for me :P')
        # if we want to go back to original resolution
        if resolution == self._ori_resolution:
            self.experiments = self._ori_experiments
            self.resolution = self._ori_resolution
            return
        # if we already changed resolution before
        if self.resolution == self._ori_resolution:
            self._ori_experiments = self.experiments
        self.resolution = resolution
        self.experiments = {}
        fact = self.resolution / self._ori_resolution
        # super for!
        for x_name, x_p in self._ori_experiments.iteritems():
            size = x_p['size']
            self.experiments[x_name] = {'hi-c':[[]], 'size': size/fact}
            for i in xrange(0, size - fact/2, fact):
                for j in xrange(0, size - fact/2, fact):
                    val = 0
                    for k in xrange(fact):
                        for l in  xrange(fact):
                            val += x_p['hi-c'][0][(i + k) * size + j + l]
                    self.experiments[x_name]['hi-c'][0].append(val)
        if not keep_original:
            del(self._ori_experiments)


    def align_experiments(self, names=None, verbose=False, randomize=False,
                          rnd_method='interpolate', **kwargs):
        """
        Align prediction of boundaries of two different experiments

        :param None names: list of names of experiments to align. If None
            align all.
        :param experiment1: name of the first experiment to align
        :param experiment2: name of the second experiment to align
        :param -0.1 penalty: penalty of inserting a gap in the alignment
        :param 500000 max_dist: Maximum distance between 2 boundaries allowing
            match
        :param False verbose: print somethings
        :param False randomize: check alignment quality by comparing
            randomization of boundaries over Slices of same size. This will
            return a extra value, the p-value of accepting that observed
            alignment is not better than random alignment
        :param interpolate rnd_method: when radomizinf use interpolation of TAD
           distribution. Alternative is 'shuffle' where TADs are simply shuffled
        """
        xpers = names or self.experiments.keys()
        tads = []
        for e in xpers:
            if not self.experiments[e]['tads']:
                raise Exception('No TADs defined, use find_tad function.\n')
            tads.append(self.experiments[e]['brks'])
        aligneds, score = align(tads, bin_size=self.resolution,
                                max_dist=self.max_tad_size, **kwargs)
        for e, ali in zip(xpers, aligneds):
            self.experiments[e]['align'] = ali
            self.experiments[e]['align'] = ali
        if verbose:
            self.print_alignment(xpers)
        if not randomize:
            return self.get_alignment(names), score
        #mean, std = self._get_tads_mean_std(xpers)
        #print 'mean', mean, 'std', std, self.r_size, self.r_size/mean
        #p_value = randomization_test(len(xpers), mean, std, score,
        #                             self.r_size, self.resolution,
        #                             verbose=verbose, **kwargs)
        p_value = self.randomization_test(xpers, score=score, method=rnd_method,
                                          verbose=verbose, **kwargs)
        return score, p_value


    def print_alignment(self, names=None, string=False):
        """
        print alignment
        
        :param None names: if None print all experiments
        :param False string: return string instead of printing
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
        
        :param None names: if None print all experiments
        :returns: alignment as :py:class:`dict`
        """
        names = names or self.experiments.keys()
        return dict([(e, self.experiments[e]['align']) \
                     for e in names if 'align' in self.experiments[e]])
                    

    def add_experiment(self, f_name, name, force=False, parser=None):
        """
        Add Hi-C experiment to Slice
        
        :param f_name: path to tsv file
        :param name: name of the experiment
        :param False force: overwrite experiments loaded under the same name
        :param None parser: a parser function that returns a tuple of lists
           representing the data matrix, and the length of a row/column, with
           this file example.tsv:

           ::
           
             chrT_001	chrT_002	chrT_003	chrT_004
             chrT_001	629	164	88	105
             chrT_002	86	612	175	110
             chrT_003	159	216	437	105
             chrT_004	100	111	146	278
           
           the output of parser('example.tsv') might be:
           ``[([629, 86, 159, 100, 164, 612, 216, 111, 88, 175, 437, 146, 105,
           110, 105, 278]), 4]``
        
        """
        nums, size = read_matrix(f_name, parser=parser)
        if not name:
            name = ''.join([letters[int(random() * len(letters))] \
                            for _ in xrange(5)])
            warn('No name provided, random name generated: {}\n'.format(name))
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
                 max_tad_size="auto", no_heuristic=False, batch_mode=False):
        """
        Call :func:`pytadbit.tadbit.tadbit` function to calculate the position
        of Topologically associated domains
        
        :param experiment: A square matrix of interaction counts in hi-C
            data or a list of such matrices for replicated experiments. The
            counts must be evenly sampled and not normalized. 'experiment'
            might be either a list of list, a path to a file or a file handler
        :param None n_cpus: The number of CPUs to allocate to tadbit. The
            value default is the total number of CPUs minus 1.
        :param auto max_tad_size: an integer defining maximum size of TAD.
            Default (auto) defines it to the number of rows/columns.
        :param False no_heuristic: whether to use or not some heuristics
        :param False batch_mode: if True, all experiments will be concatenated
            into one for the search of TADs. The resulting TADs found are stored
            under the name 'batch' plus a concatenation of the experiment names
            passed (i.e.: if experiments=['exp1', 'exp2'], the name would be:
            'batch_exp1_exp2').
        
        """
        if batch_mode:
            matrix = []
            name = 'batch'
            for xpr in experiments:
                matrix.append(self.experiments[xpr]['hi-c'][0])
                name += '_' + xpr
            result, weights = tadbit(matrix,
                                     n_cpus=n_cpus, verbose=verbose,
                                     max_tad_size=max_tad_size,
                                     no_heuristic=no_heuristic,
                                     get_weights=True)
            self.add_tad_def(result, name=name, weights=weights)
            return
        for experiment in experiments:
            result, weights = tadbit(self.experiments[experiment]['hi-c'],
                                     n_cpus=n_cpus, verbose=verbose,
                                     max_tad_size=max_tad_size,
                                     no_heuristic=no_heuristic,
                                     get_weights=True)
            self.add_tad_def(result, name=experiment, weights=weights)
            # search for centromere
            self._search_centromere(experiment)
        

    def visualize(self, name, tad=None, paint_tads=False, ax=None, show=False):
        """
        Visualize the matrix of Hi-C interactions

        :param name: name of the experiment to visualize
        
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
            matrix = [[self.experiments[name]['hi-c'][0][i+size*j]\
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
            if not tad['brk']:
                for j in xrange(int(tad['start']), int(tad['end']), 4):
                    ax.hlines(j, tad['start'], tad['end'], colors='k')
        if show:
            plt.show()


    def add_tad_def(self, f_name, name=None, weights=None):
        """
        Add Topologically Associated Domains defintinion detection to Slice

        :param f_name: path to file
        :param None name: name of the experiment, if None f_name will be
            used
        :param None weights: Store information about the weights, corresponding
            to the normalization of the Hi-C data (see tadbit function
            documentation)
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


    def _interpolation(self, experiments):
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
        

    def randomization_test(self, xpers, score=None, num=1000, verbose=False,
                           method='interpolate'):
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
        tads = []
        for e in xpers:
            if not self.experiments[e]['tads']:
                raise Exception('No TADs defined, use find_tad function.\n')
            tads.append([t['end']-t['start']for t in self.experiments[e]['tads'].values()])
        rnd_distr = []
        rnd_len = []
        distr = self._interpolation(xpers) if method is 'interpolate' else None
        num_sequences = len(tads)
        rnd_exp = lambda : tads[int(random() * num_sequences)]
        for val in xrange(num):
            if verbose:
                val = float(val)
                if not val / num * 100 % 5:
                    stdout.write('\r' + ' ' * 10 + 
                                 ' randomizing: '
                                 '{:.2%} completed'.format(val/num))
                    stdout.flush()
            if method is 'interpolate':
                rnd_tads = [generate_rnd_tads(self.r_size, distr,
                                              self.resolution)\
                            for _ in xrange(num_sequences)]
                rnd_len.append(float(sum([len(r) for r in rnd_tads])) \
                               / len(rnd_tads))
            else:
                rnd_tads = [generate_shuffle_tads(rnd_exp()) \
                            for _ in xrange(num_sequences)]
                rnd_len.append(num_sequences)
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
        
        :param tad: a given TAD -> :py:class:`dict`
        :param x_name: name og the experiment
        :param True normed: if Hi-C data has to be normalized
        
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
        
        :param x_name: name of the experiment
        :param True normed: normalize Hi-C data returned
        
        :yields: Hi-C data corresponding to each TAD
        """
        if not self.experiments[x_name]['hi-c']:
            raise Exception('No Hi-c data for {} experiment\n'.format(x_name))
        for name, ref in self.experiments[x_name]['tads'].iteritems():
            if not ref['brk']:
                continue
            yield name, self.get_tad_hic(ref, x_name, normed=normed)


    def set_max_tad_size(self, value):
        """
        Change maximum size allowed for TADs. Also apply it to computed
        experiments.

        :param 3000000 value: an int
        """
        self.max_tad_size = value
        for xpr in self.experiments:
            tads = self.experiments[xpr]['tads']
            for tad in tads:
                if (tads[tad]['end'] - tads[tad]['start']) \
                   * self.resolution < self.max_tad_size:
                    tads[tad]['brk'] = tads[tad]['end']
                else:
                    tads[tad]['brk'] = None
            self.experiments[xpr]['brks'] = []
            for tad in tads:
                if tads[tad]['brk']:
                    self.experiments[xpr]['brks'].append(tads[tad]['brk'])
            

    def _search_centromere(self, x_name):
        """
        Search for centromere in slice, presuposing that slice corresponds to a
        chromosome.
        Add a boundary to all experiments where the centromere is.
         * A centromere is defined as the area where the rows/columns of the
           Hi-C matrix is empty.
        """
        beg = 0
        end = 0
        size = self.experiments[x_name]['size']
        hic = self.experiments[x_name]['hi-c'][0]
        for pos, raw in enumerate(xrange(0, size*size, size)):
            if sum(hic[raw:raw+size]) == 0 and not beg:
                beg = pos
            if sum(hic[raw:raw+size]) != 0 and beg:
                end = pos
                break
        if not beg or not end:
            return
        tads = self.experiments[x_name]['tads']
        # in case we already have a centromere defined, check it can be reduced
        if self._centromere:
            if beg > self._centromere[0]:
                for tad in tads:
                    if tads[tad]['end'] == self._centromere[0]:
                        tads[tad]['end'] = beg
                        self._centromere[0] = beg
            if end < self._centromere[1]:
                for tad in tads:
                    if tads[tad]['beg'] == self._centromere[1]:
                        tads[tad]['start'] = end
                        self._centromere[1] = end
        else:
            self._centromere = [beg, end]
            tad = len(tads)
            plus = 0
            while tad + plus >= 1:
                start = tads[tad - 1 + plus]['start']
                final = tads[tad - 1 + plus]['end']
                # centromere found?
                if start < beg < final and start < end < final:
                    tads[tad]     = copy(tads[tad - 1])
                    tads[tad]['start'] = end
                    if (tads[tad]['end'] - tads[tad]['start']) \
                           * self.resolution > self.max_tad_size:
                        tads[tad]['brk'] = None
                    else:
                        tads[tad]['brk'] = tads[tad]['end']
                    tad -= 1
                    tads[tad] = copy(tads[tad])
                    tads[tad]['end'] = beg
                    if (tads[tad]['end'] - tads[tad]['start']) \
                           * self.resolution > self.max_tad_size:
                        tads[tad]['brk'] = None
                    else:
                        tads[tad]['brk'] = tads[tad]['end']
                    plus = 1
                else:
                    tads[tad] = copy(tads[tad - 1 + plus])
                tad -= 1
        self.experiments[x_name]['brks'] = []
        for tad in tads:
            if tads[tad]['brk']:
                self.experiments[x_name]['brks'].append(tads[tad]['brk'])


    def _get_tads_mean_std(self, experiments):
        """
        returns mean and standard deviation of TAD lengths. Value is for both
        TADs.

        ..Note: no longer used in core functions
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


def generate_rnd_tads(slice_len, distr, bin_size, start=0):
    """
    Generates random TADs over a slice of a given size according to a given
    distribution of lengths of TADs.
    
    :param slice_len: length of the slice
    :param distr: function that returns a TAD length depending on a p value
    :param bin_size: size of the bin of the Hi-C experiment
    :param 0 start: starting position in the slice
    
    :returns: list of TADs
    """
    pos = start
    tads = []
    while True:
        pos += distr(random())
        if pos > slice_len:
            break
        tads.append(float(int(pos / bin_size + .5)))
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
    for t in rnd_tads:
        if tads:
            t += tads[-1]
        tads.append(t)
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

