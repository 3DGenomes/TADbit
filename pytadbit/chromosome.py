"""
26 Nov 2012


"""

from math import sqrt, log, exp
from sys import stdout
from pytadbit.parsers.tad_parser import parse_tads
from pytadbit.tads_aligner.aligner import align

from scipy.stats import lognorm


class Chromosome():
    """
    Chromosome object designed to deal with Topologically associated domains
    predictions from different experiments, in different cell types and compare
    them.
    """
    def __init__(self, name, resolution, experiment_files=None, experiment_names=None,
                 max_tad_size=3000000, chr_len=None):
        """
        :argument name: name of the chromosome
        :argument resolution: resolution of the experiments. All experiments may have
        the same resolution
        :argument None experiment_files: list of paths to files containing the
        definition of TADs corresponding to different experiments
        :argument None experiment_names: list of names for each experiment
        :argument 3000000 max_tad_size: maximum size of TAD allowed. TADs longer than
        this will not be considered, and relative chromosome size will be reduced
        accordingly
        :argument None chr_len: size of the chromosome in bp. By default it will be
        inferred from the distribution of TADs.

        :return: Chromosome object
        """
        self.name = name
        self.max_tad_size = max_tad_size
        self.size = self.r_size = chr_len
        self.resolution = resolution
        self.forbidden = {}
        self.experiments = {}
        for i, f_name in enumerate(experiment_files):
            name = experiment_names[i] if experiment_names else None
            self.add_experiment(f_name, name)


    def align_experiments(self, names=None, verbose=False, randomize=False,
                          **kwargs):
        """
        Align prediction of boundaries of two different experiments

        :argument None names: list of names of experiments to align. If None
        align all.
        :argument experiment1: name of the first experiment to align
        :argument experiment2: name of the second experiment to align
        :argument -0.1 penalty: penalty of inserting a gap in the alignment
        :argument 500000 max_dist: Maximum distance between 2 boundaries allowing match
        :argument False verbose: print somethings
        :argument False randomize: check alignment quality by comparing randomization
        of boundaries over chromosomes of same size. This will return a extra value,
        the p-value of accepting that observed alignment is not better than random
        alignment
        """
        experiments = names or self.experiments.keys()
        tads = []
        for exp in experiments:
            tads.append(self.experiments[exp]['brks'])
        aligneds, score = align(tads, bin_size=self.resolution,
                                chr_len=self.r_size, **kwargs)
        for exp, ali in zip(experiments, aligneds):
            self.experiments[exp]['align'] = ali
            self.experiments[exp]['align'] = ali
        if verbose:
            self.print_alignment(experiments)
        if not randomize:
            return self.get_alignment(names), score
        mean, std = self.get_tads_mean_std(experiments)
        print 'mean', mean, 'std', std, self.r_size, self.r_size/mean
        p_value = randomization_test(len(experiments), mean, std, score,
                                     self.r_size, self.resolution,
                                     verbose=verbose, **kwargs)
        return score, p_value


    def print_alignment(self, names=None, string=False):
        names = names or self.experiments.keys()
        length = max([len(n) for n in names])
        out = 'Alignment (%s TADs)\n' % (len(names))
        for exp in names:
            if not 'align' in self.experiments[exp]:
                continue
            out += '{1:{0}}:'.format(length, exp)
            out += '|'.join(['%5s' % (str(x)[:-2] if x!='-' else '-' * 4)\
                             for x in self.experiments[exp]['align']]) + '\n'
        if string:
            return out
        print out


    def get_alignment(self, names=None):
        names = names or self.experiments.keys()
        return dict([(exp, self.experiments[exp]['align']) \
                     for exp in names if 'align' in self.experiments[exp]])
                    

    def add_experiment(self, f_name, name=None):
        """
        Add experiment of Topologically Associated Domains detection to chromosome

        :argument f_name: path to file
        :argument None name: name of the experiment, if None f_name will be used:
        """
        name = name or f_name
        self.experiments[name] = {}
        tads, forbidden = parse_tads(f_name, max_size=self.max_tad_size,
                                     bin_size=self.resolution)
        brks = [t['brk'] for t in tads.values() if t['brk']]
        self.experiments[name] = {'tads': tads,
                                  'brks': brks}
        if not self.forbidden:
            self.forbidden = dict([(f, None) for f in forbidden])
        else:
            self.forbidden = dict([(f, None) for f in \
                                   forbidden.intersection(self.forbidden)])
        if not self.size:
            self.size = tads[max(tads)]['end'] * self.resolution
        self.r_size = self.size - len(self.forbidden) * self.resolution


    def get_tads_mean_std(self, experiments):
        """
        returns mean and standard deviation of TAD lengths. Value is for both TADs
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


def generate_random_tads(chr_len, mean, sigma, bin_size,
                         max_tad_size=3000000, start=0):
    pos = start
    tads = []
    dist = lognorm(sigma, scale=exp(mean))
    while True:
        val = dist.rvs()
        if val > max_tad_size:
            continue
        pos += val
        if pos > chr_len:
            tads.append(float(int(chr_len/bin_size)))
            break
        tads.append(float(int(pos/bin_size+.5)))
    return tads


def randomization_test(num_sequences, mean, std, score, chr_len, bin_size,
                       num=1000, verbose=False):
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
        random_tads = [generate_random_tads(chr_len, mean,
                                            std, bin_size) \
                       for _ in xrange(num_sequences)]
        rand_len.append(float(sum([len(r) for r in random_tads]))/len(random_tads))
        rand_distr.append(align(random_tads,
                                bin_size=bin_size, chr_len=chr_len,
                                verbose=False)[1])
        if rand_distr[-1] > best[0]:
            best = rand_distr[-1], random_tads
    align(random_tads, bin_size=bin_size, chr_len=chr_len, verbose=True)
    p_value = float(len([n for n in rand_distr if n > score]))/len(rand_distr)
    if verbose:
        stdout.write('\n {} randomizations finished.'.format(num))
        stdout.flush()
        print 'Observed alignment score: {}'.format(score)
        print '  Randomized scores between {} and {}'.format(min(rand_distr),
                                                             max(rand_distr))
        print 'p-value: {}'.format(p_value if p_value else '<{}'.format(1./num))
        print sum(rand_len)/len(rand_len)
    return p_value    



def main():
    """
    main function
    """

    from pytadbit import Chromosome
    from itertools import combinations
    path = '/home/fransua/db/hi-c/dixon_hsap-mmus/hsap/20Kb/'
    sample1 = path + 'chr21/chr21.H1.hESC.rep1.20000.txt'
    sample2 = path + 'chr21/chr21.IMR90.rep1.20000.txt'
    sample3 = path + 'chr21/chr21.H1.hESC.rep2.20000.txt'
    sample4 = path + 'chr21/chr21.IMR90.rep2.20000.txt'
    chr_len  = 48000000
    bin_size = 20000
    # given that mean TAD size around 1Mb we set:
    max_size = 3000000
    max_dist = 500000
    
    chr21 = Chromosome('chr21', 20000,
                       experiment_files = [sample1, sample2, sample3, sample4],
                       experiment_names = ['hESC.rep1','IMR90.rep1',
                                           'hESC.rep2','IMR90.rep2'],
                       max_tad_size=max_size, chr_len=chr_len)

    names = ['hESC.rep1','IMR90.rep1', 'hESC.rep2','IMR90.rep2']
    for experiments in combinations(names, 2):
        print [len(chr21.experiments[n]['brks']) for n in experiments]
        print '-'.join(experiments)
        print '', chr21.align_experiments(experiments, verbose=True,
                                          randomize=True, num=10000)
        print float(sum([len(chr21.experiments[n]['brks']) for n in experiments]))/2
        print '='*80

    for experiments in combinations(names, 3):
        print '-'.join(experiments)
        print '', chr21.align_experiments(experiments, verbose=True,
                                          randomize=True, num=10000)
        print '='*80

    for experiments in combinations(names, 4):
        print '-'.join(experiments)
        print '', chr21.align_experiments(experiments, verbose=True,
                                          randomize=True, num=10000)
        print '='*80

    chr21.align_experiments(method='global', bin_size=20000, chr_len=chr_len,
                            penalty=-0.1, max_dist=500000, verbose=True)
    
if __name__ == "__main__":
    exit(main())

