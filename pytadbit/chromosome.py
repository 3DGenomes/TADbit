"""
26 Nov 2012


"""

from random import gauss
from math import sqrt
from sys import stdout
from pytadbit.parsers.tad_parser import parse_tads
from pytadbit.tads_aligner.aligner import align

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


    def align_experiments(self, experiment1, experiment2, penalty=-0.1,
                          max_dist=500000, verbose=False, randomize=False):
        """
        Align prediction of boundaries of two different experiments

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
        tads1 = self.experiments[experiment1]['brks']
        tads2 = self.experiments[experiment2]['brks']
        align1, align2, score = align(tads1, tads2, bin_size=self.resolution,
                                      chr_len=self.r_size, penalty=penalty,
                                      max_dist=max_dist, verbose=verbose)
        self.experiments[experiment1]['align'] = align1
        self.experiments[experiment2]['align'] = align2
        if not randomize:
            return score
        mean, std = get_tads_mean_std(tads1, tads2, self.resolution)
        p_value = randomization_test(mean, std, score, self.r_size,
                                     self.resolution, num=1000, verbose=verbose)
        return score, p_value


    def add_experiment(self, f_name, name=None):
        """
        Add experiment of Topologically Associated Domains detection to chromosome

        :argument f_name: path to file
        :argument None name: name of the experiment, if None f_name will be used:
        """
        name = name or f_name
        self.experiments[name] = {}
        tads, forbidden = parse_tads(f_name, self.max_tad_size, self.resolution)
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


    def __str__(self):
        """
        none sense now
        """
        string = 'Topologically Associated Domains\' boundaries:\n'
        string += 'FOUND: '
        string += '|'.join(['%4s' % str(self.tads[x]['end'])[:-2] \
                            for x in self.tads])
        string += '\n'
        string += 'VALID: '
        string += '|'.join(['%4s' % (str(self.tads[x]['brk'])[:-2] \
                                     if self.tads[x]['brk'] else '-'*3) \
                            for x in self.tads])
        return string


def get_tads_mean_std(tads1, tads2, bin_size):
    """
    returns mean and standard deviation of TAD lengths. Value is for both TADs
    """
    norm_tads = [(tads1[t] - tads1[t-1]) * bin_size for t in xrange(1, len(tads1))] +\
                [(tads2[t] - tads2[t-1]) * bin_size for t in xrange(1, len(tads2))]
    length = len(norm_tads)
    mean   = sum(norm_tads)/length
    std    = sqrt(sum([(t-mean)**2 for t in norm_tads])/length)
    return mean, std


def my_gauss(mean, sigma, minimum=0):
    """
    returns random gaussian value according to mean and sigma
    """
    while True:
        val = gauss(mean, sigma)
        if val > minimum:
            return val


def generate_random_tads(chr_len, mean, sigma, bin_size, start=1):
    """
    """
    pos = start
    tads = [1.0]
    while True:
        pos += my_gauss(mean, sigma, minimum=bin_size)
        if pos > chr_len: break
        tads.append(float(int(pos/bin_size)))
    return tads


def randomization_test(mean, std, score, chr_len, bin_size, num=1000,
                       max_dist=500000, verbose=False):
    rand_distr = []
    for n in xrange(num):
        if verbose:
            n = float(n)
            if not n / num * 100 % 5:
                stdout.write('\r' + ' ' * 10 + \
                             ' randomizing: {:.2%} completed'.format(n/num))
                stdout.flush()
        rand_distr.append(needleman_wunsch(generate_random_tads(chr_len, mean,
                                                                std, bin_size),
                                           generate_random_tads(chr_len, mean,
                                                                std, bin_size),
                                           bin_size=bin_size, chr_len=chr_len,
                                           verbose=False)[2])
    if verbose:
        print '\n Randomization Finished.'
    print score, min(rand_distr), max(rand_distr)
    p_value = float(len([n for n in rand_distr if n > score]))/len(rand_distr)
    return p_value    



def main():
    """
    main function

    sample1 = 'chr21/chr21.H1.hESC.rep1.20000.tsv'
    sample2 = 'chr21/chr21.H1.hESC.rep2.20000.tsv'
    
    sample1 = 'chr21/chr21.H1.hESC.rep1.20000.tsv'
    sample2 = 'chr21/chr21.IMR90.rep1.20000.tsv'
    
    sample1 = '/home/fransua/scratch/helping_marc/chr19/chr19.hindIII.rep2.20000.tsv'
    sample2 = '/home/fransua/scratch/helping_marc/chr19/chr19.cortex.rep2.20000.tsv'
    
    sample1 = '/home/fransua/scratch/helping_marc/chr19/chr19.hindIII.rep2.20000.tsv'
    sample2 = '/home/fransua/scratch/helping_marc/chr19/chr19.J1_mESC_NcoI.20000.tsv'
    
    sample1 = '/home/fransua/scratch/helping_marc/chr1/chr1.J1.mESC.NcoI.20000.tsv'
    sample2 = '/home/fransua/scratch/helping_marc/chr1/chr1.hindIII.rep2.20000.tsv'

    sample1 = 'chr21/chr21.H1.hESC.test.20000.tsv'

    """

    sample1 = '/home/francisco/scratch/pytadbit_tries/chr21/chr21.H1.hESC.rep1.20000.tsv'
    sample2 = '/home/francisco/scratch/pytadbit_tries/chr21/chr21.IMR90.rep1.20000.tsv'
    sample3 = '/home/francisco/scratch/pytadbit_tries/chr21/chr21.H1.hESC.rep2.20000.tsv'
    sample4 = '/home/francisco/scratch/pytadbit_tries/chr21/chr21.IMR90.rep2.20000.tsv'
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

    chr21.align_experiments('hESC.rep1', 'IMR90.rep1', verbose=True)

    sequences = []
    sequences.append(chr21.experiments['hESC.rep1' ]['brks'][:])
    sequences.append(chr21.experiments['hESC.rep2' ]['brks'][:])
    sequences.append(chr21.experiments['IMR90.rep1']['brks'][:])
    sequences.append(chr21.experiments['IMR90.rep2']['brks'][:])
    names = ['hESC.rep1' , 'hESC.rep2', 'IMR90.rep1', 'IMR90.rep2']

    align(sequences, method='global', bin_size=20000, chr_len=chr_len, penalty=-0.1,
          max_dist=500000, verbose=False)
    

    for ali, name in zip(aligneds, names):
        print '%10s :' % (name) +'|'.join(['%4s' % (str(x)[:-2] if x!='-' else '-'*3) for x in ali])
        
    print 'TADS 1: '+'|'.join(['%4s' % (str(x)[:-2] if x!='-' else '-'*3) for x in reference])
    print 'TADS 1: '+'|'.join(['%4s' % (str(x)[:-2] if x!='-' else '-'*3) for x in align1])
    
if __name__ == "__main__":
    exit(main())

