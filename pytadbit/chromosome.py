"""
26 Nov 2012


"""

from os.path                           import exists
from pytadbit.boundary_aligner.aligner import align
from pytadbit                          import tadbit
from pytadbit.experiment               import Experiment
from string                            import ascii_lowercase as letters
from warnings                          import warn
from copy                              import deepcopy as copy
from cPickle                           import load, dump
from pytadbit.alignment                import Alignment, randomization_test
from numpy                             import log2
from random                            import random

try:
    from matplotlib import pyplot as plt
except ImportError:
    warn('matplotlib not found\n')




def load_chromosome(in_f, fast=2):
    """
    Load a Chromosome object from a file. A Chromosome object can be saved with
    the :func:`Chromosome.save_chromosome` function. 
    
    :param in_f: path to a saved Chromosome object file
    :param 2 fast: if fast=2 do not load the Hi-C data (in the case that they 
       were saved in a separate file see :func:`Chromosome.save_chromosome`).
       If fast is equal to 1, the weights will be skipped from load to save 
       memory. Finally if fast=0, both the weights and Hi-C data will be loaded
    
    :returns: a Chromosome object

    TODO: remove first try/except type error... this is loading old experiments
    """
    dico = load(open(in_f))
    name = ''
    crm = Chromosome(dico['name'])
    for name in dico['experiments']:
        xpr = Experiment(name, dico['experiments'][name]['resolution'], 
                         no_warn=True)
        xpr.tads       = dico['experiments'][name]['tads']
        xpr.wght       = dico['experiments'][name]['wght']
        xpr.hic_data   = dico['experiments'][name]['hi-c']
        xpr.conditions = dico['experiments'][name]['cond']
        xpr.size       = dico['experiments'][name]['size']
        try:
            crm.experiments.append(xpr)
        except TypeError:
            continue
    crm.size            = dico['size']
    crm.r_size          = dico['r_size']
    crm.max_tad_size    = dico['max_tad_size']
    crm.forbidden       = dico['forbidden']
    crm._centromere     = dico['_centromere']
    if type(dico['experiments'][name]['hi-c']) == str and fast!= int(2):
        try:
            dicp = load(open(in_f + '_hic'))
        except IOError:
            raise Exception('ERROR: file %s not found\n' % (
                dico['experiments'][name]['hi-c']))
        for name in dico['experiments']:
            crm.get_experiment(name).hic_data = dicp[name]['hi-c']
            if fast != 1:
                crm.get_experiment(name).wght = dicp[name]['wght']
    elif not fast:
        warn('WARNING: data not saved correctly for fast loading.\n')
    return crm


class Chromosome(object):
    """
    A Chromosome object designed to deal with Topologically Associating Domains
    predictions from different experiments, in different cell types for a given
    chromosome of DNA, and to compare them.

    :param name: name of the chromosome (might be a chromosome name for example)
    :param None resolutions: list of resolutions corresponding to a list of
       experiments passed.
    :param None experiment_hic_data: :py:func:`list` of paths to files
       containing the Hi-C count matrices corresponding to different experiments
    :param None experiment_tads: :py:func:`list` of paths to files
       containing the definition of TADs corresponding to different experiments
    :param None experiment_names: :py:func:`list` of the names of each 
        experiment
    :param 5000000 max_tad_size: maximum TAD size allowed. TADs longer than
        this value will not be considered, and size of the corresponding
        chromosome size will be reduced accordingly
    :param 0 chr_len: size of the DNA chromosome in bp. By default it will be
        inferred from the distribution of TADs
    :param None parser: a parser function that returns a tuple of lists 
       representing the data matrix and the length of a row/column. With
       the file example.tsv:

       ::
       
         chrT_001	chrT_002	chrT_003	chrT_004
         chrT_001	629	164	88	105
         chrT_002	164	612	175	110
         chrT_003	88	175	437	100
         chrT_004	105	110	100	278

       the output of parser('example.tsv') would be be:
       ``[([629, 164, 88, 105, 164, 612, 175, 110, 88, 175, 437, 100, 105,
       110, 100, 278]), 4]``

    :return: Chromosome object


    """
    def __init__(self, name, experiment_resolutions=None, experiment_tads=None,
                 experiment_hic_data=None, experiment_names=None,
                 max_tad_size=5000000, chr_len=0, parser=None):
        self.name             = name
        self.max_tad_size     = max_tad_size
        self.size             = self._given_size = self.r_size = chr_len
        self.size             = ChromosomeSize(self.size)
        self.r_size           = RelativeChromosomeSize(self.size)
        self.forbidden        = {}
        self.experiments      = ExperimentList([], self)
        self._centromere      = None
        self.alignment        = AlignmentDict()
        if experiment_tads:
            for i, handler in enumerate(experiment_tads or []):
                name = experiment_names[i] if experiment_names else None
                self.add_experiment(name, experiment_resolutions[i],
                                    tad_def=handler, parser=parser)
        if experiment_hic_data:
            for i, handler in enumerate(experiment_hic_data or []):
                name = experiment_names[i] if experiment_names else None
                try:
                    xpr = self.get_experiment(name)
                    xpr.load_experiment(handler)
                    continue
                except:
                    pass
                if type(handler) == Experiment:
                    name = name or handler.name
                    self.experiments.append(handler)
                else:
                    self.add_experiment(name, experiment_resolutions[i],
                                        hic_data=handler, parser=parser)


    def _get_forbidden_region(self, xpr):
        """
        Find the regions for which there is no information in any of the
        experiments. This is used to infer the relative chromosome size.
        """
        if not xpr.tads:
            return
        forbidden = []
        for pos in xpr.tads:
            start = float(xpr.tads[pos]['start'])
            end   = float(xpr.tads[pos]['end'])
            diff  = end - start
            if diff * xpr.resolution > self.max_tad_size:
                forbidden += range(int(start), int(end+1))
                xpr.tads[pos]['score'] = -abs(xpr.tads[pos]['score'])
        if not self.forbidden:
            self.forbidden = dict([(f, None) for f in forbidden])
        else:
            self.forbidden = dict([(f, None) for f in 
                                   set(forbidden).intersection(self.forbidden)])
        # search for centromere:
        self._search_centromere(xpr)
        # add centromere as forbidden region:
        if self._centromere:
            for pos in xrange(int(self._centromere[0]),
                              int(self._centromere[1])):
                self.forbidden[pos] = 'Centromere'
        self.__update_size(xpr)


    def get_experiment(self, name):
        """
        This can also be done directly with Chromosome.experiments[name].
        
        :param name: name of the experiment to select
        :returns: :class:`pytadbit.Experiment`
        """
        for exp in self.experiments:
            if exp.name == name:
                return exp
        raise Exception('ERROR: experiment ' +
                        '%s not found\n' % (name))
                

    def save_chromosome(self, out_f, fast=True, divide=True, force=False):
        """
        Save a Chromosome object to a file (it uses :py:func:`pickle.load` from
        the :py:mod:`cPickle`). Once saved, the object can be loaded with
        :func:`load_chromosome`.

        :param out_f: path to the file where to store the :py:mod:`cPickle`
           object
        :param True fast: if True, skip Hi-C data and weights
        :param True divide: if True writes two pickles, one with what would
           result by using the fast option, and the second with the Hi-C and 
           weights data. The second file name will be extended by '_hic' (ie:
           with out_f='chromosome12.pik' we would obtain chromosome12.pik and
           chromosome12.pik_hic). When loaded :func:`load_chromosome` will
           automatically search for both files
        :param False force: overwrite the existing file

        """
        while exists(out_f) and not force:
            out_f += '_'
        dico = {}
        dico['experiments'] = {}
        if divide:
            dicp = {}
        for xpr in self.experiments:
            dico['experiments'][xpr.name] = {
                'size'      : xpr.size,
                'cond'      : xpr.conditions,
                'tads'      : xpr.tads,
                'resolution': xpr.resolution,
                'hi-c'      : None,
                'wght'      : None}
            if fast:
                continue
            if divide:
                dicp[xpr.name] = {
                    'wght': xpr.wght,
                    'hi-c': xpr.hic_data}
                dico['experiments'][xpr.name]['wght'] = None
                dico['experiments'][xpr.name]['hi-c'] = None
            else:
                dico['experiments'][xpr.name]['wght'] = xpr.wght
                dico['experiments'][xpr.name]['hi-c'] = xpr.hic_data
        dico['name']            = self.name
        dico['size']            = self.size
        dico['r_size']          = self.r_size
        dico['max_tad_size']    = self.max_tad_size
        dico['forbidden']       = self.forbidden
        dico['_centromere']     = self._centromere
        out = open(out_f, 'w')
        dump(dico, out)
        out.close()
        if not fast:
            out = open(out_f + '_hic', 'w')
            dump(dicp, out)
            out.close()


    def align_experiments(self, names=None, verbose=False, randomize=False,
                          rnd_method='interpolate', rnd_num=1000, **kwargs):
        """
        Align the predicted boundaries of two different experiments. The 
        resulting alignment will be stored in the self.experiment list.
        
        :param None names: list of names of the experiments to align. If None,
            align all
        :param experiment1: name of the first experiment to align
        :param experiment2: name of the second experiment to align
        :param -0.1 penalty: penalty for inserting a gap in the alignment
        :param 100000 max_dist: maximum distance between two boundaries
            allowing match (100Kb seems fair with HUMAN chromosomes)
        :param False verbose: if True, print some information about the 
            alignments
        :param False randomize: check the alignment quality by comparing
            randomized boundaries over Chromosomes of the same size. This will
            return a extra value, the p-value of accepting that the observed
            alignment is not better than a random alignment
        :param interpolate rnd_method: by default uses the interpolation of TAD
           distribution. The alternative method is 'shuffle', where TADs are
           simply shuffled
        :param 1000 rnd_num: number of randomizations to do
        :param reciprocal method: if global, Needleman-Wunsch is used to align
            (see :func:`pytadbit.boundary_aligner.globally.needleman_wunsch`);
            if reciprocal, a method based on reciprocal closest boundaries is
            used (see :func:`pytadbit.boundary_aligner.reciprocally.reciprocal`)

        :returns: the alignment and the score of the alignment (by default)
        """
        if names:
            xpers = ExperimentList([self.get_experiment(n) for n in names],
                                   self)
        else:
            xpers = self.experiments
        tads = []
        for xpr in xpers:
            if not xpr.tads:
                raise Exception('No TADs defined, use find_tad function.\n')
            tads.append([xpr.tads[x]['brk'] * xpr.resolution for x in xpr.tads])
        # new
        aligneds, score = align(tads, verbose=verbose, **kwargs)
        name = tuple(sorted([x.name for x in xpers]))
        ali = Alignment(name, aligneds, xpers, score=score)
        self.alignment[name] = ali
        if verbose:
            print self.alignment[name]
        # old
        # self.alignment[name] = {}
        # for xpr, ali in zip(xpers, aligneds):
        #     self.alignment[name][xpr.name] = ali
        # if verbose:
        #     self.print_alignment(xpers=xpers)
        if not randomize:
            # return self.get_alignment(name), score
            return ali
        p_value = randomization_test(xpers, score=score, rnd_method=rnd_method,
                                     verbose=verbose, r_size=self.r_size,
                                     num=rnd_num, **kwargs)
        return score, p_value


    def add_experiment(self, name, resolution=None, tad_def=None,
                       hic_data=None, replace=False, parser=None,
                       conditions=None, **kwargs):
        """
        Add a Hi-C experiment to Chromosome
        
        :param name: name of the experiment or of the Experiment object
        :param resolution: resolution of the experiment (needed if name is not
           an Experiment object)
        :param None hic_data: whether a file or a list of lists corresponding to
           the Hi-C data
        :param None tad_def: a file or a dict with precomputed TADs for this
           experiment
        :param False replace: overwrite the experiments loaded under the same
           name
        :param None parser: a parser function that returns a tuple of lists 
           representing the data matrix and the length of a row/column. With 
           a file example.tsv containing:

           ::
           
             chrT_001	chrT_002	chrT_003	chrT_004
             chrT_001	629	164	88	105
             chrT_002	164	612	175	110
             chrT_003	88	175	437	100
             chrT_004	105	110	100	278
           
           the output of parser('example.tsv') would be:
           ``[([629, 164, 88, 105, 164, 612, 175, 110, 88, 175, 437, 100, 105,
           110, 100, 278]), 4]``
        
        """
        if not name:
            name = ''.join([letters[int(random() * len(letters))] \
                            for _ in xrange(5)])
            warn('No name provided, random name generated: %s\n' % (name))
        if name in self.experiments:
            if 'hi-c' in self.get_experiment(name) and not replace:
                warn('''Hi-C data already loaded under the name: %s.
                This experiment will be kept under %s.\n''' % (name,
                                                               name + '_'))
                name += '_'
        if type(name) == Experiment:
            self.experiments.append(name)
        elif resolution:
            self.experiments.append(Experiment(name, resolution, hic_data,
                                               tad_def, parser=parser,
                                               conditions=conditions, **kwargs))
        else:
            raise Exception('resolution param is needed\n')


    def find_tad(self, experiments, name=None, n_cpus=1, verbose=True,
                 max_tad_size="auto", no_heuristic=False, batch_mode=False):
        """
        Call the :func:`pytadbit.tadbit.tadbit` function to calculate the
        position of Topologically Associated Domains
        
        :param experiment: A square matrix of interaction counts of Hi-C
           data or a list of such matrices for replicated experiments. The
           counts must be evenly sampled and not normalized. 'experiment'
           can be either a list of lists, a path to a file or a file handler
        :param 1 n_cpus: The number of CPUs to allocate to TADBit. If
           n_cpus='max' the total number of CPUs will be used
        :param auto max_tad_size: an integer defining the maximum size of a 
           TAD. Default (auto) defines it as the number of rows/columns
        :param False no_heuristic: whether to use or not some heuristics
        :param False batch_mode: if True, all the experiments will be 
           concatenated into one for the search of TADs. The resulting TADs 
           found are stored under the name 'batch' plus a concatenation of the
           experiment names passed (e.g.: if experiments=['exp1', 'exp2'], the
           name would be: 'batch_exp1_exp2').

        TODO: check option -> name for batch mode... some dirty changes....

        """
        if batch_mode:
            matrix = []
            if not name:
                name = 'batch'
            experiments = experiments or self.experiments
            xprs = []
            for xpr in experiments:
                if not type(xpr) == Experiment:
                    xprs.append(self.get_experiment(xpr))
                else:
                    xprs.append(xpr)
            resolution = xprs[0].resolution
            for xpr in sorted(xprs, key=lambda x: x.name):
                if xpr.resolution != resolution:
                    raise Exception('All Experiments might have the same ' +
                                    'resolution\n')
                matrix.append(xpr.hic_data[0])
                if name.startswith('batch'):
                    name += '_' + xpr.name
            result, weights = tadbit(matrix,
                                     n_cpus=n_cpus, verbose=verbose,
                                     max_tad_size=max_tad_size,
                                     no_heuristic=no_heuristic,
                                     get_weights=True)
            experiment = Experiment(name, resolution, hic_data=matrix,
                                    tad_def=result, weights=weights)
            self.add_experiment(experiment)
            return
        if type(experiments) is not list:
            experiments = [experiments]
        for experiment in experiments:
            if not type(experiment) == Experiment:
                xpr = self.get_experiment(experiment)
            result, weights = tadbit(xpr.hic_data,
                                     n_cpus=n_cpus, verbose=verbose,
                                     max_tad_size=max_tad_size,
                                     no_heuristic=no_heuristic,
                                     get_weights=True)
            xpr.load_tad_def(result, weights=weights)
            self._get_forbidden_region(xpr)


    def __update_size(self, xpr):
        """
        Update the chromosome size and relative size after loading new Hi-C
        experiments, unless the Chromosome size was defined by hand.
        
        """
        if not self._given_size:
            self.size = max(xpr.tads[max(xpr.tads)]['end'] * xpr.resolution,
                            self.size)
            self.size   = ChromosomeSize(self.size)
        self.r_size = self.size - len(self.forbidden) * xpr.resolution
        self.r_size = RelativeChromosomeSize(self.size)


    def visualize(self, name, tad=None, focus=None, paint_tads=False, axe=None,
                  show=True, logarithm=True, normalized=False, relative=True):
        """
        Visualize the matrix of Hi-C interactions.

        :param name: name of the experiment to visualize
        :param None tad: a given TAD in the form:
           ::
           
             {'start': start,
              'end'  : end,
              'brk'  : end,
              'score': score}
              
           **Alternatively** a list of the TADs can be passed (all the TADs
           between the first and last one passed will be showed. Thus, passing
           more than two TADs might be superfluous)
        :param None focus: a tuple with the start and end positions of the 
           region to visualize
        :param False paint_tads: draw a box around the TADs defined for this
           experiment
        :param None axe: an axe object from matplotlib can be passed in order
           to customize the picture
        :param True show: either to pop-up matplotlib image or not
        :param True logarithm: show the logarithm values
        :param True normalized: show the normalized data (weights might have
           been calculated previously). *Note: white rows/columns may appear in
           the matrix displayed; these rows correspond to filtered rows (see*
           :func:`pytadbit.utils.hic_filtering.hic_filtering_for_modelling` *)*
        :param True relative: color scale is relative to the whole matrix of
           data, not only to the region displayed
        """
        xper = self.get_experiment(name)
        if logarithm:
            fun = log2
        else:
            fun = lambda x: x
        size = xper.size
        if normalized and not xper.wght:
            raise Exception('ERROR: weights not calculated for this ' +
                            'experiment. Run Experiment.normalize_hic\n')
        if tad and focus:
            raise Exception('ERROR: only one of "tad" or "focus" might be set')
        start = end = None
        if focus:
            start, end = focus
            if start == 0:
                warn('Hi-C matrix starts at 1, setting starting point to 1.\n')
                start = 1
        elif type(tad) == dict:
            start = int(tad['start'])
            end   = int(tad['end'])
        elif type(tad) == list:
            if type(tad[0]) == dict:
                start = int(sorted(tad,
                                   key=lambda x: int(x['start']))[0 ]['start'])
                end   = int(sorted(tad,
                                   key=lambda x: int(x['end'  ]))[-1]['end'  ])
        if relative:
            if normalized:
                # find minimum, if value is non-zero... for logarithm
                mini = min([i for i in xper.wght[0] if i])
                if mini == int(mini):
                    vmin = min(xper.wght[0])
                else:
                    vmin = mini
                vmin = fun(vmin or (1 if logarithm else 0))
                vmax = fun(max(xper.wght[0]))
            else:
                vmin = fun(min(xper.hic_data[0]) or (1 if logarithm else 0))
                vmax = fun(max(xper.hic_data[0]))
        plt.figure(figsize=(8, 6))
        if not axe:
            axe = plt.subplot(111)
        if tad or focus:
            if start > -1:
                if normalized:
                    matrix = [
                        [xper.wght[0][i+size*j]
                         if (not i in xper._zeros
                             and not j in xper._zeros) else vmin
                         for i in xrange(start - 1, end)]
                        for j in xrange(start - 1, end)]
                else:
                    matrix = [
                        [xper.hic_data[0][i+size*j]
                         for i in xrange(start - 1, end)]
                        for j in xrange(start - 1, end)]
            elif type(tad) is list:
                if normalized:
                    warn('List passed, not going to be normalized.')
                matrix = tad
            else:
                # TODO: something... matrix not declared...
                pass
        else:
            if normalized:
                matrix = [[xper.wght[0][i+size*j]
                           if (not i in xper._zeros
                               and not j in xper._zeros) else vmin
                           for i in xrange(size)]
                          for j in xrange(size)]
            else:
                matrix = [[xper.hic_data[0][i+size*j]\
                           for i in xrange(size)] \
                          for j in xrange(size)]
        if relative:
            img = axe.imshow(fun(matrix), origin='lower', vmin=vmin, vmax=vmax,
                             interpolation="nearest",
                             extent=(int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5,
                                     int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5))
        else:
            img = axe.imshow(fun(matrix), origin='lower',
                             interpolation="nearest",
                             extent=(int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5,
                                     int(start or 1) - 0.5,
                                     int(start or 1) + len(matrix) - 0.5))
        cbar = axe.figure.colorbar(img)
        cbar.ax.set_ylabel('Log2 Hi-C interactions count')
        axe.set_title(('Chromosome %s experiment %s' +
                       ' %s') % (self.name, xper.name,
                                 'focus: %s-%s' % (start, end) if tad else ''))
        axe.set_xlabel('Genomic bin (resolution: %s)' % (xper.resolution))
        axe.set_ylabel('Genomic bin (resolution: %s)' % (xper.resolution))
        if not paint_tads:            
            axe.set_ylim(int(start or 1) - 0.5,
                         int(start or 1) + len(matrix) - 0.5)
            axe.set_xlim(int(start or 1) - 0.5,
                         int(start or 1) + len(matrix) - 0.5)
            if show:
                plt.show()
            return img
        for i, tad in xper.tads.iteritems():
            if start:
                print int(tad['start']) + 1, start
                print int(tad['end']) + 1, end
                if int(tad['start']) + 1 < start:
                    continue
                if int(tad['end']) + 1 > end:
                    continue
                t_start = int(tad['start']) + 1.5
                t_end   = int(tad['end'])   + 2.5
            else:
                t_start = int(tad['start']) + .5
                t_end   = int(tad['end']) + 1.5
            axe.hlines(t_start, t_start, t_end, colors='k')
            axe.hlines(t_end, t_start, t_end, colors='k')
            axe.vlines(t_start, t_start, t_end, colors='k')
            axe.vlines(t_end, t_start, t_end, colors='k')
            if not i % (len(xper.tads) / 10):
                if i % 2:
                    axe.text(t_start + abs(t_start-t_end)/2,
                             t_end + 1, str(i), va='bottom', ha='center')
                else:
                    axe.text(t_start + abs(t_start-t_end)/2,
                             t_start - 1, str(i), va='top', ha='center')
            if tad['score'] < 0:
                for j in xrange(0, int(t_end) - int(t_start), 2):
                    axe.plot((t_start    , t_start + j),
                             (t_end   - j, t_end      ), color='k')
                    axe.plot((t_end      , t_end   - j),
                             (t_start + j, t_start    ), color='k')
        axe.set_ylim(int(start or 1) - 0.5,
                     int(start or 1) + len(matrix) - 0.5)
        axe.set_xlim(int(start or 1) - 0.5,
                     int(start or 1) + len(matrix) - 0.5)
        if show:
            plt.show()


    def get_tad_hic(self, tad, x_name, normed=True, matrix_num=0):
        """
        Retrieve the Hi-C data matrix corresponding to a given TAD.
        
        :param tad: a given TAD (:py:class:`dict`)
        :param x_name: name of the experiment
        :param True normed: if True, normalize the Hi-C data
        
        :returns: Hi-C data matrix for the given TAD
        """
        beg, end = int(tad['start']), int(tad['end'])
        xpr = self.get_experiment(x_name)
        size = xpr.size
        matrix = [[0 for _ in xrange(beg, end)]\
                  for _ in xrange(beg, end)]
        for i, tadi in enumerate(xrange(beg, end)):
            tadi = tadi * size
            for j, tadj in enumerate(xrange(beg, end)):
                if normed:
                    matrix[j][i] = xpr.hic_data[matrix_num][tadi + tadj]
                else:
                    matrix[j][i] = xpr.wght[0][tadi + tadj]
        return matrix


    def iter_tads(self, x_name, normed=True):
        """
        Iterate over the TADs corresponding to the given experiment.
        
        :param x_name: name of the experiment
        :param True normed: normalize Hi-C data returned
        
        :yields: Hi-C data corresponding to each TAD
        """
        if not self.get_experiment(x_name).hic_data:
            raise Exception('No Hi-c data for %s experiment\n' % (x_name))
        for name, ref in self.get_experiment(x_name).tads.iteritems():
            yield name, self.get_tad_hic(ref, x_name, normed=normed)


    def set_max_tad_size(self, value):
        """
        Change the maximum size allowed for TADs. It also applies to the
        computed experiments.

        :param value: an int value (default is 5000000)
        """
        self.max_tad_size = value
        for xpr in self.experiments:
            for tad in xpr.tads:
                xpr.tads[tad]['brk'] = xpr.tads[tad]['end']
                if ((xpr.tads[tad]['end'] - xpr.tads[tad]['start']) 
                    * xpr.resolution) > self.max_tad_size:
                    xpr.tads[tad]['score'] = -abs(xpr.tads[tad]['score'])
            

    def _search_centromere(self, xpr):
        """
        Search for the centromere in a chromosome, assuming that
        :class:`Chromosome` corresponds to a real chromosome.
        Add a boundary to all the experiments where the centromere is.
         * A centromere is defined as the largest area where the rows/columns
           of the Hi-C matrix are empty.
        """
        beg = end = 0
        size = xpr.size
        try:
            hic = xpr.hic_data[0]
        except TypeError:
            return
        # search for largest empty region of the chromosome
        best = (0, 0, 0)
        pos = 0
        for pos, raw in enumerate(xrange(0, size * size, size)):
            if sum(hic[raw:raw + size]) == 0 and not beg:
                beg = float(pos)
            if sum(hic[raw:raw + size]) != 0 and beg:
                end = float(pos)
                if (end - beg) > best[0]:
                    best = ((end - beg), beg, end)
                beg = end = 0
        # this is for weared cases where centromere is at the end of Hi-C data
        if beg and not end:
            end = float(pos)
            if (end - beg) > best[0]:
                best = ((end - beg), beg, end)
        beg, end = best[1:]
        if not beg or not end:
            return
        tads = xpr.tads
        # if we already have a centromere defined, check if it can be reduced
        if self._centromere:
            if beg > self._centromere[0]:
                # readjust TADs that have been split around the centromere
                for tad in tads:
                    if tads[tad]['end'] == self._centromere[0]:
                        tads[tad]['end'] = beg
                self._centromere[0] = beg
            if end < self._centromere[1]:
                # readjust TADs that have been split around the centromere
                for tad in tads:
                    if tads[tad]['start'] == self._centromere[1]:
                        tads[tad]['start'] = end
                self._centromere[1] = end
        else:
            self._centromere = [beg, end]
        # split TADs overlapping  with the centromere
        if [True for t in tads.values() \
            if t['start'] < beg < t['end'] \
            and t['start'] < end < t['end']]:
            tad  = len(tads) + 1
            plus = 0
            while tad + plus > 1:
                start = tads[tad - 1 + plus]['start']
                final = tads[tad - 1 + plus]['end']
                # centromere found?
                if start < beg < final and start < end < final:
                    tads[tad] = copy(tads[tad - 1])
                    tads[tad]['start'] = end
                    tads[tad]['score'] = -abs(tads[tad]['score'])
                    if (tads[tad]['end'] - tads[tad]['start']) \
                           * xpr.resolution > self.max_tad_size:
                        xpr.tads[tad]['score'] = -abs(xpr.tads[tad]['score'])
                    tads[tad]['brk'] = tads[tad]['end']
                    tad -= 1
                    tads[tad] = copy(tads[tad])
                    tads[tad]['score'] = -abs(tads[tad]['score'])
                    tads[tad]['end'] = beg
                    if (tads[tad]['end'] - tads[tad]['start']) \
                           * xpr.resolution > self.max_tad_size:
                        xpr.tads[tad]['score'] = -abs(xpr.tads[tad]['score'])
                    tads[tad]['brk'] = tads[tad]['end']
                    plus = 1
                else:
                    tads[tad] = copy(tads[tad - 1 + plus])
                tad -= 1


class ExperimentList(list):
    """
    Inherited from python built in :py:func:`list`, modified for tadbit
    :class:`pytadbit.Experiment`.
    
    Mainly, `getitem`, `setitem`, and `append` were modified in order to
    be able to search for experiments by index or by name, and to add 
    experiments simply using Chromosome.experiments.append(Experiment).

    The whole ExperimentList object is linked to a Chromosome instance
    (:class:`pytadbit.Chromosome`).

    """
    def __init__(self, thing, crm):
        super(ExperimentList, self).__init__(thing)
        self.crm = crm
        

    def __getitem__(self, i):
        try:
            return super(ExperimentList, self).__getitem__(i)
        except TypeError:
            for nam in self:
                if nam.name == i:
                    return nam
            raise KeyError('Experiment %s not found\n' % (i))


    def __setitem__(self, i, exp):
        try:
            super(ExperimentList, self).__setitem__(i, exp)
            exp.crm = self.crm
            self.crm._get_forbidden_region(exp)
        except TypeError:
            for j, nam in enumerate(self):
                if nam.name == i:
                    exp.crm = self.crm
                    self[j] = exp
                    self.crm._get_forbidden_region(exp)
                    break
            else:
                exp.crm = self.crm
                self.append(exp)
                self.crm._get_forbidden_region(exp)


    def __delitem__(self, i):
        try:
            super(ExperimentList, self).__delitem__(i)
        except TypeError:
            for j, nam in enumerate(self):
                if nam.name == i:
                    exp = self.pop(j)
                    del(exp)
                    break
            else:
                raise KeyError('Experiment %s not found\n' % (i))


    def append(self, exp):
        if exp.name in [e.name for e in self]:
            self[exp.name] = exp
            self.crm._get_forbidden_region(exp)
        else:
            super(ExperimentList, self).append(exp)
            self.crm._get_forbidden_region(exp)
            exp.crm = self.crm


class AlignmentDict(dict):
    """
    :py:func:`dict` of :class:`pytadbit.Alignment`
    
    Modified getitem, setitem, and append in order to be able to search
    alignments by index or by name.

    linked to a :class:`pytadbit.Chromosome`
    """

    def __getitem__(self, nam):
        try:
            return super(AlignmentDict, self).__getitem__(nam)
        except KeyError:
            for i, key in enumerate(self):
                if nam == i:
                    return self[key]
            raise TypeError('Alignment %s not found\n' % (i))


class ChromosomeSize(int):
    """
    This is an integer.
    
    Chromosome size in base pairs
    """
    def __init__(self, thing):
        super(ChromosomeSize, self).__init__(thing)


class RelativeChromosomeSize(int):
    """
    This is an integer.
    
    Relative Chromosome size in base pairs. Equal to Chromosome size minus
    forbidden regions (eg: the centromere)
    """
    def __init__(self, thing):
        super(RelativeChromosomeSize, self).__init__(thing)
