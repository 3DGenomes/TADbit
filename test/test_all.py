"""
30 Oct 2012

unittest for pytadbit functions
"""

import matplotlib
matplotlib.use('Agg')

import unittest
from pytadbit                             import Chromosome, load_chromosome
from pytadbit                             import tadbit, batch_tadbit
from pytadbit.tad_clustering.tad_cmo      import optimal_cmo
from pytadbit.modelling.structuralmodels        import load_structuralmodels
from pytadbit.modelling.impmodel                import load_impmodel_from_cmm
from pytadbit.eqv_rms_drms                import rmsdRMSD_wrapper
from pytadbit.parsers.genome_parser       import parse_fasta
from pytadbit.mapping.restriction_enzymes import map_re_sites, RESTRICTION_ENZYMES
from pytadbit.parsers.hic_parser          import load_hic_data_from_reads, read_matrix
from pytadbit.mapping.analyze             import hic_map, plot_distance_vs_interactions
from pytadbit.mapping.analyze             import insert_sizes, plot_iterative_mapping
from pytadbit.mapping.analyze             import correlate_matrices, eig_correlate_matrices
from pytadbit.mapping.filter              import filter_reads, apply_filter

from random                               import random, seed
from os                                   import system, path, chdir
from re                                   import finditer
from warnings                             import warn, catch_warnings, simplefilter
from distutils.spawn                      import find_executable

import sys


PATH = path.abspath(path.split(path.realpath(__file__))[0])

ONLY = None#'10'

def check_hic(hic, size):
    """
    check if hi-c data is symmetric
    """
    for i in xrange(size):
        for j in xrange(i + 1, size):
            if not hic[i * size + j] == hic[j * size + i]:
                raise AttributeError('ERROR: matrix should be symmetric.\n')
    return True


class TestTadbit(unittest.TestCase):
    """
    test main tadbit functions
    """
    
    def test_01_tadbit(self):

        print 'PYTHON SIDE'
        print '-----------'

        # if ONLY and ONLY != '01':
        #     return
        
        if CHKTIME:
            t0 = time()

        
        global exp1, exp2, exp3, exp4
        exp1 = tadbit(PATH + '/40Kb/chrT/chrT_A.tsv', max_tad_size="max",
                      verbose=False, no_heuristic=False, n_cpus='max')
        exp2 = tadbit(PATH + '/20Kb/chrT/chrT_B.tsv', max_tad_size="max",
                      verbose=False, no_heuristic=False, n_cpus='max')
        exp3 = tadbit(PATH + '/20Kb/chrT/chrT_C.tsv', max_tad_size="max",
                      verbose=False, no_heuristic=False, n_cpus='max')
        exp4 = tadbit(PATH + '/20Kb/chrT/chrT_D.tsv', max_tad_size="max",
                      n_cpus='max',
                      verbose=False, no_heuristic=False, get_weights=True)

        # Breaks and scores with square root normalization.
        #breaks = [0, 4, 10, 15, 23, 29, 38, 45]
        #scores = [7.0, 7.0, 5.0, 7.0, 4.0, 6.0, 8.0, None]
        breaks = [0, 4, 10, 15, 20, 25, 31, 36, 45]
        scores = [7.0, 7.0, 4.0, 4.0, 4.0, 4.0, 4.0, 7.0, None]
        self.assertEqual(exp1['start'], breaks)
        self.assertEqual(exp1['score'], scores)

        if CHKTIME:
            print '1', time() - t0


    def test_02_batch_tadbit(self):
        if ONLY and ONLY != '02':
            return
        if CHKTIME:
            t0 = time()

        global batch_exp
        batch_exp = batch_tadbit(PATH + '/20Kb/chrT/', max_tad_size=20, 
                                 verbose=False, no_heuristic=True)
        # Breaks and scores with square root normalization.
        breaks = [0, 4, 14, 19, 34, 39, 44, 50, 62, 67, 72, 90, 95]
        scores = [4.0, 6.0, 5.0, 5.0, 4.0, 8.0, 5.0, 4.0, 6.0, 5.0,
                  6.0, 6.0, None]
        self.assertEqual(batch_exp['start'], breaks)
        self.assertEqual(batch_exp['score'], scores)
        if CHKTIME:
            print '2', time() - t0


    def test_03_tad_multi_aligner(self):

        if ONLY and ONLY != '03':
            return
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', centromere_search=True,
                              experiment_tads=[exp1, exp2, exp3, exp4],
                              experiment_hic_data=[
                                  PATH + '/40Kb/chrT/chrT_A.tsv',
                                  PATH + '/20Kb/chrT/chrT_B.tsv',
                                  PATH + '/20Kb/chrT/chrT_C.tsv',
                                  PATH + '/20Kb/chrT/chrT_D.tsv'],
                              experiment_names=['exp1', 'exp2', 'exp3', 'exp4'],
                              experiment_resolutions=[40000,20000,20000,20000],
                              silent=True)
        for exp in test_chr.experiments:
            exp.normalize_hic(silent=True, factor=None)

        test_chr.align_experiments(verbose=False, randomize=False,
                                   method='global')
        _, (score1, pval1,
            perc1, perc2) = test_chr.align_experiments(verbose=False,
                                                       method='global',
                                                       randomize=True, rnd_num=100)
        _, (_, pval2,
            perc1, perc2) = test_chr.align_experiments(verbose=False, randomize=True,
                                                       rnd_method='shuffle', rnd_num=100)
        # Values with alignments obtained with square root normalization.
        #self.assertEqual(round(-26.095, 3), round(score1, 3))
        #self.assertEqual(round(0.001, 1), round(pval1, 1))
        #self.assertTrue(abs(0.175 - pval2) < 0.2)
        self.assertEqual(round(-11.002, 3), round(score1, 3))
        self.assertEqual(round(0.001, 1), round(pval1, 1))
        self.assertTrue(abs(0.04 - pval2) < 0.1)
        if CHKTIME:
            print '3', time() - t0

                              
    def test_04_chromosome_batch(self):
        if ONLY and ONLY != '04':
            return
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome',
                              experiment_resolutions=[20000]*3,
                              experiment_hic_data=[
                                  PATH + '/20Kb/chrT/chrT_A.tsv',
                                  PATH + '/20Kb/chrT/chrT_D.tsv',
                                  PATH + '/20Kb/chrT/chrT_C.tsv'],
                              experiment_names=['exp1', 'exp2', 'exp3'],
                              silent=True)
        test_chr.find_tad(['exp1', 'exp2', 'exp3'], batch_mode=True,
                          verbose=False, silent=True)
        tads = test_chr.get_experiment('batch_exp1_exp2_exp3').tads
        found = [tads[t]['end'] for t in tads if tads[t]['score'] > 0]
        # Values obtained with square root normalization.
        #self.assertEqual([3.0, 8.0, 16.0, 21.0, 28.0, 35.0, 43.0,
        #                  49.0, 61.0, 66.0, 75.0, 89.0, 94.0, 99.0], found)
        self.assertEqual([3.0, 14.0, 19.0, 33.0, 43.0, 49.0, 61.0, 66.0,
                           71.0, 89.0, 94.0, 99.0], found)
        
        if CHKTIME:
            print '4', time() - t0


    def test_05_save_load(self):
        if ONLY and ONLY != '05':
            return
        if CHKTIME:
            t0 = time()

        test_chr1 = Chromosome(name='Test Chromosome',
                               experiment_tads=[exp1, exp2],
                               experiment_names=['exp1', 'exp2'],
                               experiment_resolutions=[20000,20000],
                               silent=True)
        test_chr1.save_chromosome('lolo', force=True)
        test_chr2 = load_chromosome('lolo')
        system('rm -f lolo')
        system('rm -f lolo_hic')
        self.assertEqual(str(test_chr1.__dict__), str(test_chr2.__dict__))
        if CHKTIME:
            print '5', time() - t0


    def test_06_tad_clustering(self):
        if ONLY and ONLY != '06':
            return
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome',
                              experiment_tads=[exp4],
                              experiment_names=['exp1'],
                              experiment_hic_data=[
                                  PATH + '/20Kb/chrT/chrT_D.tsv'],
                              experiment_resolutions=[20000,20000],
                              silent=True)
        all_tads = []
        for _, tad in test_chr.iter_tads('exp1', normed=False):
            all_tads.append(tad)
        #align1, align2, _ = optimal_cmo(all_tads[7], all_tads[10], 7,
        #                                method='score')
        align1, align2, _ = optimal_cmo(all_tads[1], all_tads[3], 7,
                                        method='score')
        # Values with square root normalization.
        #self.assertEqual(align1, [0, 1, '-', 2, 3, '-', 4, 5, 6, 7, 8, 9, 10])
        #self.assertEqual(align2,[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        self.assertEqual(align1, [0, 1, 2, '-', '-', 3, 4, 5, 6, 7, 8, '-', 9])
        self.assertEqual(align2, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        if CHKTIME:
            print '6', time() - t0
        

    def test_07_forbidden_regions(self):
        if ONLY and ONLY != '07':
            return
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000,
                              centromere_search=True,)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv',
                                silent=True)
        # Values with square root normalization.
        #brks = [2.0, 7.0, 12.0, 18.0, 38.0, 43.0, 49.0,
        #        61.0, 66.0, 75.0, 89.0, 94.0, 99.0]
        brks = [3.0, 14.0, 19.0, 33.0, 38.0, 43.0, 49.0, 61.0,
                  66.0, 71.0, 83.0, 89.0, 94.0, 99.0]
        tads = test_chr.experiments['exp1'].tads
        found = [tads[t]['end'] for t in tads if tads[t]['score'] > 0]
        self.assertEqual(brks, found)
        items1 = test_chr.forbidden.keys(), test_chr.forbidden.values()
        test_chr.add_experiment('exp2', 20000, tad_def=exp3,
                                hic_data=PATH + '/20Kb/chrT/chrT_C.tsv',
                                silent=True)
        items2 = test_chr.forbidden.keys(), test_chr.forbidden.values()
        know1 = ([38, 39], ['Centromere', 'Centromere'])
        #know1 = ([32, 33, 34, 38, 39, 19, 20, 21, 22,
        #          23, 24, 25, 26, 27, 28, 29, 30, 31],
        #         [None, None, None, 'Centromere', 'Centromere',
        #          None, None, None, None, None, None, None,
        #          None, None, None, None, None, None])
        know2 = ([38], ['Centromere'])
        self.assertEqual(items1, know1)
        self.assertEqual(items2, know2)
        if CHKTIME:
            print '7', time() - t0


    def test_08_changing_resolution(self):
        if ONLY and ONLY != '08':
            return
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv',
                                silent=True)
        exp = test_chr.experiments['exp1']
        sum20 = sum(exp.hic_data[0].values())
        exp.set_resolution(80000)
        sum80 = sum(exp.hic_data[0].values())
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(160000)
        sum160 = sum(exp.hic_data[0].values())
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(360000)
        sum360 = sum(exp.hic_data[0].values())
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(2400000)
        sum2400 = sum(exp.hic_data[0].values())
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(40000)
        sum40 = sum(exp.hic_data[0].values())
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(20000)
        sum21 = sum(exp.hic_data[0].values())
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(40000)
        sum41 = sum(exp.hic_data[0].values())
        check_hic(exp.hic_data[0], exp.size)
        self.assertTrue(sum20 == sum80 == sum160 == sum360 == sum40 \
                        == sum21 == sum2400 == sum41)
        if CHKTIME:
            print '8', time() - t0


    def test_09_hic_normalization(self):
        """
        writes interaction pair file.
        """
        if ONLY and ONLY != '09':
            return
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv',
                                silent=True)
        exp = test_chr.experiments[0]
        exp.load_hic_data(PATH + '/20Kb/chrT/chrT_A.tsv', silent=True)
        exp.normalize_hic(silent=True)
        exp.get_hic_zscores()
        exp.get_hic_zscores(zscored=False)
        sumz = sum([exp._zscores[k1][k2] for k1 in exp._zscores.keys()
                    for k2 in exp._zscores[k1]])
        self.assertEqual(round(sumz, 4), round(4059.2877, 4))
        if CHKTIME:
            print '9', time() - t0


    def test_10_compartments(self):
        """
        """
        if ONLY and ONLY != '10':
            return
        if CHKTIME:
            t0 = time()
        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv',
                                silent=True)
        exp = test_chr.experiments[0]
        exp.load_hic_data(PATH + '/20Kb/chrT/chrT_A.tsv', silent=True)
        hic_data = exp.hic_data[0]
        hic_data.find_compartments(label_compartments='cluster')
        self.assertEqual(len(hic_data.compartments[None]), 39)
        # self.assertEqual(round(hic_data.compartments[None][24]['dens'], 5),
        #                  0.75434)
        if CHKTIME:
            print '10', time() - t0
        

    # def test_10_generate_weights(self):
    #     """
    #     """
    #     if CHKTIME:
    #         t0 = time()

    #     test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
    #     test_chr.add_experiment('exp1', 20000, tad_def=exp4,
    #                             hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
    #     exp = test_chr.experiments[0]
    #     B = iterative(exp.hic_data[0], bads=[False] * exp.size)
    #     tadbit_weights = [float(exp.hic_data[0][i, j]) / B[i] / B[j] * exp.size
    #                       for i in B for j in B]
    #     exp.normalize_hic(factor=None, iterations=0)
    #     self.assertEqual([round(i, 3) for i in tadbit_weights[:100]],
    #                      [round(exp.norm[0][i], 3)
    #                       for i in xrange(100)])
    #     if CHKTIME:
    #         print '10', time() - t0


    def test_11_write_interaction_pairs(self):
        if ONLY and ONLY != '11':
            return
        """
        writes interaction pair file.
        """
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
        exp = test_chr.experiments[0]
        exp.load_hic_data(PATH + '/20Kb/chrT/chrT_A.tsv', silent=True)
        exp.filter_columns(silent=True)
        exp.normalize_hic(factor=None, silent=True)
        exp.get_hic_zscores(zscored=False)
        exp.write_interaction_pairs('lala')
        lines = open('lala').readlines()
        self.assertEqual(len(lines), 4674)
        self.assertEqual(lines[25], '1\t28\t0.612332461036\n')
        self.assertEqual(lines[2000], '26\t70\t0.0738742984321\n')
        system('rm -f lala')
        if CHKTIME:
            print '11', time() - t0


    def test_12_3d_modelling_optimization(self):
        """
        quick test to generate 3D coordinates from 3? simple models???
        """
        if ONLY and ONLY != '12':
            return
        if CHKTIME:
            t0 = time()

        try:
            __import__('IMP')
        except ImportError:
            warn('IMP not found, skipping test\n')
            return
        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
        exp = test_chr.experiments[0]
        exp.load_hic_data(PATH + '/20Kb/chrT/chrT_A.tsv')
        exp.filter_columns(silent=True)
        exp.normalize_hic(silent=True, factor=None)
        result = exp.optimal_imp_parameters(50, 70, n_cpus=4,
                                            n_models=8, n_keep=2,
                                            lowfreq_range=[-0.6],
                                            upfreq_range=(0, 1.1, 1.1),#from 0 till 1.1 in step of 1.1 with ()
                                            maxdist_range=[500, 600],# it will use 500 and 600 with []
                                            verbose=False)

        # get best correlations
        config = result.get_best_parameters_dict()#dict with parameters
        wanted = {'maxdist': 600.0, 'upfreq': 0.0, 'kforce': 5,
                  'dcutoff': 2,
                  'reference': '', 'lowfreq': -0.6, 'scale': 0.01}
        self.assertEqual([round(i, 4) for i in config.values()if not type(i) is str],
                         [round(i, 4) for i in wanted.values()if not type(i) is str])
        if CHKTIME:
            print '12', time() - t0


    def test_13_3d_modelling_centroid(self):#model with no optimisation
        """
        quick test to generate 3D coordinates from 3? simple models???
        """
        if ONLY and ONLY != '13':
            return
        if CHKTIME:
            t0 = time()

        try:
            __import__('IMP')
        except ImportError:
            warn('IMP not found, skipping test\n')
            return
        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv',
                                silent=True)
        exp = test_chr.experiments[0]
        exp.load_hic_data(PATH + '/20Kb/chrT/chrT_A.tsv', silent=True)
        exp.filter_columns(silent=True)
        exp.normalize_hic(silent=True, factor=None)
        models = exp.model_region(51, 71, n_models=40, n_keep=25,
                                  n_cpus=4,
                                  config={'kforce': 5, 'maxdist': 500,
                                          'scale': 0.01,
                                          'upfreq': 1.0, 'lowfreq': -0.6})
        models.save_models('models.pick')

        avg = models.average_model()
        nmd = len(models)
        dev = rmsdRMSD_wrapper(
            [models[m]['x'] for m in xrange(nmd)] + [avg['x']],
            [models[m]['y'] for m in xrange(nmd)] + [avg['y']],
            [models[m]['z'] for m in xrange(nmd)] + [avg['z']],
            models._zeros,
            models.nloci, 200, range(len(models)+1),
            len(models)+1, int(False), 'rmsd', 0)
        centroid = models[models.centroid_model()]
        # find closest
        model = min([(k, dev[(k, nmd)] )
                     for k in range(nmd)], key=lambda x: x[1])[0]
        self.assertEqual(centroid['rand_init'], models[model]['rand_init'])
        if CHKTIME:
            print '13', time() - t0


    def test_14_3d_clustering(self):
        """
        """
        if ONLY and ONLY != '14':
            return
        if CHKTIME:
            t0 = time()

        models = load_structuralmodels('models.pick')
        if find_executable('mcl'):
            models.cluster_models(method='mcl', fact=0.9, verbose=False,
                                  dcutoff=200)
            self.assertTrue(5 <= len(models.clusters.keys()) <= 7)
        models.cluster_models(method='ward', verbose=False, dcutoff=200)
        self.assertTrue(2 <= len(models.clusters.keys()) <= 3)
        d = models.cluster_analysis_dendrogram()
        self.assertEqual(d['icoord'], [[5., 5., 15., 15.]])
        # align models
        m1, m2 = models.align_models(models=[1,2])
        nrmsd = (sum([((m1[0][i] - m2[0][i])**2 + (m1[1][i] - m2[1][i])**2 + (m1[2][i] - m2[2][i])**2)**.5
                      for i in xrange(len(m1[0]))]) / (len(m1[0])))
        self.assertTrue(nrmsd < 160)
        # fetching models
        models.define_best_models(5)
        m = models.fetch_model_by_rand_init('1', all_models=True)
        self.assertEqual(m, 6)
        models.define_best_models(25)
        m = models.fetch_model_by_rand_init('1', all_models=False)
        self.assertEqual(m, 6)
        if CHKTIME:
            print '14', time() - t0


    def test_15_3d_modelling(self):
        """
        """
        if ONLY and ONLY != '15':
            return
        if CHKTIME:
            t0 = time()

        models = load_structuralmodels('models.pick') 
        models.cluster_models(method='ward', verbose=False)
        # density
        models.density_plot(savedata='lala', plot=False)
        lines = open('lala').readlines()
        self.assertEqual(len(lines), 22)
        self.assertEqual([round(float(i), 1) if i != 'nan' else i for i in lines[1].split('\t')[:3]],
                         [1.0, 'nan', 'nan'])
        self.assertEqual([round(float(i), 1) for i in lines[15].split('\t')[:3]],
                         [15, 100.0, 0.0])
        # contacts
        cmap = models.get_contact_matrix(cutoff=300)
        self.assertEqual(round(
            round(sum([i if i >=0 else 0 for i in
                       reduce(lambda x, y: x+y, cmap)])/10, 0),
            3), 8)
        # define best models
        models.define_best_models(10)
        self.assertEqual(len(models), 10)
        m1 = models[9]
        models.define_best_models(25)
        self.assertEqual(len(models), 25)
        self.assertEqual(m1, models[9])
        # correlation
        corr, pval = models.correlate_with_real_data(cutoff=300)
        self.assertTrue(0.6 <= round(corr, 1) <= 0.7)
        self.assertEqual(round(pval, 4), round(0, 4))
        # consistency
        models.model_consistency(cutoffs=(50, 100, 150, 200), plot=False,
                                 savedata='lala')
        lines = open('lala').readlines()
        self.assertEqual(len(lines), 22)
        self.assertEqual([round(float(i)/15, 0) for i in lines[1].split('\t')],
                         [0, 2, 3, 3, 3])
        self.assertEqual([round(float(i)/15, 0) for i in lines[15].split('\t')],
                         [1, 5, 6, 7, 7])
        # measure angle
        self.assertTrue(13 <= round(models.angle_between_3_particles(2,8,15)/10,
                                    0) <= 14)
        self.assertEqual(round(models.angle_between_3_particles(19,20,21), 0),
                         60)
        self.assertEqual(round(models.angle_between_3_particles(15,14,11)/5, 0),
                         13)
        # coordinates
        self.assertEqual([round(x, 2) for x in models.particle_coordinates(15)],
                         [2098.32, 1565.63, -4319.62])
        # dihedral_angle
        self.assertTrue (round(models.dihedral_angle(2 ,  8, 15,  8, 16, [0])[0], 2), -13.44)
        self.assertEqual(round(models.dihedral_angle(15, 19, 20, 19, 21, [0])[0], 2),  76.24)
        self.assertEqual(round(models.dihedral_angle(15, 14, 11, 14, 12, [0])[0], 2),   0.07)
        # median distance
        self.assertEqual(round(models.median_3d_dist(3, 20, plot=False)/100, 0),
                         15)
        self.assertEqual(round(models.median_3d_dist(3, 20, cluster=1,
                                                     plot=False)/200, 0), 8)
        self.assertEqual(round(models.median_3d_dist(7, 10, models=range(5),
                                                     plot=False), 0), 250)
        # accessibility
        models.accessibility(radius=75, nump=10, plot=False, savedata='model.acc')
        vals = [l.split() for l in open('model.acc').readlines()[1:]]
        self.assertEqual(vals[0][1:3], ['0.68', '0.933'])
        self.assertEqual(vals[20][1:3], ['1.0', '0.0'])
        # contact map
        models.contact_map(savedata='model.contacts')
        vals = [l.split() for l in open('model.contacts').readlines()[1:]]
        self.assertEqual(vals[0], ['0', '1', '1.0'])
        self.assertEqual(vals[1], ['0', '2', '0.96'])
        self.assertEqual(vals[192], ['14', '18', '0.12'])
        # interactions
        models.interactions(plot=False, savedata='model.inter')
        vals = [[float(i) for i in l.split()] for l in open('model.inter').readlines()[1:]]
        self.assertEqual(vals[2], [3.0, 4.92, 1.12, 3.88, 0.65, 4.69, 0.82, 4.01, 0.62, 4.81, 0.5])
        # walking angle
        models.walking_angle(savedata='model.walkang')
        vals = [[round(float(i), 2) if i != 'None' else i for i in l.split()] for l in open('model.walkang').readlines()[1:]]
        self.assertEqual(vals[17], [18.0, -45.42, 100.0, -14.14, 137.0],)
        self.assertEqual(vals[3],  [4.0, 125.36, 273.0, 1.96, 253.0],)
        self.assertEqual(vals[16], [17.0, -70.14, 200.0, -2.5, 84.0])
        self.assertEqual(vals[15], [16.0, -134.88, 287.0, -20.48, 121.0])
        # write cmm
        models.write_cmm('.', model_num=2)
        models.write_cmm('.', models=range(5))
        models.write_cmm('.', cluster=2)
        # write xyz
        models.write_xyz('.', model_num=2)
        models.write_xyz('.', models=range(5))
        models.write_xyz('.', cluster=2)
        # write json
        models.write_json('model.json', model_num=2)
        models.write_json('model.json', models=range(5))
        models.write_json('model.json', cluster=2)
        # clean
        system('rm -f model.*')
        system('rm -rf lala*')
        if CHKTIME:
            print '15', time() - t0


    def test_16_models_stats(self):
        if ONLY and ONLY != '16':
            return
        if CHKTIME:
            t0 = time()

        models = load_structuralmodels('models.pick') 
        # write cmm
        models.write_cmm('.', model_num=2)
        model = load_impmodel_from_cmm('model.%s.cmm' % models[2]['rand_init'])
        # clean
        system('rm -f model.*')
        # stats
        self.assertEqual(200, round(model.distance(2, 3), 0))
        self.assertTrue(9 <= round(model.distance(8, 20)/100, 0) <= 10)
        self.assertEqual(round(30, 0),
                         round(model.radius_of_gyration() / 20, 0))
        self.assertEqual(400, round(model.contour()/10, 0))
        self.assertTrue(21 <= round((model.shortest_axe() +
                                     model.longest_axe()) / 100,
                                    0) <= 22)
        self.assertEqual([11, 15, 16], model.inaccessible_particles(1000))

        acc, num, acc_area, tot_area, bypt = model.accessible_surface(
            150, superradius=200, nump=150)
        self.assertTrue(210 <= acc <= 240)
        self.assertTrue(500 <= num<= 600)
        self.assertEqual(0.4, round(acc_area, 1))
        self.assertEqual(4, round(tot_area, 0))
        self.assertEqual(101, len(bypt))
        self.assertTrue(19 <= bypt[100][0] <= 22 and
                         8 <= bypt[100][1] <= 38 and
                         8 <= bypt[100][2] <= 23)
        if CHKTIME:
            print '16', time() - t0


    def test_17_map_re_sites(self):
        """
        test fasta parsing and mapping re sites
        """
        if ONLY and ONLY != '17':
            return
        if CHKTIME:
            t0 = time()
        ref_genome = parse_fasta(PATH + '/ref_genome/chr2L_chr4_dm3.bz2',
                                 verbose=False)
        self.assertEqual(len(ref_genome['chr4']), 1351857)
        frags = map_re_sites('dpnIi', ref_genome)
        self.assertEqual(len(frags['chr2L']), 231)
        self.assertEqual(len(frags['chr2L'][230]), 16)
        self.assertEqual(frags['chr4'][10][50], 1018069)
        frags = map_re_sites('hindiii', ref_genome)
        self.assertEqual(len(frags['chr2L']), 231)
        self.assertEqual(len(frags['chr2L'][230]), 3)
        self.assertEqual(frags['chr4'][10][5], 1017223)
        if CHKTIME:
            self.assertEqual(True, True)
            print '17', time() - t0

    def test_18_filter_reads(self):
        if ONLY and ONLY != '18':
            return
        if CHKTIME:
            t0 = time()
        for ali in ['map', 'sam']:
            seed(1)
            if 13436 == int(random()*100000):
                same_seed = True
                genome = generate_random_ali(ali)
                genome_bis = parse_fasta('test.fa~', verbose=False)
                self.assertEqual(genome, genome_bis)
            else:
                same_seed = False
                genome = parse_fasta('test.fa~')
            # PARSE SAM
            if ali == 'map':
                from pytadbit.parsers.map_parser import parse_map as parser
            else:
                try:
                    from pytadbit.parsers.sam_parser import parse_sam as parser
                except ImportError:
                    print 'ERROR: PYSAM not found, skipping test\n'
                    continue

            parser(['test_read1.%s~' % (ali)], ['test_read2.%s~' % (ali)],
                   './lala1-%s~' % (ali), './lala2-%s~' % (ali), genome,
                   re_name='DPNII', mapper='GEM')

            # GET INTERSECTION
            from pytadbit.mapping import get_intersection
            get_intersection('lala1-%s~' % (ali), 'lala2-%s~' % (ali),
                             'lala-%s~' % (ali))
            # FILTER
            masked = filter_reads('lala-%s~' % (ali), verbose=False,
                                  fast=(ali=='map'))
            self.assertEqual(masked[1]['reads'], 1000)
            self.assertEqual(masked[2]['reads'], 1000)
            self.assertEqual(masked[3]['reads'], 1000)
            self.assertEqual(masked[4]['reads'], 1000)
            if same_seed:
                self.assertEqual(masked[5]['reads'], 1110)
                self.assertEqual(masked[6]['reads'], 2332)
                self.assertEqual(masked[7]['reads'], 0)
                self.assertEqual(masked[8]['reads'], 141)
                self.assertEqual(masked[10]['reads'], 1)
            else:
                self.assertTrue (masked[5]['reads'] > 1000)
            self.assertEqual(masked[9]['reads'], 1000)
        apply_filter('lala-map~', 'lala-map-filt~', masked, filters=[1],
                     reverse=True, verbose=False)
        self.assertEqual(len([True for l in open('lala-map-filt~')
                              if not l.startswith('#')]), 1000)
        d = plot_iterative_mapping('lala1-map~', 'lala2-map~')
        self.assertEqual(d[0][1], 6000)

        if CHKTIME:
            self.assertEqual(True, True)
            print '18', time() - t0

    def test_19_matrix_manip(self):
        if ONLY and ONLY != '19':
            return
        if CHKTIME:
            t0 = time()
        hic_data1 = load_hic_data_from_reads('lala-map~', resolution=10000)
        hic_map(hic_data1, savedata='lala-map.tsv~', savefig='lala.pdf~')
        hic_map(hic_data1, by_chrom='intra', savedata='lala-maps~', savefig='lalalo~')
        hic_map(hic_data1, by_chrom='inter', savedata='lala-maps~', savefig='lalala~')
        # slowest part of the all test:
        hic_data2 = read_matrix('lala-map.tsv~', resolution=10000)
        self.assertEqual(hic_data1, hic_data2)
        vals = plot_distance_vs_interactions(hic_data1)
        
        self.assertEqual([round(i, 2) if str(i)!='nan' else 0.0 for i in
                          reduce(lambda x, y: x + y, vals)],
                         [-1.68, -2.08, 0.02, 2.76, -8.99, 0.0, 0.82, -6.8, 0.0])
        
        a, b = insert_sizes('lala-map~')
        self.assertEqual([int(a),int(b)], [43, 1033])

        hic_data1 = read_matrix('20Kb/chrT/chrT_A.tsv', resolution=20000)
        hic_data2 = read_matrix('20Kb/chrT/chrT_B.tsv', resolution=20000)
        
        corr = correlate_matrices(hic_data1, hic_data2)
        corr =  [round(i,3) for i in corr[0]]
        self.assertEqual(corr, [0.755, 0.729, 0.804, 0.761, 0.789, 0.776, 0.828,
                                0.757, 0.797, 0.832])
        
        ecorr = eig_correlate_matrices(hic_data1, hic_data2)
        ecorr = [round(i,3) for i in reduce(lambda x, y:x+y, ecorr)]
        self.assertEqual(ecorr, [0.997, 0.322, 0.442, 0.017, 0.243, 0.014,
                                 0.321, 0.999, 0.01, 0.006, 0.0, 0.007, 0.451,
                                 0.012, 0.996, 0.031, 0.013, 0.004, 0.002,
                                 0.006, 0.029, 0.974, 0.076, 0.03, 0.219, 0.013,
                                 0.031, 0.08, 0.974, 0.018, 0.028, 0.004, 0.0,
                                 0.028, 0.034, 0.89])
        system('rm -rf lala*')
        if CHKTIME:
            self.assertEqual(True, True)
            print '19', time() - t0

    def test_20_tadbit_c(self):
        """
        Runs tests written in c, around the detection of TADs
        """
        if ONLY and ONLY != '20':
            return
        if CHKTIME:
            t0 = time()

        print '\n\nC SIDE'
        print '------'
        chdir(PATH + '/../src/test/')
        system('make clean')
        system('make')
        return_code = system('make test')
        chdir(PATH)
        if return_code != 0:
            print 'ERROR problem with C test'
            self.assertEqual(True, False)
        if CHKTIME:
            self.assertEqual(True, True)
            print '20', time() - t0


def generate_random_ali(ali='map'):
    # VARIABLES
    num_crms      = 9
    mean_crm_size = 1000000

    r_enz         = 'DpnII'

    selfcircle    = 1000
    dangling      = 1000
    error         = 1000
    extradangling = 1000
    duplicates    = 1000

    re_seq = RESTRICTION_ENZYMES[r_enz].replace('|', '')
    enz_cut = RESTRICTION_ENZYMES[r_enz].index('|')

    nts = ('ATGC' * 100) +'N'

    # RANDOM GENOME GENERATION
    genome = {}
    for crm in xrange(1, num_crms + 1):
        crm_len = int(mean_crm_size * random())
        genome['chr' + str(crm)] = ''.join([nts[int(401 * random())]
                                            for _ in xrange(crm_len)])

    out = open('test.fa~', 'w')
    for crm in xrange(1, num_crms + 1):
        out.write('>chr%d\n' % crm)
        crm = 'chr' + str(crm)
        for p in xrange(0, len(genome[crm]), 60):
            out.write(genome[crm][p:p+60] + '\n')
    out.close()

    # RE FRAGMENTS
    frags = {}
    for crm in genome:
        frags[crm] = {}
        beg = 0
        for pos in finditer(re_seq, genome[crm]):
            end = pos.start() + 1 + enz_cut
            if beg == end:
                continue
            frags[crm][beg] = [beg, end]
            beg = end
        if beg != end:
            frags[crm][beg] = [beg, len(genome[crm])]

    # Prepare files
    sam_crm = '@SQ\tSN:%s\tLN:%d\n'
    sam_head = """@HD\tVN:1.3\n%s@RG\tID:0\tPG:GEM\tPL:ILLUMINA\tSM:0\n@PG\tID:GEM\tPN:gem-2-sam\tVN:1.847\n"""
    crm_heads = ''
    for crm in genome:
        crm_heads += sam_crm % (crm, len(genome[crm]))

    sam_head = sam_head % crm_heads

    if ali=='sam':
        read = '{id}\t{flag}\t{crm}\t{pos}\t254\t3M\t*\t0\t0\tAAA\tHHH\tRG:Z:0\tNH:i:1\tNM:i:0\tXT:A:U\tmd:Z:3\n'
    else:
        read = "{id}\tAAA\tHHH\t1\t{crm}:{flag}:{pos}:3\n"

    out1 = open('test_read1.%s~' % (ali), 'w')
    out2 = open('test_read2.%s~' % (ali), 'w')
    if ali == 'sam':
        out1.write(sam_head)
        out2.write(sam_head)

    flags = ['+', '-'] if ali=='map' else [66 , 82 ]

    # SELF-CIRCLES
    for i in xrange(selfcircle):
        crm  = genome.keys()    [int(random() * len(genome))]
        while True:
            frag = frags[crm].keys()[int(random() * len(frags[crm]))]
            if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
                break
        while True:
            pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                       + frags[crm][frag][0])
            pos2 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                       + frags[crm][frag][0])
            if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
                continue
            if pos2 > pos1:
                sd1 = 1
                sd2 = 0
                pos1 -= 2 if ali=='map' else 3
                break
            elif pos2 < pos1:
                sd1 = 0
                sd2 = 1
                pos2 -= 2 if ali=='map' else 3
                break
        read1 = {'crm': crm, 'pos': pos1, 'flag': flags[sd1], 'id': 'lala01.%012d' % (i)}
        read2 = {'crm': crm, 'pos': pos2, 'flag': flags[sd2], 'id': 'lala01.%012d' % (i)}
        out1.write(read.format(**read1))
        out2.write(read.format(**read2))

    # DANGLING-ENDS
    for i in xrange(dangling):
        crm  = genome.keys()    [int(random() * len(genome))]
        while True:
            frag = frags[crm].keys()[int(random() * len(frags[crm]))]
            if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
                break
        while True:
            pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                       + frags[crm][frag][0])
            pos2 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                       + frags[crm][frag][0])
            if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
                continue
            if pos2 > pos1:
                sd1 = 0
                sd2 = 1
                pos2 -= 2 if ali=='map' else 3
                break
            elif pos2 < pos1:
                sd1 = 1
                sd2 = 0
                pos1 -= 2 if ali=='map' else 3
                break
        read1 = {'crm': crm, 'pos': pos1, 'flag': flags[sd1], 'id': 'lala02.%012d' % (i)}
        read2 = {'crm': crm, 'pos': pos2, 'flag': flags[sd2], 'id': 'lala02.%012d' % (i)}
        out1.write(read.format(**read1))
        out2.write(read.format(**read2))

    # ERRORS
    for i in xrange(error):
        crm  = genome.keys()    [int(random() * len(genome))]
        while True:
            frag = frags[crm].keys()[int(random() * len(frags[crm]))]
            if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
                break
        while True:
            pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                       + frags[crm][frag][0])
            pos2 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                       + frags[crm][frag][0])
            sd1 = sd2 = int(random()*2)
            if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
                continue
            if pos1 != pos2:
                if sd1 == 1:
                    pos1 -= 2 if ali=='map' else 3
                    pos2 -= 2 if ali=='map' else 3
                break
        read1 = {'crm': crm, 'pos': pos1, 'flag': flags[sd1], 'id': 'lala03.%012d' % (i)}
        read2 = {'crm': crm, 'pos': pos2, 'flag': flags[sd2], 'id': 'lala03.%012d' % (i)}
        out1.write(read.format(**read1))
        out2.write(read.format(**read2))

    # EXTRA-DANGLING-ENDS
    for i in xrange(extradangling):
        crm  = genome.keys()    [int(random() * len(genome))]
        while True:
            frag = frags[crm].keys()[int(random() * len(frags[crm]))]
            if (frags[crm][frag][1] - frags[crm][frag][0]) > 6:
                break
        while True:
            while True:
                pos1 = int(random() * (frags[crm][frag][1] - frags[crm][frag][0])
                           + frags[crm][frag][0])
                if frags[crm][frag][1] - pos1 < 490 or pos1 - frags[crm][frag][0] < 490:
                    break
            while True:
                pos2 = int(random() * 999) - 499 + pos1
                if not frags[crm][frag][0] <= pos2 <= frags[crm][frag][1]:
                    break
            if not (3 < pos1 < len(genome[crm]) - 3 and 3 < pos2 < len(genome[crm]) - 3):
                continue
            if pos2 > pos1:
                sd1 = 0
                sd2 = 1
                pos2 -= 2 if ali=='map' else 3
                break
            elif pos2 < pos1:
                sd1 = 1
                sd2 = 0
                pos1 -= 2 if ali=='map' else 3
                break
        read1 = {'crm': crm, 'pos': pos1, 'flag': flags[sd1], 'id': 'lala04.%012d' % (i)}
        read2 = {'crm': crm, 'pos': pos2, 'flag': flags[sd2], 'id': 'lala04.%012d' % (i)}
        out1.write(read.format(**read1))
        out2.write(read.format(**read2))

    # TOO CLOSE FROM RE

    # VAID PAIRS

    # DUPPLICATES
    i = 0
    while i < duplicates * 2:
        crm1  = genome.keys()    [int(random() * len(genome))]
        crm2  = genome.keys()    [int(random() * len(genome))]
        while True:
            frag1 = frags[crm1].keys()[int(random() * len(frags[crm1]))]
            frag2 = frags[crm2].keys()[int(random() * len(frags[crm2]))]
            if frags[crm2][frag2][1] - frags[crm2][frag2][0] < 6:
                continue
            if frags[crm1][frag1][1] - frags[crm1][frag1][0] < 6:
                continue
            if frag1 != frag2:
                break
        while True:
            pos1 = int(random() * (frags[crm1][frag1][1] - frags[crm1][frag1][0])
                       + frags[crm1][frag1][0])
            pos2 = int(random() * (frags[crm2][frag2][1] - frags[crm2][frag2][0])
                       + frags[crm2][frag2][0])
            if not (3 < pos1 < (len(genome[crm1]) - 3) and 3 < pos2 < (len(genome[crm2]) - 3)):
                continue
            break
        if pos2 > pos1:
            sd1 = 0
            sd2 = 1
            pos2 -= 2 if ali=='map' else 3
        elif pos2 < pos1:
            sd1 = 1
            sd2 = 0
            pos1 -= 2 if ali=='map' else 3
        read1 = {'crm': crm1, 'pos': pos1, 'flag': flags[sd1], 'id': 'lala05.%012d' % (i)}
        read2 = {'crm': crm2, 'pos': pos2, 'flag': flags[sd2], 'id': 'lala05.%012d' % (i)}
        out1.write(read.format(**read1))
        out2.write(read.format(**read2))
        i += 1
        if random() > 0.5:
            read2 = {'crm': crm1, 'pos': pos1, 'flag': flags[sd1], 'id': 'lala05.1%011d' % (i)}
            read1 = {'crm': crm2, 'pos': pos2, 'flag': flags[sd2], 'id': 'lala05.1%011d' % (i)}
            out1.write(read.format(**read1))
            out2.write(read.format(**read2))
        else:
            read2 = {'crm': crm1, 'pos': pos1, 'flag': flags[sd1], 'id': 'lala05.1%011d' % (i)}
            read1 = {'crm': crm2, 'pos': pos2, 'flag': flags[sd2], 'id': 'lala05.1%011d' % (i)}
            out1.write(read.format(**read1))
            out2.write(read.format(**read2))
        i += 1

    out1.close()
    out2.close()
    return genome


if __name__ == "__main__":
    if len(sys.argv) > 1:
        CHKTIME = bool(int(sys.argv.pop()))
        from time import time
    else:
        CHKTIME = False

    with catch_warnings():
        simplefilter("ignore")
        unittest.main()
    
