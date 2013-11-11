"""
30 Oct 2012

unittest for pytadbit functions
"""

import unittest
from pytadbit                        import Chromosome, load_chromosome
from pytadbit                        import tadbit, batch_tadbit
from pytadbit.tad_clustering.tad_cmo import optimal_cmo
from pytadbit.imp.structuralmodels   import load_structuralmodels
from pytadbit.imp.impmodel           import load_impmodel_from_cmm
from os                              import system, path, chdir
from warnings                        import warn
from distutils.spawn                 import find_executable

CHKTIME = False

if CHKTIME:
    from time import time

PATH = path.abspath(path.split(path.realpath(__file__))[0])


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
        
        if CHKTIME:
            t0 = time()

        
        global exp1, exp2, exp3, exp4
        exp1 = tadbit(PATH + '/40Kb/chrT/chrT_A.tsv', max_tad_size="auto",
                      verbose=False, no_heuristic=False, n_cpus='max')
        exp2 = tadbit(PATH + '/20Kb/chrT/chrT_B.tsv', max_tad_size="auto",
                      verbose=False, no_heuristic=False, n_cpus='max')
        exp3 = tadbit(PATH + '/20Kb/chrT/chrT_C.tsv', max_tad_size="auto",
                      verbose=False, no_heuristic=False, n_cpus='max')
        exp4 = tadbit(PATH + '/20Kb/chrT/chrT_D.tsv', max_tad_size="auto",
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
        if CHKTIME:
            t0 = time()

        global batch_exp
        batch_exp = batch_tadbit('20Kb/chrT/', max_tad_size=20, verbose=False,
                                 no_heuristic=True)
        # Breaks and scores with square root normalization.
        #breaks = [0, 4, 9, 15, 20, 29, 36, 44, 50, 62, 67, 76, 90, 95]
        #scores = [4.0, 7.0, 4.0, 8.0, 4.0, 4.0, 7.0, 7.0, 10.0, 10.0, 9.0, 8.0, 7.0, None]
        breaks = [0, 4, 14, 19, 34, 44, 50, 62, 67, 72, 90, 95]
        scores = [4.0, 6.0, 6.0, 6.0, 6.0, 6.0, 5.0, 6.0, 4.0, 6.0, 5.0, None]
        self.assertEqual(batch_exp['start'], breaks)
        self.assertEqual(batch_exp['score'], scores)
        if CHKTIME:
            print '2', time() - t0


    def test_03_tad_multi_aligner(self):

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
            exp.normalize_hic(silent=True)

        test_chr.align_experiments(verbose=False, randomize=False,
                                   method='global')
        score1, pval1 = test_chr.align_experiments(verbose=False,
                                                   method='global',
                                                   randomize=True, rnd_num=100)
        _, pval2 = test_chr.align_experiments(verbose=False, randomize=True,
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
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome',
                              experiment_resolutions=[20000]*3,
                              experiment_hic_data=[
                                  PATH + '/20Kb/chrT/chrT_A.tsv',
                                  PATH + '/20Kb/chrT/chrT_D.tsv',
                                  PATH + '/20Kb/chrT/chrT_C.tsv'],
                              experiment_names=['exp1', 'exp2', 'exp3'])
        test_chr.find_tad(['exp1', 'exp2', 'exp3'], batch_mode=True,
                          verbose=False)
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
        if CHKTIME:
            t0 = time()

        test_chr1 = Chromosome(name='Test Chromosome',
                              experiment_tads=[exp1, exp2],
                              experiment_names=['exp1', 'exp2'],
                              experiment_resolutions=[20000,20000])
        test_chr1.save_chromosome('lolo', force=True)
        test_chr2 = load_chromosome('lolo')
        system('rm -f lolo')
        system('rm -f lolo_hic')
        self.assertEqual(str(test_chr1.__dict__), str(test_chr2.__dict__))
        if CHKTIME:
            print '5', time() - t0


    def test_06_tad_clustering(self):
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome',
                              experiment_tads=[exp4],
                              experiment_names=['exp1'],
                              experiment_hic_data=[
                                  PATH + '/20Kb/chrT/chrT_D.tsv'],
                              experiment_resolutions=[20000,20000])
        all_tads = []
        for _, tad in test_chr.iter_tads('exp1'):
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
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000,
                              centromere_search=True,)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
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
                                hic_data=PATH + '/20Kb/chrT/chrT_C.tsv')
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
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
        exp = test_chr.experiments['exp1']
        sum20 = sum(exp.hic_data[0])
        exp.set_resolution(80000)
        sum80 = sum(exp.hic_data[0])
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(160000)
        sum160 = sum(exp.hic_data[0])
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(360000)
        sum360 = sum(exp.hic_data[0])
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(2400000)
        sum2400 = sum(exp.hic_data[0])
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(40000)
        sum40 = sum(exp.hic_data[0])
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(20000)
        sum21 = sum(exp.hic_data[0])
        check_hic(exp.hic_data[0], exp.size)
        exp.set_resolution(40000)
        sum41 = sum(exp.hic_data[0])
        check_hic(exp.hic_data[0], exp.size)
        self.assertTrue(sum20 == sum80 == sum160 == sum360 == sum40 \
                        == sum21 == sum2400 == sum41)
        if CHKTIME:
            print '8', time() - t0


    def test_09_hic_normalization(self):
        """
        writes interaction pair file.
        """
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
        exp = test_chr.experiments[0]
        exp.load_hic_data('20Kb/chrT/chrT_A.tsv', silent=True)
        exp.get_hic_zscores()
        exp.get_hic_zscores(zscored=False)
        sumz = sum([exp._zscores[k1][k2] for k1 in exp._zscores.keys()
                    for k2 in exp._zscores[k1]])
        self.assertEqual(round(sumz, 4), round(3983.2639, 4))
        if CHKTIME:
            print '9', time() - t0


    def test_10_generate_weights(self):
        """
        """
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
        exp = test_chr.experiments[0]
        tadbit_weights = exp.norm[:]
        exp.norm = None
        exp.normalize_hic()
        self.assertEqual([round(i, 3) for i in tadbit_weights[0][:100]],
                         [round(i, 3) for i in exp.norm[0][:100]])
        if CHKTIME:
            print '10', time() - t0


    def test_11_write_interaction_pairs(self):
        """
        writes interaction pair file.
        """
        if CHKTIME:
            t0 = time()

        test_chr = Chromosome(name='Test Chromosome', max_tad_size=260000)
        test_chr.add_experiment('exp1', 20000, tad_def=exp4,
                                hic_data=PATH + '/20Kb/chrT/chrT_D.tsv')
        exp = test_chr.experiments[0]
        exp.load_hic_data('20Kb/chrT/chrT_A.tsv', silent=True)
        exp.get_hic_zscores(zscored=False)
        exp.write_interaction_pairs('lala')
        lines = open('lala').readlines()
        self.assertEqual(len(lines), 4674)
        self.assertEqual(lines[25], '1\t28\t0.933380667098\n')
        self.assertEqual(lines[2000], '26\t70\t0.115108273108\n')
        system('rm -f lala')
        if CHKTIME:
            print '11', time() - t0


    def test_12_3d_modelling_optimization(self):
        """
        quick test to generate 3D coordinates from 3? simple models???
        """
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
        exp.load_hic_data('20Kb/chrT/chrT_A.tsv', silent=True)
        exp.normalize_hic(silent=True)
        result = exp.optimal_imp_parameters(50,70, n_cpus=2,
                                            n_models=10, n_keep=2,
                                            lowfreq_range=(-0.6, 0, 0.5),
                                            upfreq_range=(0, 1.1, 1),
                                            maxdist_range=(500, 1000, 500),
                                            verbose=False)

        # get best correlations
        config = result.get_best_parameters_dict()
        wanted = {'maxdist': 500.0, 'upfreq': 0.0, 'kforce': 5,
                  'reference': '', 'lowfreq': -0.1, 'scale': 0.005}
        self.assertEqual([round(i, 4) for i in config.values()if not type(i) is str],
                         [round(i, 4) for i in wanted.values()if not type(i) is str])
        if CHKTIME:
            print '12', time() - t0


    def test_13_3d_modelling_centroid(self):
        """
        quick test to generate 3D coordinates from 3? simple models???
        """
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
        exp.load_hic_data('20Kb/chrT/chrT_A.tsv', silent=True)
        exp.normalize_hic(silent=True)
        models = exp.model_region(50, 70, n_models=110, n_keep=25,
                                  n_cpus=4,
                                  config={'kforce': 5, 'maxdist': 500,
                                          'scale': 0.005,
                                          'upfreq': 1.0, 'lowfreq': -0.6})
        models.save_models('models.pick')
        
        avg = models.average_model()
        xis = [691.894, 576.03, 461.998, 350.654, 228.832, 249.183, 184.247,
               66.02, 131.339, 105.844, 65.856, -46.289, -167.628, -254.734,
               -196.974, -305.433, -391.167, -502.536, -438.747, -424.81,
               -383.576]
        yis = [-491.148, -424.611, -344.702, -256.15, -180.741, -354.946,
               -334.481, -272.075, -398.058, -221.028, -68.147, 57.421,
               31.059, 111.967, 239.04, 368.515, 408.923, 527.132, 470.932,
               558.696, 572.405]
        zis = [-751.145, -618.017, -489.136, -364.219, -231.151, -243.215,
               -149.716, -15.725, -65.75, -64.688, -57.14, 49.595, 185.598,
               291.651, 191.279, 292.817, 426.3, 520.773, 419.838, 367.156,
               304.889]
        self.assertEqual([round(x, 3) for x in avg['x']],
                         [round(x, 3) for x in xis])
        self.assertEqual([round(x, 3) for x in avg['y']],
                         [round(x, 3) for x in yis])
        self.assertEqual([round(x, 3) for x in avg['z']],
                         [round(x, 3) for x in zis])

        centroid = models.centroid_model()
        self.assertEqual(centroid['rand_init'], 97)
        if CHKTIME:
            print '13', time() - t0


    def test_14_3d_clustering(self):
        """
        """
        if CHKTIME:
            t0 = time()

        models = load_structuralmodels('models.pick')
        wnt = {1: [103, 69, 80, 40, 32, 30, 95, 101, 72, 77, 52],
               2: [78, 57, 97, 34, 96, 14, 88, 71, 56],
               3: [106, 27, 13, 41, 43]}
        if find_executable('mcl'):
            models.cluster_models(method='mcl', verbose=False)
            self.assertEqual(models.clusters, wnt)
        wnt = {1: [78, 57, 97, 34, 106, 27, 96, 13, 41, 14, 88, 43, 77],
               2: [103, 69, 80, 40, 32, 30, 95, 101, 72, 71, 52, 56]}
        models.cluster_models(method='ward', verbose=False)
        self.assertEqual(models.clusters, wnt)
        if CHKTIME:
            print '14', time() - t0


    def test_15_3d_modelling(self):
        """
        """
        if CHKTIME:
            t0 = time()

        models = load_structuralmodels('models.pick') 
        models.cluster_models(method='ward', verbose=False)
        # density
        models.density_plot(savedata='lala', plot=False)
        lines = open('lala').readlines()
        self.assertEqual(len(lines), 22)
        self.assertEqual(lines[1], '1\t99.979\t100.009\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\n')
        self.assertEqual(lines[15], '15\t99.937\t99.972\t99.966\t99.988\t99.965\t99.984\t99.972\t99.987\t99.972\t99.988\n')
        # contacts
        cmap = models.get_contact_matrix(cutoff=300)
        self.assertEqual(round(
            sum([i if i >=0 else 0 for i in reduce(lambda x, y: x+y, cmap)]),
            3), 81.44)
        # define best models
        models.define_best_models(10)
        self.assertEqual(len(models), 10)
        m1 = models[9]
        models.define_best_models(25)
        self.assertEqual(len(models), 25)
        self.assertEqual(m1, models[9])
        # correlation
        corr, pval = models.correlate_with_real_data(cutoff=300)
        self.assertEqual(round(corr, 4), 0.6482)
        self.assertEqual(round(pval, 4), round(0, 4))
        # consistency
        models.model_consistency(plot=False, savedata='lala')
        lines = open('lala').readlines()
        self.assertEqual(len(lines), 22)
        self.assertEqual(lines[1], '1\t15.0\t28.333\t42.333\t53.0\n')
        self.assertEqual(lines[15], '15\t83.667\t100.0\t100.0\t100.0\n')
        # measure angle
        self.assertEqual(round(models.angle_between_3_particles(2,8,15), 3),
                         129.571)
        self.assertEqual(round(models.angle_between_3_particles(19,20,21), 3),
                         60.176)
        self.assertEqual(round(models.angle_between_3_particles(15,14,11), 3),
                         65.447)
        # coordinates
        self.assertEqual([round(x, 3) for x in models.particle_coordinates(15)],
                         [6545.505, -5979.032, 1252.336])
        # dihedral_angle
        self.assertEqual(round(models.dihedral_angle(2,8,15, 16), 3), -59.15)
        self.assertEqual(round(models.dihedral_angle(15,19,20,21), 3), -30.979)
        self.assertEqual(round(models.dihedral_angle(15,14,11, 12), 3), -1.224)
        # median distance
        self.assertEqual(round(models.median_3d_dist(3, 20, plot=False), 3),
                         1563.867)
        self.assertEqual(round(models.median_3d_dist(3, 20, cluster=1,
                                                     plot=False), 3), 1574.936)
        self.assertEqual(round(models.median_3d_dist(7, 10, models=range(5),
                                                     plot=False), 3), 266.02)
        # write cmm
        models.write_cmm('.', model_num=2)
        models.write_cmm('.', models=range(5))
        models.write_cmm('.', cluster=2)
        # write xyz
        models.write_xyz('.', model_num=2)
        models.write_xyz('.', models=range(5))
        models.write_xyz('.', cluster=2)
        # clean
        system('rm -f model.*')
        system('rm -f lala')
        if CHKTIME:
            print '15', time() - t0


    def test_16_models_stats(self):
        if CHKTIME:
            t0 = time()

        models = load_structuralmodels('models.pick') 
        # write cmm
        models.write_cmm('.', model_num=2)
        model = load_impmodel_from_cmm('model.78.cmm')
        # clean
        system('rm -f model.*')
        # stats
        self.assertEqual(199.968, round(model.distance(2, 3), 3))
        self.assertEqual(1053.752, round(model.distance(8, 20), 3))
        self.assertEqual(622.85,
                         round(model.radius_of_gyration(), 3))
        self.assertEqual(4000.523, round(model.contour(), 3))
        self.assertEqual(2273.677,
                         round(model.shortest_axe()+model.longest_axe(), 3))
        self.assertEqual([15, 16], model.inaccessible_particles(1000))
        acc_inacc = [(87, 23), (24, 49), (22, 50), (21, 48), (17, 70), (27, 62),
                     (27, 46), (29, 60), (41, 49), (27, 44), (15, 56), (28, 41),
                     (24, 46), (34, 52), (26, 56), (6, 64), (23, 48), (53, 36),
                     (25, 60), (40, 48)]
        self.assertEqual([round(i, 3) for i in 830.0, 1947.0, 6.258, 8.672] + [acc_inacc],
                         [i if type(i) is list else round(i, 3)
                          for i in model.accessible_surface(300, nump=150)])
        if CHKTIME:
            print '16', time() - t0


    def test_17_tadbit_c(self):
        """
        Runs tests written in c, around the detection of TADs
        """
        if CHKTIME:
            t0 = time()

        print '\n\nC SIDE'
        print '------'
        chdir(PATH + '/../src/test/')
        system('make clean')
        system('make')
        return_code = system('make test')
        chdir(PATH)
        self.assertEqual(return_code, 0)
        if CHKTIME:
            print '17', time() - t0


if __name__ == "__main__":
    unittest.main()
    
