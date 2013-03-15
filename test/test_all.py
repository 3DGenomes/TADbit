"""
30 Oct 2012

unittest for pytadbit functions
"""

import unittest
from pytadbit import tadbit, batch_tadbit, Chromosome, load_chromosome
from pytadbit.tad_clustering.tad_cmo import optimal_cmo
from os import system


class TestTadbit(unittest.TestCase):
    """
    test main tadbit functions
    """
    def test_tadbit(self):
        out = tadbit('20Kb/chrT/chrT_A.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        breaks = [0, 3, 9, 17, 22, 29, 36, 44, 50, 62, 67, 76, 90, 95]
        scores = [5.0, 10.0, 4.0, 8.0, 9.0, 5.0, 9.0, 9.0, 9.0, 10.0, 6.0, 5.0,
                  3.0, None]
        self.assertEqual(out['start'], breaks)
        self.assertEqual(out['score'], scores)


    def test_batch_tadbit(self):
        out = batch_tadbit('20Kb/chrT/', max_tad_size=20, verbose=False,
                           no_heuristic=True)
        breaks = [0, 4, 9, 15, 20, 33, 38, 44, 50, 62, 67, 76, 90, 95]
        scores = [4.0, 7.0, 3.0, 7.0, 3.0, 3.0, 10.0, 7.0, 10.0, 10.0, 9.0, 8.0,
                  5.0, None]
        self.assertEqual(out['start'], breaks)
        self.assertEqual(out['score'], scores)


    def test_tad_multi_aligner(self):

        exp1 = tadbit('40Kb/chrT/chrT_A.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        exp2 = tadbit('20Kb/chrT/chrT_B.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        exp3 = tadbit('20Kb/chrT/chrT_C.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        exp4 = tadbit('20Kb/chrT/chrT_D.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        test_chr = Chromosome(name='Test Chromosome',
                              tad_handlers=[exp1, exp2, exp3, exp4],
                              experiment_names=['exp1', 'exp2', 'exp3', 'exp4'],
                              experiment_resolutions=[40000,20000,20000,20000])

        test_chr.align_experiments(verbose=False, randomize=False)
        score1, pval1 = test_chr.align_experiments(verbose=False, randomize=True)
        _, pval2 = test_chr.align_experiments(verbose=False, randomize=True,
                                              rnd_method='shuffle')
        self.assertEqual(round(-26.095, 3), round(score1, 3))
        self.assertEqual(round(0.001, 1), round(pval1, 1))
        self.assertEqual(round(0.175, 1), round(pval2, 1))


    def test_chromosome_batch(self):
        test_chr = Chromosome(name='Test Chromosome',
                              experiment_resolutions=[20000]*4,
                         experiment_handlers=['20Kb/chrT/chrT_A.tsv',
                                              '20Kb/chrT/chrT_D.tsv',
                                              '20Kb/chrT/chrT_C.tsv'],
                         experiment_names=['exp1', 'exp2', 'exp3'])
        test_chr.find_tad(['exp1', 'exp2', 'exp3'], batch_mode=True,
                          verbose=False)
        self.assertEqual([3.0, 8.0, 16.0, 21.0, 28.0, 33.0, 38.0, 43.0, 49.0,
                          61.0, 66.0, 75.0, 89.0, 99.0],
                         test_chr.get_experiment('batch_exp1_exp2_exp3').brks)


    def test_save_load(self):
        test_chr = Chromosome(name='Test Chromosome',
                              experiment_resolutions=[20000]*2,
                              experiment_handlers=['20Kb/chrT/chrT_A.tsv',
                                                   '20Kb/chrT/chrT_C.tsv'],
                              experiment_names=['exp1', 'exp2'])
        test_chr.find_tad(['exp1', 'exp2'], verbose=False)
        test_chr.save_chromosome('lolo')
        test_chr = load_chromosome('lolo')
        system('rm -f lolo')
        system('rm -f lolo_hic')


    def test_tad_clustering(self):
        test_chr = Chromosome(name='Test Chromosome')
        test_chr.add_experiment('exp1', 20000, xp_handler='20Kb/chrT/chrT_D.tsv')
        test_chr.find_tad(['exp1'], verbose=False)
        all_tads = []
        for _, tad in test_chr.iter_tads('exp1'):
            all_tads.append(tad)
        align1, align2, _ = optimal_cmo(all_tads[4], all_tads[11], 10)
        self.assertEqual(align1,
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, '-'])
        self.assertEqual(align2,
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        
        

        
if __name__ == "__main__":
    unittest.main()
    
