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
        out = tadbit('chrT/chrT_A.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        breaks = [0, 3, 9, 14, 20, 30, 38, 44, 50, 67, 72, 77, 82, 89, 94]
        scores = [10.0, 10.0, 8.0, 10.0, 10.0, 6.0, 8.0, 5.0,
                  10.0, 10.0, 10.0, 10.0, 10.0, 10.0, None]
        self.assertEqual(out['start'], breaks)
        self.assertEqual(out['score'], scores)


    def test_batch_tadbit(self):
        out = batch_tadbit('chrT/', max_tad_size=20, verbose=False,
                           no_heuristic=True)
        breaks = [0, 3, 9, 14, 20, 36, 41, 46, 51, 56, 62, 67, 74, 79, 84]
        scores = [8.0,  10.0,  7.0,  7.0,  5.0,  5.0,  7.0,  7.0,  9.0,  10.0,
                  10.0,  5.0,  3.0,  6.0,  None]
        self.assertEqual(out['start'], breaks)
        self.assertEqual(out['score'], scores)


    def test_tad_multi_aligner(self):
        exp1 = tadbit('chrT/chrT_A.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        exp2 = tadbit('chrT/chrT_B.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        exp3 = tadbit('chrT/chrT_C.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        exp4 = tadbit('chrT/chrT_D.tsv', max_tad_size="auto",
                     verbose=False, no_heuristic=False)
        test_chr = Chromosome(name='Test Chromosome', resolution=20000,
                         tad_handlers=[exp1, exp2, exp3, exp4],
                         experiment_names=['exp1', 'exp2', 'exp3', 'exp4'])
        test_chr.align_experiments(verbose=False, randomize=False)
        score1, pval1 = test_chr.align_experiments(verbose=False, randomize=True)
        _, pval2 = test_chr.align_experiments(verbose=False, randomize=True,
                                                 rnd_method='shuffle')
        self.assertEqual(round(9.569707498, 3), round(score1, 3))
        self.assertEqual(round(0.066, 1), round(pval1, 1))
        self.assertEqual(round(0.175, 1), round(pval2, 1))


    def test_chromosome_batch(self):
        test_chr = Chromosome(name='Test Chromosome', resolution=20000,
                         experiment_handlers=['chrT/chrT_A.tsv',
                                              'chrT/chrT_D.tsv',
                                              'chrT/chrT_C.tsv'],
                         experiment_names=['exp1', 'exp2', 'exp3'])
        test_chr.find_tad(['exp1', 'exp2', 'exp3'], batch_mode=True, verbose=False)
        self.assertEqual([2.0, 8.0, 19.0, 35.0, 40.0, 45.0, 50.0, 55.0, 61.0,
                          66.0, 73.0, 78.0, 83.0, 88.0, 93.0, 99.0],
                         test_chr.experiments['batch_exp3_exp2_exp1'].brks)


    def test_save_load(self):
        test_chr = Chromosome(name='Test Chromosome', resolution=20000,
                         experiment_handlers=['chrT/chrT_A.tsv',
                                              'chrT/chrT_C.tsv'],
                         experiment_names=['exp1', 'exp2'])
        test_chr.find_tad(['exp1', 'exp2'], verbose=False)
        test_chr.save_chromosome('lolo')
        test_chr = load_chromosome('lolo')
        system('rm -f lolo')


    def test_tad_clustering(self):
        test_chr = Chromosome(name='Test Chromosome')
        test_chr.add_experiment('exp1', 20000, xp_handler='chrT/chrT_D.tsv')
        test_chr.find_tad(['exp1'], verbose=False)
        all_tads = []
        for _, tad in test_chr.iter_tads('exp1'):
            all_tads.append(tad)
        align1, align2, _ = optimal_cmo(all_tads[9], all_tads[10], 11)
        self.assertEqual(align1,
                         [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        self.assertEqual(align2,
                         [0, 1, '-', 2, '-', 3, 4, 5, 6, 7, 8, 9, 10])
        
        

        
if __name__ == "__main__":
    unittest.main()
    
