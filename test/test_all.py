"""
30 Oct 2012

unittest for pytadbit functions
"""

import unittest
from pytadbit import tadbit, batch_tadbit, Chromosome
from pytadbit.tad_clustering.tad_cmo import optimal_cmo


class TestTadbit(unittest.TestCase):
    """
    test main tadbit funtion
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
                              experiments=[exp1, exp2, exp3, exp4],
                              experiment_names=['exp1', 'exp2', 'exp3', 'exp4'])
        score, pval = test_chr.align_experiments(verbose=False, randomize=True)
        self.assertEqual(round(19.555803, 3), round(score,3))
        self.assertEqual(round(0.4, 1), round(pval,1))


    def tad_clustering(self):
        test_chr = Chromosome(name='Test Chromosome', resolution=20000)
        test_chr.add_experiment('chrT/chrT_A.tsv', name='exp1')
        test_chr.find_TAD(['exp1'])
        all_tads = list(test_chr.iter_tads('exp1'))
        align1, align2 = optimal_cmo(all_tads[4], all_tads[8], 9)
        self.assertEqual(align1,
                         [1, 2, '-', '-', '-', '-', 3, '-', 4, 5, 6, '-', 7, 8])
        self.assertEqual(align2,
                         [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
        
        

        
if __name__ == "__main__":
    unittest.main()
    
