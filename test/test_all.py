"""
30 Oct 2012

unittest for pytadbit functions
"""

import unittest
from pytadbit import tadbit, batch_tadbit


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

        
if __name__ == "__main__":
    unittest.main()
    
