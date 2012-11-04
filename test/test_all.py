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
        out = tadbit('chrT/chrT_A.tsv', max_tad_size=4, verbose=False, heuristic=True)
        breaks = [7, 13, 18, 33, 38, 43, 49, 61, 66, 71, 89]
        scores = [3.0, 6.0, 7.0, 4.0, 4.0, 8.0, 9.0, 6.0, 8.0, 7.0, 9.0]
        self.assertEqual(out[0], breaks)
        self.assertEqual(out[1], scores)

    def test_batch_tadbit(self):
        out = batch_tadbit('chrT/', max_tad_size=4, verbose=False, heuristic=True)
        breaks = [4, 13, 18, 33, 38, 46, 52, 61, 66, 71, 89, 94]
        scores = [3.0, 6.0, 8.0, 8.0, 8.0, 5.0, 4.0, 5.0, 5.0, 7.0, 7.0, 7.0]
        self.assertEqual(out[0], breaks)
        self.assertEqual(out[1], scores)

        
if __name__ == "__main__":
    unittest.main()
