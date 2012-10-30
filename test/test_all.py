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
        out = tadbit('chrT/chrT_A.tsv', speed=4, verbose=False, heuristic=True)
        breaks = [7, 13, 18, 32, 37, 43, 49, 61, 66, 71, 76, 89]
        mlik = [float('nan'), 99.4382, float('nan'), float('nan'), 117.1356,
                 251.0043, float('nan'), float('nan'), 203.3079, 75.4968,
                 float('nan'), float('nan')]
        mlik = ['nan' if str(m)=='nan' else round(m, 3) for m in mlik]
        self.assertEqual(out[0], breaks)
        self.assertEqual(['nan' if str(m)=='nan' else round(m, 3) for m in out[1]], mlik)

    def test_batch_tadbit(self):
        out = batch_tadbit('chrT/', speed=4, verbose=False, heuristic=True)
        breaks = [4, 13, 18, 32, 38, 46, 52, 61, 66, 74, 80, 87, 93]
        mlik = [float('nan'), float('nan'), float('nan'), float('nan'),
                float('nan'), float('nan'), float('nan'), float('nan'),
                500.322110096, float('nan'), 331.156014405, 92.8811097334,
                555.467753751]
        mlik = ['nan' if str(m)=='nan' else round(m, 3) for m in mlik]
        self.assertEqual(out[0], breaks)
        self.assertEqual(['nan' if str(m)=='nan' else round(m, 3) for m in out[1]], mlik)

        
if __name__ == "__main__":
    unittest.main()
