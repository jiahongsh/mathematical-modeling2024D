import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))


class X4PythonPortTest(unittest.TestCase):
    def test_x4_matches_matlab_result_file(self):
        import x4

        matlab = loadmat(ROOT / "matlab" / "x4_result.mat")
        d, I, peak, peak_d = x4.x4()
        np.testing.assert_allclose(d, matlab["d"].ravel(), rtol=0, atol=1e-12)
        np.testing.assert_allclose(I, matlab["I"].ravel(), rtol=0, atol=1e-12)
        self.assertAlmostEqual(peak, float(matlab["peak"].ravel()[0]), places=12)
        self.assertAlmostEqual(peak_d, float(matlab["peak_d"].ravel()[0]), places=12)


if __name__ == "__main__":
    unittest.main()
