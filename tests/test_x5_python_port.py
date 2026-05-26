import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))


class X5PythonPortTest(unittest.TestCase):
    def test_x5_matches_matlab_result_file(self):
        import x5

        matlab = loadmat(ROOT / "matlab" / "x5_result.mat")
        d, I, peak, peak_d = x5.x5()
        np.testing.assert_allclose(d, matlab["d"].ravel(), rtol=0, atol=1e-12)
        np.testing.assert_allclose(I, matlab["I"].ravel(), rtol=0, atol=1e-12)
        self.assertAlmostEqual(peak, float(matlab["peak"].ravel()[0]), places=12)
        self.assertAlmostEqual(peak_d, float(matlab["peak_d"].ravel()[0]), places=12)


if __name__ == "__main__":
    unittest.main()
