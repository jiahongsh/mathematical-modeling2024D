import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = ROOT / "python"
sys.path.insert(0, str(PYTHON_DIR))


class X2PythonPortTest(unittest.TestCase):
    def test_x2_module_exists_and_returns_matlab_like_vectors(self):
        spec = importlib.util.find_spec("x2")
        self.assertIsNotNone(spec)

        import x2

        d, I, peak, peak_d = x2.x2()
        self.assertEqual(len(d), 28)  # 152.5 to 180.0 inclusive is 28 values
        self.assertEqual(len(I), len(d))
        self.assertAlmostEqual(float(d[0]), 152.5)
        self.assertAlmostEqual(float(d[-1]), 179.5)
        self.assertAlmostEqual(float(peak), float(max(I)))
        self.assertIn(float(peak_d), [float(v) for v in d])

    def test_x2_matches_matlab_result_file(self):
        import x2

        result_file = ROOT / "matlab" / "x2_result.mat"
        self.assertTrue(result_file.exists(), f"MATLAB result file not found: {result_file}")
        matlab = loadmat(result_file)

        d, I, peak, peak_d = x2.x2()
        np.testing.assert_allclose(d, matlab["d"].ravel(), rtol=0, atol=1e-12)
        np.testing.assert_allclose(I, matlab["I"].ravel(), rtol=0, atol=1e-12)
        self.assertAlmostEqual(peak, float(matlab["peak"].ravel()[0]), places=12)
        self.assertAlmostEqual(peak_d, float(matlab["peak_d"].ravel()[0]), places=12)


if __name__ == "__main__":
    unittest.main()
