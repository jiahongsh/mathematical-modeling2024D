import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = ROOT / "python"
sys.path.insert(0, str(PYTHON_DIR))


class X7PythonPortTest(unittest.TestCase):
    def test_x7_module_exists_and_returns_matlab_like_vectors(self):
        spec = importlib.util.find_spec("x7")
        self.assertIsNotNone(spec)

        import x7

        d, I, peak, peak_d = x7.x7()
        self.assertEqual(len(d), 13)  # 152.5 to 165.0 inclusive with 1.0 step is 13 values
        self.assertEqual(len(I), len(d))
        self.assertAlmostEqual(float(d[0]), 152.5)
        self.assertAlmostEqual(float(d[-1]), 164.5)
        self.assertAlmostEqual(float(peak), float(max(I)))
        self.assertIn(float(peak_d), [float(v) for v in d])

    def test_x7_matches_matlab_result_file(self):
        import x7

        result_file = ROOT / "matlab" / "x7_result.mat"
        self.assertTrue(result_file.exists(), f"MATLAB result file not found: {result_file}")
        matlab = loadmat(result_file)

        d, I, peak, peak_d = x7.x7()
        np.testing.assert_allclose(d, matlab["d"].ravel(), rtol=0, atol=1e-12)
        np.testing.assert_allclose(I, matlab["I"].ravel(), rtol=0, atol=1e-12)
        self.assertAlmostEqual(peak, float(matlab["peak"].ravel()[0]), places=12)
        self.assertAlmostEqual(peak_d, float(matlab["peak_d"].ravel()[0]), places=12)


if __name__ == "__main__":
    unittest.main()
