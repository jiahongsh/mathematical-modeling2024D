import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = ROOT / "python"
sys.path.insert(0, str(PYTHON_DIR))


class X1PythonPortTest(unittest.TestCase):
    def test_x1_module_exists_and_returns_correct_probability(self):
        spec = importlib.util.find_spec("appendix_1_problem1_probability_p00")
        self.assertIsNotNone(spec)

        import appendix_1_problem1_probability_p00

        poe = appendix_1_problem1_probability_p00.appendix_1_problem1_probability_p00()
        self.assertTrue(0.0 < poe < 1.0)

    def test_x1_matches_matlab_result_file(self):
        import appendix_1_problem1_probability_p00

        result_file = ROOT / "matlab" / "appendix_1_result.mat"
        self.assertTrue(result_file.exists(), f"MATLAB result file not found: {result_file}")
        matlab = loadmat(result_file)

        poe = appendix_1_problem1_probability_p00.appendix_1_problem1_probability_p00()
        matlab_poe = float(matlab["poe"].ravel()[0])

        np.testing.assert_allclose(poe, matlab_poe, rtol=0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
