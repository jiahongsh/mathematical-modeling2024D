# D033 Mathematical Modeling Code

This repository contains Python ports of the appendix code for the D033 mathematical modeling problem, together with the original MATLAB reference files and lightweight teaching notes.

## Contents

- `python/`: runnable Python implementations and appendix entry scripts.
- `matlab/`: original/reference MATLAB scripts and `.mat` result files used for comparison.
- `tests/`: regression tests comparing Python output with MATLAB reference results.
- `video_scripts/`: two combined Markdown handouts:
  - `classroom_display_all.md`: concise classroom display version.
  - `read_aloud_all.md`: teacher read-aloud version.

The raw OCR materials and temporary working files are intentionally not included in this clean release.

## Setup

Create and activate a Python environment, then install dependencies:

```powershell
python -m pip install -r requirements.txt
```

Required packages:

- `numpy`
- `scipy`
- `pytest`

## Run Examples

Run the first problem:

```powershell
python python/appendix_1_problem1_probability_p00.py
```

Run problem 2 fine search:

```powershell
python python/appendix_2_5_problem2_fine_depth_156_to_160.py
```

Run problem 3 nine-charge array:

```powershell
python python/appendix_3_problem3_nine_charges.py
```

Some scripts perform three-dimensional numerical integration and may take a long time to finish.

## Tests

Quick syntax check:

```powershell
python -m compileall -q python tests
```

Run all tests:

```powershell
python -m pytest tests
```

The full tests may be slow because they execute the numerical integration code.

## File Index

The main Python entry points are:

- `python/appendix_1_problem1_probability_p00.py`
- `python/appendix_2_1_problem2_depth_152_5_to_180.py`
- `python/appendix_2_2_problem2_depth_140_to_152_5.py`
- `python/appendix_2_3_problem2_depth_100_to_140.py`
- `python/appendix_2_4_problem2_depth_87_5_to_100.py`
- `python/appendix_2_5_problem2_fine_depth_156_to_160.py`
- `python/appendix_3_problem3_nine_charges.py`

See `python/CODE_INDEX.md` for the detailed code index.

