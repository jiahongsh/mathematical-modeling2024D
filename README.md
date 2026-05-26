#2024D数学建模附录代码整理

本仓库整理了 D033 数学建模题目附录中的 MATLAB 代码，并提供对应的 Python 移植版本、MATLAB 参考结果、回归测试和课堂讲义。

项目目标是让附录代码更容易运行、验证和讲解。

## 目录结构

```text
.
├── python/          # Python 移植代码和各附录入口脚本
├── matlab/          # MATLAB 原始/参考代码，以及 .mat 结果文件
├── tests/           # Python 与 MATLAB 结果对照的测试
├── video_scripts/   # 课堂讲义合集
├── requirements.txt # Python 依赖
└── README.md        # 项目说明
```

## 主要内容

`python/` 中的主要入口文件：

- `appendix_1_problem1_probability_p00.py`：问题一，单枚深弹命中概率。
- `appendix_2_1_problem2_depth_152_5_to_180.py`：问题二，定深 152.5 到 180 米。
- `appendix_2_2_problem2_depth_140_to_152_5.py`：问题二，定深 140 到 152.5 米。
- `appendix_2_3_problem2_depth_100_to_140.py`：问题二，定深 100 到 140 米。
- `appendix_2_4_problem2_depth_87_5_to_100.py`：问题二，定深 87.5 到 100 米。
- `appendix_2_5_problem2_fine_depth_156_to_160.py`：问题二，156 到 160 米精细搜索。
- `appendix_3_problem3_nine_charges.py`：问题三，九枚深弹阵列。

`python/x2.py` 到 `python/x7.py` 是从 MATLAB 代码移植过来的核心计算文件。

`matlab/` 中保留 MATLAB 版本和 `.mat` 结果文件，用于和 Python 结果进行对照。

`video_scripts/` 中保留两份讲义：

- `classroom_display_all.md`：课堂展示版，适合录屏展示。

## 安装依赖

建议先创建虚拟环境，然后安装依赖：

```powershell
python -m pip install -r requirements.txt
```

依赖包括：

- `numpy`
- `scipy`
- `pytest`

如果运行时报 `ModuleNotFoundError: No module named 'scipy'`，说明当前 Python 环境没有安装依赖，需要先执行上面的安装命令。

## 运行示例

运行问题一：

```powershell
python python/appendix_1_problem1_probability_p00.py
```

运行问题二精细搜索：

```powershell
python python/appendix_2_5_problem2_fine_depth_156_to_160.py
```

运行问题三九枚深弹阵列：

```powershell
python python/appendix_3_problem3_nine_charges.py
```

## 测试

只检查 Python 文件语法：

```powershell
python -m compileall -q python tests
```

运行全部测试：

```powershell
python -m pytest tests
```

测试会调用部分数值积分程序，运行时间可能较长。

## 运行时间说明

部分脚本包含三维数值积分，运行时间会比较长，尤其是：

- `appendix_2_5_problem2_fine_depth_156_to_160.py`
- `appendix_3_problem3_nine_charges.py`
- 对应的 `x6.py`、`x7.py`


## 说明

本仓库中的 Python 代码主要用于复现 MATLAB 附录代码的计算逻辑。为了便于运行，部分边界处的平方根计算加入了安全处理，避免由于浮点误差导致开方负数。

原始 OCR 文本、临时缓存和压缩包没有放入本开源整理版。
