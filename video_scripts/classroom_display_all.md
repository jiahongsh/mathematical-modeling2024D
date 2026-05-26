## 问题一代码讲解，单枚深弹命中概率

## 本节展示

源码文件：`python/appendix_1_problem1_probability_p00.py`

### 1. 参数和二维正态密度

```python
def appendix_1_problem1_probability_p00() -> float:
    sigma = 120.0  # 投弹落点偏差的标准差（误差散布半径指标）
    L = 100.0      # 潜艇本体的有效长度
    R = 20.0       # 深弹的破坏/杀伤半径
    W = 20.0       # 潜艇本体的有效宽度

    # dblquad 要求被积函数参数顺序为 (y, x)。
    def f(y: float, x: float) -> float:
        return math.exp(-(x**2 + y**2) / (2.0 * sigma**2))
```

**关键要点**

- 这里先定义第一问的主函数。
- `sigma` 是水平投弹误差的标准差，`L` 是潜艇长度，`R` 是深弹杀伤半径，`W` 是潜艇宽度。

### 2. 左侧半圆区域

```python
def ymin1(x: float) -> float:
    val = R**2 - (x + L/2.0) ** 2
    return -W / 2.0 - math.sqrt(max(0.0, val))

def ymax1(x: float) -> float:
    val = R**2 - (x + L/2.0) ** 2
    return W / 2.0 + math.sqrt(max(0.0, val))

I1, _ = dblquad(f, -R - L / 2.0, -L / 2.0, ymin1, ymax1)
```

**关键要点**

- 第一块区域是潜艇左端向外扩展出来的半圆区域。
- `ymin1` 和 `ymax1` 分别给出这个区域的下边界和上边界。

### 3. 中间矩形区域和右侧半圆区域

```python
I2, _ = dblquad(
    f,
    -L / 2.0,
    L / 2.0,
    lambda _: -W / 2.0 - R,
    lambda _: W / 2.0 + R,
)

def ymin3(x: float) -> float:
    val = R**2 - (x - L/2.0) ** 2
    return -W / 2.0 - math.sqrt(max(0.0, val))

def ymax3(x: float) -> float:
    val = R**2 - (x - L/2.0) ** 2
    return W / 2.0 + math.sqrt(max(0.0, val))

I3, _ = dblquad(f, L / 2.0, R + L / 2.0, ymin3, ymax3)
```

**关键要点**

- `I2` 是中间矩形区域。
- x 从潜艇左端到右端，y 从潜艇宽度下边界再向下扩展一个杀伤半径，到上边界再向上扩展一个杀伤半径。

### 4. 合成概率并输出

```python
poe = (1.0 / (2.0 * math.pi * sigma**2)) * (I1 + I2 + I3)
return poe

if __name__ == "__main__":
    poe = appendix_1_problem1_probability_p00()
    print(f"第一问单枚深弹命中概率 P(0,0) = {poe:.6f} (即约 {poe*100:.4f}%)")
```

**关键要点**

- 最后把三个积分结果相加，再乘上二维正态分布的归一化系数，就得到单枚深弹命中概率。
- 这里可以强调：这一段代码的核心不是求面积，而是在命中区域上累加概率密度。

---

## 问题二代码框架，三维概率模型

## 本节展示

源码文件：`python/x2.py`

### 1. 导入库和辅助函数

```python
import math
from typing import Callable

import numpy as np
from scipy.integrate import quad, tplquad
from scipy.stats import norm

def _matlab_range(start: float, step: float, stop: float) -> np.ndarray:
    if not (math.isfinite(start) and math.isfinite(stop)):
        return np.array([], dtype=float)
    eps = abs(step) * 1e-9
    if start > stop + eps:
        return np.array([], dtype=float)
    num = int(math.floor((stop - start + eps) / step)) + 1
    return np.linspace(start, start + (num - 1) * step, num)

def _safe_sqrt(value: float) -> float:
    if value < 0:
        return 0.0
    return math.sqrt(value)
```

**关键要点**

- 这里先导入数值计算需要的库。
- `numpy` 用来处理数组，`quad` 和 `tplquad` 用来做积分，`norm` 用来计算正态分布函数。

### 2. 三维黎曼积分工具

```python
def _riemann_3d_vectorized(
    total: float,
    zmin: float,
    zmax: float,
    dz: float,
    xmin: Callable[[float], float],
    xmax: Callable[[float], float],
    dx: float,
    ymin: Callable[[float, float], float],
    ymax: Callable[[float, float], float],
    dy: float,
    fun: Callable[[float, np.ndarray, float], np.ndarray],
) -> float:
    z_arr = _matlab_range(zmin, dz, zmax)
    if len(z_arr) == 0:
        return total
    for z in z_arr:
        x_arr = _matlab_range(xmin(z), dx, xmax(z))
        for x in x_arr:
            y_arr = _matlab_range(ymin(x, z), dy, ymax(x, z))
            if len(y_arr) > 0:
                total += np.sum(fun(x, y_arr, z)) * dx * dy * dz
    return total
```

**关键要点**

- 这段是三维数值积分工具。
- 它按 z、x、y 三层循环，把三维空间切成很多小网格。

### 3. 问题二参数和概率密度

```python
def x2() -> tuple[np.ndarray, np.ndarray, float, float]:
    sigma = 120.0  # 水平投弹误差标准差
    L = 100.0      # 潜艇本体长度
    R = 20.0       # 深弹杀伤半径
    W = 20.0       # 潜艇本体宽度
    H = 25.0       # 潜艇本体高度
    sigma_z = 40.0 # 潜艇深度分布的标准差
    l1 = 120.0     # 截断深度（潜艇只在 120 米以下活动）
    h0 = 150.0     # 潜艇深度分布的均值

    def f(x: float, y: np.ndarray | float) -> np.ndarray | float:
        return (1.0 / (2.0 * math.pi * sigma**2)) * np.exp(
            -((x**2 + y**2) / (2.0 * sigma**2))
        )

    Phi = lambda x: norm.cdf(x, 0.0, 1.0)
    dm = 1.0 / (1.0 - Phi((l1 - h0) / sigma_z))

    def g_z(z: float) -> float:
        return (
            (1.0 / sigma_z)
            * dm
            * (1.0 / math.sqrt(2.0 * math.pi))
            * math.exp(-((z - h0) ** 2) / (2.0 * sigma_z**2))
        )

    def fun(x: float, y: np.ndarray | float, z: float) -> np.ndarray | float:
        return f(x, y) * g_z(z)
```

**关键要点**

- 这里进入问题二主函数。
- 相比第一问，多了潜艇高度 `H`、深度标准差 `sigma_z`、截断深度 `l1` 和平均深度 `h0`。

### 4. 搜索深度

```python
d = _matlab_range(152.5, 1.0, 180.0)
```

**关键要点**

- 这行代码表示本节对应的深度区间是 152.5 米到 180 米，步长是 1 米。
- 这里可以这样总结：`x2.py` 的前半部分负责把问题二从二维概率升级成三维概率，后半部分才开始按不同空间区域计算命中概率。

---

## 问题二代码讲解，深水区间 152.5 到 180 米

## 本节展示

源码文件：`python/x2.py`

### 1. 直接覆盖项 I1

```python
prob_x = norm.cdf(L / (2.0 * sigma)) - norm.cdf(-L / (2.0 * sigma))
prob_y = norm.cdf(W / (2.0 * sigma)) - norm.cdf(-W / (2.0 * sigma))
area_prob = prob_x * prob_y
I1 = []
for di in d:
    val = area_prob * quad(g_z, l1, di - R - H / 2.0)[0]
    I1.append(val)
I1 = np.array(I1)
```

**关键要点**

- `I1` 是直接覆盖项。
- 水平面上先算潜艇主体矩形范围内的概率，x 方向和 y 方向可以分开用正态分布 CDF 计算。

### 2. 中间覆盖项 I5

```python
I5 = np.array([quad(g_z, di - H / 2.0, di + H / 2.0)[0] for di in d])
I5 = 0.083734 * I5
```

**关键要点**

- `I5` 用的是第一问得到的水平覆盖概率常数 `0.083734`。
- 这表示水平面上落入扩展命中区域的概率。

### 3. 过渡区域的三维积分

```python
for di in d:
    dx = 0.5
    dy = 0.5
    dz = 0.5

    total = 0.0
    zmin = di - R - 0.5 * H
    zmax = di - 0.5 * H

    def xmin(z: float) -> float:
        return -L / 2.0 - _safe_sqrt(R**2 - (di - z - H / 2.0) ** 2)

    def xmax(z: float) -> float:
        return -L / 2.0

    def ymin(x: float, z: float) -> float:
        return -W / 2.0 - _safe_sqrt(
            R**2 - (di - z - H / 2.0) ** 2 - (x + L / 2.0) ** 2
        )

    def ymax(x: float, z: float) -> float:
        return W / 2.0 + _safe_sqrt(
            R**2 - (di - z - H / 2.0) ** 2 - (x + L / 2.0) ** 2
        )

    total = _riemann_3d_vectorized(
        total, zmin, zmax, dz, xmin, xmax, dx, ymin, ymax, dy, fun
    )
    I2.append(total)
```

**关键要点**

- 这段展示的是 `I2` 的计算方式，后面的 `I3`、`I4`、`I6`、`I7`、`I8` 也是类似思路。
- 先固定当前投放深度 `di`，再设置网格步长 `dx`、`dy`、`dz`。

### 4. 合成总概率并找最优深度

```python
I = (
    I1
    + np.array(I2)
    + np.array(I3)
    + np.array(I4)
    + np.array(I5)
    + np.array(I6)
    + np.array(I7)
    + np.array(I8)
)
idx = int(np.argmax(I))
peak = float(I[idx])
peak_d = float(d[idx])
return d, I, peak, peak_d
```

**关键要点**

- 最后把八个区域的概率全部相加，得到每一个投放深度对应的总命中概率。
- `np.argmax(I)` 找到最大概率的位置，`peak` 是最大命中概率，`peak_d` 是对应的最优定深。

---

## 问题二代码讲解，深度 140 到 152.5 米

## 本节展示

源码文件：`python/x3.py`

### 1. 深度区间和保留的积分项

```python
d = _matlab_range(140.0, 1.0, 152.5)

I2 = []
I3 = []
I4 = []
I5 = np.array([quad(g_z, di - H / 2.0, di + H / 2.0)[0] for di in d])
I5 = 0.083734 * I5
I6 = []
I7 = []
I8 = []
```

**关键要点**

- 这里的搜索深度从 140 米到 152.5 米。
- 和上一节的 `x2.py` 相比，这段没有单独计算 `I1`，而是重点计算 `I2` 到 `I8`。

### 2. 深度下界改成 l1

```python
for di in d:
    dx = 0.5
    dy = 0.5
    dz = 0.5

    total = 0.0
    zmin = l1
    zmax = di - 0.5 * H
```

**关键要点**

- 这是 `x3.py` 最需要讲清楚的地方。
- 在 `x2.py` 的深水区间里，某些区域的 z 下界由爆炸半径决定。

### 3. 边缘区域积分

```python
def xmin(z: float) -> float:
    return -L / 2.0 - _safe_sqrt(R**2 - (di - z - H / 2.0) ** 2)

def xmax(z: float) -> float:
    return -L / 2.0

def ymin(x: float, z: float) -> float:
    return -W / 2.0 - _safe_sqrt(
        R**2 - (di - z - H / 2.0) ** 2 - (x + L / 2.0) ** 2
    )

def ymax(x: float, z: float) -> float:
    return W / 2.0 + _safe_sqrt(
        R**2 - (di - z - H / 2.0) ** 2 - (x + L / 2.0) ** 2
    )

total = _riemann_3d(total, zmin, zmax, dz, xmin, xmax, dx, ymin, ymax, dy, fun)
I2.append(total)
```

**关键要点**

- 这一段和上一讲的结构相似，仍然是用球面边界确定 x 和 y 的范围。
- 区别在于，这里调用的是 `_riemann_3d`，也就是普通三层循环版本；上一节 `x2.py` 用的是向量化版本。

### 4. 合成概率

```python
I = (
    np.array(I2)
    + np.array(I3)
    + np.array(I4)
    + np.array(I5)
    + np.array(I6)
    + np.array(I7)
    + np.array(I8)
)
idx = int(np.argmax(I))
peak = float(I[idx])
peak_d = float(d[idx])
return d, I, peak, peak_d
```

**关键要点**

- 最后仍然是把有效积分项加起来，找到最大命中概率和对应深度。
- 这一课的代码重点是：同一个三维模型在不同深度区间里，积分项和积分边界会发生变化。

---

## 问题二代码讲解，深度 100 到 140 米

## 本节展示

源码文件：`python/x4.py`

### 1. 导入复用工具

```python
from x3 import _matlab_range, _riemann_3d, _safe_sqrt
```

**关键要点**

- `x4.py` 没有重新写 `_matlab_range`、`_riemann_3d` 和 `_safe_sqrt`，而是直接从 `x3.py` 复用。
- 这说明从这一段开始，模型结构已经固定，主要变化是深度区间和保留的积分区域。

### 2. 深度区间和 I5

```python
d = _matlab_range(100.0, 1.0, 140.0)
I5 = np.array([quad(g_z, l1, di + H / 2.0)[0] for di in d])
I5 = 0.083734 * I5
I6 = []
I7 = []
I8 = []
```

**关键要点**

- 这里搜索深度从 100 米到 140 米。
- `I5` 的深度积分从 `l1` 到 `di + H/2`。

### 3. I6 的三维积分

```python
for di in d:
    dx = 0.5
    dy = 0.5
    dz = 0.5

    total = 0.0
    zmin = l1
    zmax = di + R + 0.5 * H

    def xmin6(z: float) -> float:
        return -L / 2.0 - _safe_sqrt(R**2 - (di - z + H / 2.0) ** 2)

    def xmax6(z: float) -> float:
        return -L / 2.0

    def ymin6(x: float, z: float) -> float:
        return -W / 2.0 - _safe_sqrt(
            R**2 - (di - z + H / 2.0) ** 2 - (x + L / 2.0) ** 2
        )

    def ymax6(x: float, z: float) -> float:
        return W / 2.0 + _safe_sqrt(
            R**2 - (di - z + H / 2.0) ** 2 - (x + L / 2.0) ** 2
        )

    total = _riemann_3d(total, zmin, zmax, dz, xmin6, xmax6, dx, ymin6, ymax6, dy, fun)
    I6.append(total)
```

**关键要点**

- 这一段是 `I6` 区域。
- 它从 `zmin = l1` 开始，到 `di + R + H/2` 结束。

### 4. 总概率

```python
I = I5 + np.array(I6) + np.array(I7) + np.array(I8)
idx = int(np.argmax(I))
return d, I, float(I[idx]), float(d[idx])
```

**关键要点**

- 这个区间的总概率只由 `I5`、`I6`、`I7`、`I8` 构成。
- 这里可以强调：深度越浅，直接覆盖项越少，代码里保留下来的积分项也越少。

---

## 问题二代码讲解，深度 87.5 到 100 米

## 本节展示

源码文件：`python/x5.py`

### 1. 深度区间和积分项

```python
d = _matlab_range(87.5, 1.0, 100.0)
I6 = []
I7 = []
I8 = []
```

**关键要点**

- 这里的投放深度范围是 87.5 米到 100 米。
- 和前面的文件相比，这里连 `I5` 都不再保留，只剩 `I6`、`I7`、`I8`。

### 2. 循环每一个定深

```python
for di in d:
    dx = 0.5
    dy = 0.5
    dz = 0.5
    zmin = l1
    zmax = di + R + 0.5 * H
```

**关键要点**

- 每一个候选定深 `di` 都要单独计算一次概率。
- `dx`、`dy`、`dz` 是三维网格积分步长。

### 3. I6 区域边界

```python
def xmin6(z: float) -> float:
    return -L / 2.0 - _safe_sqrt(R**2 - (di - z + H / 2.0) ** 2)

def xmax6(z: float) -> float:
    return -L / 2.0

def ymin6(x: float, z: float) -> float:
    return -W / 2.0 - _safe_sqrt(
        R**2 - (di - z + H / 2.0) ** 2 - (x + L / 2.0) ** 2
    )

def ymax6(x: float, z: float) -> float:
    return W / 2.0 + _safe_sqrt(
        R**2 - (di - z + H / 2.0) ** 2 - (x + L / 2.0) ** 2
    )

total = _riemann_3d(total, zmin, zmax, dz, xmin6, xmax6, dx, ymin6, ymax6, dy, fun)
I6.append(total)
```

**关键要点**

- 这里仍然使用球面边界。
- 平方根里的内容可以理解为：爆炸半径平方，减去深度方向已经消耗的距离平方，再减去 x 方向距离平方，剩下的才是 y 方向可覆盖范围。

### 4. 只合成 I6、I7、I8

```python
I = np.array(I6) + np.array(I7) + np.array(I8)
idx = int(np.argmax(I))
return d, I, float(I[idx]), float(d[idx])
```

**关键要点**

- 最后只把三项相加。
- 这个文件非常适合讲“分段建模”的意义：随着定深区间变化，有效命中区域也在减少。

---

## 问题二代码讲解，156 到 160 米精细搜索

## 本节展示

源码文件：`python/appendix_2_5_problem2_fine_depth_156_to_160.py` / `python/x6.py`

### 1. 包装文件入口

```python
from x6 import x6

def run():
    """Appendix 2.5: problem 2, d = 156:0.01:160."""
    return x6()

if __name__ == "__main__":
    d, I, peak, peak_d = run()
    print("d =", d)
    print("I =", I)
    print("peak =", peak)
    print("peak_d =", peak_d)
```

**关键要点**

- 这个 appendix 文件本身很短，它只是调用 `x6()`。
- 这里可以先讲清楚：真正的计算逻辑在 `x6.py`，这个文件只是方便直接运行和打印结果。

### 2. 精细搜索深度

```python
def x6() -> tuple[np.ndarray, np.ndarray, float, float]:
    sigma = 120.0
    L = 100.0
    R = 20.0
    W = 20.0
    H = 25.0
    sigma_z = 40.0
    l1 = 120.0
    h0 = 150.0

    d = _matlab_range(156.0, 0.01, 160.0)
```

**关键要点**

- `x6.py` 的参数和 `x2.py` 基本一致，变化最大的是这行深度搜索。
- `d = _matlab_range(156.0, 0.01, 160.0)` 表示从 156 米到 160 米，每隔 0.01 米取一个定深。

### 3. 概率计算结构仍然是 I1 到 I8

```python
I1 = []
for di in d:
    val = area_prob * quad(g_z, l1, di - R - H / 2.0)[0]
    I1.append(val)
I1 = np.array(I1)

I2 = []
I3 = []
I4 = []
I5 = np.array([quad(g_z, di - H / 2.0, di + H / 2.0)[0] for di in d])
I5 = 0.083734 * I5
I6 = []
I7 = []
I8 = []
```

**关键要点**

- 这一段说明 `x6.py` 沿用 `x2.py` 的八个区域结构。
- 不同的是，`d` 里面的点更多，所以同样一套积分要重复更多次，运行时间自然更长。

### 4. 找最大概率

```python
I = (
    I1
    + np.array(I2)
    + np.array(I3)
    + np.array(I4)
    + np.array(I5)
    + np.array(I6)
    + np.array(I7)
    + np.array(I8)
)
idx = int(np.argmax(I))
peak = float(I[idx])
peak_d = float(d[idx])
return d, I, peak, peak_d
```

**关键要点**

- 最后仍然是合成概率数组，找到最大值和对应定深。
- 这里重点讲：精细搜索提高的是最优深度的精度，不改变前面的概率模型。

---

## 问题三代码讲解，九枚深弹阵列

## 本节展示

源码文件：`python/x7.py`

### 1. 单枚深弹函数入口

```python
def compute_single_depth_charge(
    x_offset: float,
    y_offset: float,
    d: np.ndarray,
    sigma: float,
    L: float,
    R: float,
    W: float,
    H: float,
    sigma_z: float,
    l1: float,
    h0: float,
) -> np.ndarray:
```

**关键要点**

- 这个函数计算“一枚深弹在某个偏移位置上”的命中概率数组。
- `x_offset` 和 `y_offset` 是这一枚深弹相对于中心深弹的水平偏移。

### 2. 带偏移的水平概率密度

```python
def f(x: float, y: np.ndarray | float) -> np.ndarray | float:
    return (1.0 / (2.0 * math.pi * sigma**2)) * np.exp(
        -((x + x_offset) ** 2 + (y + y_offset) ** 2) / (2.0 * sigma**2)
    )
```

**关键要点**

- 这是问题三最关键的代码。
- 第一问和第二问里，水平误差是围绕中心点计算的。

### 3. 水平覆盖常数

```python
area_prob, _ = dblquad(
    f_y_x,
    -L / 2.0,
    L / 2.0,
    lambda _: -W / 2.0,
    lambda _: W / 2.0,
)

C_k1, _ = dblquad(f_y_x, -R - L / 2.0, -L / 2.0, ymin_c1, ymax_c1)
C_k2, _ = dblquad(f_y_x, -L / 2.0, L / 2.0, lambda _: -W / 2.0 - R, lambda _: W / 2.0 + R)
C_k3, _ = dblquad(f_y_x, L / 2.0, R + L / 2.0, ymin_c3, ymax_c3)
C_k = C_k1 + C_k2 + C_k3

I5 = C_k * I5
```

**关键要点**

- 这里计算的是带偏移后的水平覆盖概率。
- 因为深弹位置变了，水平概率密度也变了，所以不能继续直接使用问题一那个固定常数。

### 4. 九枚深弹阵列合成

```python
def x7() -> tuple[np.ndarray, np.ndarray, float, float]:
    a = L + 2.0 * R
    b = W + 2.0 * R

    d = _matlab_range(152.5, 1.0, 165.0)

    I_center = compute_single_depth_charge(0.0, 0.0, d, sigma, L, R, W, H, sigma_z, l1, h0)
    I_corner = compute_single_depth_charge(a, b, d, sigma, L, R, W, H, sigma_z, l1, h0)
    I_updown = compute_single_depth_charge(0.0, b, d, sigma, L, R, W, H, sigma_z, l1, h0)
    I_leftright = compute_single_depth_charge(a, 0.0, d, sigma, L, R, W, H, sigma_z, l1, h0)

    I = I_center + 4.0 * I_corner + 2.0 * I_updown + 2.0 * I_leftright
```

**关键要点**

- 这里是九枚深弹阵列的合成。
- `a = L + 2R` 是 x 方向的阵列间距，`b = W + 2R` 是 y 方向的阵列间距。

### 5. 输出最优结果

```python
idx = int(np.argmax(I))
peak = float(I[idx])
peak_d = float(d[idx])
return d, I, peak, peak_d
```

**关键要点**

- 最后和问题二一样，找到最大概率和对应定深。
- 这一节代码的核心就是一句话：把单枚深弹模型加上水平偏移，再利用阵列对称性合成九枚深弹结果。
