from __future__ import annotations

import math
from typing import Callable

import numpy as np
from scipy.integrate import quad, dblquad
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
        xl = xmin(z)
        xu = xmax(z)
        x_arr = _matlab_range(xl, dx, xu)
        if len(x_arr) == 0:
            continue
        for x in x_arr:
            yl = ymin(x, z)
            yu = ymax(x, z)
            y_arr = _matlab_range(yl, dy, yu)
            if len(y_arr) == 0:
                continue
            total += np.sum(fun(x, y_arr, z)) * dx * dy * dz
    return total


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
    def f(x: float, y: np.ndarray | float) -> np.ndarray | float:
        return (1.0 / (2.0 * math.pi * sigma**2)) * np.exp(
            -((x + x_offset) ** 2 + (y + y_offset) ** 2) / (2.0 * sigma**2)
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

    def f_y_x(y: float, x: float) -> float:
        return (1.0 / (2.0 * math.pi * sigma**2)) * math.exp(
            -((x + x_offset) ** 2 + (y + y_offset) ** 2) / (2.0 * sigma**2)
        )

    area_prob, _ = dblquad(f_y_x, -L / 2.0, L / 2.0, lambda _: -W / 2.0, lambda _: W / 2.0)

    I1 = []
    for di in d:
        val = area_prob * quad(g_z, l1, di - R - H / 2.0)[0]
        I1.append(val)
    I1 = np.array(I1)

    I2 = []
    I3 = []
    I4 = []
    I5 = np.array([quad(g_z, di - H / 2.0, di + H / 2.0)[0] for di in d])

    # 水平面覆盖积分概率常数。
    def ymin_c1(x: float) -> float:
        val = R**2 - (x + L / 2.0) ** 2
        return -W / 2.0 - math.sqrt(max(0.0, val))

    def ymax_c1(x: float) -> float:
        val = R**2 - (x + L / 2.0) ** 2
        return W / 2.0 + math.sqrt(max(0.0, val))

    def ymin_c3(x: float) -> float:
        val = R**2 - (x - L / 2.0) ** 2
        return -W / 2.0 - math.sqrt(max(0.0, val))

    def ymax_c3(x: float) -> float:
        val = R**2 - (x - L / 2.0) ** 2
        return W / 2.0 + math.sqrt(max(0.0, val))

    C_k1, _ = dblquad(f_y_x, -R - L / 2.0, -L / 2.0, ymin_c1, ymax_c1)
    C_k2, _ = dblquad(f_y_x, -L / 2.0, L / 2.0, lambda _: -W / 2.0 - R, lambda _: W / 2.0 + R)
    C_k3, _ = dblquad(f_y_x, L / 2.0, R + L / 2.0, ymin_c3, ymax_c3)
    C_k = C_k1 + C_k2 + C_k3

    I5 = C_k * I5
    I6 = []
    I7 = []
    I8 = []

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

        total = 0.0
        xmin_const = -0.5 * L
        xmax_const = 0.5 * L

        def ymin_z(z: float) -> float:
            return -W / 2.0 - _safe_sqrt(R**2 - (di - z - H / 2.0) ** 2)

        def ymax_z(z: float) -> float:
            return W / 2.0 + _safe_sqrt(R**2 - (di - z - H / 2.0) ** 2)

        for z in _matlab_range(zmin, dz, zmax):
            for x in _matlab_range(xmin_const, dx, xmax_const):
                y_arr = _matlab_range(ymin_z(z), dy, ymax_z(z))
                if len(y_arr) > 0:
                    total += np.sum(fun(x, y_arr, z)) * dx * dy * dz
        I3.append(total)

        total = 0.0

        def xmin4(z: float) -> float:
            return L / 2.0

        def xmax4(z: float) -> float:
            return L / 2.0 + _safe_sqrt(R**2 - (di - z - H / 2.0) ** 2)

        def ymin4(x: float, z: float) -> float:
            return -W / 2.0 - _safe_sqrt(
                R**2 - (di - z - H / 2.0) ** 2 - (x - L / 2.0) ** 2
            )

        def ymax4(x: float, z: float) -> float:
            return W / 2.0 + _safe_sqrt(
                R**2 - (di - z - H / 2.0) ** 2 - (x - L / 2.0) ** 2
            )

        total = _riemann_3d_vectorized(
            total, zmin, zmax, dz, xmin4, xmax4, dx, ymin4, ymax4, dy, fun
        )
        I4.append(total)

        total = 0.0
        zmin = di + 0.5 * H
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

        total = _riemann_3d_vectorized(
            total, zmin, zmax, dz, xmin6, xmax6, dx, ymin6, ymax6, dy, fun
        )
        I6.append(total)

        total = 0.0

        def xmin7(z: float) -> float:
            return -L / 2.0

        def xmax7(z: float) -> float:
            return L / 2.0

        def ymin7(x: float, z: float) -> float:
            return -W / 2.0 - _safe_sqrt(R**2 - (di - z + H / 2.0) ** 2)

        def ymax7(x: float, z: float) -> float:
            return W / 2.0 + _safe_sqrt(R**2 - (di - z + H / 2.0) ** 2)

        total = _riemann_3d_vectorized(
            total, zmin, zmax, dz, xmin7, xmax7, dx, ymin7, ymax7, dy, fun
        )
        I7.append(total)

        total = 0.0

        def xmin8(z: float) -> float:
            return L / 2.0

        def xmax8(z: float) -> float:
            return L / 2.0 + R

        def ymin8(x: float, z: float) -> float:
            return -W / 2.0 - _safe_sqrt(
                R**2 - (di - z + H / 2.0) ** 2 - (x - L / 2.0) ** 2
            )

        def ymax8(x: float, z: float) -> float:
            return W / 2.0 + _safe_sqrt(
                R**2 - (di - z + H / 2.0) ** 2 - (x - L / 2.0) ** 2
            )

        total = _riemann_3d_vectorized(
            total, zmin, zmax, dz, xmin8, xmax8, dx, ymin8, ymax8, dy, fun
        )
        I8.append(total)

    I_out = (
        I1
        + np.array(I2)
        + np.array(I3)
        + np.array(I4)
        + np.array(I5)
        + np.array(I6)
        + np.array(I7)
        + np.array(I8)
    )
    return I_out


def x7() -> tuple[np.ndarray, np.ndarray, float, float]:
    sigma = 120.0  # 水平投弹误差标准差
    L = 100.0      # 潜艇长度
    R = 20.0       # 深弹杀伤半径
    W = 20.0       # 潜艇宽度
    H = 25.0       # 潜艇高度
    sigma_z = 40.0 # 深度Z方向的标准差
    l1 = 120.0     # 截断深度
    h0 = 150.0     # 均值深度
    
    a = L + 2.0 * R
    b = W + 2.0 * R

    d = _matlab_range(152.5, 1.0, 165.0)

    # 3x3 阵列关于 X/Y 轴对称，只需计算中心、角点、上下、左右四类位置。
    I_center = compute_single_depth_charge(0.0, 0.0, d, sigma, L, R, W, H, sigma_z, l1, h0)
    I_corner = compute_single_depth_charge(a, b, d, sigma, L, R, W, H, sigma_z, l1, h0)
    I_updown = compute_single_depth_charge(0.0, b, d, sigma, L, R, W, H, sigma_z, l1, h0)
    I_leftright = compute_single_depth_charge(a, 0.0, d, sigma, L, R, W, H, sigma_z, l1, h0)

    # 按近似互斥事件叠加 9 枚深弹的命中概率。
    I = I_center + 4.0 * I_corner + 2.0 * I_updown + 2.0 * I_leftright
    
    idx = int(np.argmax(I))  # 寻找命中率最大值对应的索引
    peak = float(I[idx])     # 最大命中概率
    peak_d = float(d[idx])   # 对应的最优定深
    return d, I, peak, peak_d


if __name__ == "__main__":
    d, I, peak, peak_d = x7()
    print("==================================================")
    print("【第三问】 9 枚深弹阵列定深优化搜索完成！")
    print(f"最优投掷定深 d* = {peak_d:.2f} 米")
    print(f"最大联合命中概率 P_max = {peak:.6f} (即 {peak*100:.4f}%)")
    print("==================================================")
