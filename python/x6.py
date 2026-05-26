from __future__ import annotations

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


def x6() -> tuple[np.ndarray, np.ndarray, float, float]:
    sigma = 120.0  # 标准差
    L = 100.0  # 潜艇长度
    R = 20.0  # 杀伤半径
    W = 20.0  # 潜艇宽度
    H = 25.0  # 高度
    sigma_z = 40.0  # Z的标准差
    l1 = 120.0
    h0 = 150.0

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

    d = _matlab_range(156.0, 0.01, 160.0)

    # I1 的水平矩形区域可分离，直接用 CDF 差值计算。
    prob_x = norm.cdf(L / (2.0 * sigma)) - norm.cdf(-L / (2.0 * sigma))
    prob_y = norm.cdf(W / (2.0 * sigma)) - norm.cdf(-W / (2.0 * sigma))
    area_prob = prob_x * prob_y
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
        zmin = di - R - 0.5 * H
        zmax = di - 0.5 * H
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
        zmin = di - R - 0.5 * H
        zmax = di - 0.5 * H

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
        zmin = di + 0.5 * H
        zmax = di + R + 0.5 * H

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
        zmin = di + 0.5 * H
        zmax = di + R + 0.5 * H

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


if __name__ == "__main__":
    d, I, peak, peak_d = x6()
    print("d =", d)
    print("I =", I)
    print("peak =", peak)
    print("peak_d =", peak_d)
