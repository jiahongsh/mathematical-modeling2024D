from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm

from x3 import _matlab_range, _riemann_3d, _safe_sqrt


def x5():
    sigma = 120.0
    L = 100.0
    R = 20.0
    W = 20.0
    H = 25.0
    sigma_z = 40.0
    l1 = 120.0
    h0 = 150.0

    def f(x: float, y: float) -> float:
        return (1.0 / (2.0 * math.pi * sigma**2)) * math.exp(
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

    def fun(x: float, y: float, z: float) -> float:
        return f(x, y) * g_z(z)

    d = _matlab_range(87.5, 1.0, 100.0)
    I6 = []
    I7 = []
    I8 = []

    for di in d:
        dx = 0.5
        dy = 0.5
        dz = 0.5
        zmin = l1
        zmax = di + R + 0.5 * H

        total = 0.0

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

        total = 0.0

        def xmin7(z: float) -> float:
            return -L / 2.0

        def xmax7(z: float) -> float:
            return L / 2.0

        def ymin7(x: float, z: float) -> float:
            return -W / 2.0 - _safe_sqrt(R**2 - (di - z + H / 2.0) ** 2)

        def ymax7(x: float, z: float) -> float:
            return W / 2.0 + _safe_sqrt(R**2 - (di - z + H / 2.0) ** 2)

        total = _riemann_3d(total, zmin, zmax, dz, xmin7, xmax7, dx, ymin7, ymax7, dy, fun)
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

        total = _riemann_3d(total, zmin, zmax, dz, xmin8, xmax8, dx, ymin8, ymax8, dy, fun)
        I8.append(total)

    I = np.array(I6) + np.array(I7) + np.array(I8)
    idx = int(np.argmax(I))
    return d, I, float(I[idx]), float(d[idx])


if __name__ == "__main__":
    d, I, peak, peak_d = x5()
    print("d =", d)
    print("I =", I)
    print("peak =", peak)
    print("peak_d =", peak_d)
