from __future__ import annotations

import math
from scipy.integrate import dblquad


def appendix_1_problem1_probability_p00() -> float:
    sigma = 120.0  # 投弹落点偏差的标准差（误差散布半径指标）
    L = 100.0      # 潜艇本体的有效长度
    R = 20.0       # 深弹的破坏/杀伤半径
    W = 20.0       # 潜艇本体的有效宽度

    # dblquad 要求被积函数参数顺序为 (y, x)。
    def f(y: float, x: float) -> float:
        return math.exp(-(x**2 + y**2) / (2.0 * sigma**2))

    def ymin1(x: float) -> float:
        # 浮点误差可能让边界处出现极小负数。
        val = R**2 - (x + L/2.0) ** 2
        return -W / 2.0 - math.sqrt(max(0.0, val))

    def ymax1(x: float) -> float:
        val = R**2 - (x + L/2.0) ** 2
        return W / 2.0 + math.sqrt(max(0.0, val))

    I1, _ = dblquad(f, -R - L / 2.0, -L / 2.0, ymin1, ymax1)

    I2, _ = dblquad(f, -L / 2.0, L / 2.0, lambda _: -W / 2.0 - R, lambda _: W / 2.0 + R)

    def ymin3(x: float) -> float:
        val = R**2 - (x - L/2.0) ** 2
        return -W / 2.0 - math.sqrt(max(0.0, val))

    def ymax3(x: float) -> float:
        val = R**2 - (x - L/2.0) ** 2
        return W / 2.0 + math.sqrt(max(0.0, val))

    I3, _ = dblquad(f, L / 2.0, R + L / 2.0, ymin3, ymax3)

    poe = (1.0 / (2.0 * math.pi * sigma**2)) * (I1 + I2 + I3)
    return poe


if __name__ == "__main__":
    poe = appendix_1_problem1_probability_p00()
    print(f"第一问单枚深弹命中概率 P(0,0) = {poe:.6f} (即约 {poe*100:.4f}%)")

    print("The probability p(0,0) is:", poe)
