from __future__ import annotations

from x4 import x4


def run():
    """Appendix 2.3: problem 2, d = 100:1:140."""
    return x4()


if __name__ == "__main__":
    d, I, peak, peak_d = run()
    print("d =", d)
    print("I =", I)
    print("peak =", peak)
    print("peak_d =", peak_d)
