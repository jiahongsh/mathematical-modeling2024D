from __future__ import annotations

from x2 import x2


def run():
    """Appendix 2.1: problem 2, d = 152.5:1:180."""
    return x2()


if __name__ == "__main__":
    d, I, peak, peak_d = run()
    print("d =", d)
    print("I =", I)
    print("peak =", peak)
    print("peak_d =", peak_d)
