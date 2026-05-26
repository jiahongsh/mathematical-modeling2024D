from __future__ import annotations

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
