from __future__ import annotations

from x3 import x3


def run():
    """Appendix 2.2: problem 2, d = 140:1:152.5."""
    return x3()


if __name__ == "__main__":
    d, I, peak, peak_d = run()
    print("d =", d)
    print("I =", I)
    print("peak =", peak)
    print("peak_d =", peak_d)
