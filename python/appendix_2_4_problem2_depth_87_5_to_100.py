from __future__ import annotations

from x5 import x5


def run():
    """Appendix 2.4: problem 2, d = 87.5:1:100."""
    return x5()


if __name__ == "__main__":
    d, I, peak, peak_d = run()
    print("d =", d)
    print("I =", I)
    print("peak =", peak)
    print("peak_d =", peak_d)
