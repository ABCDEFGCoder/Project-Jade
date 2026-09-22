#!/usr/bin/env python3
"""Exact certificate for the finite boundary signs in the integer staircase.

This file is part of the proof only for the finitely many cases with n >= 6
listed in `finite_cases()`.  The pentagon row n=5 is proved analytically and is
intentionally outside this certificate.  Every theorem-level decision uses
fractions only.
A rational enclosure for pi is derived internally from Machin's identity
    pi = 16 arctan(1/5) - 4 arctan(1/239)
using alternating-series remainders.  Sine values are then enclosed by
alternating Taylor sums on [0, pi/2].

No floating-point number is used to decide a sign.
"""
from fractions import Fraction
from functools import lru_cache


def c_q(q: int) -> int:
    return q * (q - 2) // 2


def arctan_inv_interval(q: int, last_index: int):
    """Enclose arctan(1/q) by an alternating partial sum.

    Sum terms k=0,...,last_index.  If the last term is positive (even k),
    the partial sum is an upper bound; if negative, it is a lower bound.
    The adjacent partial sum supplies the opposite bound.
    """
    x = Fraction(1, q)
    terms = [((-1) ** k) * x ** (2 * k + 1) / (2 * k + 1)
             for k in range(last_index + 2)]
    s0 = sum(terms[:last_index + 1], Fraction(0))
    s1 = s0 + terms[last_index + 1]
    return (s1, s0) if last_index % 2 == 0 else (s0, s1)


def pi_interval():
    # More than enough for all subsequent sign separations.
    a5 = arctan_inv_interval(5, 18)
    a239 = arctan_inv_interval(239, 4)
    lo = 16 * a5[0] - 4 * a239[1]
    hi = 16 * a5[1] - 4 * a239[0]
    assert lo < hi
    # Independent coarse sanity checks, still exact.
    assert Fraction(333, 106) < lo
    assert hi < Fraction(355, 113)
    return lo, hi


PI_LO, PI_HI = pi_interval()


def sin_point_interval(y: Fraction, depth: int = 12):
    """Rigorous enclosure for sin(y), 0 <= y <= PI_HI/2."""
    assert 0 <= y and 2 * y <= PI_HI
    # On [0,pi/2] the sine Taylor terms decrease in magnitude from the
    # first step onward, since y^2/6 < 1; hence the alternating partial
    # sums below are rigorous lower/upper bounds.
    assert y * y < 6
    total = Fraction(0)
    term = y
    upper = lower = None
    for k in range(2 * depth + 2):
        if k == 0:
            term = y
        else:
            term = term * y * y / ((2 * k) * (2 * k + 1))
        total += term if k % 2 == 0 else -term
        if k == 2 * depth:
            upper = total
        elif k == 2 * depth + 1:
            lower = total
    assert lower <= upper
    return lower, upper


@lru_cache(None)
def sin_coeff_pi_interval(num: int, den: int):
    coeff = Fraction(num, den)
    if coeff < 0:
        lo, hi = sin_coeff_pi_interval(-num, den)
        return -hi, -lo
    x_lo = coeff * PI_LO
    x_hi = coeff * PI_HI
    assert 2 * x_hi <= PI_HI
    return sin_point_interval(x_lo)[0], sin_point_interval(x_hi)[1]


def add_interval(A, B):
    return A[0] + B[0], A[1] + B[1]


def scale_interval(c, A):
    c = Fraction(c)
    return (c * A[0], c * A[1]) if c >= 0 else (c * A[1], c * A[0])


def centered_sine_terms(d: int, q: int, z: int):
    c = c_q(q)
    n = c * d + z
    m = n + d
    terms = []
    if q % 2 == 0:
        r = q // 2
        s = q - 1
        terms.append((1, Fraction(s, m)))
        for k in range(1, r):
            terms.append((2, Fraction(s, m) - Fraction(2 * k, n)))
            terms.append((-2, Fraction(2 * k, m)))
    else:
        r = (q - 1) // 2
        s = q - 1
        terms.append((1, Fraction(s, m)))
        terms.append((-1, Fraction(s * d, n * m)))
        for j in range(1, r):
            terms.append((2, Fraction(2 * j, m)
                          - Fraction(2 * (r - j) * d, n * m)))
        for j in range(1, r + 1):
            terms.append((-1, Fraction(2 * j - 1, m)
                          - Fraction(d, n * m)))
        for j in range(1, r):
            terms.append((-1, Fraction(2 * j + 1, m)
                          + Fraction(d, n * m)))
    return n, m, terms


def centered_S_interval(d: int, q: int, z: int):
    n, m, terms = centered_sine_terms(d, q, z)
    out = (Fraction(0), Fraction(0))
    for coeff, theta_coeff in terms:
        I = sin_coeff_pi_interval(theta_coeff.numerator, theta_coeff.denominator)
        out = add_interval(out, scale_interval(coeff, I))
    return n, m, out


def certify(d: int, q: int, z: int, positive: bool):
    n, m, I = centered_S_interval(d, q, z)
    if positive:
        assert I[0] > 0, (d, q, z, n, m, I)
    else:
        assert I[1] < 0, (d, q, z, n, m, I)
    return d, q, z, n, m, I


def finite_cases():
    """Return every finite sign needed by the n>=6 integer-threshold proof.

    The analytic part proves, for d=1, positivity for z>=4 and (for q>=17)
    negativity for z<=3.  Hence for 4<=q<=16 every still-undecided z in
    {1,2,3} must be checked separately (subject to n>=6).  When the asserted
    threshold is z=4 we also retain the z=4 positive check as an independent
    boundary audit.  In particular, for (d,q)=(1,4) both z=2 and z=3 are
    separately certified positive; no monotonicity in z is used.
    """
    cases = []

    # d=1, 4 <= q <= 16.
    zmin = {
        4: 2,
        5: 3, 6: 3, 7: 3, 8: 3,
        9: 4,
        10: 3,
        11: 4, 12: 4, 13: 4, 14: 4, 15: 4, 16: 4,
    }
    for q in range(4, 17):
        z0 = zmin[q]
        # z>=4 is covered analytically.  Thus z=1,2,3 are precisely the
        # potentially unresolved layers; exclude only n=5, which lies outside
        # the n>=6 staircase theorem.
        for z in range(1, 4):
            if c_q(q) + z >= 6:
                cases.append((1, q, z, z >= z0))
        # For rows whose threshold itself is z=4, keep an exact positive
        # boundary check as a redundant audit of the analytic sufficiency.
        if z0 == 4:
            cases.append((1, q, 4, True))

    # d=2: only q=4,z=1 is not decided by the analytic tail.
    cases.append((2, 4, 1, False))

    # d=3: q=4,z=1 is positive; q=5,6,z=1 are negative.  For q>=7 the
    # analytic lower shadow handles z=1, while z>=2 is analytically feasible.
    cases.extend([
        (3, 4, 1, True),
        (3, 5, 1, False),
        (3, 6, 1, False),
    ])
    return cases


def audit_case_coverage(cases):
    """Mechanically verify that the finite proof ledger has no omitted layer."""
    zmin = {
        4: 2,
        5: 3, 6: 3, 7: 3, 8: 3,
        9: 4,
        10: 3,
        11: 4, 12: 4, 13: 4, 14: 4, 15: 4, 16: 4,
    }
    case_map = {(d, q, z): positive for d, q, z, positive in cases}
    assert len(case_map) == len(cases), "duplicate finite certificate case"

    # Every d=1 unresolved layer z=1,2,3 with n>=6 must occur, with the sign
    # prescribed by the claimed threshold zmin[q].
    for q in range(4, 17):
        for z in range(1, 4):
            n = c_q(q) + z
            if n >= 6:
                key = (1, q, z)
                assert key in case_map, ("missing finite layer", key, n, n + 1)
                assert case_map[key] == (z >= zmin[q]), (
                    "wrong expected sign", key, case_map[key], zmin[q])
        if zmin[q] == 4:
            assert case_map.get((1, q, 4)) is True

    # The exact small-d boundary ledger.
    assert case_map.get((2, 4, 1)) is False
    assert case_map.get((3, 4, 1)) is True
    assert case_map.get((3, 5, 1)) is False
    assert case_map.get((3, 6, 1)) is False

    # The previously omitted point must be present explicitly.
    assert case_map.get((1, 4, 3)) is True
    assert len(cases) == 49, len(cases)



def audit_analytic_tail_boundaries():
    """Redundant exact checks at the first points covered analytically.

    These cases are not needed for finite-case completeness.  They guard
    transcription of the analytic tail interfaces in the paper.
    """
    sentinels = [
        (2, 5, 1, False),  # first odd d=2 tail
        (2, 6, 1, False),  # first even d=2 tail
        (3, 7, 1, False),  # first odd d=3 tail
        (3, 8, 1, False),  # first even d=3 tail
        (1, 17, 3, False), # first odd d=1 infinite tail
        (1, 18, 3, False), # first even d=1 infinite tail
    ]
    for case in sentinels:
        certify(*case)
    return len(sentinels)

def main():
    cases = finite_cases()
    audit_case_coverage(cases)
    tail_sentinels = audit_analytic_tail_boundaries()
    assert all(c_q(q) * d + z >= 6 for d, q, z, _ in cases)
    rows = [certify(*case) for case in cases]
    margins = []
    for d, q, z, n, m, I in rows:
        gap = I[0] if I[0] > 0 else -I[1]
        margins.append((gap, d, q, z, n, m, I))
    margins.sort(key=lambda x: x[0])
    g, d, q, z, n, m, I = margins[0]
    print("integer staircase certificate: PASS")
    print("derived pi enclosure:", PI_LO, "< pi <", PI_HI)
    print("certified finite signs:", len(rows))
    print("finite-case ledger coverage audit: PASS")
    print("analytic-tail boundary cross-audit: PASS", f"({tail_sentinels} sentinels)")
    s_d, s_q, s_z, s_n, s_m, s_I = certify(1, 4, 3, True)
    print("coverage sentinel d=1,q=4,z=3, (n,m)=(7,8): positive")
    print("smallest exact distance from zero occurs at",
          f"d={d}, q={q}, z={z}, (n,m)=({n},{m})")
    sign = "positive" if I[0] > 0 else "negative"
    print("certified interval sign:", sign)
    print("exact gap numerator bits:", g.numerator.bit_length())
    print("exact gap denominator bits:", g.denominator.bit_length())
    print("All theorem-level computations used exact Fraction arithmetic.")


if __name__ == "__main__":
    main()
