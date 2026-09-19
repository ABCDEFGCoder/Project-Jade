#!/usr/bin/env python3
"""Exact audit for the consolidated f(m) paper.

All checks use integer arithmetic.  The program verifies:
  * explicit tensor certificates for m=1,...,6,13,14;
  * the theta-path parameters and exact P-colour covers for m=7,...,12;
  * the stated support sizes and all three marginal families.
"""

from __future__ import annotations

import json

from f5_audit import F5_TENSOR
from round2_audit import audit_m14
from round3_audit import audit_m13


SMALL = {
    1: [(1, 1, 2, 2), (1, 1, 3, 1), (1, 2, 1, 2),
        (1, 2, 3, 1)],
    2: [(1, 1, 2, 2), (1, 1, 4, 6), (1, 3, 2, 4),
        (2, 2, 1, 6), (2, 2, 3, 2), (2, 3, 3, 4)],
    3: [(1, 1, 5, 3), (1, 2, 1, 10), (1, 4, 4, 7),
        (2, 1, 3, 12), (2, 4, 1, 2), (2, 4, 2, 6),
        (3, 2, 4, 5), (3, 3, 2, 6), (3, 3, 5, 9)],
    4: [(1, 1, 5, 12), (1, 2, 2, 18), (2, 1, 6, 12),
        (2, 3, 3, 16), (2, 5, 6, 2), (3, 2, 6, 6),
        (3, 4, 3, 4), (3, 4, 4, 20), (4, 3, 5, 8),
        (4, 5, 1, 20), (4, 5, 2, 2)],
    5: F5_TENSOR,
    6: [(1, 3, 6, 8), (1, 5, 3, 6), (1, 5, 8, 42),
        (2, 1, 5, 30), (2, 4, 5, 12), (2, 7, 7, 14),
        (3, 1, 1, 2), (3, 2, 7, 28), (3, 6, 2, 26),
        (4, 2, 4, 20), (4, 4, 3, 36), (5, 1, 2, 16),
        (5, 3, 1, 40), (6, 6, 4, 22), (6, 7, 6, 34)],
}


THETA = {
    7: {
        "r": (2, 3, 4), "X": (29, 25, 9),
        "groups": [
            [(2, 1), (3, 2)], [(2, 3), (3, 4)],
            [(2, 5), (3, 6)], [(1, 2), (1, 3), (3, 1)],
            [(2, 2), (2, 4), (2, 6)],
            [(1, 1), (1, 4), (3, 5)],
            [(3, 3), (3, 7), (3, 8)],
        ],
    },
    8: {
        "r": (2, 2, 6), "X": (30, 39, 11),
        "groups": [
            [(2, 3), (3, 9)], [(3, 2), (3, 10)],
            [(3, 4), (3, 8)], [(1, 2), (3, 5), (3, 12)],
            [(1, 1), (2, 4), (3, 7)],
            [(1, 3), (2, 2), (3, 3)],
            [(1, 4), (3, 1), (3, 6)],
            [(2, 1), (3, 11)],
        ],
    },
    9: {
        "r": (2, 4, 5), "X": (68, 10, 21),
        "groups": [
            [(1, 3), (3, 10)], [(2, 2), (3, 3)],
            [(2, 4), (3, 5)], [(2, 6), (3, 7)],
            [(1, 1), (3, 8)], [(2, 8), (3, 9)],
            [(1, 2), (2, 3), (3, 2)],
            [(1, 4), (2, 7), (3, 4)],
            [(2, 1), (2, 5), (3, 1), (3, 6)],
        ],
    },
    10: {
        "r": (2, 3, 7), "X": (26, 58, 36),
        "groups": [
            [(1, 2), (3, 10), (3, 14)], [(2, 3), (3, 4)],
            [(2, 5), (3, 6)], [(2, 4), (3, 7), (3, 12)],
            [(3, 3), (3, 11)], [(3, 5), (3, 9)],
            [(1, 1), (1, 4), (2, 6)], [(2, 1), (3, 2)],
            [(1, 3), (3, 13)], [(2, 2), (3, 1), (3, 8)],
        ],
    },
    11: {
        "r": (1, 6, 6), "X": (29, 69, 45),
        "groups": [
            [(1, 2), (2, 12), (3, 1)], [(2, 1), (3, 2)],
            [(2, 2), (2, 4), (2, 6)], [(2, 3), (3, 4)],
            [(2, 5), (3, 6)], [(2, 7), (3, 8)],
            [(1, 1), (2, 8), (2, 10), (3, 7)],
            [(2, 9), (3, 10)], [(2, 11), (3, 12)],
            [(3, 3), (3, 11)], [(3, 5), (3, 9)],
        ],
    },
    12: {
        "r": (2, 6, 6), "X": (62, 66, 40),
        "groups": [
            [(1, 2), (3, 9)], [(1, 4), (3, 11)],
            [(2, 7), (3, 8)], [(2, 9), (3, 10)],
            [(2, 11), (3, 12)], [(2, 1), (3, 2)],
            [(2, 8), (3, 3), (3, 7)],
            [(1, 3), (2, 6), (2, 10)],
            [(1, 1), (2, 2), (2, 12)],
            [(2, 5), (3, 6)], [(2, 3), (3, 4)],
            [(2, 4), (3, 1), (3, 5)],
        ],
    },
}


def audit_tensor(m: int, rows):
    p = [0] * m
    q = [0] * (m + 1)
    r = [0] * (m + 2)
    triples = set()
    for i, j, k, weight in rows:
        assert weight > 0
        assert (i, j, k) not in triples
        triples.add((i, j, k))
        p[i - 1] += weight
        q[j - 1] += weight
        r[k - 1] += weight
    assert p == [(m + 1) * (m + 2)] * m
    assert q == [m * (m + 2)] * (m + 1)
    assert r == [m * (m + 1)] * (m + 2)
    return {"support": len(rows), "P": p[0], "Q": q[0], "R": r[0]}


def theta_weights(m: int, radii, starts):
    weights = {}
    for ell, (radius, start) in enumerate(zip(radii, starts), 1):
        for s in range(radius):
            weights[(ell, 2 * s + 1)] = start + m * s
            weights[(ell, 2 * s + 2)] = m * (m + 1) - start - m * s
    return weights


def audit_theta(m: int, data):
    radii = data["r"]
    starts = data["X"]
    groups = data["groups"]
    assert sum(radii) == m + 2
    assert sum(starts) == m * (m + 2)
    assert all(0 < starts[i] < m * (m + 2 - radii[i])
               for i in range(3))
    weights = theta_weights(m, radii, starts)
    target = (m + 1) * (m + 2)
    assert len(weights) == 2 * m + 4
    assert len(groups) == m
    flattened = [edge for group in groups for edge in group]
    assert len(flattened) == len(set(flattened)) == len(weights)
    assert set(flattened) == set(weights)
    assert all(sum(weights[edge] for edge in group) == target
               for group in groups)
    return {
        "support": len(weights),
        "r": radii,
        "X": starts,
        "P_target": target,
    }


def main():
    report = {"small_tensor_certificates": {}, "theta_certificates": {}}
    for m, rows in SMALL.items():
        report["small_tensor_certificates"][m] = audit_tensor(m, rows)
    for m, data in THETA.items():
        report["theta_certificates"][m] = audit_theta(m, data)
    report["m13"] = audit_m13()
    report["m14"] = audit_m14()
    report["status"] = "all consolidated certificates exactly verified"
    print(json.dumps(report, ensure_ascii=True))


if __name__ == "__main__":
    main()
