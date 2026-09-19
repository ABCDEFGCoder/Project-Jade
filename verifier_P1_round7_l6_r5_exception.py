#!/usr/bin/env python3
"""
Exact finite verifier for Stage P1, Round 7:
the exceptional point (ell,r)=(6,5).

Claim verified:
There is no 45-edge convex lattice polygon whose edge directions are
pairwise distinct and which fits in a 45 x 45 axis-parallel box.

Reduction used by the research note:
- A closed lattice polygon with 45 distinct edge directions has total
  weighted l1 edge cost C = 2(W+H) <= 180.
- The 45 cheapest distinct primitive oriented directions have primitive
  cost 178: all directions from shells H_1,...,H_5, plus five from H_6.
- Hence only primitive costs 178,179,180 can occur. Closure makes C even,
  so the only weighted totals to test are C=178 and C=180.
- The complete shell-pattern list at primitive cost <=180 is enumerated
  below. Multiplicity surplus is also exhausted exactly.
No floating point arithmetic, MILP, random search, or external solver is used.
"""

from itertools import combinations
from math import gcd


def shell(L):
    out = []
    for x in range(-L, L + 1):
        ay = L - abs(x)
        ys = [0] if ay == 0 else [ay, -ay]
        for y in ys:
            if (x, y) == (0, 0):
                continue
            if gcd(abs(x), abs(y)) == 1:
                out.append((x, y))
    # deterministic de-duplication
    return list(dict.fromkeys(out))


H = {L: shell(L) for L in range(1, 9)}

assert {L: len(H[L]) for L in range(1, 9)} == {
    1: 4, 2: 4, 3: 8, 4: 8, 5: 16, 6: 8, 7: 24, 8: 16
}

LOW = H[1] + H[2] + H[3] + H[4] + H[5]
assert len(LOW) == 40
assert sum((abs(x) + abs(y)) for x, y in LOW) == 148


def evaluate(selected, extra=None):
    if extra is None:
        extra = {}
    sx = sy = W = Hh = cost = 0
    for v in selected:
        lam = 1 + extra.get(v, 0)
        x, y = v
        sx += lam * x
        sy += lam * y
        W += lam * max(x, 0)
        Hh += lam * max(y, 0)
        cost += lam * (abs(x) + abs(y))
    return sx, sy, W, Hh, cost


assert evaluate(LOW)[:2] == (0, 0)

tested = 0
feasible = []


def check(tag, selected, extra=None):
    global tested
    tested += 1
    ev = evaluate(selected, extra)
    if ev[0] == 0 and ev[1] == 0 and ev[2] <= 45 and ev[3] <= 45:
        feasible.append((tag, ev))


# C = 178: necessarily all low shells plus five H_6 directions,
# all multiplicities equal to one.
for S6 in combinations(H[6], 5):
    check("C178", LOW + list(S6))

# C = 180, primitive cost 178: multiplicity surplus exactly 2.
for S6 in combinations(H[6], 5):
    selected = LOW + list(S6)

    # +2 on one H_1 direction
    for v in H[1]:
        check("P178+2H1", selected, {v: 2})

    # +1 on two distinct H_1 directions
    for v, w in combinations(H[1], 2):
        check("P178+H1+H1", selected, {v: 1, w: 1})

    # +1 on one H_2 direction
    for v in H[2]:
        check("P178+H2", selected, {v: 1})

# C = 180, primitive cost 179, multiplicity surplus exactly 1.
# Pattern A: all low shells + four H_6 + one H_7.
for S6 in combinations(H[6], 4):
    for v7 in H[7]:
        selected = LOW + list(S6) + [v7]
        for e in H[1]:
            check("P179-A+H1", selected, {e: 1})

# Pattern B: omit one H_5, then take six H_6.
for miss in H[5]:
    low = [v for v in LOW if v != miss]
    for S6 in combinations(H[6], 6):
        selected = low + list(S6)
        for e in H[1]:
            check("P179-B+H1", selected, {e: 1})

# C = 180, primitive cost 180, no multiplicity surplus.

# Pattern A: all low shells + three H_6 + two H_7.
for S6 in combinations(H[6], 3):
    for S7 in combinations(H[7], 2):
        check("P180-A", LOW + list(S6) + list(S7))

# Pattern B: all low shells + four H_6 + one H_8.
for S6 in combinations(H[6], 4):
    for v8 in H[8]:
        check("P180-B", LOW + list(S6) + [v8])

# Pattern C: omit one H_5, then take five H_6 + one H_7.
for miss in H[5]:
    low = [v for v in LOW if v != miss]
    for S6 in combinations(H[6], 5):
        for v7 in H[7]:
            check("P180-C", low + list(S6) + [v7])

# Pattern D: omit one H_4, then take six H_6.
for miss in H[4]:
    low = [v for v in LOW if v != miss]
    for S6 in combinations(H[6], 6):
        check("P180-D", low + list(S6))

# Pattern E: omit two H_5, then take seven H_6.
for miss2 in combinations(H[5], 2):
    miss = set(miss2)
    low = [v for v in LOW if v not in miss]
    for S6 in combinations(H[6], 7):
        check("P180-E", low + list(S6))


assert tested == 48616, tested
assert feasible == [], feasible

print("PASS")
print("tested =", tested)
print("feasible =", len(feasible))
