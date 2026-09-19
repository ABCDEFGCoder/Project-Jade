#!/usr/bin/env python3
"""Search beta=2 theta-graph certificates with exact integer arithmetic."""

from __future__ import annotations

import argparse
import itertools
import json
import random
import time
from functools import lru_cache


def path_weights(m: int, r: int, x: int):
    ans = []
    for s in range(r):
        ans.extend((x + m * s, m * (m + 1) - x - m * s))
    return ans


def exact_partition(weights, target, groups, forbidden_pairs=None):
    n = len(weights)
    containing = [[] for _ in range(n)]
    # Positivity and target > max(weights) make singleton blocks impossible.
    # If every other block has at least two elements, the largest possible
    # block has n-2(groups-1) elements.  Computing this value keeps the
    # routine valid for the 2m+5 edge templates used in later rounds.
    max_group_size = n - 2 * (groups - 1)
    order = sorted(range(n), key=weights.__getitem__)
    sw = [weights[i] for i in order]
    prefix = [0]
    for w in sw:
        prefix.append(prefix[-1] + w)
    forbidden_pairs = set() if forbidden_pairs is None else {
        frozenset(pair) for pair in forbidden_pairs
    }

    def subsets_of_size(size):
        chosen = []
        def rec(start, left, total):
            if left == 0:
                if total == target:
                    candidate = tuple(order[i] for i in chosen)
                    if not any(frozenset((candidate[i], candidate[j]))
                               in forbidden_pairs
                               for i in range(len(candidate))
                               for j in range(i + 1, len(candidate))):
                        yield candidate
                return
            if len(sw) - start < left:
                return
            if total + prefix[start + left] - prefix[start] > target:
                return
            if total + prefix[len(sw)] - prefix[len(sw) - left] < target:
                return
            for pos in range(start, len(sw) - left + 1):
                nt = total + sw[pos]
                if nt >= target and left > 1:
                    break
                chosen.append(pos)
                yield from rec(pos + 1, left - 1, nt)
                chosen.pop()
        yield from rec(0, size, 0)

    for size in range(2, min(max_group_size, n) + 1):
        for comb in subsets_of_size(size):
            mask = sum(1 << i for i in comb)
            for i in comb:
                containing[i].append(mask)
    if any(not x for x in containing):
        return None

    full = (1 << n) - 1

    @lru_cache(None)
    def dfs(rem, left):
        if not rem:
            return () if left == 0 else None
        if (left <= 0 or rem.bit_count() < 2 * left
                or rem.bit_count() > max_group_size * left):
            return None
        # Branch on the uncovered item with the fewest currently applicable
        # target-sum subsets.
        candidates = []
        rr = rem
        while rr:
            bit = rr & -rr
            i = bit.bit_length() - 1
            opts = [s for s in containing[i] if s & rem == s]
            if not opts:
                return None
            candidates.append((len(opts), i, opts))
            rr ^= bit
        _, _, opts = min(candidates)
        for s in opts:
            tail = dfs(rem ^ s, left - 1)
            if tail is not None:
                return (s,) + tail
        return None

    masks = dfs(full, groups)
    if masks is None:
        return None
    return [[i for i in range(n) if mask >> i & 1] for mask in masks]


def compositions3(total):
    for a in range(1, total - 1):
        for b in range(a, total - a):
            c = total - a - b
            if b <= c:
                yield (a, b, c)


def search(m, seconds, seed):
    rng = random.Random(seed)
    rs = list(compositions3(m + 2))
    deadline = time.monotonic() + seconds
    trials = 0
    while time.monotonic() < deadline:
        r = rng.choice(rs)
        caps = [m * (m + 2 - a) - 1 for a in r]
        # Positive x_l with prescribed sum m(m+2).
        x1 = rng.randint(1, caps[0])
        lo = max(1, m * (m + 2) - x1 - caps[2])
        hi = min(caps[1], m * (m + 2) - x1 - 1)
        if lo > hi:
            continue
        x2 = rng.randint(lo, hi)
        x = (x1, x2, m * (m + 2) - x1 - x2)
        if x[2] <= 0 or x[2] > caps[2]:
            continue
        weights = sum((path_weights(m, a, b) for a, b in zip(r, x)), [])
        trials += 1
        part = exact_partition(weights, (m + 1) * (m + 2), m)
        if part is not None:
            return {"m": m, "r": r, "x": x, "weights": weights,
                    "partition": part, "trials": trials,
                    "target": (m + 1) * (m + 2)}
    return {"m": m, "trials": trials, "found": False}


def systematic(m):
    """Exhaust all positive integral theta parameters, up to path ordering."""
    total = m * (m + 2)
    target = (m + 1) * (m + 2)
    trials = 0
    for r in compositions3(m + 2):
        caps = [m * (m + 2 - a) - 1 for a in r]
        for x1 in range(1, caps[0] + 1):
            lo = max(1, total - x1 - caps[2])
            hi = min(caps[1], total - x1 - 1)
            for x2 in range(lo, hi + 1):
                x = (x1, x2, total - x1 - x2)
                if not (1 <= x[2] <= caps[2]):
                    continue
                weights = sum((path_weights(m, a, b)
                               for a, b in zip(r, x)), [])
                trials += 1
                part = exact_partition(weights, target, m)
                if part is not None:
                    return {"m": m, "r": r, "x": x, "weights": weights,
                            "partition": part, "trials": trials,
                            "target": target, "systematic": True}
    return {"m": m, "trials": trials, "found": False,
            "systematic": True}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("m", type=int)
    ap.add_argument("--seconds", type=float, default=600)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--systematic", action="store_true")
    args = ap.parse_args()
    ans = systematic(args.m) if args.systematic else search(
        args.m, args.seconds, args.seed)
    print(json.dumps(ans, ensure_ascii=False))


if __name__ == "__main__":
    main()
