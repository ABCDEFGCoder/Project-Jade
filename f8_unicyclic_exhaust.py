#!/usr/bin/env python3
"""Exhaust the unicyclic (19-atom) case for m=8 over exact integers.

After multiplying the original masses by 360, Q-vertices have demand 40,
R-vertices demand 36, and every P-colour has demand 45.  The unicyclic
fractional-part lemma restricts the unique cycle to length 12 or 14.  This
program enumerates every bipartite forest attached to that cycle, computes
all tree weights by leaf elimination, enumerates the remaining cycle
parameter, and tests exact P-colourability.
"""

from __future__ import annotations

import argparse
import itertools
import json
import time

from theta_search import exact_partition


def q(name):
    return name[0] == "q"


def demand(name):
    return 40 if q(name) else 36


def forest_weights(off_q, off_r, roots, edges):
    """Return (tree edge weights, root residuals), or None."""
    adj = {v: [] for v in off_q + off_r + list(roots)}
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    if any(not adj[v] for v in off_q + off_r):
        return None
    seen = set()
    all_weights = []
    residual = {}
    rootset = set(roots)
    for start in off_q + off_r:
        if start in seen:
            continue
        stack = [start]
        comp = []
        comp_roots = set()
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v); comp.append(v)
            if v in rootset:
                comp_roots.add(v)
            stack.extend(w for w in adj[v] if w not in seen)
        if len(comp_roots) != 1:
            return None
        # A connected component must be a tree.  Count its internal edges.
        ecount = sum(len(adj[v]) for v in comp) // 2
        if ecount != len(comp) - 1:
            return None
        root = next(iter(comp_roots))
        parent = {root: None}
        order = [root]
        for v in order:
            for w in adj[v]:
                if w == parent[v]:
                    continue
                if w in parent:
                    return None
                parent[w] = v; order.append(w)
        edge_to_parent = {}
        for v in reversed(order[1:]):
            child_sum = sum(edge_to_parent[w] for w in adj[v]
                            if parent.get(w) == v)
            val = demand(v) - child_sum
            if val <= 0:
                return None
            edge_to_parent[v] = val
            all_weights.append(val)
        root_child_sum = sum(edge_to_parent[w] for w in adj[root]
                             if parent.get(w) == root)
        rem = demand(root) - root_child_sum
        if rem <= 0:
            return None
        residual[root] = rem
    return all_weights, residual


def cycle_weights(k, residual, t):
    # Vertex order qC0,rC0,qC1,rC1,...; edge i joins vertex i to i+1.
    verts = sum(([f"qC{i}", f"rC{i}"] for i in range(k)), [])
    ds = [residual.get(v, demand(v)) for v in verts]
    es = [t]
    for i in range(1, 2 * k):
        es.append(ds[i] - es[-1])
    if es[-1] + es[0] != ds[0] or any(x <= 0 for x in es):
        return None
    return es


def attachment_edge_sets(k):
    """Yield every labelled attached forest allowed by Q support degrees."""
    oqn = 9 - k
    orn = 10 - k
    off_q = [f"q{i}" for i in range(oqn)]
    off_r = [f"r{i}" for i in range(orn)]
    core_q = [f"qC{i}" for i in range(k)]
    core_r = [f"rC{i}" for i in range(k)]
    # Exactly 2*oqn plus one Q-incidence occurs outside the cycle.  Either
    # one core-Q attachment is present and every off-Q has degree 2, or no
    # core-Q attachment is present and q0 has degree 3.
    for core_attach in (0, 1):
        degrees = [2] * oqn
        if not core_attach:
            degrees[0] = 3
        choices = [list(itertools.combinations(off_r + core_r, d))
                   for d in degrees]
        # With a unique core-Q attachment, cycle rotation lets us fix its
        # core endpoint to qC0 without loss of generality.
        cq_options = [(None, None)] if not core_attach else [
            (core_q[0], rr) for rr in off_r]
        for nbr_tuple in itertools.product(*choices):
            base = []
            for u, nbrs in zip(off_q, nbr_tuple):
                base.extend((u, v) for v in nbrs)
            for cq, rr in cq_options:
                edges = list(base)
                if cq is not None:
                    edges.append((cq, rr))
                roots = set(v for e in edges for v in e if "C" in v)
                yield off_q, off_r, roots, edges


def run(verbose=False):
    start = time.monotonic()
    stats = {"candidate_edge_sets": 0, "positive_forests": 0,
             "cycle_parameter_trials": 0, "distinct_weight_multisets": 0}
    seen_weights = set()
    for k in (6, 7):
        for off_q, off_r, roots, edges in attachment_edge_sets(k):
            stats["candidate_edge_sets"] += 1
            fw = forest_weights(off_q, off_r, roots, edges)
            if fw is None:
                continue
            tree_w, residual = fw
            stats["positive_forests"] += 1
            for t in range(1, 40):
                cw = cycle_weights(k, residual, t)
                if cw is None:
                    continue
                stats["cycle_parameter_trials"] += 1
                weights = tuple(sorted(tree_w + cw))
                if weights in seen_weights:
                    continue
                seen_weights.add(weights)
                if sum(weights) != 360 or len(weights) != 19:
                    raise AssertionError("bad edge accounting")
                part = exact_partition(list(weights), 45, 8)
                if part is not None:
                    return {"feasible": True, "cycle_length": 2 * k,
                            "weights": weights, "partition": part,
                            "edges": edges, "t": t, "stats": stats,
                            "seconds": time.monotonic() - start}
            if verbose and stats["candidate_edge_sets"] % 100000 == 0:
                print(json.dumps(stats), flush=True)
    stats["distinct_weight_multisets"] = len(seen_weights)
    return {"feasible": False, "stats": stats,
            "seconds": time.monotonic() - start}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    print("RESULT_JSON=" + json.dumps(run(args.verbose), ensure_ascii=False))
