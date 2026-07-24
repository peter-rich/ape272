"""
Exact computation of t(N) = clique number of the graph on the 2^N - 1 nonempty
subsets of [N], with A ~ B iff A ∩ B is a nonempty arithmetic progression.

This is a from-scratch, independent check of Theorem 1.1: it does NOT use the
Section 2 construction at all, it only uses the definition of t(N).

Method: bitmask Bron-Kerbosch with pivoting for maximum clique. Python big
integers serve as bitsets, so adjacency row AND/OR/popcount are single machine
operations (via CPython's bigint arithmetic) even though N up to ~10 gives
graphs with over 1000 vertices.

This is intentionally a *reference/verification* implementation, not the
"bit-parallel branch-and-bound" / "star pruning" production solver described
in the paper for N = 11, 12; those require substantially more engineering to
finish in reasonable time. This script is practical up to about N = 8-9 within
a few minutes; see README for timing notes.

Run:  python3 exact_tN.py N
"""
import sys
from itertools import combinations


def is_ap(s):
    if len(s) <= 2:
        return True
    xs = sorted(s)
    d = xs[1] - xs[0]
    return all(xs[i + 1] - xs[i] == d for i in range(len(xs) - 1))


def all_nonempty_subsets(N):
    subsets = []
    for r in range(1, N + 1):
        for c in combinations(range(1, N + 1), r):
            subsets.append(frozenset(c))
    return subsets


def build_graph(N):
    """Returns (vertex_list, adj) where adj[i] is a bitmask (python int) of
    neighbours of vertex i."""
    verts = all_nonempty_subsets(N)
    n = len(verts)
    adj = [0] * n
    for i in range(n):
        Ai = verts[i]
        bit_i = 0
        for j in range(n):
            if i == j:
                continue
            Bj = verts[j]
            I = Ai & Bj
            if I and is_ap(I):
                bit_i |= (1 << j)
        adj[i] = bit_i
    return verts, adj


def popcount(x):
    return x.bit_count()  # Python 3.10+


def max_clique_bronkerbosch(adj, n):
    best = [0]
    full = (1 << n) - 1

    def bk(R_size, P, X):
        if P == 0 and X == 0:
            if R_size > best[0]:
                best[0] = R_size
            return
        if R_size + popcount(P) <= best[0]:
            return  # bound: can't beat current best
        # choose pivot u in P|X maximizing |P & adj[u]|
        PX = P | X
        u = -1
        max_deg = -1
        t = PX
        while t:
            v = (t & -t).bit_length() - 1
            cnt = popcount(P & adj[v])
            if cnt > max_deg:
                max_deg = cnt
                u = v
            t &= t - 1
        # candidates = P minus neighbours of pivot
        cand = P & ~adj[u] & full
        c = cand
        while c:
            vbit = c & -c
            v = vbit.bit_length() - 1
            bk(R_size + 1, P & adj[v], X & adj[v])
            P &= ~vbit
            X |= vbit
            c &= c - 1

    bk(0, full, 0)
    return best[0]


def t_N_exact(N):
    verts, adj = build_graph(N)
    n = len(verts)
    return max_clique_bronkerbosch(adj, n), n


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    t, n = t_N_exact(N)
    print(f"N={N}: graph has {n} vertices; exact t(N) = {t}")
