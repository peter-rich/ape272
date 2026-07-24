"""
Verifies the construction of Section 2 ("A construction achieving Szabo's bound").

For a given N, builds the family F described in the paper and checks:
  (a) every pairwise intersection A_i ∩ A_j is a nonempty arithmetic progression (AP);
  (b) |F| equals the claimed closed form C(N,2) + 1 + floor((N-1)/4).

Run:  python3 verify_construction.py [N_max]
Default N_max = 40, matching the paper's claim
"The construction has been machine-verified (full pairwise check) for all N <= 40."
"""
import sys
from itertools import combinations
from math import comb


def is_ap(s):
    """A finite set of integers is an AP iff it has size <= 2, or (sorted) has
    constant consecutive difference."""
    if len(s) <= 2:
        return True
    xs = sorted(s)
    d = xs[1] - xs[0]
    return all(xs[i + 1] - xs[i] == d for i in range(len(xs) - 1))


def build_family(N):
    m = -(-N // 2)  # ceil(N/2)
    k = (N - 1) // 4
    F = []

    # (i) {m} and all pairs {m, x}, x != m
    F.append(frozenset({m}))
    for x in range(1, N + 1):
        if x != m:
            F.append(frozenset({m, x}))

    # blocked triples B_d, B'_d for 1 <= d <= k
    blocked = set()
    for d in range(1, k + 1):
        Bd = frozenset({m - 2 * d, m, m + d})
        Bpd = frozenset({m - d, m, m + 2 * d})
        blocked.add(Bd)
        blocked.add(Bpd)

    # (ii) all 3-sets {m, u, v} except the blocked ones
    others = [x for x in range(1, N + 1) if x != m]
    for u, v in combinations(others, 2):
        T = frozenset({m, u, v})
        if T not in blocked:
            F.append(T)

    # (iii) for each d, three APs: P_d and its two 4-term sub-APs containing m
    for d in range(1, k + 1):
        Pd = frozenset({m - 2 * d, m - d, m, m + d, m + 2 * d})
        F.append(Pd)
        F.append(Pd - {m + 2 * d})
        F.append(Pd - {m - 2 * d})

    return F, m, k


def verify(N, verbose=False):
    F, m, k = build_family(N)

    # all members must be subsets of [N] and nonempty
    for A in F:
        assert A, "empty set found"
        assert all(1 <= x <= N for x in A), f"member {A} escapes [N] for N={N}"

    # distinctness
    assert len(F) == len(set(F)), f"duplicate members found for N={N}"

    # pairwise intersection nonempty AP
    for A, B in combinations(F, 2):
        I = A & B
        if not I:
            return False, f"empty intersection: {sorted(A)} vs {sorted(B)}"
        if not is_ap(I):
            return False, f"non-AP intersection {sorted(I)} from {sorted(A)} vs {sorted(B)}"

    expected = comb(N, 2) + 1 + k
    if len(F) != expected:
        return False, f"|F|={len(F)} but expected C(N,2)+1+floor((N-1)/4)={expected}"

    if verbose:
        print(f"N={N:3d}: |F|={len(F):6d} (= C(N,2)+1+floor((N-1)/4) = {expected}), OK")
    return True, expected


if __name__ == "__main__":
    N_max = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    all_ok = True
    for N in range(1, N_max + 1):
        ok, info = verify(N, verbose=True)
        if not ok:
            print(f"FAILED at N={N}: {info}")
            all_ok = False
            break
    if all_ok:
        print(f"\nAll pairwise-intersection and cardinality checks passed for N = 1..{N_max}.")
