"""
Verifies Theorem 5.3 (crooked-member bound): for any family Z of members
z (subset containing 0, |z|>=4, z subset of [-lambda,rho]) with pairwise
AP intersections,
    |Z| - kill(Z) <= floor(min(lambda,rho)/2)
where kill(Z) = number of distinct bad pairs contained in >= 1 member of Z.

This is formulated as an ILP:
  variables:  x_z in {0,1} for each candidate member z (|z|>=4, 0 in z)
              y_p in {0,1} for each candidate bad pair p
  maximize:   sum_z x_z  -  sum_p y_p
  subject to: x_z1 + x_z2 <= 1   whenever z1 ∩ z2 is empty or not an AP
              y_p >= x_z          whenever bad pair p is contained in z
  (y_p is forced to 1 exactly when some chosen z covers p, since the
   objective penalizes y_p and the solver will not set it higher than
   forced -- this correctly computes kill(Z) for the chosen Z)

Run: python3 verify_theorem_53.py
Checks the exact windows (lambda,rho) cited in the paper:
  (2,4), (2,5), (3,3), (3,4), (4,4)
each against the claimed bound floor(min(lambda,rho)/2).
"""
from itertools import combinations
import pulp


def is_ap(s):
    if len(s) <= 2:
        return True
    xs = sorted(s)
    d = xs[1] - xs[0]
    return all(xs[i + 1] - xs[i] == d for i in range(len(xs) - 1))


def candidate_members(lam, rho):
    universe = [x for x in range(-lam, rho + 1) if x != 0]
    members = []
    for r in range(3, len(universe) + 1):
        for combo in combinations(universe, r):
            z = frozenset(combo) | {0}
            members.append(z)
    return members


def bad_pairs_of(z):
    rest = [x for x in z if x != 0]
    return [frozenset((u, v)) for u, v in combinations(rest, 2) if not is_ap({0, u, v})]


def solve_window(lam, rho, verbose=False):
    members = candidate_members(lam, rho)
    n = len(members)

    # precompute bad pairs per member, and global pair -> covering members
    member_bad_pairs = [bad_pairs_of(z) for z in members]
    all_pairs = {}
    for i, bps in enumerate(member_bad_pairs):
        for p in bps:
            all_pairs.setdefault(p, []).append(i)

    # incompatibility: z_i, z_j incompatible if intersection empty or not AP
    incompatible = []
    for i in range(n):
        for j in range(i + 1, n):
            inter = members[i] & members[j]
            if not inter or not is_ap(inter):
                incompatible.append((i, j))

    prob = pulp.LpProblem("theorem_53", pulp.LpMaximize)
    x = [pulp.LpVariable(f"x_{i}", cat="Binary") for i in range(n)]
    y = {p: pulp.LpVariable(f"y_{k}", cat="Binary") for k, p in enumerate(all_pairs)}

    prob += pulp.lpSum(x) - pulp.lpSum(y.values())

    for (i, j) in incompatible:
        prob += x[i] + x[j] <= 1

    for p, covering in all_pairs.items():
        for i in covering:
            prob += y[p] >= x[i]

    prob.solve(pulp.PULP_CBC_CMD(msg=0))
    obj = pulp.value(prob.objective)

    if verbose:
        chosen = [members[i] for i in range(n) if x[i].value() > 0.5]
        print(f"  chosen family size: {len(chosen)}")

    return obj, n


if __name__ == "__main__":
    windows = [(2, 4), (2, 5), (3, 3), (3, 4), (4, 4)]
    for lam, rho in windows:
        bound = min(lam, rho) // 2
        obj, n_candidates = solve_window(lam, rho)
        obj_int = round(obj)
        status = "OK" if obj_int <= bound else "VIOLATION"
        match = "exact match" if obj_int == bound else f"(claimed exact max = {bound})"
        print(f"(lambda,rho)=({lam},{rho}): {n_candidates:3d} candidate members, "
              f"max(|Z|-kill(Z)) = {obj_int}, claimed bound = floor(min/2) = {bound}  "
              f"{status} {match}")
