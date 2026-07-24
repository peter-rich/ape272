r""" (Defect-one counting inequality) computationally.

Definitions (matching the paper, Section 4 exactly):
  A family S consists of integer intervals [-l, r] containing 0, with
  l + r >= 3, 0 <= l <= Lambda, 0 <= r <= P. The extremal S (maximizing
  D - P_count) is the down-set of an antichain, described by a nonincreasing
  staircase profile L(0) >= L(1) >= ... >= L(P) >= 0 with L(0) <= Lambda;
  members are all cells (l, r) with 0 <= l <= L(r) and l + r >= 3.

    D       = number of demand cells (l, r) in the down-set
    P_count = number of PRIMITIVE bad pairs {a,b} subset of [-Lambda,P]\{0}
              (gcd(|a|,|b|)=1, b not in {-a, 2a, a/2}) contained in >= 1 member

  Lemma 3.3 claims: D <= P_count + eps, eps = 1 if min(Lambda,P)>=2 else 0,
  i.e. max(D - P_count) over all profiles is 0 (if min(Lambda,P)<=1) or 1
  (if min(Lambda,P)>=2), the latter attained at the (2,2) window.

Three pair types, by direct inspection of the containment condition:
  - RIGHT  pair {a,b}, 0 < a < b: contained in [-l,r] iff r >= b (any l>=0
    works). Always covered once b <= P (and the demand cell exists, which
    needs b>=3 or l>=3-b; since l can be 0, need b>=3). Independent of the
    chosen profile once P is fixed -- constant contribution.
  - LEFT   pair {-a,-b}, 0 < a < b: contained in [-l,r] iff l >= b. Always
    covered once L(0) = Lambda >= b (we always choose L(0)=Lambda, which is
    optimal since it only helps). Independent of the r-profile -- constant.
  - MIXED  pair {-a,b}, a,b >= 1, a != b: contained in [-l,r] iff l>=a AND
    r>=b, i.e. (at column r=b) L(b) >= a. This is the ONLY part that
    depends on the chosen profile L(1..P).

So max(D - P_count) = [D(column 0) - 0] + max_profile sum_r [D(col r,L(r))
  - mixed_covered(r,L(r))]  -  RIGHT_count(P)  -  LEFT_count(Lambda)

with D(col 0) = #{l in [0,Lambda] : l>=3} = max(Lambda-2,0).

Run: python3 verify_lemma_B.py [LAMBDA] [P]        -> runs the fast DP
     python3 verify_lemma_B.py --check K           -> brute-force cross-check
                                                       for all Lambda,P<=K
"""
import sys
from math import gcd


_phi_cache = {}


def _phi_sieve(n):
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:  # p is prime
            for m in range(p, n + 1, p):
                phi[m] -= phi[m] // p
    return phi


def right_or_left_count(M):
    """#{(a,b): 3<=b<=M, 1<=a<b, gcd(a,b)=1, b != 2a}.
    Equals sum_{b=3}^{M} phi(b): the exclusion "b != 2a" is vacuous for
    b>=3, since gcd(a,2a)=a=1 forces (a,b)=(1,2), outside this range.
    (Verified against the direct O(M^2) definition for M<=200.)"""
    if M < 3:
        return 0
    key = "phi"
    if key not in _phi_cache or len(_phi_cache[key]) <= M:
        _phi_cache[key] = _phi_sieve(max(M, 512))
    phi = _phi_cache[key]
    return sum(phi[3:M + 1])


def demand_count(l, r):
    """#{ll in [0,l] : ll + r >= 3}, i.e. size of column r up to level l."""
    lo = max(0, 3 - r)
    if l < lo:
        return 0
    return l - lo + 1


def mixed_covered_count(l, r):
    """#{a in [1,l] : gcd(a,r)=1, a != r} -- primitive mixed pairs {-a,r}
    covered once L(r) = l."""
    cnt = 0
    for a in range(1, l + 1):
        if gcd(a, r) == 1 and a != r:
            cnt += 1
    return cnt


def col_score(r, l):
    return demand_count(l, r) - mixed_covered_count(l, r)


def fast_dp(Lambda, P):
    """Exact DP (not an approximation): returns max(D - P_count).
    Runs in O(Lambda*P) time: for each column r, the per-l mixed-coverage
    count and demand count are both accumulated incrementally rather than
    recomputed from scratch."""
    if P == 0:
        D0 = max(Lambda - 2, 0)
        left = right_or_left_count(Lambda)
        return D0 - left

    def col_scores(r):
        """Returns array cs[0..Lambda] = col_score(r, l) for all l at once,
        in O(Lambda) time via incremental accumulation."""
        cs = [0] * (Lambda + 1)
        mixed_running = 0
        for l in range(0, Lambda + 1):
            if l >= 1 and gcd(l, r) == 1 and l != r:
                mixed_running += 1
            cs[l] = demand_count(l, r) - mixed_running
        return cs

    dp = col_scores(P)
    for r in range(P - 1, 0, -1):
        cs = col_scores(r)
        new_dp = [0] * (Lambda + 1)
        running_max = -10**9
        for l in range(Lambda + 1):
            running_max = max(running_max, dp[l])
            new_dp[l] = cs[l] + running_max
        dp = new_dp
    best_cols = max(dp)
    D0 = max(Lambda - 2, 0)
    right = right_or_left_count(P)
    left = right_or_left_count(Lambda)
    return D0 + best_cols - right - left


def brute_force(Lambda, P):
    """Exhaustive search over ALL nonincreasing profiles L(1)>=...>=L(P)>=0,
    L(1)<=Lambda. Exponential -- only for small Lambda, P (sanity check)."""
    best = [-10**9]

    def rec(r, cap, acc):
        if r > P:
            best[0] = max(best[0], acc)
            return
        for l in range(cap, -1, -1):
            rec(r + 1, l, acc + col_score(r, l))

    rec(1, Lambda, 0)
    D0 = max(Lambda - 2, 0)
    right = right_or_left_count(P)
    left = right_or_left_count(Lambda)
    return D0 + best[0] - right - left


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--check":
        K = int(sys.argv[2]) if len(sys.argv) > 2 else 8
        print(f"Cross-checking fast_dp vs brute_force for all Lambda,P in [0,{K}]:")
        mismatches = 0
        for Lam in range(0, K + 1):
            for P in range(0, K + 1):
                a = fast_dp(Lam, P)
                b = brute_force(Lam, P)
                status = "OK" if a == b else "MISMATCH"
                if a != b:
                    mismatches += 1
                    print(f"  Lambda={Lam:2d} P={P:2d}: fast_dp={a:3d} brute={b:3d}  {status}")
        if mismatches == 0:
            print("All matched. fast_dp is validated against brute force on this range.")
        else:
            print(f"{mismatches} mismatches found -- fast_dp has a bug, do not trust larger runs.")
    else:
        Lambda = int(sys.argv[1]) if len(sys.argv) > 1 else 500
        P = int(sys.argv[2]) if len(sys.argv) > 2 else 500
        val = fast_dp(Lambda, P)
        eps = 1 if min(Lambda, P) >= 2 else 0
        print(f"Lambda={Lambda}, P={P}: max(D - P_count) = {val}  "
              f"(Lemma 3.3 claims <= {eps}; {'OK' if val <= eps else 'VIOLATION'})")
