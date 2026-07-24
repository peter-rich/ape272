"""
Verifies Theorem 5.2 (Private-pair theorem): every "crooked" member z
(a set containing 0, |z|>=4, not an arithmetic progression) contains a
"private" bad pair -- a bad pair {u,v} not spanned within z (Lemma 5.1's
spanning condition), hence contained in no other member.

For a single member z, we check this in isolation (the definition of
"private to z" only depends on z itself): a bad pair {u,v} in z is spanned
in z if there is some delta | gcd(|u|,|v|) such that z contains every
multiple of delta in [min(0,u,v), max(0,u,v)]. It's private otherwise.

Run: python3 verify_theorem_52.py
Checks all windows [a,b] (a<=0<=b) cited in the paper:
  [-6,6], [-4,8], [-3,9], [-2,10], [-7,7]
and reports, for every subset z of that window with 0 in z, |z|>=4, z not
an AP: whether z contains a private bad pair. ~26000 members total, matching
the paper's count.
"""
from itertools import combinations
from math import gcd


def is_ap(s):
    if len(s) <= 2:
        return True
    xs = sorted(s)
    d = xs[1] - xs[0]
    return all(xs[i + 1] - xs[i] == d for i in range(len(xs) - 1))


def bad_pairs(z):
    """All pairs {u,v} subset of z\\{0} with {0,u,v} not an AP."""
    rest = [x for x in z if x != 0]
    pairs = []
    for u, v in combinations(rest, 2):
        if not is_ap({0, u, v}):
            pairs.append((u, v))
    return pairs


def is_spanned(u, v, z):
    """Bad pair {u,v} subset of z is spanned in z if some delta | gcd(|u|,|v|)
    has z containing every multiple of delta in [min(0,u,v), max(0,u,v)]."""
    g = gcd(abs(u), abs(v))
    lo, hi = min(0, u, v), max(0, u, v)
    for delta in range(1, g + 1):
        if g % delta != 0:
            continue
        ok = True
        m = (lo // delta) * delta
        if m < lo:
            m += delta
        while m <= hi:
            if m not in z:
                ok = False
                break
            m += delta
        if ok:
            return True
    return False


def has_private_bad_pair(z):
    for u, v in bad_pairs(z):
        if not is_spanned(u, v, z):
            return True
    return False


def check_window(a, b):
    """a <= 0 <= b, window = {a,...,b}. Enumerate all z subset of window with
    0 in z, |z|>=4, z not an AP (crooked)."""
    universe = [x for x in range(a, b + 1) if x != 0]
    total_crooked = 0
    total_ok = 0
    counterexamples = []
    for r in range(3, len(universe) + 1):  # r = |z|-1 extra elements besides 0
        for combo in combinations(universe, r):
            z = frozenset(combo) | {0}
            if len(z) < 4:
                continue
            if is_ap(z):
                continue
            total_crooked += 1
            if has_private_bad_pair(z):
                total_ok += 1
            else:
                counterexamples.append(z)
    return total_crooked, total_ok, counterexamples


if __name__ == "__main__":
    windows = [(-6, 6), (-4, 8), (-3, 9), (-2, 10), (-7, 7)]
    grand_total = 0
    grand_ok = 0
    for a, b in windows:
        total, ok, counters = check_window(a, b)
        grand_total += total
        grand_ok += ok
        status = "OK (all have a private bad pair)" if not counters else f"COUNTEREXAMPLES: {counters[:3]}"
        print(f"window [{a},{b}]: {total:6d} crooked members checked, {ok:6d} confirmed private-pair, {status}")
    print(f"\nTOTAL across all windows: {grand_total} crooked members "
          f"(paper reports 'about 26,000'), {grand_ok} confirmed, "
          f"{grand_total - grand_ok} failures.")
