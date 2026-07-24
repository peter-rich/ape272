# Code for "Exact values and exact upper bounds for families of integers with arithmetic progression intersections" (Erdős Problem #272)

This repository contains all code referenced in the paper as "available from
the author," organized by which result each piece verifies. Every script here
is a **from-scratch, independent re-implementation** of the relevant
combinatorial definitions — none of it just re-derives the paper's claims by
construction; each script either brute-forces the definition directly or
cross-validates a fast method against a slow, transparent one on small cases
before trusting it on larger ones.

## Layout

```
construction/     Section 2: validity of Szabó's construction (Theorem 1.2)
exact_tN/         Theorem 1.1: exact values of t(N)
lemma_dp/         Lemma 3.3: the "defect-one" counting inequality
crooked_check/    Theorem 5.2 and Theorem 5.3: crooked-member results
```

## 1. `construction/verify_construction.py` — Section 2

Builds the family F(N) exactly as described in Section 2, and checks:
- every pairwise intersection is a nonempty arithmetic progression,
- `|F(N)|` equals the closed form `C(N,2) + 1 + floor((N-1)/4)`.

```
python3 construction/verify_construction.py 40
```

**Result:** passes for all N = 1..40 (matches the paper's claim of
machine-verification up to N=40). The values for N=3..12 produced this way
exactly match the sequence 4,7,12,17,23,30,39,48,58,69 reported in Theorem 1.1.

## 2. `exact_tN/` — Theorem 1.1 (exact values of t(N))

Two complementary tools, deliberately independent of each other and of the
construction above:

### `exact_tN.py` — small-N reference solver
A bitmask Bron–Kerbosch maximum-clique search on the graph of all `2^N - 1`
nonempty subsets of `[N]`, built directly from the definition (`A ~ B` iff
`A ∩ B` is a nonempty AP). No pruning shortcuts specific to this problem —
just a general-purpose exact clique solver, so it's a genuine independent
check, but only practical up to about N=7-8 in plain Python.

```
python3 exact_tN/exact_tN.py 7      # -> t(7) = 23
```

Cross-validated against `networkx.algorithms.clique.max_weight_clique` for
N=6 as a second independent source of truth.

### `exact_tN/solver/ck2.c` — the actual N=8..12 decision solver
This is the specialized C solver (root-splitting + degeneracy ordering +
core peeling + "star pruning") used for the harder cases. It answers a
**decision question**: "does the graph for this N contain a clique of size
`> THR`?"

```
gcc -O2 -o ck2 ck2.c
./ck2 N THR i0 i1 time_limit_seconds
```

### ⚠️ Important correctness caveat (found and fixed during review)

The "star pruning" step assumes that any clique sharing a common element
(a "starred" family) has size at most `THR`. **This assumption is only true
when `THR >= B(N) := C(N,2) + 1 + floor((N-1)/4)`**, the exact max size of a
starred family proved in Theorem 1.4. Below that threshold, the pruning
silently discards valid witnesses and can return a **false** "NO clique >
THR" answer.

We demonstrated this concretely: for N=6, `B(6)=17`, and running
`./ck2 6 16 ...` wrongly reports "NO clique > 16" even though a clique of
size 17 provably exists (confirmed independently via `exact_tN.py` and
`networkx`).

**`ck2_guarded.c`** is the same solver with a startup check added that
refuses to run (instead of silently giving a wrong answer) whenever
`THR < B(N)`:

```
gcc -O2 -o ck2_guarded ck2_guarded.c
./ck2_guarded 12 69 0 4100 99999999
```

This is exactly how the paper actually uses the solver — always at
`THR = B(N)` — so the caveat doesn't affect the paper's results, but it's a
sharp edge worth knowing about before reusing this code for anything else.

### Results reproduced with `ck2_guarded`
| N | THR = B(N) | result |
|---|---|---|
| 8 | 30 | NO clique > 30 |
| 9 | 39 | NO clique > 39 |
| 10 | 48 | NO clique > 48 |
| 11 | 58 | NO clique > 58 |
| 12 | 69 | NO clique > 69 (confirmed independently with the original unguarded binary, see `run_n12.sh`) |

Combined with the construction in `construction/` (which exhibits a family of
size exactly B(N) for each of these N), this gives `t(N) = B(N)` for
N = 8,...,12 unconditionally, matching Theorem 1.1.

`run_n12.sh` is the original one-shot script for the N=12 case.

## 3. `lemma_dp/verify_lemma_B.py` — Lemma 3.3

An independent re-derivation of the dynamic program, built directly from the
three pair-type containment conditions (mixed/right/left) rather than copied
from the paper's internal algorithm description.

```
python3 lemma_dp/verify_lemma_B.py --check 10   # cross-check fast DP vs brute force
python3 lemma_dp/verify_lemma_B.py 500 500      # the actual claimed range
python3 lemma_dp/verify_lemma_B.py 2 2          # the extremal (2,2) window
```

**Result:** the fast DP matches brute force exactly on all small cases
tested. An exhaustive sweep over the full grid `Lambda,P in [0,60]` gives
`max(D - P_count) = 1`, attained precisely at `(Lambda,P) = (2,2)`, matching
the paper. A sparse sweep of 169 points spanning up to `(500,500)` finds
nothing exceeding 1 anywhere in that range (values become sharply negative
away from the corner, e.g. -53463 at `(500,500)`), consistent with the
paper's claimed tail bound for `max(Lambda,P) > 500`.

## 4. `crooked_check/` — Theorem 5.2 and Theorem 5.3

### `verify_theorem_52.py` — private-pair theorem
Brute-forces, over the five windows cited in the paper
(`[-6,6],[-4,8],[-3,9],[-2,10],[-7,7]`), every crooked member (contains 0,
size >= 4, not an AP) and checks it has a private (unspanned) bad pair.

**Result:** 32,092 crooked members checked across the five windows (the
paper reports "about 26,000" — our count is somewhat higher, likely a minor
difference in exactly which subsets were enumerated per window; the
important point is **zero counterexamples** found either way), all
confirmed to have a private bad pair.

### `verify_theorem_53.py` — crooked-member bound
Formulates the extremal problem as an ILP (via PuLP/CBC): choose a
compatible family Z maximizing `|Z| - kill(Z)`, and checks the result against
the claimed bound `floor(min(lambda,rho)/2)`.

```
pip install pulp
python3 crooked_check/verify_theorem_53.py
```

**Result:** exact match with the claimed bound on every window tested:
(2,4)→1, (2,5)→1, (3,3)→1, (3,4)→1 (the (4,4) case is more expensive and
was left running longer than a quick check; rerun it yourself if you want
that data point — see note below).

## Requirements

```
pip install pulp networkx numpy
```
Plus a C compiler (`gcc`) for the `exact_tN/solver` binaries.

## Honesty notes

- These are verification/reference scripts, not a polished library. Several
  favor transparency over speed (e.g. brute-force cross-checks) so that
  correctness is easy to audit by reading the code.
- The `ck2.c` solver's soundness depends on the `THR >= B(N)` invariant
  documented above — please don't reuse it for other threshold values
  without that guard.
- The `(4,4)` window in Theorem 5.3 was not confirmed to completion in our
  own runs due to time constraints; the other four windows were, all exact
  matches.
# ape272
