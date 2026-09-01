# Tie-break sensitivity

Barcodes that share a read count have no natural order. The abundance-ranked
greedy pass has to visit them in *some* order, and whichever comes first is the
one allowed to seed a cluster, so the choice moves a handful of barcodes between
"recovered" and "absorbed". `super_cluster2()` fixes the order deterministically
(`tie_break = "sequence"`) so a result never depends on how the input rows
happened to be sorted — but the rule it settles on is still one arbitrary choice
among many equally defensible ones.

This directory measures how much a result rests on that choice, which is the
error bar the small between-method differences should be read against.

## Running it

```sh
Rscript run_tiebreak_sensitivity.R ../indel_experiment/results 15
Rscript plot_tiebreak_sensitivity.R
```

The first script clusters each benchmark condition once with the shipped default
and once per alternative tie order (`tie_break = "hash"`, `tie_seed = 1..n`),
scoring each run with the same FN / FP / WS definitions as the Python benchmark
harness, and writes `tiebreak_seeds.csv`. The second draws
`tiebreak_sensitivity.png` from that file plus `reference_methods.csv`, which
holds one run of each competing method on the same inputs.

Every tie order is itself fully reproducible: a given `tie_seed` always yields
the same clustering, and none of them depend on input row order.

## What it shows

False negatives across 15 alternative tie orders, against the shipped default:

| condition | default | alternative orders (min–median–max) |
|---|---|---|
| Random 20 bp, no indels | 54 | 45 – 48 – 53 |
| Random 20 bp, 0.5% indels | 170 | 146 – 158 – 169 |
| Anchored 28 bp, no indels | 11 | 12 – 14 – 17 |
| Anchored 28 bp, 0.5% indels | 45 | 42 – 47 – 50 |

Two things follow.

The spread is wider than the gaps between methods. On the random 20 bp designs
every competing method lands inside it — Starcode's 158 false negatives under
0.5% indels is exactly barbac's median across tie orders. Differences of that
size are not evidence that one algorithm recovers more barcodes than another.

The default ordering is not systematically good or bad. It is the worst of the
16 orders on random 20 bp barcodes and the best of them on the anchored design,
which is what luck looks like: sorting ties by the barcode string is neither
better nor worse than any other rule, it simply lands somewhere different in each
dataset. Picking the ordering that scores best on a benchmark would be fitting to
that benchmark's answer key, so the default stays fixed and the spread is
reported instead.

Order sensitivity is a separate property from this one, and worth keeping
distinct: permuting the input rows leaves barbac, Starcode and Bartender
unchanged, while Shepherd's false negatives move (44 → 47 on random 20 bp with no
indels). barbac behaved the same way until the tie-break was fixed in July 2026.
