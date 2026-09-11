# LV optimization and an expected-error indel exception

The reference experiment reduces LV with support ordering from **192.0s to
43.4s end to end**, and improves FN/FP from **471/85 to 471/83**. Incorrectly
assigned reads fall from **231 to 131**, with all 24,996,128 reads retained.
The new indel option is experimental and opt-in. The exact search optimization
applies to LV automatically.

```r
result <- super_cluster2(
  input,
  method = "lv",
  tie_break = "support",
  error_rate = 0.005,
  indel_model = "poisson"
)
```

`error_rate` remains a configured approximate per-base rate. No truth labels
are used to fit it. `indel_model = "none"` remains the default, retaining the
previous merge decisions while using the faster search.

## Reference result and ablation

| LV configuration | FN | FP | Incorrect reads | Core seconds | Process seconds |
|---|---:|---:|---:|---:|---:|
| Previous v12, support | 471 | 85 | 231 | 187.442 | 192.019 |
| v13 search only, support | 471 | 85 | 231 | 39.997 | 44.320 |
| v13 search + Poisson, support | 471 | 83 | 131 | 38.980 | 43.399 |
| v13 search only, sequence | 472 | 86 | 234 | 27.852 | 31.900 |
| v13 search + Poisson, sequence | 472 | 84 | 134 | 28.539 | 32.567 |

The reference inputs and the v12 measurements are from the immediately
preceding [reference comparison](../reference_comparison/README.md), unchanged
and hash-verified. The two v13 search-only runs reproduce **all centroids,
counts, and sequence-to-centroid assignments** of v12. Model-enabled runs are
scored independently against the supplied labels. Old benchmark results are
preserved, rather than overwritten.

Core time includes CSV reading, ordering, clustering, and result construction.
Process time also includes R/package startup and both centroid/member exports.
These are single serial runs on the same development machine, with lightweight
interactive work during some measurements; small timing differences are not
established speed effects. No competitor rerun was necessary because its input
and configuration did not change. The previous report retains those results:
Shepherd FN 470 / FP 85, 217 wrong or unassigned reads, 111.0s workflow time;
Hamming with support FN 469 / FP 87, 228 wrong reads, 18.4s process time.
LV with the Poisson option now narrowly leads their centroid F1 on this one
reference dataset, but Hamming still has fewer false negatives and is faster.

The truth table includes 409 zero-read barcodes, 439 exact true sequences absent
from the input, and four inconsistent parent totals. Their interpretation and
the label-sensitivity calculation remain as documented in the reference audit.
Neither the new option nor the scoring rewrites those source files.

## Why the search is faster and still exact

The one-edit index is queried first without abundance filtering. Once an
absorbing parent is available, every still-unseen candidate must be at least
two edits away. Its score is bounded above by its abundance, a distance of two,
and the query length. The likelihood decreases with distance and comparison
length, so this is conservative even when actual candidate lengths differ.

Integer bisection finds the smallest parent count that can reach the current
best score. Each partition posting list is in descending abundance under the
public API, so scanning stops at the first count below that threshold. All
possible equal-score ties remain eligible. Unsorted native inputs filter each
posting without early stopping. Unsupported strings still use the complete
fallback candidate population, subject to the same proven score bound. If no
absorbing parent exists, the wider query is unfiltered.

This changes the amount of searching, not the scoring or merging criteria.
With support ordering, full edit-distance verifications decrease from
726,901,780 in v12 to 97,393,407 with the new model-enabled search. Cached
centroid bit masks or a different distance kernel are not part of this change.

## Why the 20x guard rejected the two deletion variants

Two 19-base deletion variants of 20-base parents occurred 68 and 24 times;
their parents occurred 1,114 and 415 times. Ratios of 16.4x and 17.3x fail the
20x guard. Each deletion is inside a long repeated-C run: deleting any one of
several identical bases produces exactly the same variant.

For a single deletion with `r` equivalent positions, the model uses:

```
expected_variant_count = observed_parent_count * r * e / (1 - e)
```

For a single insertion, it divides that expectation by four to represent the
specific inserted base under a uniform-base assumption. The equivalent gap
positions are counted directly from the strings, without an alignment matrix.
Only repeated-base single indels qualify. Ordinary substitutions and other
blocked merges retain their existing rules.

For an otherwise blocked merge, let `X ~ Poisson(expected_variant_count)`.
The exception permits merging when `P(X >= observed_variant_count) >= 0.01`,
provided the parent is more abundant. Thus an observed count well above the
error expectation is still protected. The supplied per-base rate is an
upper-bound proxy for deletion error, not a learned platform-specific rate.
The tail cutoff is fixed at 0.01 for these experiments; it is not a posterior
probability, false-discovery guarantee, or proof that a variant is erroneous.

This exception removes the two false deletion clusters and allows their error
families to follow the correct parent. It does not globally lower the merge
ratio, choose centroids using true labels, infer a mandatory barcode length,
or change the parent likelihood function.

## Validation and limitations

`results_all.csv` records the reference ablation and four original conditions
with seeds 42, 43, and 44. Every run conserves member counts. Search-only
centroid/count outputs match v12 on all twelve condition/seed combinations.
The indel option changes no FN and removes five FP in total across those cases.

`holdout_results.csv` adds a new seed (20260909) for the original four
conditions and a separate homopolymer-enriched stress simulation. The stress
test is intentionally favorable to the proposed error mechanism; it is not
representative of every barcode library. Its independently generated truth
and counts test the mechanism beyond the two reference families.

On the fresh seed, all four conditions retain the same FN; one FP is removed
from random low-indel data, and the other three cases have unchanged FN/FP.
In the 1,000-truth, one-million-read homopolymer stress test, FN remains zero
and FP falls from 1,042 to 913. F1 improves from 65.75% to 68.66%, but remains
poor: this single-indel exception does not solve the many remaining error
clusters. The exact generators, settings, and hashes are saved in the holdout
manifest. No parameters were adjusted after viewing the holdout results.

The reference assignment audit independently confirms that exactly nine
observed sequences (100 reads) change parents, and all nine changes agree with
their supplied parent labels. Errors among the 1,001 off-length reads drop
from 100 to zero. See `reference_changed_assignments.csv` and
`reference_indel_audit.json`.

Meaningful tests cover indexed/full-scan equivalence, equal-score candidates,
unsorted native inputs, both model settings, error-rate sensitivity, insertion
versus deletion, repeats at sequence boundaries, nonrepeated gaps, protection
of abundant length variants, API validation, and file/directory forwarding.

A genuine length variant with the same sequence and count pattern as an error
cannot be distinguished from counts alone. In particular, this option can
merge real variable-length barcodes. Overdispersion, early PCR errors,
nonuniform inserted bases, and a platform-specific indel rate different from
`error_rate` can invalidate the Poisson approximation. That is why the option
remains off by default and requires independent controls for new libraries.
Faster exact search does not remove these statistical ambiguities.

## Reproduce

Build v12 (revision `f2b6468`, native v12) and this source into separate R
libraries. The experiment uses the current package's native v13 build marker.

```sh
R CMD INSTALL --library=/private/tmp/barbac-lv-v13 .
python3 benchmark/lv_optimization/run_experiment.py \
  --library /private/tmp/barbac-lv-v13 \
  --baseline-library /private/tmp/barbac-exact-final \
  --source /path/to/barbac-benchmark
python3 benchmark/lv_optimization/run_holdout.py
```

The main runner reuses hash-verified local simulations under
`benchmark/four_condition_comparison/generated/revision-accuracy/`; generate
these with the existing revision-comparison runner if absent. Large inputs,
memberships, and checkpoints stay in ignored `generated/lv-v13*` directories.
`run_holdout.py` generates its own inputs and uses the temporary library path
shown above. Manifests record source, installed-library, and input hashes.
