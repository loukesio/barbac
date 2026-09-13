# Completing the Shepherd comparison

This is a separate, post hoc sensitivity analysis. The original publication campaign, its ten failed Shepherd calls, its statistical analysis and its PDFs remain unchanged in `benchmark/publication_final` in the original publication worktree. This directory must not be presented as the original frozen configuration.

## Outcome

All **121/121** new Shepherd calls completed successfully. Barbac LV retains the highest mean centroid F1 and lowest mean FP in both simulated designs, and the highest F1 / lowest FP on fixed Milo. Barbac Hamming is fastest. Starcode MP retains fewer FN in the simulated designs; no all-metric dominance is claimed.

The supplementary LV-minus-Shepherd mean F1 differences are +0.016172 percentage points (random) and +0.014868 (anchored). Both paired-t and bootstrap lower bounds at alpha 0.05/8 are above zero. These remain post hoc sensitivity findings. Milo's tiny numerical difference is descriptive, not statistical superiority.

[Complete results](RESULTS.md) · [Updated PDF](benchmark_table.pdf) · [Independent verification](independent_validation.json). Every new centroid score was independently recomputed; complete read mappings were checked on three selected inputs, and paired-t calculations were reproduced in base R. No sleep/wake event overlapped a timed call. The full one-page PDF was visually inspected.

## The problem and correction

The original benchmark explicitly used `-eps 3 -bft 4` and automatic substitution-rate estimation. Shepherd failed on ten of 120 simulated inputs. Reproducing its installed estimator shows estimates of approximately 18–66%, exceeding its 10% rejection threshold, although the generating substitution rate was 0.4%. Repeated contributions from more abundant one-substitution neighbours account for approximately 98–100% of the estimator numerator in these failed inputs. The input files themselves were not rejected as malformed.

[Shepherd's author documentation](https://github.com/Nik-Tavakolian/Shepherd) recommends supplying `-e` when automatic estimates are unreliable, defines this as the per-base substitution rate, and documents `-bft -4` as the default. The original benchmark's explicit `+4` is preserved and disclosed. No provenance establishing that it was the author's recommended value was found in the local archived simulation notebook.

The supplementary configuration uses `-eps 3 -bft -4 -e 0.004` on **every** registered simulated input. The rate comes directly from the simulator; it is not optimized against accuracy outcomes. This supplies Shepherd with the generating substitution rate, which is an advantage unavailable to an automatic estimator on unknown real data. Barbac keeps its originally configured 0.005 error proxy; these parameters do not represent identical statistical models. Fixed Milo retains automatic rate estimation, with `-eps 3 -bft -4`.

This corrects the documented threshold and handles the estimator problem consistently across all inputs. It does not isolate the individual effects of the two parameter changes. Distance remains explicitly 3; this is not an all-default Shepherd run.

## Execution and interpretation

`protocol.json`, `worker.py`, and `run.py` were fingerprinted before the first new invocation. One new Shepherd invocation is scheduled for each of the original 121 inputs, in their original dataset order. Every original outcome for Barbac, Bartender and Starcode is reused; none of those tools is rerun. No inputs, labels, simulator settings, final seeds, native libraries or original results are changed. Failures remain visible and cannot be silently replaced on resume.

Fresh-process timings use the original workflow boundary, but Shepherd is measured in a later session. They are reported descriptively; the supplementary timing comparison is not a randomized paired-session experiment. A task-owned keep-awake assertion spans execution and reporting. Some lightweight diagnostic file reads occurred during the early calls; this is another reason not to interpret small timing differences as controlled speed effects.

The accuracy comparison uses all 60 paired libraries per design when complete. It reports mean LV-minus-Shepherd F1 differences, a one-sided paired-t lower bound at alpha 0.05/8, and a paired 100,000-resample bootstrap lower bound at the same level. Keeping the original eight-comparison correction is conservative for these two supplementary contrasts. It does **not** undo the post hoc choice of competitor configuration or turn the supplementary analysis into the original confirmatory test. Milo remains a single development reference without a population-level significance claim.

Local commands, from this worktree:

```sh
python3 benchmark/shepherd_completion/run.py freeze
python3 benchmark/shepherd_completion/run.py run
python3 benchmark/shepherd_completion/diagnose.py
python3 benchmark/shepherd_completion/analyze.py
python3 benchmark/shepherd_completion/report.py
```

## Barbac recall diagnosis

The separate recall audit uses only the five **development** references, never the consumed final test seeds. On random and anchored mixed development inputs, 102 of 125 missing identities had zero reads. All remaining missed identities had only one or two reads. On Milo, 409 of 471 misses had zero reads; almost all remaining cases also had one or two reads.

A true sequence with no reads contributes no identifying evidence. A pair of one-count sequences differing by one substitution can be equally consistent with either sequence being the true parent. Retaining both trades a potential FN for an FP; choosing the correct one requires additional information or an explicit prior. Changing lexicographic ties after inspecting truth would be tuning to the answer.

One Milo development family has three distinct erroneous observations and an unobserved true sequence that a strict majority consensus could reconstruct. That is a legitimate, narrowly bounded candidate for future evaluation. It is not evidence that a general consensus or more permissive merge rule improves accuracy: it must also be tested on genuine-neighbour controls and separate development inputs. No such candidate replaces v14 in this comparison, and no new FASTQ or quality-score dependency is introduced.
