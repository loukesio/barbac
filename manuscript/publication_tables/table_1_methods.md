# Table 1 methods and statistical comparisons

Table 1 compares six configurations across two independently simulated library designs and one fixed published simulation. It reports the complete Shepherd sensitivity configuration alongside unchanged results for the other methods.

## Datasets and generation

Random libraries have 20 variable bases. Anchored libraries use NNNNNAANNNNNAANNNNNTTNNNNN (20 variable bases, 26 bases total). Each design contains 60 independent libraries with 10,000 true identities and one million expected reads. Within each seed, the two designs share variable identities and parent abundances. Abundances follow the recorded scaled Johnson exponential mixture (9,989 ordinary, 10 intermediate and one high-abundance identity).

Substitutions occur with probability 0.004 per base. Repeat-dependent indel allocation uses the unscaled archived empirical homopolymer-rate table. Repeat coordinates are updated after events. When a repeat exceeds the measured support (length 13), further indel recursion stops for those reads; reads, origin labels and subsequent substitutions are retained. This boundary affected 249 of 119,992,068 reads across 16 of 120 libraries. All registered seeds and reads remain in the analysis.

The fixed Milo / Johnson deposited simulation contains 100,000 true identities and 24,996,128 reads and was used during development. Its reference counts and origin labels differ at four parents by six reads in total absolute count difference; totals and identities reconcile. Its scores are descriptive. Simulation source: https://pmc.ncbi.nlm.nih.gov/articles/PMC10276077/ ; archived analysis: https://zenodo.org/records/7411747 .

## Configurations

All methods retain an explicit distance limit of 3 and one thread. Barbac uses the preserved v14 native implementation with support ordering, merge ratio 20, configured error proxy 0.005 and design scoring disabled. LV uses the existing Poisson option; Hamming retains its rare-indel rescue. Starcode sphere uses sphere clustering; Starcode MP uses message-passing ratio 5. Bartender uses seed length 5, step 1, cutoff 1, z = 5 and forward direction.

Shepherd uses nominal length, distance 3 and its documented default Bayes-factor threshold −4. It receives the known generating substitution rate 0.004 on every simulated library; the fixed Milo reference uses automatic rate estimation. The earlier frozen configuration used threshold +4 and automatic estimation and failed on ten simulated inputs. The complete supplementary configuration was registered after observing those outcomes. Both changes apply consistently; no accuracy-based parameter search was performed. It supplies the generating substitution rate to Shepherd, not truth identities or repeat-specific indel rates. This is a post hoc sensitivity comparison, not an all-default run or a retroactive replacement for the original confirmatory protocol. Author documentation: https://github.com/Nik-Tavakolian/Shepherd .

## Metrics and timing

Exact-centroid F1 includes all designed truth identities, including those with zero reads. For each library, F1 (%) = 200TP/(2TP + FN + FP). Table entries average library-level accuracy scores without pooling reads across libraries. Read accuracy divides correctly assigned reads by all input reads; incorrect and unassigned reads are both failures. The displayed F1 spread is the sample standard deviation, not a confidence interval.

Runtime measures a fresh Python worker including tool/package startup, required input conversion and canonical centroid/member exports. Common staging, scoring and hashing are excluded. Runs were serial. Original accepted times are retained for all methods other than the new Shepherd configuration, including the previously registered 28 sleep-affected calls repeated once with identical clustering outputs. The 121 supplementary Shepherd calls were timed in a later session; lightweight diagnostic file reads overlapped early calls. No sleep/wake event overlapped those new timed calls. Comparisons involving the new Shepherd times are descriptive. Barbac package startup uses the verified lazy BAM-loading change; its clustering implementation is unchanged.

## Paired accuracy comparisons

The original primary analysis compared Barbac LV with four external configurations in each simulated design, keeping an eight-contrast family. Paired-t one-sided lower bounds use alpha 0.05/8; paired bootstrap lower bounds resample whole libraries 100,000 times at the same level. A positive difference favours Barbac LV. The six complete original contrasts with Bartender and Starcode meet both criteria. The two original Shepherd contrasts remain recorded as unavailable.

For the completed Shepherd sensitivity configuration, the same library pairing and conservative alpha 0.05/8 are retained. These two comparisons remain explicitly post hoc. The table below combines the six original contrasts and the two supplementary Shepherd contrasts for inspection, with their scope identified. Milo is excluded from population-level inference.

| Design / comparator | LV − comparator F1 (pp) | t lower | Bootstrap lower | Scope |
|---|---:|---:|---:|---|
| Random / Shepherdᵃ | +0.016172 | +0.011819 | +0.012066 | Post hoc |
| Random / Starcode sphere | +0.044684 | +0.033221 | +0.033717 | Original |
| Random / Starcode MP | +0.584555 | +0.559418 | +0.560488 | Original |
| Random / Bartender | +0.149749 | +0.138516 | +0.139150 | Original |
| Anchored / Shepherdᵃ | +0.014868 | +0.009718 | +0.009941 | Post hoc |
| Anchored / Starcode sphere | +0.042949 | +0.030841 | +0.031382 | Original |
| Anchored / Starcode MP | +0.790526 | +0.757208 | +0.759126 | Original |
| Anchored / Bartender | +0.565086 | +0.544386 | +0.545500 | Original |

All displayed summary values were reconciled with 726 execution records. Earlier independent checks reproduced every new Shepherd centroid score, three complete read mappings and the supplementary paired-t calculations in base R. Source fingerprints and the export checks accompany the table in `table_1_provenance.json`.
