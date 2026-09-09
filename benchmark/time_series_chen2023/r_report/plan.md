# R report plan

Audience: technical; question: what do barbac's extraction, clustering and
time-series diagnostics show for the completed Chen hBFA1/YPD assay?
Delivery: one portable Quarto HTML report using the existing R package.
The user's explicit R/package-function requirement overrides the analytics
skills' generic JavaScript artifact packaging. No Sites deployment is requested.

## Function inventory and data contracts

| Existing function | Report use | Input / unit |
|---|---|---|
| `summarise_bam_stats()` | Recount mapped/unmapped alignments | Eight completed merged BAMs |
| `plot_bam_stats()` | Mapping overview | Returned counts, merged alignments |
| `barbac_xtr.stats()` | Sample/component length, abundance, entropy panels | Shared extracted distinct sequences, UMI molecule counts |
| `cluster_stats()` | Component and paired-lineage abundance summaries | Positive centroid counts; zeros excluded from cluster counting |
| `barbac_ts_area()` | Interactive LV/Hamming composition per replicate | Top four pooled LV pairs plus Other; complete grid, actual zeros |
| `theme_barbac()` | Consistent report plots | Existing theme; explicit palettes |

Scope: generations 8, 16, 24, 40; R1/R2; all eight complete libraries. Cluster
assignments are the saved first repeat. No distance changes or new clustering.
Method agreement is recalculated in R against saved published counts and checked
against the previous results. Publication counts are not truth labels.

## Reading order / technical specification mapping

1. Title and technical summary: retention, agreement, timing and its scope.
2. Mapping and extraction evidence: counts and package QC panels; cohort and
   metric definitions appear before their plots.
3. Clustering evidence: component / paired-lineage summaries and interpretation.
4. Composition and barcode exploration: four measured generations, replicate
   context, no fitted dynamics; searchable complete paired count tables.
5. Publication comparison: all six methods, agreement, coverage and timing;
   unmatched mass remains visible. No FN/FP from experimental reference.
6. Methods, checks, limitations: denominators, shared inputs, tool assumptions.
7. Next steps and open questions: investigate abundant unpublished pairs;
   distinguish publication filtering from clustering before accuracy claims.

## Chart contract

| Section | Question / form | Rows and sufficiency | Renderer and palette |
|---|---|---|---|
| Mapping | Which samples lose mapped reads? Stacked horizontal bars | 8 samples × 2 states | `plot_bam_stats`, blue / neutral, labels |
| Extraction | Length, abundance, sequence complexity? Three histograms and table | Every distinct sequence, 16 sample/component combinations | `barbac_xtr.stats`, blue; static package panels embedded with sample controls |
| Clustering | How concentrated are molecule counts? Ranked abundance shares plus exact table | 4 pooled component configurations, 16 paired sample summaries | `cluster_stats`, neutral DT tables; exact lookup takes priority |
| Composition | How do major pairs change? Stacked areas | 4 generations × 2 replicates × 2 modes × 5 groups | `barbac_ts_area(interactive='ggiraph')`, fixed blue/gold/olive/pink/gray mapping |
| Agreement | Which methods reproduce counts? Sample scatter/point evidence and exact summary | 48 method/sample comparisons | R ggplot/plotly, neutral references; agreement axis explicitly focused near 1 |
| Speed | What is the comparison cost? Horizontal bars with range | 6 methods × 3 repeats | R ggplot/plotly, blue, zero baseline, median and min/max |

Four generations are the complete assay, not a sparse query. Replicate panels
provide the meaningful additional breakdown; lines/areas interpolate visually,
without inventing observations. Composition is normalized over all assigned
molecules; Other preserves the denominator. No pseudocounts or abundance filter.

Validation: compare numeric outputs with existing results, test the extraction
length-bin defect and zero-aware report adapters, run relevant package tests,
inspect exported plots and the complete HTML, exercise interactive controls and
downloads when browser automation is available. Preserve inputs and hashes.
