# Time-series figure contract

Scientific outputs for the manuscript and local inspection: reproducible
Matplotlib PNG, PDF and SVG files, plus their source tables. No web publication.

- **Trajectories:** How do LV and Hamming barbac frequencies compare with the
  publication over the four sampled generations? Six line-and-point panels:
  the three most abundant published pairs, each in both biological replicates.
  The cohort has only four measured generations (8, 16, 24, 40); there are no
  intervening measurements to add. Showing the publication and both methods
  provides the relevant comparison; lines connect observations without fitting
  a growth model. Condition every curve on the same 2,314 published pair IDs,
  and report each method's retained mass outside that set separately. Use blue
  circles for LV, gold triangles for Hamming and gray dashed squares for the
  publication; do not rely on color alone. Export at approximately 10 by 8 inches.
- **Composition:** What proportion of assigned molecules belongs to the largest
  LV barcode pairs at each sampled generation? Two panels of stacked bars, one
  per replicate, with four largest pooled LV pairs plus Other. Use an explicit
  four-color palette plus gray, thin outlines and a barcode-ID lookup. Counts
  normalize over all LV-assigned molecules, including pairs absent from the
  publication. The y-axis spans 0–100%; x labels identify discrete generations.
- **Clustering time:** How long do the six methods take on the two pooled
  component tables? Horizontal bars from zero, medians of three serial combined
  component workflow times, with min–max ranges. Shared FASTQ preprocessing and
  downstream trajectory reconstruction are excluded and reported separately.
  Use one blue root and neutral error bars; label seconds directly.

Titles are descriptive, with explicit units, cohort and denominator subtitles.
Use white backgrounds, dark text and restrained grid lines. This is third-party
scientific work; no OpenAI logo or decorative branding. Inspect PNG exports for
overlap, clipping, marker/legend consistency, zero treatment and readable labels.
The source tables retain sample, replicate, generation, numerator, denominators,
method, pair ID, published counts and normalized frequencies.
