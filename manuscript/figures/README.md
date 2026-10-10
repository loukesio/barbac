# Workflow figure

`workflow.html` preserves the user's wide six-stage design with corrected
implementation descriptions. The original is retained locally before editing.
`workflow.svg` presents the same corrected flow at 178 mm width with readable
publication typography; PDF and PNG are exported from this vector source.

Corrections: R1-only means single-end input; paired reads must overlap; the
mapping wrapper ends at BAM/QC; flank extraction preserves observed indels;
length filtering is optional; clustering uses native abundance-aware centroid
assignment, not igraph connections; independent count tables enter without
FASTQ; trajectory bands represent frequencies; table exports are CSV.

The native SVG text and arrows remain editable. The wide design is intended for
screen viewing; reducing all six columns to a journal-width panel would make
its labels too small. The compact version retains charcoal headers, coral
function labels and restrained grey connectors.

## Pipeline figure and barcode-design schematic (10 October 2026)

`render_pipeline_figures.py` writes `pipeline.svg` and `barcode_designs.svg`
at 178 mm width and renders PNG (400 dpi) and PDF through headless Chrome.
Colours are the bundled LTC palette `alger`. The pipeline figure corrects the
root `Barbac Pipeline.dc.html` draft: clustering is native abundance-aware
centroid assignment in `barbac_cluster()` (no igraph or pairwise distance
matrix), extraction preserves observed indels, count tables enter clustering
without FASTQ, mapping statistics are CSV, and function names match the
exported API. Pass `--design` strings such as `ATTA-NNNNNN-ATCCC` to draw
other barcode designs.
