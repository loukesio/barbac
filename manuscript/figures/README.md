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
