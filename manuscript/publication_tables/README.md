# Publication Table 1

[Print-ready comparison PDF](table_1_benchmark_comparison.pdf) · [Editable Word table](table_1_benchmark_comparison.docx) · [CSV](table_1_benchmark_comparison.csv)

All six method configurations are shown for each of the three benchmarks. The table includes FN, FP, exact-centroid F1, read-assignment accuracy and workflow runtime. Bold and pale orange shading indicate the numerical leader within each benchmark and metric.

The [companion methods PDF](table_1_methods.pdf) contains full settings, simulation provenance, timing details and the original versus supplementary statistical comparisons. [Unrounded data](table_1_unrounded_data.csv) and [source receipts](table_1_provenance.json) support reuse.

Rebuild from this worktree with `python3 benchmark/shepherd_completion/publication_table.py`. This only reads saved results and renders files; it runs no clustering or simulation. The PDF is the visually verified print layout; the Word table remains editable and its serialized values are checked against the same source.
