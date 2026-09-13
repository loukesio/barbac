# Publication manuscript

- [Current manuscript source](barbac_manuscript.md)
- [Editable manuscript](barbac_manuscript_publication.docx)
- [Manuscript PDF](barbac_manuscript_publication.pdf)
- [Table 1 PDF](publication_tables/table_1_benchmark_comparison.pdf) and
  [editable Word table](publication_tables/table_1_benchmark_comparison.docx)
- [Benchmark methods and paired comparisons](publication_tables/table_1_methods.pdf)
- [Experimental-reference supplement](supplementary_chen2023.md)
- [Corrected workflow figure](figures/workflow.pdf)
- [Remaining author submission details](submission_readiness.md)
- [Cover letter draft](cover_letter_draft.md)

The current comparison uses the preserved v14 clustering engine. Original
manuscripts, figures and benchmark records remain in Git history and archived
analysis folders. `barbac_manuscript.docx` is the earlier editable source;
use the explicitly named publication DOCX above for the revised paper.

After inserting completed runtime values in the Markdown, rebuild the revised
Word and PDF files with `python3 manuscript/build_publication.py`. This uses
Quarto/Pandoc, python-docx, XeLaTeX and Times New Roman. Word receives the original
editable, orange-highlighted Table 1 rather than a screenshot. The manuscript
PDF embeds the comparison preview; the standalone Table 1 PDF remains the
vector, submission-quality table.
