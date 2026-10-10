#!/usr/bin/env python3
"""Assemble the supplementary information PDF from its reviewed sources.

Sources (all tracked):
  publication_tables/table_1_methods.md            -> Supplementary Methods, Table S1
  supplementary_chen2023.md                        -> Table S2, Figure S1
  publication_tables/table_s3_development_reference.csv -> Table S3

Writes manuscript/supplementary_information.md and .pdf. Uses Quarto's pandoc
and XeLaTeX with Times New Roman, like build_publication.py.
"""
import csv
import pathlib
import re
import subprocess
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
METHODS = HERE / 'publication_tables' / 'table_1_methods.md'
CHEN = HERE / 'supplementary_chen2023.md'
S3_CSV = HERE / 'publication_tables' / 'table_s3_development_reference.csv'

CONDITIONS = {
    'random_substitutions': 'Random N20: substitutions only',
    'random_mixed': 'Random N20: substitutions + repeat-dependent indels',
    'anchored_substitutions': 'Anchored N20 + AA/TT: substitutions only',
    'anchored_mixed': 'Anchored N20 + AA/TT: substitutions + repeat-dependent indels',
}
METHODS_LABEL = {
    'hamming': 'barbac Hamming', 'lv': 'barbac LV + Poisson', 'shepherd': 'Shepherd',
    'starcode_sphere': 'Starcode sphere', 'starcode_mp': 'Starcode MP', 'bartender': 'Bartender',
}


def table_s3():
    rows = list(csv.DictReader(S3_CSV.open()))
    build_ids = sorted({r['build_id'] for r in rows if r['method'] in ('hamming', 'lv')})
    assert build_ids == ['barbac-2026-09-10-eligible-parent-v14'], build_ids
    out = ['## Table S3. Development reference: substitution-only and indel-containing designs\n',
           '**Table S3. Single-library development reference for the random and anchored designs, '
           'with and without repeat-dependent indels.** One library per condition (10,000 true identities, '
           'one million expected reads, 0.4% substitutions per base; indel conditions use the unscaled '
           'archived homopolymer-rate table). These libraries were used during development and are '
           'distinct from the 60 final seeds per design in Table 1. Barbac uses the same native v14 '
           'implementation and settings as Table 1. Shepherd uses the original frozen configuration '
           '(Bayes-factor threshold +4, automatic substitution-rate estimation), which completed on all '
           'four inputs; Table 1 reports its completed sensitivity configuration. FN, FP and F1 are '
           'exact-centroid metrics including zero-read true identities; read accuracy counts reads '
           'assigned to their true parent. Runtime is one fresh-worker observation including startup, '
           'conversion and exports, recorded before the 0.2.1 startup change, so these times are not '
           'comparable with Table 1. Bold marks the best observed value per metric within a condition. '
           'Source: `publication_tables/table_s3_development_reference.csv`, with per-cell receipt signatures.\n']
    for cond, title in CONDITIONS.items():
        sub = [r for r in rows if r['condition'] == cond]
        assert len(sub) == 6, cond
        best = {
            'fn': min(float(r['fn']) for r in sub), 'fp': min(float(r['fp']) for r in sub),
            'f1': max(float(r['f1_percent']) for r in sub),
            'acc': max(float(r['read_assignment_accuracy_percent']) for r in sub),
            't': min(float(r['workflow_seconds']) for r in sub),
        }
        out.append(f'### {title}\n')
        out.append('| Method | FN | FP | F1 (%) | Read accuracy (%) | Workflow s |')
        out.append('|---|---:|---:|---:|---:|---:|')
        for key in METHODS_LABEL:
            r = next(x for x in sub if x['method'] == key)
            def cell(val, fmt, is_best):
                s = fmt.format(val)
                return f'**{s}**' if is_best else s
            out.append('| {} | {} | {} | {} | {} | {} |'.format(
                METHODS_LABEL[key],
                cell(float(r['fn']), '{:.0f}', float(r['fn']) == best['fn']),
                cell(float(r['fp']), '{:.0f}', float(r['fp']) == best['fp']),
                cell(float(r['f1_percent']), '{:.5f}', float(r['f1_percent']) == best['f1']),
                cell(float(r['read_assignment_accuracy_percent']), '{:.6f}',
                     float(r['read_assignment_accuracy_percent']) == best['acc']),
                cell(float(r['workflow_seconds']), '{:.2f}', float(r['workflow_seconds']) == best['t'])))
        out.append('')
    out.append('On random substitution-only barcodes Starcode sphere has the highest exact-centroid F1 '
               'in this single library, while both barbac configurations have the highest read-assignment '
               'accuracy; the 18 positive-read identities missed by Starcode sphere account for 383 '
               'misassigned reads, compared with 36 for barbac. On anchored substitution-only barcodes '
               'barbac leads F1 and read accuracy. No statistical test is applied to single-library '
               'development results.\n')
    return '\n'.join(out)


def main():
    methods = METHODS.read_text()
    methods = methods.replace('# Table 1 methods and statistical comparisons',
                              '## Supplementary Methods: Table 1 settings, provenance and statistical comparisons')
    methods = re.sub(r'^## ', '### ', methods, flags=re.M)
    methods = methods.replace('### Supplementary Methods', '## Supplementary Methods')
    methods = methods.replace('### Paired accuracy comparisons',
                              '### Paired accuracy comparisons (Table S1)')
    methods = methods.replace('The table below combines', '**Table S1.** The table below combines')
    chen = CHEN.read_text().replace('# Supplementary application: Chen 2023 time series',
                                    '## Supplementary application: Chen et al. (2023) time series (Table S2, Figure S1)')
    chen = chen.replace('![](media/time_series_trajectories.png)', '![](media/time_series_trajectories.png){width=160mm}')
    text = ('# Supplementary information\n\n'
            '**barbac: fast DNA barcode clustering from reads to lineage trajectories**\n\n'
            'Loukas Theodosiou, Andrew D. Farr and Paul B. Rainey\n\n'
            + methods + '\n\n' + chen + '\n\n' + table_s3())
    (HERE / 'supplementary_information.md').write_text(text)
    with tempfile.TemporaryDirectory(prefix='barbac-supp-') as folder:
        md = pathlib.Path(folder) / 'supp.md'
        md.write_text(text)
        subprocess.run(['quarto', 'pandoc', str(md), '--from=markdown-implicit_figures',
                        '--resource-path=' + str(HERE), '--pdf-engine=xelatex',
                        '-V', 'mainfont=Times New Roman', '-V', 'geometry:margin=20mm',
                        '-V', 'fontsize=10pt', '-V', 'colorlinks=true',
                        '-o', str(HERE / 'supplementary_information.pdf')], cwd=HERE, check=True)
    print('Rendered supplementary_information.md and .pdf')


if __name__ == '__main__':
    main()
