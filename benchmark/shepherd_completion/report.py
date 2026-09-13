"""Render the complete comparison in the established publication-table format."""
import json
from pathlib import Path
import re
import subprocess
import sys
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle

HERE = Path(__file__).resolve().parent
FINAL = Path('/Users/theodosiou/Documents/Projects/test_barbac/.codex/publication-final-2026-09-12/source/benchmark/publication_final')
sys.path.insert(0,str(FINAL))
from common import load,save,sha
from report import table_style, TITLES, BLUE, PALE, GOLD

NAMES = {'hamming':'barbac Hamming','lv':'barbac LV + Poisson','shepherd_documented':'Shepherd*',
         'starcode_sphere':'Starcode sphere','starcode_mp':'Starcode MP','bartender':'Bartender'}


def footer(canvas,doc):
    canvas.setFont('Helvetica',7)
    canvas.setFillColor(BLUE)
    canvas.drawString(16*mm,8*mm,'barbac | supplementary Shepherd comparison | 13 September 2026')
    canvas.drawRightString(194*mm,8*mm,str(doc.page))


def main():
    data=load(HERE/'summary.json')
    contrasts=load(HERE/'accuracy_contrasts.json')
    validation=load(HERE/'analysis_validation.json')
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='Note',fontName='Helvetica',fontSize=8,leading=10.5,spaceAfter=5))
    styles.add(ParagraphStyle(name='Tiny',fontName='Helvetica',fontSize=7.2,leading=9,spaceAfter=4))
    story=[Paragraph('Barcode clustering benchmark',styles['Title']),
        Paragraph('Preserved v14 + faster startup | completed Shepherd sensitivity comparison',styles['Note']),
        Paragraph('60 independent libraries per simulated design. FN / FP: mean counts. F1: mean % (SD). Time: median workflow seconds. Milo: one fixed reference.',styles['Note']),Spacer(1,5)]
    rows=[['Method','FN (lower)','FP (lower)','F1 % (higher)','Time s (lower)']]
    commands=table_style()
    displayed=[]
    for condition,title in TITLES.items():
        i=len(rows)
        rows.append([title,'','','',''])
        commands += [('SPAN',(0,i),(-1,i)),('BACKGROUND',(0,i),(-1,i),PALE),
            ('FONTNAME',(0,i),(-1,i),'Helvetica-Bold'),('TEXTCOLOR',(0,i),(-1,i),BLUE)]
        selected=[r for r in data if r['condition']==condition]
        full=[r for r in selected if r['complete_population']]
        metrics=['fn','fp','f1_percent','workflow_seconds']
        best={m:(max if m=='f1_percent' else min)(r[m] for r in full) for m in metrics}
        for r in selected:
            fixed=condition=='milos'
            if r['complete_population']:
                digits=0 if fixed else 2
                f1=f"{r['f1_percent']:.5f}"+(f" ({r['f1_percent_sd']:.4f})" if not fixed else '')
                values=[NAMES[r['method']],f"{r['fn']:,.{digits}f}",f"{r['fp']:,.{digits}f}",f1,f"{r['workflow_seconds']:.2f}"]
            else:
                values=[NAMES[r['method']]+f" [{r['n_successful']}/{r['planned_n']}]",'n/a','n/a','n/a','n/a']
            i=len(rows);rows.append(values);displayed.append(values)
            if r['method'] in ['hamming','lv']:
                commands.append(('TEXTCOLOR',(0,i),(0,i),colors.HexColor('#255D8C')))
            for col,metric in enumerate(metrics,1):
                if r[metric] is not None and abs(r[metric]-best[metric])<1e-10:
                    commands += [('FONTNAME',(col,i),(col,i),'Helvetica-Bold'),('BACKGROUND',(col,i),(col,i),GOLD)]
    commands.append(('LINEBELOW',(0,len(rows)-1),(-1,len(rows)-1),.7,BLUE))
    tab=Table(rows,colWidths=[52*mm,25*mm,25*mm,47*mm,29*mm]);tab.setStyle(TableStyle(commands))
    story += [tab,Spacer(1,10)]
    notes=[
      'Bold/shading marks numerical leaders, not statistical significance. All other methods retain their frozen outputs and accepted timing observations. No Barbac code or simulation data changed.',
      '*Shepherd supplementary configuration: distance 3, author-documented Bayes threshold -4; known generating substitution rate 0.004 on every new simulation. Milo retains automatic rate estimation. The original configuration used threshold +4 and automatic rates, failing on 10/120 new inputs. Its results and PDF remain preserved separately.',
      'This configuration was registered after the original outcomes were available. It is an explicitly post hoc sensitivity comparison, not the original confirmatory result. Supplying Shepherd the generating substitution rate favours a correctly specified model; no truth identities or repeat-specific indel rates are supplied. This is not an all-default run.',
      'New simulations: 10,000 true identities and one million expected reads per library; paired designs share variable identities and parent counts. Substitutions: 0.4% per base. Indels: unscaled archived homopolymer-specific rates. For repeats beyond calibrated length 13, further indel recursion stops while reads and subsequent substitutions are retained (249 affected reads among 119,992,068; 16/120 libraries).',
      'FN counts true identities absent from inferred centroids; FP counts inferred identities absent from truth. F1 includes zero-read true identities. Positive-read F1, read assignment accuracy and abundance metrics are in the companion results. The fixed Milo simulation contains 100,000 identities and 24,996,128 reads and was used during development; it cannot establish population-level significance.',
      'Timing includes fresh-worker/tool startup, required conversion and canonical exports; common staging and scoring are excluded. Shepherd times come from a later session and are descriptive, not a randomized paired-session speed comparison. Lightweight diagnostic file reads overlapped early calls. Original methods retain the previously documented sleep-repair observations.',
      'Other settings: distance 3, one thread. Barbac support ordering, ratio 20, configured proxy 0.005, design scoring off; LV Poisson on, Hamming rare-indel rescue retained. Starcode MP ratio 5. Bartender seed length 5, step 1, z=5, cutoff 1. Full protocol, preserved failures and supplementary paired uncertainty: README.md and RESULTS.md.'
    ]
    story += [Paragraph(escape(s),styles['Tiny']) for s in notes]
    pdf=HERE/'benchmark_table.pdf'
    SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=16*mm,rightMargin=16*mm,topMargin=11*mm,
                      bottomMargin=14*mm,title='Barbac supplementary complete comparison').build(story,onFirstPage=footer,onLaterPages=footer)
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
    assert re.search(r'^Pages:\s+1$',info,re.M),info
    extracted=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True)
    actual=[re.sub(r'\s+',' ',line).strip() for line in extracted.splitlines() if line.strip().startswith(tuple(NAMES.values()))]
    assert actual==[' '.join(r) for r in displayed],(actual,displayed)
    subprocess.run(['pdftoppm','-scale-to','1800','-png','-singlefile',str(pdf),str(HERE/'benchmark_table')],check=True)
    lines=['# Complete comparison with Shepherd sensitivity configuration','',
        'The original confirmatory benchmark is unchanged. Shepherd uses the separately registered post hoc configuration explained in [README.md](README.md). All other scores and accepted times are reused.','',
        '[PDF table](benchmark_table.pdf) · [CSV table](benchmark_table.csv) · [execution receipts](results.json)','']
    for condition,title in TITLES.items():
        lines += ['## '+title,'','| Method | FN | FP | F1 % | Workflow s | Correct reads % | Wrong reads | Unassigned |',
                  '|---|---:|---:|---:|---:|---:|---:|---:|']
        for r in [r for r in data if r['condition']==condition]:
            if not r['complete_population']:
                lines.append(f"| {NAMES[r['method']]} ({r['n_successful']}/{r['planned_n']}) | — | — | — | — | — | — | — |")
            else:
                lines.append(f"| {NAMES[r['method']]} | {r['fn']:.2f} | {r['fp']:.2f} | {r['f1_percent']:.5f} | {r['workflow_seconds']:.2f} | {r['read_assignment_accuracy_percent']:.6f} | {r['incorrect_reads']:.2f} | {r['unassigned_reads']:.2f} |")
        lines += ['','Simulations: means over all 60 libraries; time is median. Milo: one observation. Shepherd was timed in a later session.','']
    lines += ['## Supplementary paired accuracy','',
              '| Design | LV − Shepherd F1, percentage points | t lower bound | Bootstrap lower | Both > 0 |',
              '|---|---:|---:|---:|---|']
    for c in contrasts:
        if c['status']=='complete':
            lines.append(f"| {c['condition']} | {c['mean_difference']:+.6f} | {c['simultaneous_one_sided_lower']:+.6f} | {c['bootstrap_lower']:+.6f} | {'Yes' if c['superiority_supported'] else 'No'} |")
        else:
            lines.append(f"| {c['condition']} | Incomplete | — | — | — |")
    lines += ['','One-sided alpha 0.05/8 lower bounds; 100,000 whole-library bootstrap resamples. These intervals are supplementary because the Shepherd configuration was selected after the original benchmark. They do not establish universal superiority or retroactively complete the original confirmatory contrasts.',
              '',f"Shepherd completed {validation['complete_shepherd_calls']}/121 calls. Original +4/automatic results, including all ten failures, remain in the original publication worktree. No final test seed was used for Barbac tuning.",'',
              '## Development recall','',
              '| Input | Total FN | Zero-read truth | True sequence unobserved despite positive reads | Observed identity merged/unassigned |',
              '|---|---:|---:|---:|---:|']
    for r in load(HERE/'development_recall_audit.json'):
        c=r['categories']
        lines.append(f"| {r['condition']} | {r['fn']} | {c.get('zero_read_truth',0)} | {c.get('exact_sequence_unobserved',0)} | {c.get('observed_identity_merged_or_unassigned',0)} |")
    lines += ['','These are development references, not final-test decompositions. Most positive misses are one-read errors or tied pairs. More permissive merging cannot identify sequences with no evidence and can erase genuine nearby identities. The previous rejected paired-indel candidate remains excluded.','']
    (HERE/'RESULTS.md').write_text('\n'.join(lines))
    save(HERE/'report_validation.json',dict(pdf_pages=1,displayed_rows=len(displayed),
        all_pdf_rows_equal_source=True,pdf_sha256=sha(pdf),rendered_image=str(HERE/'benchmark_table.png'),
        visual_inspection_pending=True))
    print(pdf)


if __name__=='__main__':
    main()
