"""Publication tables in the established A4 format, from verified final results."""
from pathlib import Path
import re
import subprocess
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from common import HERE, load, save, sha
from run import check_freeze

NAMES = {'hamming':'barbac Hamming', 'lv':'barbac LV + Poisson', 'shepherd':'Shepherd',
         'starcode_sphere':'Starcode sphere', 'starcode_mp':'Starcode MP', 'bartender':'Bartender'}
TITLES = {'random_mixed':'Random N20: substitutions + repeat-dependent indels',
          'anchored_mixed':'Anchored N20 + AA/TT: substitutions + repeat-dependent indels',
          'milos':'Milo / Johnson published reference simulation'}
BLUE, GOLD, PALE = (colors.HexColor(s) for s in ['#2C4A63','#F3E7CB','#EDF1F4'])


def table_style():
    return [('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,1),(-1,-1),'Helvetica'),
        ('FONTSIZE',(0,0),(-1,-1),8),('LEADING',(0,0),(-1,-1),10),
        ('ALIGN',(1,0),(-1,-1),'RIGHT'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),
        ('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),
        ('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),
        ('LINEBELOW',(0,0),(-1,0),.7,BLUE)]


def display(row):
    fixed = row['condition'] == 'milos'
    if not row.get('complete_population', True):
        return [NAMES[row['method']]+f" [{row['n_successful']}/{row['planned_n']}]", 'n/a', 'n/a', 'n/a', 'n/a']
    return [NAMES[row['method']], f"{row['fn']:,.0f}" if fixed else f"{row['fn']:,.2f}",
        f"{row['fp']:,.0f}" if fixed else f"{row['fp']:,.2f}",
        f"{row['f1_percent']:.5f}" if fixed else f"{row['f1_percent']:.5f} ({row['f1_percent_sd']:.4f})",
        f"{row['workflow_seconds']:.2f}"]


def number(value, digits=5, signed=False):
    if value is None:
        return 'n/a'
    return format(value, ('+' if signed else '')+f'.{digits}f')


def footer(canvas, doc):
    canvas.setFont('Helvetica', 7)
    canvas.setFillColor(BLUE)
    canvas.drawString(16*mm, 8*mm, 'barbac | frozen publication benchmark')
    canvas.drawRightString(194*mm, 8*mm, str(doc.page))


def main():
    frozen = check_freeze()
    validation = load(HERE/'execution_validation.json')
    analysis = load(HERE/'analysis_validation.json')
    data = load(HERE/'summary.json')
    raw = pd.DataFrame(load(HERE/'results.json'))
    protocol = load(HERE/'final_protocol.json')
    generation = load(HERE/'generation_validation.json')
    contrasts = load(HERE/'accuracy_contrasts.json')
    n = protocol['n_independent_libraries_per_design']
    # Independently reconcile every displayed aggregate with execution receipts.
    for row in data:
        cells = raw[(raw.condition==row['condition']) & (raw.method==row['method'])]
        assert len(cells) == row['planned_n']
        assert int((cells.status=='complete').sum()) == row['n_successful']
        if row['complete_population']:
            for metric in ['fn','fp','f1_percent','incorrect_reads','unassigned_reads','read_assignment_accuracy_percent','abundance_total_variation']:
                assert np.isclose(row[metric], cells[metric].mean(), rtol=0, atol=1e-10)
            assert np.isclose(row['workflow_seconds'], cells.workflow_seconds.median(), rtol=0, atol=1e-10)
        else:
            assert row['f1_percent'] is None and row['workflow_seconds'] is None
        for cell in cells.to_dict('records'):
            key = 'milos' if cell['scope']=='published_reference' else f"{int(cell['seed'])}/{cell['condition']}"
            dest = HERE/'generated/results'/key/cell['method']
            for name,digest in cell['output_sha256'].items():
                assert sha(dest/name)==digest
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='Note',fontName='Helvetica',fontSize=8,leading=10.5,spaceAfter=5))
    styles.add(ParagraphStyle(name='Tiny',fontName='Helvetica',fontSize=7.2,leading=9,spaceAfter=4))
    story=[Paragraph('Barcode clustering benchmark',styles['Title']),
           Paragraph(f'Preserved v14 + faster startup | {n} independent libraries per simulated design | six configurations',styles['Note']),
           Paragraph('FN and FP: mean counts. F1: mean % (SD). Time: median workflow seconds. Milo: one fixed reference.',styles['Note']),Spacer(1,5)]
    rows=[['Method','FN (lower)','FP (lower)','F1 % (higher)','Time s (lower)']]
    commands=table_style()
    displayed=[]
    for condition,title in TITLES.items():
        i=len(rows);rows.append([title,'','','',''])
        commands += [('SPAN',(0,i),(-1,i)),('BACKGROUND',(0,i),(-1,i),PALE),('FONTNAME',(0,i),(-1,i),'Helvetica-Bold'),('TEXTCOLOR',(0,i),(-1,i),BLUE)]
        selected=[r for r in data if r['condition']==condition]
        full=[r for r in selected if r['complete_population']]
        best={metric:(max if metric=='f1_percent' else min)(r[metric] for r in full) for metric in ['fn','fp','f1_percent','workflow_seconds']} if full else {}
        for row in selected:
            i=len(rows); values=display(row); rows.append(values);displayed.append(values)
            if row['method'] in ['hamming','lv']:commands.append(('TEXTCOLOR',(0,i),(0,i),colors.HexColor('#255D8C')))
            for col,metric in enumerate(best,1):
                if row[metric] is not None and abs(row[metric]-best[metric])<1e-10:
                    commands += [('FONTNAME',(col,i),(col,i),'Helvetica-Bold'),('BACKGROUND',(col,i),(col,i),GOLD)]
    commands += [('LINEBELOW',(0,len(rows)-1),(-1,len(rows)-1),.7,BLUE)]
    tab=Table(rows,colWidths=[52*mm,25*mm,25*mm,47*mm,29*mm]);tab.setStyle(TableStyle(commands))
    story += [tab,Spacer(1,10)]
    notes=[
      f'Bold/shading identifies numerical leaders among complete rows, not significance. {validation["failed_cells"]} failed tool calls are retained. An incomplete row shows successes/planned in brackets; n/a means the full-design result is unavailable. Conditional summaries and all eight accuracy contrasts appear in the supplement.',
      'FN: true identities absent from inferred centroids. FP: inferred identities absent from truth. Exact-identity F1 includes zero-read true barcodes. Read assignment, positive-read truth and abundance metrics are reported separately.',
      f'New simulations: {n} independent seeds, 10,000 true identities and one million expected reads per library; paired designs share variable identities and parent counts. Substitutions: 0.4% per base. Indels: unscaled archived homopolymer-specific rates.',
      f'The calibrated table ends at repeat length 13. Further indel recursion is stopped for affected reads beyond that support; reads, substitutions and truth are retained. Boundary use: {generation["boundary_reads"]:,} reads across {generation["inputs_reaching_boundary"]} of {generation["inputs"]} libraries. This is a declared finite-support model, not an empirical estimate beyond the table.',
      'Milo: unchanged deposited simulation, 100,000 true barcodes and 24,996,128 reads; used during development. Its single-reference scores do not establish independent-library statistical superiority. Its timing is one observation in this final campaign.',
      'All methods use distance 3 and one thread. Barbac: support ordering, ratio 20, configured error proxy 0.005, design scoring off, LV Poisson on. Hamming retains rare-indel rescue. Shepherd: nominal length, Bayes threshold 4. Starcode MP ratio 5. Bartender: seed length 5, step 1, z=5, cutoff 1.',
      'Timing: serial fresh Python-worker wall time including worker/tool startup, required conversion and centroid/member exports; common staging, scoring and hashing excluded. One call per tool/input. Input-dependent timing distributions and paired ratios are in the supplement.',
      'The publication scope was selected after development inspection. Code, settings, seeds and test methods were frozen before final generation. A disclosed reporting addendum retains tool failures and tests only fully observed 60-pair contrasts, keeping the eight-comparison correction. All registered libraries remain. Development controls are preserved.'
    ]
    story += [Paragraph(escape(s),styles['Tiny']) for s in notes]
    pdf=HERE/'benchmark_table.pdf'
    SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=16*mm,rightMargin=16*mm,topMargin=11*mm,bottomMargin=14*mm,title='Barbac final publication benchmark').build(story,onFirstPage=footer,onLaterPages=footer)
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
    assert re.search(r'^Pages:\s+1$',info,re.M),info
    extracted=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True)
    observed=[re.sub(r'\s+',' ',line).strip() for line in extracted.splitlines() if line.strip().startswith(tuple(NAMES.values()))]
    assert observed==[' '.join(row) for row in displayed],(observed,displayed)
    supplement=[Paragraph('Statistical comparisons and full metrics',styles['Title']),
        Paragraph('Accuracy: paired differences between LV and each external competitor. Units are F1 percentage points; positive differences favour LV.',styles['Note'])]
    crows=[['Design / competitor','Mean difference','Lower bound*','Bootstrap lower','Holm p','Supported']]
    for row in contrasts:
        name=('Random' if row['condition']=='random_mixed' else 'Anchored')+' / '+NAMES[row['competitor']]
        crows.append([name,number(row['mean_difference'],signed=True),number(row['simultaneous_one_sided_lower'],signed=True),
            number(row['bootstrap_lower'],signed=True),f"{row['holm_adjusted_p']:.3g}" if row['holm_adjusted_p'] is not None else 'n/a',
            'n/a' if row['status']!='complete' else ('Yes' if row['superiority_supported'] else 'No')])
    ct=Table(crows,colWidths=[60*mm,26*mm,26*mm,28*mm,20*mm,18*mm]);ct.setStyle(TableStyle(table_style()))
    supplement += [ct,Spacer(1,10),Paragraph('*Simultaneous one-sided paired-t lower bounds use Bonferroni alpha/8 (familywise alpha 0.05). The paired bootstrap resamples whole libraries 100,000 times. A superiority label requires both lower bounds above zero. Incomplete contrasts are n/a; all eight positions remain in the correction. The failure-reporting addendum was introduced during execution. Milo is excluded from these tests.',styles['Note']),
        Paragraph('The sample size was fixed at 60 per design after a four-library development variance pilot. The resource cap was binding: the planned 0.01 percentage-point difference is not detectable with 90% power in every contrast. These uncertainty estimates rely on the stated statistical assumptions; no sample-size extension followed final results.',styles['Note']),Spacer(1,8),
        Paragraph('Read assignments and abundance',styles['Heading2'])]
    srows=[['Design / method','Correct reads %','Wrong reads','Unassigned','Abundance TV']]
    for row in data:
        design={'random_mixed':'Random','anchored_mixed':'Anchored','milos':'Milo'}[row['condition']]
        srows.append([design+' / '+NAMES[row['method']],number(row['read_assignment_accuracy_percent'],6),
            number(row['incorrect_reads'],2),number(row['unassigned_reads'],2),number(row['abundance_total_variation'],8)])
    st=Table(srows,colWidths=[68*mm,31*mm,26*mm,24*mm,29*mm]);st.setStyle(TableStyle(table_style()))
    supplement += [st,Spacer(1,8),Paragraph('Simulation entries are means over all 60 independent libraries; Milo entries are one fixed reference. Wrong reads exclude unassigned reads, which are counted separately. Read accuracy counts both as failures. Abundance TV is half the normalized absolute count difference, including unassigned reads.',styles['Tiny']),PageBreak(),
        Paragraph('Timing, simulation limits and reproducibility',styles['Title'])]
    speeds=load(HERE/'speed_contrasts.json')
    trows=[['Design / competitor','LV / competitor time','95% interval']]
    for row in speeds:
        design='Random' if row['condition']=='random_mixed' else 'Anchored'
        trows.append([design+' / '+NAMES[row['competitor']],number(row['geometric_time_ratio'],3),
            number(row['ci95_lower'],3)+' to '+number(row['ci95_upper'],3) if row['status']=='complete' else 'n/a'])
    tt=Table(trows,colWidths=[90*mm,45*mm,43*mm]);tt.setStyle(TableStyle(table_style()))
    supplement += [tt,Spacer(1,8),Paragraph('A time ratio below one favours LV. Ratios are geometric means of within-library workflow ratios. Intervals use paired log times and are secondary descriptive 95% intervals; they are not the multiplicity-adjusted primary accuracy tests.',styles['Note']),
        Paragraph('Simulation boundary',styles['Heading2'])]
    boundary=[r for r in load(HERE/'datasets.json') if r['event_counts'].get('boundary_reads',0)>0]
    if boundary:
        brows=[['Seed','Design','Boundary reads','Total reads']]+[[str(r['seed']),r['condition'],str(r['event_counts']['boundary_reads']),f"{r['input_reads']:,}"] for r in boundary]
        bt=Table(brows,colWidths=[45*mm,65*mm,32*mm,36*mm]);bt.setStyle(TableStyle(table_style()));supplement.append(bt)
    else:
        supplement.append(Paragraph('No final library reached an unmeasured repeat length. The explicit terminal-state policy was therefore unused in the final inputs.',styles['Note']))
    supplement += [Spacer(1,8),Paragraph('Previous error-rate and genuine-neighbour challenges',styles['Heading2']),
        Paragraph('The earlier 0x/1x/3x/10x pilot retains all 24 planned cells, including seven unsupported cases under the previous strict simulator. The rejected paired-indel candidate corrected some errors but lost 31 true identities in six deliberately enriched genuine-neighbour controls. That candidate is excluded here; released v14 preserved those identities. These controls characterize failure modes, not real-world prevalence.',styles['Note']),
        Paragraph('Reproducibility',styles['Heading2']),
        Paragraph(f'Code/protocol commit: {escape(frozen["source_git"])}. Frozen at: {escape(frozen["frozen_at"])}. Full parameters, sample-size calculations, seeds, environment, binary fingerprints, commands, individual metrics and output hashes accompany this table. Existing main, previous benchmark outputs and earlier candidate branches are preserved.',styles['Note']),
        Paragraph('Sources: Johnson, Venkataram and Kryazhimskiy (2023), Best Practices in Designing, Sequencing, and Identifying Random DNA Barcodes (PMC10276077); archived analysis, Zenodo 7411747. Simulation-study design guidance: Morris, White and Crowther (2019), doi:10.1002/sim.8086. Our updated-coordinate and finite-support simulation is an explicitly described adaptation.',styles['Note'])]
    failures=load(HERE/'failure_audit.json')
    if failures:
        supplement += [PageBreak(),Paragraph('Retained tool failures',styles['Title']),
            Paragraph('These registered calls did not produce a valid clustering output. Inputs, commands and logs remain preserved. No failing call was rerun with altered settings and no failed score was replaced by zero. A success-only average is conditional on the method completing and does not establish full-design accuracy or runtime.',styles['Note'])]
        frows=[['Design / method','Seed','Failure']]
        for row in failures:
            reason='Automatic error-rate estimate rejected' if 'Error rate could not be reliably estimated' in row['logs'].get('tool.log','') else 'Tool execution failed; see retained logs'
            frows.append([row['condition']+' / '+NAMES[row['method']],str(row['seed']),reason])
        ft=Table(frows,colWidths=[71*mm,36*mm,71*mm],repeatRows=1);ft.setStyle(TableStyle(table_style()))
        supplement += [ft,Spacer(1,10),Paragraph('Successful-only summaries for incomplete rows',styles['Heading2'])]
        conditional=[r for r in load(HERE/'successful_only_summary.json') if not r['complete_population']]
        urows=[['Design / method','Succeeded','FN','FP','F1 %','Time s']]
        for row in conditional:
            selected=raw[(raw.condition==row['condition'])&(raw.method==row['method'])&(raw.status=='complete')]
            for metric in ['fn','fp','f1_percent']:
                assert np.isclose(row[metric],selected[metric].mean(),rtol=0,atol=1e-10)
            urows.append([row['condition']+' / '+NAMES[row['method']],f"{row['n_successful']}/{row['planned_n']}",
                number(row['fn'],2),number(row['fp'],2),number(row['f1_percent'],5),number(row['workflow_seconds'],2)])
        ut=Table(urows,colWidths=[71*mm,24*mm,20*mm,20*mm,25*mm,18*mm]);ut.setStyle(TableStyle(table_style()))
        supplement += [ut,Spacer(1,10),Paragraph('FN, FP and F1 above are means over successful calls only; time is their median. They use fewer libraries than the main table and must not be compared as if they represented the same complete population. Full per-input results retain every failure.',styles['Note']),
            Paragraph('Failure handling: FAILURE_REPORTING.md and analyze_available.py were added after the first failure. The original frozen analyze.py is retained. The complete-pair calculations, bootstrap seeds and correction family are unchanged; unavailable contrasts have no superiority conclusion. Pawel et al. (2025), Handling Missingness, Failures, and Non-Convergence in Simulation Studies, arXiv:2409.18527v3, discusses the consequences of excluding or replacing failed outputs.',styles['Note'])]
    spdf=HERE/'benchmark_supplement.pdf'
    SimpleDocTemplate(str(spdf),pagesize=A4,leftMargin=16*mm,rightMargin=16*mm,topMargin=12*mm,bottomMargin=15*mm,title='Barbac benchmark statistical supplement').build(supplement,onFirstPage=footer,onLaterPages=footer)
    subprocess.run(['pdftoppm','-png','-singlefile','-scale-to','1800',str(pdf),str(HERE/'benchmark_table')],check=True)
    save(HERE/'report_validation.json',dict(rows=18,main_pdf_pages=1,all_displayed_aggregates_reconciled=True,
        all_execution_output_hashes_verified=True,main_pdf_text_rows_verified=True,pdf_sha256=sha(pdf),supplement_sha256=sha(spdf)))
    lines=['# Final publication benchmark','', '[One-page table](benchmark_table.pdf) · [Statistical supplement](benchmark_supplement.pdf) · [Full per-input results](all_results.csv) · [Protocol](README.md)','',
        f'{n} independent libraries per simulated design; all {validation["cells"]} registered tool cells attempted, {validation["failed_cells"]} failed. {analysis["supported_contrasts"]} of eight specified LV-versus-competitor F1 contrasts meet both lower-bound criteria; {8-analysis["available_contrasts"]} are unavailable. Only full 60-pair contrasts are tested under the disclosed [failure-reporting addendum](FAILURE_REPORTING.md). Milo is a fixed reference, not an independent replication test.','']
    for condition,title in TITLES.items():
        lines += ['## '+title,'','| Method | FN | FP | F1 % (SD) | Workflow seconds |','|---|---:|---:|---:|---:|']
        lines += ['| '+' | '.join(display(row))+' |' for row in data if row['condition']==condition]
        lines.append('')
    lines += ['## Accuracy contrasts','','| Design | Competitor | Mean F1 difference (points) | Simultaneous lower bound | Bootstrap lower bound | Supported |','|---|---|---:|---:|---:|---|']
    for row in contrasts:
        outcome='n/a' if row['status']!='complete' else ('Yes' if row['superiority_supported'] else 'No')
        lines.append(f"| {row['condition']} | {NAMES[row['competitor']]} | {number(row['mean_difference'],signed=True)} | {number(row['simultaneous_one_sided_lower'],signed=True)} | {number(row['bootstrap_lower'],signed=True)} | {outcome} |")
    lines += ['', '## Interpretation and scope','']+notes
    (HERE/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print('VERIFIED',pdf,spdf,flush=True)


if __name__=='__main__':
    main()
