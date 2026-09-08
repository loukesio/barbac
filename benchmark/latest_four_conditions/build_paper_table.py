"""Build the five-dataset paper table and update the manuscript comparison.

Includes only distance-three results. Accuracy and timing sample sizes are
explicit; read-assignment accuracy is confined to the labeled reference data.
"""
import json
from pathlib import Path
import re
import textwrap

import pandas as pd
from run_comparison import HERE, ROOT, BASE, CONDITIONS, sha

METHODS=['hamming','lv','shepherd','starcode_sphere','starcode_mp','bartender']
NAMES={'hamming':'barbac Hamming','lv':'barbac LV + Poisson','shepherd':'Shepherd',
       'starcode_sphere':'Starcode sphere','starcode_mp':'Starcode MP','bartender':'Bartender'}
DATASETS={
    'random_substitutions':'Random: substitutions',
    'random_low_indels':'Random: substitutions + indels',
    'anchored_substitutions':'Anchored: substitutions',
    'anchored_low_indels':'Anchored: substitutions + indels',
    'milos':'Milos / Johnson reference'}
CODES=dict(zip(CONDITIONS+['milos'],['R-S','R-I','A-S','A-I','Milos']))


def validate_reference_outputs():
    provenance=json.loads((ROOT/'benchmark/reference_comparison/provenance.json').read_text())
    for name in ['peer_manifest.json','completion_manifest.json']:
        current=json.loads((HERE/name).read_text())
        for file,expected in current['tool_sha256'].items():
            assert expected==provenance['tool_file_sha256'][file],file
    directory=BASE/'generated/reference_2026-09-08'
    for file,expected in provenance['input_sha256'].items(): assert sha(directory/file)==expected
    for name in ['shepherd','starcode_sphere','starcode_mp','bartender']:
        saved=json.loads((ROOT/'benchmark/reference_comparison'/f'{name}.json').read_text())
        for file,expected in saved['output_sha256'].items(): assert sha(directory/name/file)==expected
    saved=json.loads((BASE/'generated/lv-v13/reference/0/poisson/support/result.json').read_text())
    assert saved['build_id']=='barbac-2026-09-08-bounded-indels-v13'
    assert saved['input_sha256']==provenance['input_sha256']['input.csv']
    for file,expected in saved['output_sha256'].items():
        assert sha(BASE/'generated/lv-v13/reference/0/poisson/support'/file)==expected
    return saved


def main():
    assert json.loads((HERE/'completion_manifest.json').read_text())['status']=='complete'
    raw=pd.read_csv(HERE/'results.csv')
    peers=pd.read_csv(HERE/'fresh_peers.csv')
    extra=pd.read_csv(HERE/'additional_peers.csv')
    assert len(extra)==24 and len(peers)==8
    barbac=raw[(raw.origin=='fresh_v13')&(raw.distance==3)].copy()
    retained=raw[(raw.origin=='retained_peer')&(raw.seed!=42)].copy()
    four=pd.concat([barbac,retained,peers,extra],ignore_index=True)
    four['method']=four.method.replace({'starcode':'starcode_sphere'})
    assert len(four)==72 and not four.duplicated(['condition','seed','method']).any()
    assert set(four.method)==set(METHODS)
    assert (four.groupby(['condition','method']).seed.nunique()==3).all()
    assert ((four.f1-2*four.tp/(2*four.tp+four.fn+four.fp)).abs()<1e-12).all()
    lv_ref=validate_reference_outputs()
    ref=pd.read_csv(ROOT/'benchmark/reference_comparison/results.csv')
    ref_rows=[]
    ham=json.loads((HERE/'milos_hamming.json').read_text())
    assert ham['build_id']=='barbac-2026-09-08-bounded-indels-v13'
    for file,expected in ham['output_sha256'].items():
        assert sha(BASE/'generated/latest-four-v13/paper-completion/milos/barbac_hamming_support'/file)==expected
    for method in METHODS:
        if method=='hamming': row=ham
        elif method=='lv': row=lv_ref
        else: row=ref[ref.method==method].iloc[0].to_dict()
        ref_rows.append(dict(condition='milos',seed=0,method=method,distance=3,
                             origin='validated_reference',tp=row['tp'],fn=row['fn'],fp=row['fp'],
                             f1=row['f1'],core_seconds=row.get('core_seconds'),
                             workflow_seconds=row.get('pipeline_seconds',row.get('process_seconds')),
                             input_reads=row['input_reads'],output_reads=row['output_reads'],
                             unassigned_reads=row['unassigned_reads'],misassigned_reads=row['misassigned_reads'],
                             read_assignment_accuracy=row['read_assignment_accuracy']))
    reference=pd.DataFrame(ref_rows)
    assert ((reference.f1-2*reference.tp/(2*reference.tp+reference.fn+reference.fp)).abs()<1e-12).all()
    assert ((reference.read_assignment_accuracy-(1-(reference.misassigned_reads+reference.unassigned_reads)/reference.input_reads)).abs()<1e-12).all()
    all_runs=pd.concat([four,reference],ignore_index=True)
    all_runs.to_csv(HERE/'paper_source_runs.csv',index=False)
    rows=[]
    for condition in CONDITIONS+['milos']:
        for method in METHODS:
            sample=all_runs[(all_runs.condition==condition)&(all_runs.method==method)]
            assert len(sample)==(1 if condition=='milos' else 3)
            timing=sample[sample.seed==(0 if condition=='milos' else 42)].iloc[0]
            assert timing.origin!='retained_peer'
            rows.append(dict(dataset=condition,dataset_label=DATASETS[condition],dataset_code=CODES[condition],
                             method=method,method_label=NAMES[method],accuracy_n=len(sample),
                             true_barcodes=100000 if condition=='milos' else 10000,
                             fn=sample.fn.mean(),fp=sample.fp.mean(),f1=sample.f1.mean(),
                             f1_percent=100*sample.f1.mean(),f1_min_percent=100*sample.f1.min(),
                             f1_max_percent=100*sample.f1.max(),workflow_seconds=timing.workflow_seconds,
                             core_seconds=timing.core_seconds,timing_n=1,
                             timing_seed=0 if condition=='milos' else 42,
                             mean_unassigned_reads=sample.unassigned_reads.mean()))
    table=pd.DataFrame(rows)
    table.to_csv(HERE/'paper_table.csv',index=False)
    columns=['Dataset','Method','FN','FP','F1 (%)','Time (s)']
    md=['| '+' | '.join(columns)+' |','|---|---|---:|---:|---:|---:|']
    tex=['% Requires booktabs and longtable. Values are generated from paper_table.csv.',
         '\\begin{longtable}{llrrrr}',
         '\\caption{Latest barcode clustering comparison at distance three.}\\label{tab:latest-barbac}\\\\',
         '\\toprule','Dataset & Method & FN & FP & F1 (\\%) & Time (s) \\\\',
         '\\midrule','\\endfirsthead','\\toprule',
         'Dataset & Method & FN & FP & F1 (\\%) & Time (s) \\\\', '\\midrule','\\endhead']
    for row in table.itertuples():
        group=table[table.dataset==row.dataset]
        f1=f'{row.f1_percent:.5f}' if row.dataset=='milos' else f'{row.f1_percent:.3f}'
        count_format=',.0f' if row.dataset=='milos' else ',.1f'
        cells=[row.dataset_code,row.method_label,format(row.fn,count_format),format(row.fp,count_format),f1,f'{row.workflow_seconds:.2f}']
        bold_f1=abs(row.f1-group.f1.max())<1e-14
        md_cells=cells.copy()
        if bold_f1: md_cells[4]='**'+md_cells[4]+'**'
        md.append('| '+' | '.join(md_cells)+' |')
        if bold_f1: cells[4]='\\textbf{'+cells[4]+'}'
        tex.append(' & '.join(cells)+' \\\\')
    tex.extend(['\\bottomrule','\\end{longtable}'])
    caption=(
        '**Table 2. Latest barcode clustering performance at distance three.** '
        'R-S: random barcodes with substitutions; R-I: random barcodes with substitutions and indels; '
        'A-S and A-I: the corresponding anchored designs. Milos denotes the Johnson et al. (2023) reference simulation. '
        'FN and FP are absolute barcode counts; for R-S, R-I, A-S and A-I, counts and centroid F1 are means across seeds 42, 43 and 44 (10,000 true barcodes and one million reads per seed). '
        'The Milos row for each method uses one fixed dataset of 100,000 true barcodes and 24,996,128 reads. '
        'Bold F1 values identify the highest score within a dataset. F1 measures exact centroid recovery and differs from the cluster-label F1 used in Section 3.1. '
        'Time is one serial workflow observation in seconds on the same development machine: seed 42 for the four smaller simulations and the full dataset for Milos. '
        'It includes program startup, required format conversion, clustering, and centroid/member exports; shared input staging and scoring are excluded. '
        'These timings have no confidence intervals. Both barbac modes use native v13 and support ordering; LV additionally enables the experimental Poisson indel option. '
        'Starcode sphere and default message passing (MP) are reported separately. All methods use distance three, with its method-specific Hamming or Levenshtein interpretation.')
    markdown=caption+'\n\n'+'\n'.join(md)+'\n'
    (HERE/'paper_table.md').write_text(markdown)
    tex_caption=caption.replace('**Table 2. ','').replace('**','').replace('%','\\%').replace('_','\\_')
    tex[2]='\\caption{'+tex_caption+'}\\label{tab:latest-barbac}\\\\'
    (HERE/'paper_table.tex').write_text('\n'.join(tex)+'\n')
    scores=table.set_index(['dataset','method'])
    beats_all_external=all(scores.loc[(c,'lv'),'f1']>table[(table.dataset==c)&~table.method.isin(['hamming','lv'])].f1.max() for c in CONDITIONS+['milos'])
    barbac_wins=sum(table[(table.dataset==c)&table.method.isin(['hamming','lv'])].f1.max() >= table[table.dataset==c].f1.max()-1e-14 for c in CONDITIONS)
    lv=scores.loc[('milos','lv')];shep=scores.loc[('milos','shepherd')]
    headline=(f'The highest-scoring barbac configuration achieved the highest mean centroid F1 in {barbac_wins} of the four simulation conditions, and LV achieved the highest centroid F1 on Milos (Table 2). '
              'Bartender narrowly led the random substitution-only condition. ')
    methods=(
        'We compared barbac Hamming and Levenshtein (LV) clustering with Shepherd, Starcode sphere, Starcode default message passing, and Bartender using four simulation conditions and the Johnson et al. (2023) reference dataset, referred to here as Milos (Table 2). '
        'The four conditions crossed fully random 20-base barcodes with anchored 28-base barcodes containing 16 variable positions, and substitutions alone with substitutions plus indels. '
        'Each condition used 10,000 true barcodes and one million reads with lognormal abundances (sigma = 1.5), independently generated with seeds 42, 43 and 44. '
        'The substitution probability was 0.005 per base; indel conditions additionally used insertion and deletion probabilities of 0.005 per opportunity each. '
        'These simulations evaluate specified error regimes rather than representing every sequencing platform. The Milos input contains 1,544,850 unique observed sequences and 24,996,128 reads for 100,000 listed true barcode sequences.\n\n'
        'All comparisons used maximum distance three. Both barbac modes used `super_cluster2()`, native build v13, merge ratio 20, configured error rate 0.005, support ordering for equal-count sequences, and the design option disabled. '
        'LV additionally enabled `indel_model = "poisson"`, an experimental exception for repeated-base single indels whose abundance is consistent with an expected-error model. '
        'This option is disabled by default; the table explicitly evaluates the enabled configuration. The error rate was supplied, not fitted to true labels. '
        'Shepherd and Bartender used their remaining default parameters, with Shepherd estimating its substitution error rate automatically. Starcode sphere and default message passing were tested separately. Timed Starcode and Bartender runs each requested one thread. '
        'Input sequences were ordered by decreasing observed count and then sequence, without using parent labels. '
        'The v13 LV search optimization excludes candidate parents whose best possible likelihood score cannot improve the current assignment; the indexed results were checked against existing outputs and full-scan tests. '
        'Hamming refinement incorporates the binomial criterion described for Shepherd (Tavakolian et al., 2022).\n\n'
        'For this comparison, a true positive (TP) is an output centroid string present in the true barcode set; FN counts missing true strings and FP counts extra centroid strings. '
        'Centroid F1 is `2 TP / (2 TP + FN + FP)`. All listed truth sequences remain in the denominator, including those whose exact sequence does not occur in the noisy input. '
        'For the Milos simulation this includes 409 truth barcodes with zero reads. Input identity, output uniqueness and count reconciliation were checked before scoring. '
        'All barbac outputs conserved the supplied reads. Per-read parent labels were not retained for the four smaller simulations, so no read-assignment accuracy is inferred from their centroid F1 values.')
    results=(headline+
        'On random substitution-only data, Hamming and LV had identical centroid recovery, while Hamming slightly outperformed LV on anchored substitution-only data. '
        'The main distinction emerged in indel-containing data: LV produced far fewer false clusters than Hamming, Shepherd, and Bartender. '
        f"For random indel data, LV yielded mean FN {scores.loc[('random_low_indels','lv'),'fn']:.1f} and FP {scores.loc[('random_low_indels','lv'),'fp']:.1f}, with F1 {scores.loc[('random_low_indels','lv'),'f1_percent']:.3f}%. "
        f"For anchored indel data, the corresponding values were FN {scores.loc[('anchored_low_indels','lv'),'fn']:.1f}, FP {scores.loc[('anchored_low_indels','lv'),'fp']:.1f}, and F1 {scores.loc[('anchored_low_indels','lv'),'f1_percent']:.3f}%. "
        'Hundreds of false clusters remained in the latter condition, so relative superiority does not imply error-free reconstruction. '
        'The Poisson exception removed five false clusters across these 12 smaller datasets without changing FN; most of the observed difference from competitors reflects the broader LV clustering strategy and support ordering.\n\n'
        f"On Milos, LV returned FN {lv.fn:.0f} and FP {lv.fp:.0f}, giving centroid F1 {lv.f1_percent:.5f}%, compared with Shepherd's FN {shep.fn:.0f}, FP {shep.fp:.0f}, and F1 {shep.f1_percent:.5f}%. "
        'This F1 difference is very small and is reported descriptively. '
        'With the supplied read-parent labels, LV made 131 wrong assignments and left no reads unassigned; Shepherd made 129 wrong assignments and left 88 reads unassigned. '
        'Thus the combined wrong-or-unassigned count was 131 for LV and 217 for Shepherd. '
        'The reference labels contain four discrepancies in per-parent totals, spanning six reads in absolute difference; these source inconsistencies were retained and audited.')
    timing=(
        'Runtime depended on dataset size, sequence design and the reported boundary. '
        f"On Milos, the LV workflow took {lv.workflow_seconds:.1f} s, compared with {shep.workflow_seconds:.1f} s for Shepherd ({shep.workflow_seconds/lv.workflow_seconds:.2f}-fold shorter). "
        f"Hamming took {scores.loc[('milos','hamming'),'workflow_seconds']:.1f} s and was the fastest tested workflow on that dataset. "
        'On the smaller simulations, interpreter/package startup contributed several seconds to barbac workflows, and native competitors were sometimes faster end to end despite poorer recovery. '
        'The table therefore reports complete measured workflows rather than comparing barbac core time against another method\'s total runtime. '
        'Times for the four simulations are fresh serial seed-42 observations; accuracy additionally includes seeds 43 and 44. Previously recorded peer accuracy outputs were reused only after verifying identical inputs and output hashes. '
        'Earlier overlapping benchmark timings were excluded. Timing repetitions sufficient for confidence intervals were not collected.')
    limitations=(
        'The results support an accuracy advantage for the tested LV configuration in the indel-containing simulations, with Hamming useful for substitution-only inputs. '
        'They do not establish universal superiority across libraries or sequencing platforms. The Milos dataset has only 1,001 reads whose lengths differ from their supplied parents, whereas the four-condition experiment explicitly includes broader mixed indel errors. '
        'Shepherd includes separate correction of simple single insertions and deletions around a supplied barcode length; its poorer performance in mixed-error simulations should not be described as complete absence of indel support. '
        'The Poisson option can merge real length variants with error-like counts, and its configured error rate is not a learned platform-specific indel rate. '
        'Independent labeled experimental controls, additional abundance distributions and repeated timings are needed before broader accuracy or speed claims.')
    section='**3.3 Comparison with established error-correction methods**\n\n'+methods+'\n\n'+markdown+'\n'+results+'\n\n**3.4 Runtime and operational cost**\n\n'+timing+'\n\n**3.5 Scope and limitations of the comparison**\n\n'+limitations+'\n\n'
    section='\n\n'.join(block if block.startswith('|') else textwrap.fill(block,width=88,break_long_words=False,break_on_hyphens=False) for block in section.strip().split('\n\n'))+'\n\n'
    (HERE/'paper_results_section.md').write_text(section)
    discussion=(
        'The current clustering routine uses abundance- and error-rate-based scoring, with an optional Poisson consistency rule for repeated-base single indels. '
        'Its parameters remain configurable and require validation for each experimental setting.\n\n'+headline+
        'The strongest practical distinction is the substantially lower false-cluster burden of LV in the simulated indel conditions. '
        'Hamming remains useful when errors are predominantly substitutions, but its low FN count in indel-containing data is accompanied by many extra centroids. '
        'The narrow Milos F1 advantage and the remaining errors on anchored indel data argue for reporting FN, FP and F1 together. '
        'Workflow runtime also matters: barbac is substantially faster than Shepherd on Milos and the anchored benchmarks, while other native workflows can be faster on small inputs. '
        'These findings support the tested barbac configurations as competitive choices for lineage analysis, with a demonstrated advantage in defined indel regimes, rather than establishing a universal ranking.\n\n'
        'The choice between missing a rare true barcode and retaining a spurious centroid depends on the downstream analysis. '
        'False clusters can inflate inferred lineage diversity, whereas missed barcodes can obscure rare lineages. '
        'The comparison was limited to three seeds per smaller simulation and one fixed reference dataset, and it did not measure per-read assignment accuracy for the smaller simulations. '
        'Additional labeled experimental controls and broader error and abundance regimes remain necessary to assess generalization.\n\n')
    discussion='\n\n'.join(textwrap.fill(block,width=88,break_long_words=False,break_on_hyphens=False) for block in discussion.strip().split('\n\n'))+'\n\n'
    paper=ROOT/'manuscript/barbac_manuscript.md'
    text=paper.read_text()
    start=text.index('**3.3 Comparison with established error-correction methods**')
    end=text.index('**4 Discussion**',start)
    before=text[:start]; after=text[end:]
    # Change only discussion text dependent on the replaced comparisons.
    if 'Furthermore, while barbac efficiently integrates' in after:
        after,n=re.subn(r'Furthermore, while barbac efficiently integrates.*?(?=In summary, barbac)',discussion,after,flags=re.S)
        assert n==1
    elif 'The current clustering routine uses abundance-' in after:
        after,n=re.subn(r'The current clustering routine uses abundance-.*?(?=In summary, barbac)',discussion,after,flags=re.S)
        assert n==1
    else: raise AssertionError('Cannot locate comparison-specific discussion')
    updated=before+section+after
    assert updated.count('**Table 2.')==1
    paper.write_text(updated)
    summary=table.pivot(index='dataset',columns='method',values='f1_percent').reindex(CONDITIONS+['milos'])[METHODS]
    (HERE/'paper_summary_f1.csv').write_text(summary.to_csv())
    validation=dict(status='complete',source_rows=78,paper_rows=30,
                    distance=3,accuracy_seeds=[42,43,44],timing_seed=42,
                    lv_f1_exceeds_external_comparators_in_all_five=bool(beats_all_external),
                    simulation_categories_with_best_barbac_f1=int(barbac_wins),
                    simulation_accuracy='Mean of three seed-level FN/FP/F1 values; no read-parent score available.',
                    timing='One complete serial workflow observation per table cell; no timing confidence intervals.',
                    manuscript_scope='Replaces Sections 3.3–3.5 and dependent comparison discussion; other manuscript sections retained.',
                    source_sha256={str(p.relative_to(ROOT)):sha(p) for p in [HERE/'results.csv',HERE/'fresh_peers.csv',HERE/'additional_peers.csv',HERE/'milos_hamming.json',ROOT/'benchmark/reference_comparison/results.csv',Path(__file__)]})
    (HERE/'paper_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(summary.to_string())
    print(table[['dataset','method','fn','fp','f1_percent','workflow_seconds']].to_string(index=False))


if __name__=='__main__': main()
