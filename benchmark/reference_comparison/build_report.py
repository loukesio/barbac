"""Build the canonical report payload from validated benchmark results."""
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

HERE=Path(__file__).resolve().parent
NAMES={
 'barbac_hamming_sequence':'barbac Hamming / sequence',
 'barbac_hamming_support':'barbac Hamming / support',
 'barbac_lv_sequence':'barbac LV / sequence',
 'barbac_lv_support':'barbac LV / support',
 'previous_hamming_sequence':'Previous Hamming / sequence',
 'previous_lv_sequence':'Previous LV / sequence',
 'shepherd':'Shepherd', 'starcode_sphere':'Starcode / sphere',
 'starcode_mp':'Starcode / message passing', 'bartender':'Bartender'}


def main():
    raw=pd.read_csv(HERE/'results.csv')
    query=(HERE/'report_metrics.sql').read_text()
    with sqlite3.connect(':memory:') as connection:
        raw.to_sql('reference_runs',connection,index=False)
        frame=pd.read_sql_query(query,connection)
    original=raw.set_index('method')
    recalculated=frame.set_index('method').reindex(original.index)
    for field in ['f1','precision','recall','positive_truth_f1','read_assignment_accuracy']:
        assert (abs(original[field]-recalculated[field]) < 1e-14).all(), field
    assert set(frame.method)==set(NAMES), 'Wait for every method before reporting.'
    assert frame.mapping_counts_reconcile.all()
    frame['label']=frame.method.map(NAMES)
    frame['f1_display']=frame.f1.map(lambda n:f'{100*n:.5f}%')
    frame['positive_f1_display']=frame.positive_truth_f1.map(lambda n:f'{100*n:.5f}%')
    frame['assignment_display']=frame.read_assignment_accuracy.map(lambda n:f'{100*n:.6f}%')
    frame['errors']=frame.misassigned_reads+frame.unassigned_reads
    scores=frame.set_index('method')
    best=frame.loc[frame.f1.idxmax()]
    fast=frame.loc[frame.pipeline_seconds.idxmin()]
    hs=scores.loc['barbac_hamming_support']; hq=scores.loc['barbac_hamming_sequence']
    ls=scores.loc['barbac_lv_support']; lq=scores.loc['barbac_lv_sequence']
    prevh=scores.loc['previous_hamming_sequence'];prevl=scores.loc['previous_lv_sequence']
    shep=scores.loc['shepherd']
    title='barbac on the 100,000-barcode reference simulation'
    blocks=[]
    def prose(id,body,source='results'):
        blocks.append(dict(id=id,type='markdown',body=body,**({'sourceId':source} if source else {})))
    prose('title','# '+title,None)
    prose('summary',f"## The result on this dataset\n\n**{best.label} has the highest centroid F1 ({best.f1_display}).** {fast.label} has the shortest measured workflow time ({fast.pipeline_seconds:.1f} seconds). Updated Hamming with support ordering takes {hs.pipeline_seconds:.1f} seconds versus Shepherd's {shep.pipeline_seconds:.1f} seconds ({shep.pipeline_seconds/hs.pipeline_seconds:.1f} times faster in this run), with FN {int(hs.fn)} / FP {int(hs.fp)} versus Shepherd's FN {int(shep.fn)} / FP {int(shep.fp)}. Its {int(hs.errors)} wrong or unassigned reads compare with {int(shep.errors)} for Shepherd, out of nearly 25 million. **The accuracy differences are very small; this is not evidence of universal superiority.** These are descriptive results from one supplied simulation and one retained serial timing observation per configuration.")
    prose('definitions',"## Barcode recovery and read assignment answer different questions\n\n**FN** counts true barcode sequences missing from the output; **FP** counts extra output centroid sequences. **F1 = 2 TP / (2 TP + FN + FP)**, where TP is the number of recovered true sequences. The historical score uses all 100,000 truth entries. A second score uses only the 99,591 truth entries with positive read counts.\n\n**Read assignment** is the count-weighted fraction of all 24,996,128 reads whose output centroid exactly equals the supplied parent label. Wrong and unassigned reads both reduce this fraction. This measures agreement with supplied labels, which contain small inconsistencies described below.",None)
    prose('accuracy',f"## Support ordering changes only a few recovery decisions\n\nHamming changes from FN {int(hq.fn)} / FP {int(hq.fp)} with sequence ordering to FN {int(hs.fn)} / FP {int(hs.fp)} with support ordering. LV changes from FN {int(lq.fn)} / FP {int(lq.fp)} to FN {int(ls.fn)} / FP {int(ls.fp)}. These small differences should be read as counts, rather than rounded percentages that would hide them. The table retains all configurations, including the previous implementation and both Starcode modes.")
    blocks.append(dict(id='accuracy-table',type='table',tableId='accuracy'))
    prose('timing',f"## Runtime depends strongly on the distance mode\n\nUpdated Hamming with sequence ordering takes {hq.core_seconds:.1f} seconds inside barbac versus {prevh.core_seconds:.1f} seconds for the previous implementation. Updated LV takes {lq.core_seconds:.1f} seconds versus {prevl.core_seconds:.1f} seconds previously. Support ordering adds preprocessing, bringing current Hamming to {hs.core_seconds:.1f} seconds and LV to {ls.core_seconds:.1f} seconds.\n\nThe bars compare complete measured workflows in seconds; shorter is faster. They include method-specific format conversion and output of cluster membership, including read expansion for Bartender. Shared input staging and scoring are excluded. Each configuration ran once, sequentially, on the same machine; these are observations without timing confidence intervals.")
    blocks.append(dict(id='runtime-chart',type='chart',chartId='runtime'))
    prose('boundaries',"## Separate clustering time from operational time\n\nCore time is available only for barbac and includes reading the CSV, ordering, clustering, and building the returned result. Process time includes startup and native outputs. Workflow time additionally includes required method-specific input preparation and conversion of outputs to common centroid and membership tables. Compare workflow time across tools; use core time to compare barbac implementations.",None)
    blocks.append(dict(id='timing-table',type='table',tableId='timings'))
    prose('assignment',f"## Most reads agree with the supplied parent labels\n\nUpdated Hamming with support ordering assigns {hs.assignment_display} of reads to their supplied parents. Shepherd has {int(shep.misassigned_reads)} wrongly assigned reads and {int(shep.unassigned_reads)} unassigned reads. The table separates these outcomes, so a method cannot improve its reported accuracy by dropping difficult reads. Every output was checked: member assignments are unique, and their summed counts reproduce the reported centroid counts.")
    blocks.append(dict(id='assignment-table',type='table',tableId='assignments'))
    prose('quality',"## Zero-read barcodes dominate the historical false negatives\n\nThe three input files contain 1,544,850 unique observed sequences, 100,000 true barcodes, and 24,996,128 reads. All observed counts and sequences match between the input and labeled source; keys are unique and nonnull, and every parent label appears in the truth table.\n\n**409 true barcodes have zero reads.** Another 30 true sequences are absent in exact form despite having reads. A method choosing centroids from observed strings therefore cannot recover those 439 exact sequences. A consensus method could reconstruct an absent sequence. Positive-count and observed-only FN metrics distinguish these cases.\n\n**Four parent barcode totals disagree across files**, with a summed absolute difference of six reads. Grand totals still agree. The source is retained unchanged. This is a small but real limitation on label-based accuracy; the CSV includes a sensitivity score excluding all 2,116 reads assigned to those four parents. The source does not establish a hard upper bound on hidden read-label errors. No temporal analysis applies to this static simulation.",'quality')
    prose('methods',"## Fixed settings and label-blind inputs\n\nAll methods use the same observed sequences and counts, sorted by decreasing count with sequence ordering for ties. Truth labels are used only after clustering. Maximum distance is three. barbac uses merge ratio 20 and error rate 0.005, with its default design option disabled. Hamming retains its existing trace-indel rescue. Shepherd uses barcode length 20 and epsilon three. Starcode is evaluated with sphere clustering and with its default message-passing mode. Bartender uses distance three and its other defaults, including unique read identifiers after expansion.\n\nRuns are sequential; Starcode and Bartender explicitly request one thread, and the environment limits BLAS/OpenMP threads to one. Exact tool commands, platform, source hashes, installed libraries, and output hashes are saved alongside the runner. Historical timings and scores used different input ordering or measurement boundaries and are not interchangeable with this run.",None)
    prose('next',"## Keep the claim specific and investigate the LV cost\n\nThis comparison establishes recovery, label agreement, count conservation, and observed runtime for the supplied dataset and settings. It does not establish universal superiority or accuracy on real sequencing data. Preserve both barbac ordering options: support ordering trades extra time for a small change in accuracy here. Profile LV candidate enumeration on this large library before promising a general speed improvement.\n\nThe next requested stage is time-series download guidance and a configurable SLURM workflow. Real time-series data can test scaling and reproducibility, but will need independent controls or replicate agreement to assess accuracy without simulated parent labels.",None)
    prose('questions',"## What could change the conclusion?\n\nIndependent simulations with different barcode densities and error models, repeated timing runs, and labeled spike-ins or real-data controls could change the ranking. How the methods behave at high indel rates is a separate question: this dataset has only 602 off-length observed sequences. The support-ordering rule and clustering parameters were fixed before this comparison; no thresholds were tuned to improve these results.",None)
    sources=[dict(id='results',label='Fresh serial benchmark results, 8 September 2026',path='benchmark/reference_comparison/results.csv',query=dict(language='sql',engine='SQLite',sql=query,description='Independently recomputes F1, precision, recall, and read-assignment accuracy from benchmark counts produced by run_comparison.py. The in-memory reference_runs table is loaded from results.csv; process and workflow timings are measured values, unchanged by SQL.',tables_used=['reference_runs (loaded from benchmark/reference_comparison/results.csv)'],metric_definitions={'F1':'2TP/(2TP+FN+FP), all 100000 truth barcode strings','read_assignment_accuracy':'sum Count where predicted centroid equals true_BC / sum all input Count','pipeline_seconds':'serial tool process plus method-specific input and output conversion; excludes common staging/scoring'})),dict(id='quality',label='Validated input profile',path='benchmark/reference_comparison/data_quality.json')]
    def table(id,title,columns,sort,direction='asc'):
        return dict(id=id,title=title,dataset='results',sourceId='results',defaultSort=dict(field=sort,direction=direction),columns=[dict(field=f,label=l) for f,l in columns])
    tables=[table('accuracy','Centroid recovery',[('label','Method'),('fn','FN'),('fp','FP'),('f1_display','F1, all truth'),('positive_truth_fn','FN, positive truth'),('positive_f1_display','F1, positive truth')],'f1_display','desc'),table('timings','Timing components, seconds',[('label','Method'),('core_seconds','barbac core'),('process_seconds','Process'),('input_preparation_seconds','Input preparation'),('pipeline_seconds','Workflow')],'pipeline_seconds'),table('assignments','Read assignment and conservation',[('label','Method'),('assignment_display','Correct label'),('misassigned_reads','Wrong'),('unassigned_reads','Unassigned'),('output_reads','Output reads')],'misassigned_reads')]
    chart=dict(id='runtime',title='Workflow runtime by method',subtitle='Seconds; one serial run per configuration, including required format conversion',showDescription=True,type='bar',dataset='results',sourceId='results',encodings=dict(x=dict(field='label',type='nominal',label='Method'),y=dict(field='pipeline_seconds',type='quantitative',label='Seconds')),options=dict(orientation='horizontal'))
    chart['source'] = sources[0]
    artifact=dict(surface='report',manifest=dict(version=1,surface='report',title=title,generatedAt=datetime.now(timezone.utc).isoformat(),blocks=blocks,tables=tables,charts=[chart],sources=sources),snapshot=dict(version=1,status='ready',datasets={'results':json.loads(frame.sort_values('pipeline_seconds').to_json(orient='records'))}),sources=sources)
    (HERE/'artifact.json').write_text(json.dumps(artifact,indent=2)+'\n')

if __name__=='__main__':main()
