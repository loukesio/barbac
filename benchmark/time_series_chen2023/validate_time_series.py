"""Independently reconcile complete time-series outputs before reporting results."""
import argparse
import json
from pathlib import Path
import re
import zipfile

import numpy as np
import pandas as pd

from extract_barcodes import sha256
from run_sample import checked

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    args=parser.parse_args();work=args.work_dir.resolve();out=HERE/'results'
    provenance=json.loads((out/'provenance.json').read_text())
    assert provenance['status']=='complete'
    for path,expected in provenance['files'].items():assert sha256(ROOT/path)==expected,path
    source=HERE/'generated/sources/hBFA1_all_freqs_tidy.csv'
    assert sha256(source)==provenance['publication_sha256']
    publication=pd.read_csv(source)
    publication=publication[publication.Test_Environment.eq('YPD') & publication.Time.isin([8,16,24,40])]
    assert len(publication)==18512
    samples=pd.read_csv(HERE/'samples.tsv',sep='\t')
    extraction=pd.read_csv(out/'extraction_summary.csv')
    assert len(extraction)==8 and extraction['sample'].is_unique
    for row in samples.to_dict('records'):
        raw=work/'raw'/row['run'];totals=[]
        for mate in (1,2):
            assert checked(raw/f'{row["run"]}_{mate}.fastq.gz',row,mate)
            qc=work/'mapping'/row['sample']/'fastQC'/f'{row["run"]}_{mate}_fastqc.zip'
            with zipfile.ZipFile(qc) as z:
                text=z.read(next(n for n in z.namelist() if n.endswith('/fastqc_data.txt'))).decode()
            totals.append(int(re.search(r'^Total Sequences\t(\d+)',text,re.M).group(1)))
        actual=extraction[extraction['sample'].eq(row['sample'])].iloc[0]
        assert totals[0]==totals[1]==actual.total_pairs
        receipt=json.loads((work/'extracted'/row['sample']/'extraction.json').read_text())
        assert receipt['inline1_top']==[[row['inline1'],int(actual.total_pairs)]]
        assert receipt['inline2_top']==[[row['inline2'],int(actual.total_pairs)]]
        assert actual.total_pairs==sum(actual[k] for k in ['short_reads','quality_failed',
            'no_complete_BAM_barcode_pair','umi_duplicates','retained_molecules'])
        assert actual.barbac_input_molecules==actual.retained_molecules-actual.non_acgt_molecules
        path=work/'extracted'/row['sample']/'pairs_barbac.csv'
        assert sha256(path)==provenance['extraction_inputs'][row['sample']]
        pairs=pd.read_csv(path)
        assert pairs.counts.sum()==actual.barbac_input_molecules
        assert pairs.diverse_barcode.str.len().between(24,28).all()
        assert pairs.environment_barcode.str.len().between(24,28).all()
    agreement=pd.read_csv(out/'sample_method_agreement.csv')
    assert len(agreement)==48 and not agreement.duplicated(['method','sample']).any()
    assert agreement.spearman_published_pairs.between(-1,1).all()
    assert agreement.published_pair_mass_percent.between(0,100).all()
    assert agreement.total_variation_published_set.between(0,1).all()
    series_rows={}
    for method in provenance['methods']:
        x=pd.read_csv(out/f'time_series_{method}.csv.gz')
        assert not x.duplicated(['Barcode','sample']).any()
        assert x.groupby('Barcode')['sample'].nunique().eq(8).all()
        assert not x[['Barcode','counts','sample','replicate','generation','method']].isna().any().any()
        assert (x.counts>=0).all() and np.equal(x.counts,np.floor(x.counts)).all()
        for name,frame in x.groupby('sample'):
            a=agreement[agreement.method.eq(method)&agreement['sample'].eq(name)].iloc[0]
            assert frame.counts.sum()==a.assigned_molecules
            assert a.input_molecules==a.assigned_molecules+a.unassigned_molecules
            assert np.isclose(frame.frequency_assigned.sum(),1)
            assert np.allclose(frame.frequency_assigned,frame.counts/a.assigned_molecules)
            reference=frame[frame.is_published_pair]
            assert len(reference)==2314 and not reference.published_counts.isna().any()
            published_sample=publication[publication.Replicate.eq(frame.replicate.iloc[0]) &
                                         publication.Time.eq(frame.generation.iloc[0])].set_index('Barcode').Count
            assert np.array_equal(reference.published_counts,reference.Barcode.map(published_sample))
            assert reference.counts.sum()==a.published_pair_molecules
            assert np.isclose(reference.frequency_published_set.sum(),1)
            assert np.isclose(reference.published_frequency.sum(),1)
            # Independent rank/correlation formulation, rather than scipy's call.
            rho=reference.counts.rank().corr(reference.published_counts.rank())
            assert np.isclose(rho,a.spearman_published_pairs,atol=1e-12)
            tv=np.abs(reference.counts/reference.counts.sum()-reference.published_counts/reference.published_counts.sum()).sum()/2
            assert np.isclose(tv,a.total_variation_published_set,atol=1e-12)
        series_rows[method]=len(x)
    timings=pd.read_csv(out/'clustering_runs.csv')
    assert len(timings)==6*2*provenance['repeats']
    assert not timings.duplicated(['method','component','repeat']).any()
    assert (timings.workflow_seconds>0).all() and (timings.max_rss_bytes>0).all()
    summary=pd.read_csv(out/'method_summary.csv').set_index('method')
    recomputed=timings.groupby(['method','repeat']).workflow_seconds.sum().groupby('method').median()
    assert np.allclose(summary.loc[recomputed.index].clustering_seconds_median,recomputed)
    stability=pd.read_csv(out/'repeat_stability.csv')
    assert len(stability)==36 and stability.identical_membership.all()
    categories=pd.read_csv(out/'published_identity_categories.csv')
    for (method,sample),frame in categories.groupby(['method','sample']):
        a=agreement[agreement.method.eq(method)&agreement['sample'].eq(sample)].iloc[0]
        assert len(frame)==5 and frame.molecules.sum()==a.assigned_molecules
        assert frame.loc[frame.category.eq('published_pair'),'molecules'].sum()==a.published_pair_molecules
        assert np.allclose(frame.percent_input,100*frame.molecules/a.input_molecules)
    diagnosis=json.loads((out/'unmatched_diagnosis.json').read_text())
    for filename,expected in diagnosis['inputs_sha256'].items():assert sha256(out/filename)==expected
    assert sha256(HERE/'diagnose_unmatched.py')==diagnosis['script_sha256']
    result=dict(status='numerical_checks_passed',assessment='Share with stated experimental-reference caveats',
        input_pairs=int(extraction.total_pairs.sum()),input_molecules=int(extraction.barbac_input_molecules.sum()),
        full_FASTQ_files_MD5_verified=16,FastQC_record_counts_reconciled=True,
        inline_indexes_all_expected=True,published_counts_checked_against_pinned_source=True,
        unmatched_component_categories_reconciled=True,
        sample_method_comparisons=48,series_rows=series_rows,
        counts_and_normalization_reconciled=True,rank_correlations_independently_recomputed=True,
        three_repeat_timing_summary_recomputed=True,
        repeat_membership_all_identical=bool(stability.identical_membership.all()),
        visual_inspection='pending manual inspection of exported PNG figures',
        caveats=['Publication is a differently processed reference, not true barcode labels.',
                 'Pooled retrospective clustering uses all four generations.',
                 'No fitness inference; four measured timepoints per biological replicate.',
                 'First repeat supplies calls; timing medians use all three repeats.'],
        artifact_sha256={p.name:sha256(p) for p in out.iterdir() if p.is_file() and p.name!='validation.json'})
    (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='artifact_sha256'},indent=2))


if __name__=='__main__':main()
