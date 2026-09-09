"""Pool each component, compare six methods serially, and reconstruct each sample.

Distance 3 and the existing benchmark settings are frozen before inspecting
publication agreement. Three fresh timing repeats; first repeat supplies calls.
"""
import argparse
from collections import Counter
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from extract_barcodes import sha256

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
METHODS=['barbac_lv','barbac_hamming','shepherd','starcode_sphere','starcode_mp','bartender']
COMPONENTS=['diverse','environment']


def measured(command, work, log, extra_env=None):
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
         'VECLIB_MAXIMUM_THREADS':'1','PYTHONHASHSEED':'0',**(extra_env or {})}
    started=time.perf_counter()
    with open(log,'w') as stream:
        with subprocess.Popen(list(map(str,command)),cwd=work,stdout=stream,stderr=subprocess.STDOUT,env=env) as child:
            _,status,usage=os.wait4(child.pid,0)
            child.returncode=os.waitstatus_to_exitcode(status)
            if child.returncode: raise subprocess.CalledProcessError(child.returncode,command)
    result=dict(process_seconds=time.perf_counter()-started,
        max_rss_bytes=int(usage.ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
        user_seconds=usage.ru_utime,system_seconds=usage.ru_stime,
        memory_definition='wait4 reported maximum child RSS; not summed concurrent process-tree RSS',
        command=list(map(str,command)))
    Path(log).with_suffix('.process.json').write_text(json.dumps(result,indent=2)+'\n')
    return result['process_seconds']


def check_members(dest, inputs):
    centroids=pd.read_csv(dest/'centroids.csv')
    members=pd.read_csv(dest/'members.csv')
    assert centroids.central_barcode.is_unique and members.member.is_unique
    assert not centroids.isna().any().any() and not members.isna().any().any()
    assert members.member.isin(inputs.barcode).all()
    assert members.central_barcode.isin(centroids.central_barcode).all()
    counts=inputs.set_index('barcode').counts
    totals=members.assign(counts=members.member.map(counts)).groupby('central_barcode').counts.sum().sort_index()
    pd.testing.assert_series_equal(totals,centroids.set_index('central_barcode').sum_counts.sort_index(),
                                   check_names=False,check_dtype=False)
    return dict(zip(members.member,members.central_barcode))


def run_one(method, args, work, inputs, peer):
    dest=work/method;dest.mkdir(exist_ok=True)
    receipt=dest/'run.json'
    if receipt.exists():
        saved=json.loads(receipt.read_text())
        for name,expected in saved['output_sha256'].items(): assert sha256(dest/name)==expected
        check_members(dest,inputs)
        return saved
    started=time.perf_counter()
    core=build=None
    if method.startswith('barbac_'):
        metric='lv' if method=='barbac_lv' else 'hamming'
        command=['Rscript',ROOT/'benchmark/latest_four_conditions/run_barbac.R',args.library,
            work/'input.csv',dest,metric,'3','poisson' if metric=='lv' else 'none']
        wall=measured(command,dest,dest/'run.log')
        for line in (dest/'run.log').read_text().splitlines():
            if line.startswith('BARBAC_ALGO_SECONDS='):core=float(line.split('=')[1])
            if line.startswith('BARBAC_BUILD_ID='):build=line.split('=')[1]
        assert build=='barbac-2026-09-08-bounded-indels-v13' and core is not None
    elif method=='shepherd':
        command=[sys.executable,args.tools/'Shepherd/shepherd_t0.py','-f',work/'input.tsv','-l','26','-eps','3','-e','0.005']
        wall=measured(command,dest,dest/'run.log')
        c=pd.read_csv(work/'input_pb_freq.csv');c.columns=['central_barcode','sum_counts']
        raw=pd.read_csv(work/'input_seq_clust.csv')
        roots=raw[raw.sequence.isin(c.central_barcode)]
        assert roots.cluster.is_unique and len(roots)==len(c)
        members=pd.DataFrame({'member':raw.sequence,
            'central_barcode':raw.cluster.map(roots.set_index('cluster').sequence)}).dropna()
        c.to_csv(dest/'centroids.csv',index=False);members.to_csv(dest/'members.csv',index=False)
    else:
        peer.execute=measured
        result=peer.run_method(method,SimpleNamespace(tools=args.tools,resume=True),work,inputs)
        result['max_rss_bytes']=json.loads((dest/'run.process.json').read_text())['max_rss_bytes']
        receipt.write_text(json.dumps(result,indent=2)+'\n')
        check_members(dest,inputs)
        return result
    result=dict(method=method,process_seconds=wall,core_seconds=core,build_id=build,
        pipeline_seconds=time.perf_counter()-started,
        max_rss_bytes=json.loads((dest/'run.process.json').read_text())['max_rss_bytes'],
        command=list(map(str,command)),output_sha256={name:sha256(dest/name) for name in ['centroids.csv','members.csv']})
    check_members(dest,inputs)
    receipt.write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--tools',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=3)
    args=parser.parse_args()
    args.work_dir=args.work_dir.resolve();args.library=args.library.resolve();args.tools=args.tools.resolve()
    assert args.repeats>=1
    work=args.work_dir;out=HERE/'results';out.mkdir(exist_ok=True)
    manifest=pd.read_csv(HERE/'samples.tsv',sep='\t')
    samples={};receipts={};pools={c:Counter() for c in COMPONENTS}
    for row in manifest.to_dict('records'):
        path=work/'extracted'/row['sample']
        receipt=json.loads((path/'extraction.json').read_text())
        assert receipt['status']=='complete' and receipt['profile']=='chen2023_mapped_PEAR_flanks_24_28_Q30_first_UMI'
        for name,expected in receipt['outputs'].items(): assert sha256(path/name)==expected
        pairs=pd.read_csv(path/'pairs_barbac.csv')
        assert not pairs.duplicated(['diverse_barcode','environment_barcode']).any()
        assert (pairs.counts>0).all() and int(pairs.counts.sum())==receipt['barbac_input_molecules']
        samples[row['sample']]=pairs;receipts[row['sample']]=receipt
        for component in COMPONENTS:
            pools[component].update(pairs.groupby(component+'_barcode').counts.sum().to_dict())
    inputs={c:pd.DataFrame(sorted(pools[c].items(),key=lambda x:(-x[1],x[0])),columns=['barcode','counts']) for c in COMPONENTS}
    spec=importlib.util.spec_from_file_location('reference_peer_runner',ROOT/'benchmark/reference_comparison/run_comparison.py')
    peer=importlib.util.module_from_spec(spec);spec.loader.exec_module(peer)
    timings=[];first_maps={};repeat_equal=[]
    for repeat in range(1,args.repeats+1):
        order=METHODS[(repeat-1)*2:]+METHODS[:(repeat-1)*2]
        for component in COMPONENTS:
            dest=work/'clustering'/f'repeat_{repeat}'/component;dest.mkdir(parents=True,exist_ok=True)
            for name,options in [('input.csv',{}),('input.tsv',dict(sep='\t',header=False))]:
                file=dest/name
                if file.exists():
                    previous=pd.read_csv(file,**(dict(sep='\t',header=None,names=['barcode','counts']) if name.endswith('tsv') else {}))
                    pd.testing.assert_frame_equal(previous,inputs[component])
                else:inputs[component].to_csv(file,index=False,**options)
            for method in order:
                print('START',repeat,component,method,flush=True)
                result=run_one(method,args,dest,inputs[component],peer)
                mapping=check_members(dest/method,inputs[component])
                key=(method,component)
                if repeat==1:first_maps[key]=mapping
                same=mapping==first_maps[key]
                repeat_equal.append(dict(repeat=repeat,component=component,method=method,identical_membership=same))
                timings.append(dict(repeat=repeat,component=component,method=method,
                    workflow_seconds=result['pipeline_seconds'],process_seconds=result['process_seconds'],
                    core_seconds=result['core_seconds'],max_rss_bytes=result['max_rss_bytes'],
                    centroids=len(set(mapping.values())),input_unique=len(inputs[component]),
                    assigned_molecules=sum(pools[component][x] for x in mapping)))
                pd.DataFrame(timings).to_csv(out/'clustering_runs.csv',index=False)
                print('DONE',round(result['pipeline_seconds'],3),'seconds',len(set(mapping.values())),'centroids',flush=True)
    source=HERE/'generated/sources/hBFA1_all_freqs_tidy.csv'
    expected=next(x['sha256'] for x in json.loads((HERE/'sources.json').read_text())['files'] if x['filename']==source.name)
    assert sha256(source)==expected
    published=pd.read_csv(source)
    published=published[published.Test_Environment.eq('YPD') & published.Time.isin([8,16,24,40])].copy()
    assert len(published)==18512 and published.Barcode.nunique()==2314
    assert not published.duplicated(['Barcode','Replicate','Time']).any()
    assert (published.Count>=0).all()
    universe=sorted(published.Barcode.unique());universe_set=set(universe)
    method_universes={method:set(universe) for method in METHODS}
    for method in METHODS:
        for pairs in samples.values():
            d=pairs.diverse_barcode.map(first_maps[(method,'diverse')])
            e=pairs.environment_barcode.map(first_maps[(method,'environment')])
            valid=d.notna() & e.notna()
            method_universes[method].update(d[valid]+'_'+e[valid])
    rows=[];series=[];extraction_rows=[]
    for sample in manifest.to_dict('records'):
        name=sample['sample'];r=receipts[name];pairs=samples[name]
        pub=published[published.Replicate.eq(sample['replicate']) & published.Time.eq(sample['generation'])].set_index('Barcode').Count.reindex(universe)
        assert len(pub)==2314 and not pub.isna().any() and pub.sum()>0
        extraction_rows.append(dict(sample=name,replicate=sample['replicate'],generation=sample['generation'],
            **r['stats'],barbac_input_molecules=r['barbac_input_molecules'],
            bam_component_seconds=r['bam_extraction']['seconds'],paired_join_seconds=r['seconds']))
        for method in METHODS:
            d=pairs.diverse_barcode.map(first_maps[(method,'diverse')]);e=pairs.environment_barcode.map(first_maps[(method,'environment')])
            valid=d.notna() & e.notna()
            grouped=pairs.loc[valid].assign(Barcode=d[valid]+'_'+e[valid]).groupby('Barcode').counts.sum()
            assigned=int(grouped.sum());total=int(pairs.counts.sum());unassigned=total-assigned
            assert assigned+unassigned==total
            reference_counts=grouped.reindex(universe,fill_value=0)
            shared_mass=int(reference_counts.sum())
            assert shared_mass>0 and assigned>0
            positive=pub>0
            rows.append(dict(sample=name,replicate=sample['replicate'],generation=sample['generation'],method=method,
                input_molecules=total,assigned_molecules=assigned,unassigned_molecules=unassigned,
                observed_pair_centroids=len(grouped),published_pair_molecules=shared_mass,
                published_pair_mass_percent=100*shared_mass/total,
                unmatched_assigned_molecules=assigned-shared_mass,
                published_positive_pairs=int(positive.sum()),
                published_positive_pairs_observed=int(((reference_counts>0)&positive).sum()),
                spearman_published_pairs=float(spearmanr(reference_counts,pub).statistic),
                published_count_total=int(pub.sum()),
                total_variation_published_set=float(0.5*np.abs(reference_counts/shared_mass-pub/pub.sum()).sum())))
            all_pairs=sorted(method_universes[method])
            counts=grouped.reindex(all_pairs,fill_value=0).astype('int64')
            frame=pd.DataFrame(dict(Barcode=all_pairs,counts=counts.values))
            frame['method']=method;frame['sample']=name;frame['replicate']=sample['replicate'];frame['generation']=sample['generation']
            frame['input_molecules']=total;frame['assigned_molecules']=assigned
            frame['frequency_assigned']=frame.counts/assigned
            frame['is_published_pair']=frame.Barcode.isin(universe_set)
            frame['published_set_molecules']=shared_mass
            frame['frequency_published_set']=np.where(frame.is_published_pair,frame.counts/shared_mass,np.nan)
            frame['published_counts']=frame.Barcode.map(pub)
            frame['published_frequency']=frame.published_counts/pub.sum()
            series.append(frame)
    agreement=pd.DataFrame(rows);agreement.to_csv(out/'sample_method_agreement.csv',index=False)
    pd.DataFrame(extraction_rows).to_csv(out/'extraction_summary.csv',index=False)
    trajectories=pd.concat(series,ignore_index=True)
    for method in METHODS:
        sub=trajectories[trajectories.method.eq(method)]
        # Explicit zero observations for pairs seen at another time in this method.
        # Published pairs already have a complete eight-sample grid above.
        sub.to_csv(out/f'time_series_{method}.csv.gz',index=False,compression=dict(method='gzip',mtime=0))
    total_times=pd.DataFrame(timings).groupby(['method','repeat']).agg(
        workflow_seconds=('workflow_seconds','sum'),max_rss_bytes=('max_rss_bytes','max')).reset_index()
    total_times.to_csv(out/'combined_clustering_times.csv',index=False)
    summary=agreement.groupby('method').agg(spearman_median=('spearman_published_pairs','median'),
        spearman_min=('spearman_published_pairs','min'),spearman_max=('spearman_published_pairs','max'),
        published_mass_percent_median=('published_pair_mass_percent','median'),
        total_variation_median=('total_variation_published_set','median'),
        unassigned_molecules=('unassigned_molecules','sum'))
    summary=summary.join(total_times.groupby('method').agg(clustering_seconds_median=('workflow_seconds','median'),
        clustering_seconds_min=('workflow_seconds','min'),clustering_seconds_max=('workflow_seconds','max'),
        reported_peak_RSS_MiB=('max_rss_bytes',lambda x:x.max()/1024**2)))
    summary.to_csv(out/'method_summary.csv')
    top=published.groupby('Barcode').Count.sum().sort_values(ascending=False).head(12).index
    selected=trajectories[trajectories.Barcode.isin(top)].copy()
    selected['display_id']=selected.Barcode.map({bc:f'L{i+1:02}' for i,bc in enumerate(top)})
    selected.to_csv(out/'selected_trajectories.csv',index=False)
    pd.DataFrame(repeat_equal).to_csv(out/'repeat_stability.csv',index=False)
    provenance=dict(status='complete',samples=8,generations=[8,16,24,40],replicates=['R1','R2'],
        input_molecules=sum(r['barbac_input_molecules'] for r in receipts.values()),
        methods=METHODS,distance=3,repeats=args.repeats,barbac_tie_break='support',barbac_merge_ratio=20,
        barbac_error_rate=0.005,barbac_LV_indel_model='poisson',barbac_Hamming_indel_model='none',barbac_use_design=False,
        shepherd_error_rate=0.005,
        shepherd_configuration_note='Automatic error-rate estimation failed on BC1. Both components use the documented -e 0.005 option, matching the pre-existing barbac configured rate; no publication agreement was used to choose it. Automatic attempts are retained outside the timed comparison.',
        publication_used_for_clustering=False,publication_sha256=sha256(source),
        membership_source='First timing repeat; subsequent repeats measure stability and runtime',
        timing_scope='Sum of two separately pooled component workflows: process startup, clustering, required format conversion and centroid/member exports; excludes shared extraction, staging, trajectory reconstruction and evaluation. Serial after preprocessing.',
        normalization='frequency_assigned uses all assigned molecules; frequency_published_set conditions on the same 2314 published barcode pairs in every method/sample. Unmatched and unassigned mass reported separately.',
        interpretation='Publication agreement, not known-truth accuracy; pooled retrospective clustering uses all four generations.',
        files={str(p.relative_to(ROOT)):sha256(p) for p in [Path(__file__),HERE/'extract_mapped_pairs.py',
            HERE/'extract_mapped_components.R',HERE/'extract_barcodes.py',HERE/'map_sample.R',
            HERE/'mapped_workflow.py',HERE/'samples.tsv',HERE/'reference/chen2023_masked_amplicon.fasta',
            ROOT/'R/10_barbac_xtr.R',ROOT/'R/10a_extract_flanked_bam.R',ROOT/'R/09_run_cli_pipeline.R',
            ROOT/'R/11_super_cluster2.R',ROOT/'src/clustering.cpp',
            ROOT/'benchmark/latest_four_conditions/run_barbac.R',ROOT/'benchmark/reference_comparison/run_comparison.py']},
        extraction_inputs={name:sha256(work/'extracted'/name/'pairs_barbac.csv') for name in samples},
        tool_sha256={name:sha256(args.tools/name) for name in ['Shepherd/shepherd_t0.py','starcode/starcode','bartender-1.1/bartender_single_com','bartender-1.1/bartender_single']},
        python=sys.version,platform=platform.platform())
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(summary.to_string(),flush=True)


if __name__=='__main__':main()
