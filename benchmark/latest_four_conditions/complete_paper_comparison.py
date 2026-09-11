"""Complete the distance-three paper table with both Starcode modes/Bartender.

Uses existing validated barbac, Shepherd and Starcode sphere observations.
Fresh Hamming on the large reference ensures both barbac modes report v13.
"""
import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pandas as pd
from run_comparison import HERE, ROOT, BASE, CONDITIONS, sha, score


def main():
    spec=importlib.util.spec_from_file_location('reference_runner',ROOT/'benchmark/reference_comparison/run_comparison.py')
    ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
    args=SimpleNamespace(tools=Path('/Users/theodosiou/Documents/Projects/Barcodes/barbac-benchmark/tools'),
                         library=Path('/private/tmp/barbac-lv-v13'),rscript='Rscript',resume=True)
    work=BASE/'generated/latest-four-v13/paper-completion'
    rows=[];receipts={}
    for seed in [42,43,44]:
        for condition in CONDITIONS:
            source=BASE/'generated/revision-accuracy'/str(seed)/condition
            cfg=json.loads((source/'simulation.json').read_text())
            for file,expected in cfg['sha256'].items(): assert sha(source/file)==expected
            dest=work/str(seed)/condition;dest.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source/'input.csv',dest/'input.csv')
            shutil.copyfile(source/'shepherd_input.txt',dest/'input.tsv')
            inputs,truth=pd.read_csv(dest/'input.csv'),pd.read_csv(source/'true_counts.csv')
            for method in ['starcode_mp','bartender']:
                print('START',seed,condition,method,flush=True)
                run=ref.run_method(method,args,dest,inputs)
                c=pd.read_csv(dest/method/'centroids.csv');m=pd.read_csv(dest/method/'members.csv')
                assert m.member.is_unique and not m.isna().any().any()
                assert m.member.isin(inputs.barcode).all() and m.central_barcode.isin(c.central_barcode).all()
                m['count']=m.member.map(inputs.set_index('barcode').counts)
                pd.testing.assert_series_equal(m.groupby('central_barcode')['count'].sum().sort_index(),
                                              c.set_index('central_barcode').sum_counts.sort_index(),check_names=False,check_dtype=False)
                metrics=score(c,inputs,truth)
                assert metrics['unassigned_reads']==0
                rows.append(dict(seed=seed,condition=condition,method=method,distance=3,
                                 origin='fresh_peer',workflow_seconds=run['pipeline_seconds'],
                                 core_seconds=run['core_seconds'],output_sha256=run['output_sha256']['centroids.csv'],**metrics))
                receipts[f'{seed}/{condition}/{method}']=run
                pd.DataFrame(rows).to_csv(HERE/'additional_peers.csv',index=False)
                print(f"DONE FN={metrics['fn']} FP={metrics['fp']} F1={100*metrics['f1']:.4f}% workflow={run['pipeline_seconds']:.2f}s",flush=True)
    source=BASE/'generated/reference_2026-09-08'
    original=json.loads((ROOT/'benchmark/reference_comparison/provenance.json').read_text())
    for file,expected in original['input_sha256'].items(): assert sha(source/file)==expected
    dest=work/'milos';dest.mkdir(parents=True,exist_ok=True)
    for file in ['input.csv','input.tsv']: shutil.copyfile(source/file,dest/file)
    inputs=pd.read_csv(dest/'input.csv')
    _,labels,truth,quality=ref.profile(Path('/Users/theodosiou/Documents/Projects/Barcodes/barbac-benchmark'))
    print('START Milos latest barbac Hamming',flush=True)
    run=ref.run_method('barbac_hamming_support',args,dest,inputs)
    assert run['build_id']=='barbac-2026-09-08-bounded-indels-v13'
    output=dest/'barbac_hamming_support'
    metrics=ref.evaluate(pd.read_csv(output/'centroids.csv'),pd.read_csv(output/'members.csv'),labels,truth)
    receipt=dict(**run,**metrics)
    (HERE/'milos_hamming.json').write_text(json.dumps(receipt,indent=2)+'\n')
    receipts['milos/barbac_hamming_support']=run
    print(f"DONE Milos Hamming FN={metrics['fn']} FP={metrics['fp']} F1={100*metrics['f1']:.5f}% workflow={run['pipeline_seconds']:.2f}s",flush=True)
    provenance=dict(status='complete',source_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'benchmark/reference_comparison/run_comparison.py',ROOT/'benchmark/reference_comparison/run_barbac.R']},
                    tool_sha256={name:sha(args.tools/name) for name in ['starcode/starcode','bartender-1.1/bartender_single_com','bartender-1.1/bartender_single']},
                    runs=receipts,protocol='Serial distance-three runs; one thread for Starcode/Bartender; workflow includes required format conversion and both centroid/member exports. Same original simulation inputs, no labels used in clustering.')
    (HERE/'completion_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__=='__main__': main()
