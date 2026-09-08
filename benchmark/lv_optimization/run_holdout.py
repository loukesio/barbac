"""Fresh-seed validation plus a deliberately homopolymer-enriched stress test."""
import json
from pathlib import Path
import sys
import pandas as pd
from run_experiment import run, validate_counts, sha, ROOT, HERE
sys.path.insert(0, str(ROOT/'benchmark/indel_experiment'))
from simulate import simulate


def main():
    work=ROOT/'benchmark/four_condition_comparison/generated/lv-v13-holdout'
    library=Path('/private/tmp/barbac-lv-v13')
    conditions=[
        ('random_substitutions',None,0.0,10000),
        ('random_low_indels',None,0.005,10000),
        ('anchored_substitutions','NNNNNNNNATGCNNNNNNNNATCGTTAA',0.0,10000),
        ('anchored_low_indels','NNNNNNNNATGCNNNNNNNNATCGTTAA',0.005,10000),
        ('homopolymer_stress','NNNNNNNNNNCCCCCCCCCCNNNNNNNNNN',0.005,1000),
    ]
    rows=[]; metadata=[]
    for name,template,indel,n in conditions:
        source=work/name
        config=dict(n_barcodes=n,barcode_length=len(template) if template else 20,
                    n_reads=1000000,sub_rate=0.005,ins_rate=indel,del_rate=indel,
                    sigma=1.5,seed=20260909,template=template,abundance='lognormal',tie_order='sequence')
        print('SIMULATE',name,flush=True)
        simulate(source,**config)
        metadata.append(dict(condition=name,configuration=config,input_sha256=sha(source/'input.csv'),truth_sha256=sha(source/'true_counts.csv')))
        truth=set(pd.read_csv(source/'true_counts.csv').iloc[:,0])
        for model in ['none','poisson']:
            output=source/model
            print('START',name,model,flush=True)
            result=run(library,source/'input.csv',output,'support',model)
            c=validate_counts(source/'input.csv',output)
            found=set(c.central_barcode);tp=len(truth&found);fn=len(truth-found);fp=len(found-truth)
            row=dict(condition=name,seed=20260909,model=model,fn=fn,fp=fp,f1=2*tp/(2*tp+fn+fp),
                     core_seconds=result['core_seconds'],process_seconds=result['process_seconds'],output_reads=int(c.sum_counts.sum()))
            rows.append(row)
            (output/'result.json').write_text(json.dumps({**result,**row},indent=2)+'\n')
            pd.DataFrame(rows).to_csv(HERE/'holdout_results.csv',index=False)
            print(json.dumps(row),flush=True)
    metadata=dict(datasets=metadata,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in [HERE/'run_holdout.py',ROOT/'benchmark/indel_experiment/simulate.py',ROOT/'src/clustering.cpp',ROOT/'R/11_super_cluster2.R']},
                  library_sha256={str(p.relative_to(library)):sha(p) for p in (library/'barbac').rglob('*') if p.is_file()},
                  scope='Four original categories with a fresh seed fixed before results, plus a synthetic homopolymer-enriched stress test; the latter is not representative of every dataset.')
    (HERE/'holdout_manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')

if __name__=='__main__': main()
