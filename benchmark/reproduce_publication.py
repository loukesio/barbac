#!/usr/bin/env python3
"""Portable reproduction entry point; never rewrites archived campaign receipts."""
import argparse, datetime, importlib.util, json, os, pathlib, shutil, subprocess, sys, time

ROOT=pathlib.Path(__file__).resolve().parent.parent
FINAL=ROOT/'benchmark/publication_final'
METHODS=('hamming','lv','shepherd','starcode_sphere','starcode_default','bartender')
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
def worker(args):
    if args.method!='shepherd':
        module('archived_worker',FINAL/'worker.py').run(args.method,args.source,args.output,
            args.tools,args.library,args.length)
        return
    # Same completed Shepherd configuration and export rules as its archived worker.
    import pandas as pd
    command=[sys.executable,str(args.tools/'Shepherd/shepherd_t0.py'),'-f',
        str(args.output/'input.tsv'),'-l',str(args.length),'-eps','3','-bft','-4']
    if args.condition!='milos':command+=['-e','0.004']
    start=time.perf_counter()
    with open(args.output/'tool.log','w') as log:
        subprocess.run(command,cwd=args.output,stdout=log,stderr=subprocess.STDOUT,
            check=True,timeout=900)
    seconds=time.perf_counter()-start
    counts=pd.read_csv(args.source/'input.csv').set_index('barcode').counts
    centroids=pd.read_csv(args.output/'input_pb_freq.csv')
    centroids.columns=['central_barcode','sum_counts']
    raw=pd.read_csv(args.output/'input_seq_clust.csv')
    own=raw[raw.sequence.isin(centroids.central_barcode)]
    assert own.cluster.is_unique and len(own)==len(centroids)
    roots=own.set_index('cluster').sequence
    members=pd.DataFrame({'member':raw.sequence,'central_barcode':raw.cluster.map(roots)}).dropna()
    members['member_count']=members.member.map(counts)
    assert not members.isna().any().any()
    centroids.to_csv(args.output/'centroids.csv',index=False)
    members.to_csv(args.output/'members.csv',index=False)
    (args.output/'worker.json').write_text(json.dumps(dict(command=command,tool_seconds=seconds),indent=2)+'\n')
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=pathlib.Path,required=True,help='New directory; existing paths are refused')
    p.add_argument('--seed',type=int,default=2026091202)
    p.add_argument('--condition',choices=('random_mixed','anchored_mixed','milos'),default='random_mixed')
    p.add_argument('--source',type=pathlib.Path,help='Existing canonical input.csv, truth.csv and labels.csv (required for Milo)')
    p.add_argument('--tools',type=pathlib.Path,default=pathlib.Path('tools'))
    p.add_argument('--library',type=pathlib.Path,default=None)
    p.add_argument('--methods',nargs='*',choices=METHODS,default=[],help='Omit to generate/validate inputs only')
    p.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--method',choices=METHODS,help=argparse.SUPPRESS)
    p.add_argument('--length',type=int,help=argparse.SUPPRESS)
    args=p.parse_args()
    if args.library is None:
        active=ROOT/'app/.runtime/active-library.txt'
        args.library=pathlib.Path(active.read_text().strip()) if active.exists() else ROOT/'app/.runtime/library'
    for key in ('output','source','tools','library'):
        val=getattr(args,key)
        if val is not None:setattr(args,key,val.resolve())
    if args.worker:return worker(args)
    if args.output.exists():p.error('Choose a new --output directory to preserve existing results.')
    if args.condition=='milos' and args.source is None:p.error('Milo requires --source with the deposited data in canonical format.')
    if any(m in ('lv','hamming') for m in args.methods) and not (args.library/'barbac').is_dir():
        p.error('Install the source release into --library before running Barbac.')
    args.output.mkdir(parents=True)
    simulation=module('publication_simulator',FINAL/'simulate.py')
    if args.source is None:
        args.source=args.output/'inputs'
        config=json.loads((FINAL/'protocol.json').read_text())
        simulation.generate_condition(config,args.seed,args.condition,args.source)
    import pandas as pd
    metric=module('publication_metrics',FINAL/'metrics.py')
    inputs=pd.read_csv(args.source/'input.csv');truth=pd.read_csv(args.source/'truth.csv')
    labels=pd.read_csv(args.source/'labels.csv')
    input_validation=metric.validate_data(inputs,truth,labels)
    length=26 if args.condition=='anchored_mixed' else 20
    env=os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[name]='1'
    env['PYTHONHASHSEED']='0'
    records=[]
    for method in args.methods:
        dest=args.output/method;dest.mkdir()
        inputs.to_csv(dest/'input.tsv',sep='\t',header=False,index=False)
        command=[sys.executable,str(pathlib.Path(__file__).resolve()),'--worker','--method',method,
          '--source',str(args.source),'--output',str(dest),'--tools',str(args.tools),
          '--library',str(args.library),'--length',str(length),'--condition',args.condition]
        start=time.perf_counter()
        with open(dest/'execution.log','w') as log:
            result=subprocess.run(command,env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        record=dict(method=method,exit_status=result.returncode,elapsed_seconds=time.perf_counter()-start)
        if result.returncode==0:
            record['metrics']=metric.score(pd.read_csv(dest/'centroids.csv'),pd.read_csv(dest/'members.csv'),inputs,truth,labels)
        else:record['status']='failed; retained without retry'
        records.append(record)
    receipt=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
      condition=args.condition,seed=args.seed if args.condition!='milos' else None,
      source=str(args.source),source_sha256={n:simulation.sha(args.source/n) for n in ('input.csv','truth.csv','labels.csv')},
      source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      input_validation=input_validation,results=records,
      scope='Independent reproduction; new machine/session/worker timings are not replacements for Table 1.')
    (args.output/'reproduction.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(args.output/'reproduction.json')
if __name__=='__main__':main()
