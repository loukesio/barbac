#!/usr/bin/env python3
"""Verify preselected local inputs, then measure one fresh R subprocess."""
import argparse, csv, datetime, hashlib, json, os, pathlib, platform, subprocess, time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
def sha(path):
    return hashlib.file_digest(open(path,'rb'),'sha256').hexdigest()
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input-root',type=pathlib.Path,required=True)
    parser.add_argument('--output',type=pathlib.Path,required=True)
    parser.add_argument('--library',type=pathlib.Path)
    args=parser.parse_args()
    if args.library is None:
        active=ROOT/'app/.runtime/active-library.txt'
        args.library=pathlib.Path(active.read_text().strip()) if active.exists() else ROOT/'app/.runtime/library'
    args.library=args.library.resolve()
    out=args.output.resolve(); out.mkdir(parents=True,exist_ok=False)
    rows=list(csv.DictReader(open(args.input_root/'samples.tsv'),delimiter='\t'))
    rows=[r for r in rows if r['well']=='A3' and int(r['passage']) in (2,4,6)]
    assert len(rows)==12
    inputs=[]
    for row in rows:
        path=args.input_root/'generated/fastq'/f"{row['run_accession']}.fastq.gz"
        assert hashlib.file_digest(open(path,'rb'),'md5').hexdigest()==row['fastq_md5']
        inputs.append(dict(accession=row['run_accession'],md5=row['fastq_md5'],
            sha256=sha(path),bytes=path.stat().st_size,reads=int(row['read_count']),
            sample=f"A3_p{int(row['passage']):02d}",generation=int(row['passage'])*6))
    sources=['DESCRIPTION','NAMESPACE','R/09_run_cli_pipeline.R','R/10_barbac_xtr.R',
             'R/11_super_cluster2.R','R/12_barbac_ts_area.R','src/clustering.cpp',
             'benchmark/workflow_runtime/run.R','benchmark/workflow_runtime/run.py']
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[key]='1'
    tool_dir=pathlib.Path.home()/'Library/r-miniconda/envs/barbac_env/bin'
    env['PATH']=str(tool_dir)+os.pathsep+env['PATH']
    cmd=['Rscript','--vanilla',str(HERE/'run.R'),str(args.input_root.resolve()),str(out),
         str(args.library)]
    receipt=dict(protocol='protocol.md',source_commit=subprocess.check_output(
        ['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),inputs=inputs,
        source_sha256={p:sha(ROOT/p) for p in sources},
        reference_sha256=sha(args.input_root/'reference/cassette.fasta'),
        hardware=dict(system=platform.platform(),cpu=subprocess.check_output(
            ['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),
            memory_bytes=int(subprocess.check_output(['sysctl','-n','hw.memsize'])),
            logical_cpus=os.cpu_count()),command=cmd,
        threads='Single-thread native clustering and BLAS; minimap2 sr default 3 threads; samtools and FastQC defaults; no simultaneous analysis jobs.',
        tool_sha256={p:sha(tool_dir/p) for p in ('minimap2','samtools','fastqc','multiqc')})
    receipt['started_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    start=time.perf_counter()
    with open(out/'run.log','w') as log:
        result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
    receipt.update(outer_elapsed_seconds=time.perf_counter()-start,exit_status=result.returncode,
        finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    receipt['artifacts_sha256']={p.name:sha(p) for p in out.iterdir()
        if p.is_file() and p.name!='receipt.json'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('outer_elapsed_seconds','exit_status','finished_utc')}))
    raise SystemExit(result.returncode)
if __name__=='__main__': main()
