"""Extract the eight mapped samples locally with bounded concurrency."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--wait-for-mapping',action='store_true')
    args=parser.parse_args();work=args.work_dir.resolve()
    assert 1<=args.jobs<=4
    with (HERE/'samples.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
    pending=set(range(len(rows)));active={};done=[];started=time.monotonic()
    (work/'logs').mkdir(parents=True,exist_ok=True)
    def task(index):
        with (work/'logs'/f'extract_{index}.log').open('w') as log:
            subprocess.run([sys.executable,str(HERE/'run_sample.py'),'extract','--index',str(index),
                            '--work-dir',str(work)],stdout=log,stderr=subprocess.STDOUT,check=True)
        return rows[index]['sample']
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        while pending or active:
            for future,index in list(active.items()):
                if future.done():
                    sample=future.result();done.append(sample);del active[future]
                    print('COMPLETE',sample,flush=True)
            for index in sorted(pending):
                if len(active)>=args.jobs:break
                marker=work/'mapping'/rows[index]['sample']/'mapping.json'
                if marker.exists():
                    print('EXTRACT',rows[index]['sample'],flush=True)
                    active[pool.submit(task,index)]=index;pending.remove(index)
                elif not args.wait_for_mapping:raise FileNotFoundError(marker)
            if time.monotonic()-started>7200:raise TimeoutError('Mapping/extraction did not finish within two hours')
            if pending or active:time.sleep(5)
    (work/'extraction_completion.json').write_text(json.dumps(dict(status='complete',samples=sorted(done)),indent=2)+'\n')


if __name__=='__main__':main()
