"""Retain original evidence and repair only an objectively selected sleep interval."""
import argparse
from datetime import datetime, timedelta, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from common import HERE, LIBRARY, TOOLS, ENVIRONMENT, load, save, sha, verify_outputs

START=datetime(2026,9,12,20,44,52,tzinfo=timezone.utc)-timedelta(seconds=1)
END=datetime(2026,9,12,21,8,13,tzinfo=timezone.utc)+timedelta(seconds=1)


def now():return datetime.now(timezone.utc).isoformat()


def key(row):
    return ('milos' if row['condition']=='milos' else f"{row['seed']}/{row['condition']}")+'/'+row['method']


def register():
    path=HERE/'timing_repair_protocol.json'
    assert not path.exists(), 'Existing timing repair registration must be retained'
    rows=load(HERE/'results.json')
    selected=[r for r in rows if r['status']=='complete' and
        datetime.fromisoformat(r['started_at'])<END and datetime.fromisoformat(r['finished_at'])>START]
    assert len(selected)==28
    schedule=load(HERE/'schedule.json')
    sources={d['data']['condition'] if d['data']['condition']=='milos' else f"{d['data']['seed']}/{d['data']['condition']}":d['data'] for d in schedule}
    power_log=subprocess.check_output(['pmset','-g','log'],text=True)
    relevant=[line for line in power_log.splitlines() if line.startswith('2026-09-12 ') and
        '22:44:52'<=line[11:19]<='23:08:13' and any(x in line for x in [' Sleep ',' DarkWake ',' Wake ',' ThermalEvent '])]
    assert any('Entering Sleep' in x for x in relevant) and any('Wake from Deep Idle' in x for x in relevant)
    cells=[]
    for row in selected:
        original=HERE/'generated/results'/key(row)
        data_key=key(row).rsplit('/',1)[0]
        cells.append(dict(key=key(row),condition=row['condition'],seed=row['seed'],method=row['method'],
            source=sources[data_key],original_receipt_sha256=sha(original/'result.json'),
            original_workflow_seconds=row['workflow_seconds'],original_output_sha256=row['output_sha256']))
    save(path,dict(registered_at=now(),reason='Sleep and partial wake; objective timestamp overlap, not result-based selection',
        window_start_utc=START.isoformat(),window_end_utc=END.isoformat(),timestamp_margin_seconds=1,
        original_campaign_attempted_at_registration=len(rows),cells=cells,n_cells=len(cells),
        order='Original relative order',maximum_new_invocations_per_selected_cell=1,
        source_sha256={name:sha(HERE/name) for name in ['timing_repair.py','TIMING_REPAIR.md','worker.py','freeze.json']},
        power_log_excerpt=relevant,accuracy_replaced=False,faster_measurement_selection=False))
    print('REGISTERED TIMING REPAIR',len(cells),flush=True)


def run():
    from run import check_freeze
    check_freeze()
    protocol=load(HERE/'timing_repair_protocol.json')
    for name,digest in protocol['source_sha256'].items():assert sha(HERE/name)==digest
    assert load(HERE/'execution_validation.json')['cells']==726
    settings=load(HERE/'final_protocol.json')['execution']
    (HERE/'generated/timing_repair').mkdir(exist_ok=True)
    with (HERE/'generated/campaign.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        records=[]
        for i,cell in enumerate(protocol['cells'],1):
            original=HERE/'generated/results'/cell['key']
            assert sha(original/'result.json')==cell['original_receipt_sha256']
            verify_outputs(original,cell['original_output_sha256'])
            source=Path(cell['source']['path'])
            verify_outputs(source,cell['source']['sha256'])
            dest=HERE/'generated/timing_repair'/cell['key']
            receipt=dest/'result.json'
            if receipt.exists():
                record=load(receipt)
                verify_outputs(dest,record.get('output_sha256',{}))
                assert record['original_receipt_sha256']==cell['original_receipt_sha256']
            else:
                assert shutil.disk_usage(HERE).free>settings['minimum_free_disk_gib']*2**30
                dest.mkdir(parents=True,exist_ok=False)
                shutil.copyfile(source/'input.tsv',dest/'input.tsv')
                command=[sys.executable,str(HERE/'worker.py'),cell['method'],str(source),str(dest),str(TOOLS),str(LIBRARY),str(cell['source']['nominal_length'])]
                print('TIMING REPAIR',i,'/',len(protocol['cells']),cell['key'],flush=True)
                before=now();begin=time.perf_counter()
                try:
                    with (dest/'worker.log').open('w') as stream:
                        subprocess.run(command,env={**os.environ,**ENVIRONMENT},stdout=stream,stderr=subprocess.STDOUT,
                            check=True,timeout=settings['timeout_seconds']+60)
                    elapsed=time.perf_counter()-begin
                    digests={name:sha(dest/name) for name in ['centroids.csv','members.csv','worker.json']}
                    assert all(digests[name]==cell['original_output_sha256'][name] for name in ['centroids.csv','members.csv'])
                    record=dict(status='complete',key=cell['key'],condition=cell['condition'],seed=cell['seed'],method=cell['method'],
                        started_at=before,finished_at=now(),workflow_seconds=elapsed,worker_command=command,
                        original_receipt_sha256=cell['original_receipt_sha256'],original_workflow_seconds=cell['original_workflow_seconds'],
                        accuracy_mapping_byte_equal=True,output_sha256=digests,**load(dest/'worker.json'))
                except (subprocess.CalledProcessError,subprocess.TimeoutExpired,AssertionError) as exc:
                    record=dict(status='failed',key=cell['key'],condition=cell['condition'],seed=cell['seed'],method=cell['method'],
                        started_at=before,finished_at=now(),error=str(exc),original_receipt_sha256=cell['original_receipt_sha256'],output_sha256={})
                save(receipt,record)
                if record['status']=='complete' and (dest/'expanded.csv').exists():(dest/'expanded.csv').unlink()
            records.append(record)
            save(HERE/'timing_repair_results.json',records)
        check_freeze()
        save(HERE/'timing_repair_validation.json',dict(cells=len(records),successful=sum(r['status']=='complete' for r in records),
            failed=sum(r['status']!='complete' for r in records),original_accuracy_and_receipts_preserved=True,
            original_source_fingerprints_preserved=True,completed_at=now()))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['register','run'])
    args=parser.parse_args()
    register() if args.action=='register' else run()
