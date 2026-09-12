"""Recorded launcher repair: supply run.py's omitted standard-library json import.

The original frozen script is retained unchanged. This adapter is fingerprinted
before the first tool measurement; it does not alter generation, tool settings,
scoring, inference, scheduling or existing data.
"""
import datetime
import fcntl
import json
import os

import run
from common import HERE, sha, save, load


def main():
    record_path = HERE / 'execution_adapter.json'
    record = dict(adapter_sha256=sha(__file__), original_run_sha256=sha(HERE/'run.py'),
        original_freeze_sha256=sha(HERE/'freeze.json'),
        generation_validation_sha256=sha(HERE/'generation_validation.json'),
        correction='Supply missing json module binding before execute(); original frozen source retained unchanged.')
    if record_path.exists():
        previous = load(record_path)
        assert all(previous[k] == v for k, v in record.items())
    else:
        assert not (HERE/'generated/results').exists(), 'Adapter must be registered before any tool cells'
        save(record_path, dict(registered_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), **record))
    run.json = json
    with (HERE/'generated/campaign.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        lock.write(str(os.getpid()))
        lock.flush()
        run.execute()
    assert sha(__file__) == record['adapter_sha256']


if __name__ == '__main__':
    main()
