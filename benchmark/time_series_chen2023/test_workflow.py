"""Check scheduler arguments using a local sbatch stand-in; submits no jobs."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent


class SubmissionTests(unittest.TestCase):
    def check_submission(self, configured):
        with tempfile.TemporaryDirectory(prefix='barbac-slurm-test-') as tmp:
            base = Path(tmp)
            stub = base/'sbatch'
            stub.write_text('#!/usr/bin/env python3\nimport json, os, sys\n'
                'with open(os.environ["BARBAC_TEST_LOG"], "a") as h: '
                'h.write(json.dumps(sys.argv[1:])+"\\n")\nprint("12345;local")\n')
            stub.chmod(0o755)
            config = base/'config.sh'
            config.write_text('BARBAC_REPO='+shlex.quote(str(HERE.parents[1]))+'\n'
                'BARBAC_WORK='+shlex.quote(str(base/'work with spaces'))+'\n'
                'BARBAC_ARRAY_LIMIT=2\n'+(
                'BARBAC_ACCOUNT=test\nBARBAC_PARTITION=compute\n'
                'BARBAC_DOWNLOAD_PARTITION=transfer\nBARBAC_QOS=normal\n' if configured else ''))
            logfile = base/'calls.jsonl'
            env = dict(os.environ, PATH=str(base)+os.pathsep+os.environ['PATH'], BARBAC_TEST_LOG=str(logfile))
            for mode in ('pilot', 'full'):
                subprocess.run(['bash', str(HERE/'submit.sh'), str(config), mode], env=env,
                               check=True, capture_output=True, text=True)
            calls = [json.loads(line) for line in logfile.read_text().splitlines()]
            self.assertEqual(len(calls), 4)
            self.assertIn('--array=0', calls[0])
            self.assertIn('--array=0-7%2', calls[2])
            self.assertIn('--dependency=afterok:12345', calls[1])
            self.assertIn('--dependency=afterok:12345', calls[3])
            self.assertEqual(calls[0][-2:], ['download', '100000'])
            self.assertEqual(calls[3][-2:], ['extract', '0'])
            if configured:
                self.assertIn('--partition=transfer', calls[0])
                self.assertIn('--partition=compute', calls[1])
                for call in calls:
                    self.assertIn('--account=test', call)
                    self.assertIn('--qos=normal', call)

    def test_explicit_scheduler_settings(self):
        self.check_submission(True)

    def test_site_defaults(self):
        self.check_submission(False)


if __name__ == '__main__':
    unittest.main()
