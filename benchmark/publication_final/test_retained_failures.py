"""Missing tool outputs must never become successful scores or favourable tests."""
import unittest

import pandas as pd

from analyze import paired_summary
from analyze_available import compare, summarize, METRICS


class RetainedFailureTests(unittest.TestCase):
    def setUp(self):
        self.protocol=dict(n_independent_libraries_per_design=4,final_seeds=[1,2,3,4],
            designs=['random_mixed','anchored_mixed'],
            methods=['hamming','lv','shepherd','starcode_sphere','starcode_mp','bartender'])
        self.data=pd.DataFrame([dict(seed=s,condition=c,method=m,status='complete',
            f1_percent=99+(s*.001 if m=='lv' else 0),workflow_seconds=1+s/10)
            for s in self.protocol['final_seeds'] for c in self.protocol['designs']
            for m in self.protocol['methods']])

    def test_failure_does_not_change_other_full_contrasts(self):
        before,_=compare(self.data,self.protocol,resamples=128)
        target=(self.data.seed==2)&(self.data.condition=='random_mixed')&(self.data.method=='shepherd')
        self.data.loc[target,'status']='failed'
        self.data.loc[target,['f1_percent','workflow_seconds']]=float('nan')
        after,speeds=compare(self.data,self.protocol,resamples=128)
        self.assertEqual(after[0]['n'],3)
        self.assertEqual(after[0]['missing_seeds'],[2])
        self.assertIsNone(after[0]['mean_difference'])
        self.assertIsNone(after[0]['holm_adjusted_p'])
        self.assertFalse(after[0]['superiority_supported'])
        self.assertIsNone(speeds[0]['geometric_time_ratio'])
        for old,new in zip(before[1:],after[1:]):
            for key in ['mean_difference','simultaneous_one_sided_lower','bootstrap_lower','one_sided_p']:
                self.assertEqual(old[key],new[key])
            self.assertGreaterEqual(new['holm_adjusted_p'],old['holm_adjusted_p'])

    def test_failed_lv_disables_all_four_affected_contrasts(self):
        target=(self.data.seed==2)&(self.data.condition=='random_mixed')&(self.data.method=='lv')
        self.data.loc[target,'status']='failed'
        contrasts,_=compare(self.data,self.protocol,resamples=128)
        self.assertTrue(all(r['status']=='unavailable_missing_pairs' for r in contrasts[:4]))
        self.assertTrue(all(r['status']=='complete' for r in contrasts[4:]))

    def test_complete_contrast_matches_frozen_function(self):
        contrasts,_=compare(self.data,self.protocol,resamples=128)
        selected=self.data[self.data.condition=='random_mixed'].pivot(index='seed',columns='method',values='f1_percent')
        expected=paired_summary(selected.lv-selected.shepherd,20260912401,family_size=8,resamples=128)
        for key,value in expected.items():
            self.assertEqual(contrasts[0][key],value)

    def test_conditional_summary_is_never_full_population_mean(self):
        cells=pd.DataFrame([dict(status='complete',workflow_seconds=2,**{m:8 for m in METRICS}),
                            dict(status='failed',workflow_seconds=None,**{m:None for m in METRICS})])
        whole,conditional=summarize(cells,'random_mixed','shepherd',2)
        self.assertEqual(whole['n_failed'],1)
        self.assertIsNone(whole['f1_percent'])
        self.assertIsNone(whole['workflow_seconds'])
        self.assertEqual(conditional['f1_percent'],8)
        self.assertFalse(conditional['complete_population'])


if __name__=='__main__':
    unittest.main()
