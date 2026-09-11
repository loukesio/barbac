"""Small independent examples for denominators, missing reads, and mapping QA."""
import unittest
import pandas as pd
from run_comparison import evaluate


class MetricTests(unittest.TestCase):
    def test_missing_and_wrong_reads_remain_in_denominator(self):
        truth=pd.DataFrame({'barcode':['AAAA','CCCC','GGGG'], 'true_count':[12,5,0]})
        labels=pd.DataFrame({'BC':['AAAA','AAAT','CCCC'], 'Count':[10,2,5], 'true_BC':['AAAA','AAAA','CCCC']})
        centroids=pd.DataFrame({'central_barcode':['AAAA','TTTT'], 'sum_counts':[10,2]})
        members=pd.DataFrame({'member':['AAAA','AAAT'], 'central_barcode':['AAAA','TTTT']})
        r=evaluate(centroids,members,labels,truth)
        self.assertEqual((r['tp'],r['fn'],r['fp']),(1,2,1))
        self.assertEqual((r['correct_reads'],r['misassigned_reads'],r['unassigned_reads']),(10,2,5))
        self.assertAlmostEqual(r['read_assignment_accuracy'],10/17)
        self.assertAlmostEqual(r['f1'],2/5)
        self.assertAlmostEqual(r['positive_truth_f1'],1/2)
        self.assertEqual(r['observed_fn'],1)

    def test_invented_output_counts_are_rejected(self):
        truth=pd.DataFrame({'barcode':['AAAA'], 'true_count':[12]})
        labels=pd.DataFrame({'BC':['AAAA'], 'Count':[12], 'true_BC':['AAAA']})
        centroids=pd.DataFrame({'central_barcode':['AAAA'], 'sum_counts':[11]})
        members=pd.DataFrame({'member':['AAAA'], 'central_barcode':['AAAA']})
        with self.assertRaises(AssertionError): evaluate(centroids,members,labels,truth)

    def test_duplicate_assignment_is_rejected(self):
        truth=pd.DataFrame({'barcode':['AAAA'], 'true_count':[12]})
        labels=pd.DataFrame({'BC':['AAAA'], 'Count':[12], 'true_BC':['AAAA']})
        centroids=pd.DataFrame({'central_barcode':['AAAA'], 'sum_counts':[12]})
        members=pd.DataFrame({'member':['AAAA','AAAA'], 'central_barcode':['AAAA','AAAA']})
        with self.assertRaises(AssertionError): evaluate(centroids,members,labels,truth)

    def test_consensus_can_recover_an_unobserved_true_sequence(self):
        truth=pd.DataFrame({'barcode':['AAAA','CCCC'], 'true_count':[10,0]})
        labels=pd.DataFrame({'BC':['AAAT'], 'Count':[10], 'true_BC':['AAAA']})
        centroids=pd.DataFrame({'central_barcode':['AAAA'], 'sum_counts':[10]})
        members=pd.DataFrame({'member':['AAAT'], 'central_barcode':['AAAA']})
        r=evaluate(centroids,members,labels,truth)
        self.assertEqual(r['absent_truth_recovered'],1)
        self.assertEqual(r['observed_fn'],0)
        self.assertEqual(r['positive_truth_f1'],1)
        self.assertEqual(r['read_assignment_accuracy'],1)

if __name__=='__main__': unittest.main()
