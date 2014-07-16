#!/usr/bin/env python

import unittest
from mock import Mock
import io
from src.oligos import *

class TestOligos(unittest.TestCase):

    def setUp(self):
        self.ipr_file = io.StringIO(\
        "primer\tGATACA\tGATACA\tsample1\n"\
        "primer\tGATACA\tGATACA\tsample2\n")

    def test_sorted_reads_reads_by_locus_sample(self):
        reads_dict = {"locus":{"sample":{"ATGC":(["SCORES1", "SCORES2"], 2)}}}
        reads = SortedReads(reads_dict)
        self.assertEqual(list(reads.reads_by_locus_sample()), [("locus", "sample", [("ATGC", ["SCORES1", "SCORES2"], 2)])])

    def test_read_oligos(self):
        oligos = read_oligos(self.ipr_file)
        self.assertEquals(2, len(oligos["primer"]))

    def test_sort_seq_positive_match(self):
        oligos = {"primer" : {"foo":PrimerPair("GATACA", "GATACA"), "dog":PrimerPair("ATGC", "ATGC")}}
        seq = Mock()
        seq.bases = "GATACAGGGGGTGTATC"
        seq.locus = ""
        sort_seq(oligos, seq)
        self.assertEqual(seq.bases, "GGGGG")
        self.assertEquals(seq.locus, "foo")

    def test_sort_seq_positive_no_match(self):
        oligos = {"primer" : {"foo":PrimerPair("GATACA", "GATACA"), "dog":PrimerPair("ATGC", "ATGC")}}
        seq = Mock()
        seq.bases = "GAGACAGGGGGTGTATC"
        seq.locus = ""
        sort_seq(oligos, seq)
        self.assertEqual(seq.bases, "GAGACAGGGGGTGTATC")
        self.assertEquals(seq.locus, "")

    def test_sort_seq_positive_match_tolerant(self):
        oligos = {"primer" : {"foo":PrimerPair("GATACA", "GATACA"), "dog":PrimerPair("ATGC", "ATGC")}}
        seq = Mock()
        seq.bases = "GAGACAGGGGGTGTATA"
        seq.locus = ""
        sort_seq(oligos, seq, 1)
        self.assertEqual(seq.bases, "GGGGG")
        self.assertEquals(seq.locus, "foo")


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOligos))
    return suite

if __name__ == '__main__':
    unittest.main()
