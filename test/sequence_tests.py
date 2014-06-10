#!/usr/bin/env python

import unittest
from mock import Mock
from src.sequence import Sequence, add_sample_name_from_header

class TestSequence(unittest.TestCase):

    def setUp(self):
        self.seq1 = Sequence("seq1", "GATTACA", [30, 30, 30, 30, 30, 30, 30])

    def test_get_subseq(self):
        self.assertEqual("ATTA", self.seq1.get_subseq(2, 5))

    def test_to_fasta(self):
        expected = ">seq1\nGATTACA\n"
        self.assertEqual(expected, self.seq1.to_fasta())

    def test_add_sample_name_from_header(self):
        header = "AC326-A:29:000000000-A92C9:1:1101:20529:1320 1:N:0:O243 ACCACTGT|0|TATCCTCT|0"
        test_seq = Sequence(header, "GATTACA")
        self.assertFalse(test_seq.sample)
        add_sample_name_from_header(test_seq)
        self.assertTrue(test_seq.sample)

    def test_to_string(self):
        self.seq1.locus = "foo_locus"
        self.seq1.sample = "foo_sample"
        expected = "Header: seq1\n"
        expected += "Bases: GATTACA\n"
        expected += "Locus: foo_locus\n"
        expected += "Sample: foo_sample\n"
        actual = str(self.seq1)
        self.assertEqual(expected, actual)

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()
