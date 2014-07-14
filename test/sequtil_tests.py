#!/usr/bin/env python

import unittest
from src.sequtil import *

class TestTranslate(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEqual('C', reverse_complement('G'))
        self.assertEqual('CAT', reverse_complement('ATG'))

    def test_reverse_complement_with_bogus_base(self):
        self.assertEqual('CATN', reverse_complement('MATG'))

    def test_reverse_complement_longer_seq(self):
        self.assertEqual('TGTAATCTGTAATCTGTAATCTGTAATCTGTAATC', reverse_complement('GATTACAGATTACAGATTACAGATTACAGATTACA'))

    def test_compare_seqs(self):
        seq1 = "GATACA"
        seq2 = "GATCCC"
        self.assertEqual(compare_seqs(seq1, seq2), 2)

        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranslate))
    return suite

if __name__ == '__main__':
    unittest.main()
