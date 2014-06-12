#!/usr/bin/env python

import unittest
from src.consensus import call_consensus

class TestConsensus(unittest.TestCase):

    def setUp(self):
        self.seqs_dict = {}
        self.seqs_dict["ACGTACGT"] = 1
        self.seqs_dict["GATTACA"] = 8 

    def test_call_consensus_no_mins(self):
        consensus_seqs = call_consensus(self.seqs_dict)
        self.assertTrue("ACGTACGT" in consensus_seqs)
        self.assertTrue("GATTACA" in consensus_seqs)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestConsensus))
    return suite

if __name__ == '__main__':
    unittest.main()
