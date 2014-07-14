#!/usr/bin/env python

import unittest
from mock import Mock
import io
from src.debarcode import *

class TestDebarcode(unittest.TestCase):

    def setUp(self):
        self.samples_file = io.StringIO(\
        "A1\tACGAGTGCGT\t20110222_001_A_001\n"\
        "B1\tACGCTCGACA\tRX100706_001\n")

    def test_read_samples(self):
        expected = {"20110222_001_A_001" : "ACGAGTGCGT",\
                    "RX100706_001" : "ACGCTCGACA"}
        actual = read_samples(self.samples_file)
        self.assertEqual(expected, actual)

    def test_debarcode_seqs(self):
        seqs = [Mock(), Mock()]
        seqs[0].bases = "ACGTGTGCGTACACTGACGACATGGTTCTACAGTCCATGTGGTCCATGCTCTGAGTTGCATTTGATAGAATTGGAGGCCGAAGCGTACCCGAACTGGTAAATGTTGATGACCCAGATGTTTAGAAATTTGGAATTTGGTGTTATCCAATACAACAGGGAGTCGAGAAGTGGTTAGGCCCTAGCAAACCGAATTCAATTGACCTAAAGGAAAGGCGAGATAACGGATTTGGTGGGAATTGGGGTCCTTTACGGAAACGTCCTTTGGGAATCGGAGGATTCCGTTAAAAGGACCAAAAAAGGCCAAAAAATTAAAGGGTTT"
        seqs[1].bases = "ACGCTCGACATACGGTAGCAGAGACTTGGTCTTATCTGATCGCCTTCGAACCTCTAACTTTCGTTCTTGGG"

        samples = {"20110222_001_A_001" : "ACGAGTGCGT",\
                   "RX100706_001" : "ACGCTCGACA"}

        debarcoded_seqs = list(debarcode_seqs(seqs, samples, 1))
        self.assertEqual(debarcoded_seqs[0].sample, "20110222_001_A_001")
        self.assertEqual(debarcoded_seqs[1].sample, "RX100706_001")


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDebarcode))
    return suite

if __name__ == '__main__':
    unittest.main()
