#!/usr/bin/env python

import unittest
import io
from src.oligos import *

class TestOligos(unittest.TestCase):

    def setUp(self):
        self.ipr_file = io.StringIO(\
        "primer\tGATACA\tTGTATC\tsample1\n"\
        "primer\tGATACA\tTGTATC\tsample2\n")

    def test_read_oligos(self):
        oligos = read_oligos(self.ipr_file)
        self.assertEquals(2, len(oligos))


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOligos))
    return suite

if __name__ == '__main__':
    unittest.main()
