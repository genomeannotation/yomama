#!/usr/bin/env python

import unittest
from mock import Mock
import io
from src.debarcode import *

class TestDebarcode(unittest.TestCase):

    def setUp(self):
        pass


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestDebarcode))
    return suite

if __name__ == '__main__':
    unittest.main()
