#!/usr/bin/env python

# import all the lovely files
import unittest
import test.sequence_tests
import test.oligos_tests
import test.fastq_tests
import test.sequtil_tests
import test.consensus_tests
import test.debarcode_tests

# get suites from test modules
suites = [
test.sequence_tests.suite(),
test.oligos_tests.suite(),
test.fastq_tests.suite(),
test.sequtil_tests.suite(),
test.consensus_tests.suite(),
test.debarcode_tests.suite(),
]

# collect suites in a TestSuite object
suite = unittest.TestSuite()
for s in suites:
    suite.addTest(s)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)
