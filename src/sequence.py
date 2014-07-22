#!/usr/bin/env python

import sys

class Sequence:

    def __init__(self, header="", bases="", scores=None, locus="", sample=""):
        self.header = header
        self.bases = bases
        if not scores:
            scores = []
        self.scores = scores
        self.locus = locus
        self.sample = sample
        self.has_both_primers = False
        self.reverse = False

    def __str__(self):
        result = "Header: " + self.header + "\n"
        result += "Bases: " + self.bases + "\n"
        result += "Locus: " + self.locus + "\n"
        result += "Sample: " + self.sample + "\n"
        return result

    def to_fasta(self):
        result = '>' + self.header + '\n'
        result += self.bases + '\n'
        return result

    def get_subseq(self, start=1, stop=None):
        """Returns subseq from 1-based index "start" to "stop", inclusive."""
        if not stop:
            stop = len(self.bases)
        if stop > len(self.bases):
            return ""
        return self.bases[start-1:stop]

def add_sample_name_from_header(seq):
    try:
        header_fields = seq.header.strip().split()
        sample_name = header_fields[1].split(":")[3]
        seq.sample = sample_name
    except IndexError:
        sys.stderr.write("Unable to find sample name in seq with header " + seq.header + "\n")
        return

def sliding_window_filter(seq, window_size, min_avg):
    for start in range(0, len(seq.scores)-window_size+1):
        avg = sum(seq.scores[start:start+window_size])/window_size
        if avg < min_avg:
            return False
    return True
