#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

def reverse_complement(seq):
    bases = ['a', 'c', 'g', 't', 'n', 'A', 'C', 'G', 'T', 'N']
    complements = ['t', 'g', 'c', 'a', 'n', 'T', 'G', 'C', 'A', 'N']
    rev_comp_dict = dict(zip(bases, complements))
    # Convert mixed or illegal bases to 'N'
    for i, base in enumerate(seq):
        if base not in 'actgnACTGN':
            seq = seq[0:i] + 'N' + seq[i+1:]
    return ''.join([rev_comp_dict.get(base) for base in reversed(seq)])

def score_as_int(char, offset=33):
    return ord(char) - offset

def score_as_char(num, offset=33):
    return chr(num + offset)

def compare_seqs(seq1, seq2, max_mismatch):
    if seq1 == seq2:
        return True
    mismatch = 0
    for i in range(0, min(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            mismatch += 1
            if mismatch > max_mismatch:
                return False
    return True

def search_seq(seq, search, max_distance, max_mismatch):
    if len(seq) < len(search):
        return None
    if max_distance:
        search_distance = max_distance
    else:
        search_distance = len(seq)-len(search)
    for start in range(0, search_distance):
        subseq = seq[start:start+len(search)]
        if compare_seqs(subseq, search, max_mismatch):
            return start
    return None
