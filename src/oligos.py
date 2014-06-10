from collections import namedtuple
from src.sequtil import reverse_complement

PrimerPair = namedtuple("Oligo", "left right")

def read_oligos(io_buffer):
    oligos = {} # Map of locus names to primer pairs
    for line in io_buffer:
        if not line:
            continue
        columns = line.strip("\t\n ").split("\t")
        if len(columns) != 4:
            continue
        oligos[columns[3]] = PrimerPair(columns[1], columns[2])
    return oligos

def sort_seq(oligos, seq):
    for locus, primer_pair in oligos.items():
        if primer_pair.left == seq.bases[:len(primer_pair.left)] and\
           primer_pair.right == reverse_complement(seq.bases[-len(primer_pair.right):]):
            seq.locus = locus # Sort
            seq.bases = seq.bases[len(primer_pair.left):-len(primer_pair.right)] # Trim
            return
