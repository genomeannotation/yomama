import sys
from collections import namedtuple
from src.sequtil import reverse_complement, compare_seqs
from src.sequence import add_sample_name_from_header
from src.debarcode import read_samples, debarcode_seqs

PrimerPair = namedtuple("Oligo", "left right")

class SortedReads:
    def __init__(self, reads):
        self.reads = reads

    def reads_by_locus_sample(self):
        for locus, locus_dict in self.reads.items():
            for sample, sample_dict in locus_dict.items():
                yield (locus, sample, [(bases, seq_data[0], seq_data[1]) for bases, seq_data in sample_dict.items()])

def read_oligos(io_buffer):
    oligos = {} # Map of locus names to primer pairs
    for line in io_buffer:
        if not line:
            continue
        columns = line.strip("\t\n ").split("\t")
        if len(columns) != 4:
            continue
        if columns[0] not in oligos:
            oligos[columns[0]] = {}
        oligos[columns[0]][columns[3]] = PrimerPair(columns[1], columns[2])
    return oligos

def strip_primer(seq, primer1, primer2, max_mismatch):
    seq_left = seq.bases[:len(primer1)] 
    seq_right = reverse_complement(seq.bases[-len(primer2):])
    if primer1 == seq_left and\
       primer2 == seq_right:
        seq.bases = seq.bases[len(primer1):-len(primer2)] # Trim
        return True
    else: # Doesn't match perfectly, check if there's few enough mismatches
        left_match = compare_seqs(primer1, seq_left, max_mismatch)
        if not left_match:
            return False
        right_match = compare_seqs(primer2, seq_right, max_mismatch)
        if not right_match:
            return False

        # It's good enough, take it
        seq.bases = seq.bases[len(primer1):-len(primer2)] # Trim
        return True
    return False

def sort_seq(oligos, seq, max_mismatch = 0):
    for locus, primer_pair in oligos["primer"].items():
        if strip_primer(seq, primer_pair.left, primer_pair.right, max_mismatch):
            seq.locus = locus # Sort
            return

def deoligo_seqs(seqs, oligos, bdiffs, ldiffs, pdiffs):
    with open("foo.samples", "r") as samples_file:
        samples = read_samples(samples_file)
        seqs = debarcode_seqs(seqs, samples, bdiffs)
        print("Debarcoded")

    sorted_reads = {}
    for seq in seqs:
        # Get sample name for each sequence
        if not seq.sample:
            add_sample_name_from_header(seq)
        
        # Skip if no sample name found
        if not seq.sample:
            sys.stderr.write("Failed to find sample for sequence:\t" +
                             seq.header + "\t" + seq.bases + "\n")
            continue

        # Deprimer each sequence
        sort_seq(oligos, seq, max_mismatch=pdiffs)

        # Skip if deprimering didn't work
        if not seq.locus:
            sys.stderr.write("Failed to deprimer sequence:\t" + seq.header +
                             "\t" + seq.bases + "\n")
            continue

        # Build dictionary of counts of unique reads for each locus/sample
        update_reads_dict(sorted_reads, seq)
    return SortedReads(sorted_reads)

def update_reads_dict(counts_dict, seq):
    # dict maps locus to a locus_dict, which maps sample to a seq dict,
    # which maps seqs to counts. omg wtf.
    if seq.locus not in counts_dict:
        counts_dict[seq.locus] = {}
    update_locus_dict(counts_dict[seq.locus], seq)

def update_locus_dict(locus_dict, seq):
    if seq.sample not in locus_dict:
        locus_dict[seq.sample] = {}
    update_sample_dict(locus_dict[seq.sample], seq)

def update_sample_dict(sample_dict, seq):
    if seq.bases not in sample_dict:
        sample_dict[seq.bases] = [[], 0]
    sample_dict[seq.bases][0].append(seq.scores)
    sample_dict[seq.bases][1] += 1

