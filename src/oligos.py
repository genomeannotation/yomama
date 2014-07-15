import sys
from collections import namedtuple
from src.sequtil import reverse_complement, compare_seqs
from src.sequence import add_sample_name_from_header

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
        oligos[columns[3]] = PrimerPair(columns[1], columns[2])
    return oligos

def sort_seq(oligos, seq, max_mismatch = 0):
    for locus, primer_pair in oligos.items():
        seq_left = seq.bases[:len(primer_pair.left)] 
        seq_right = reverse_complement(seq.bases[-len(primer_pair.right):])
        if primer_pair.left == seq_left and\
           primer_pair.right == seq_right:
            seq.locus = locus # Sort
            seq.bases = seq.bases[len(primer_pair.left):-len(primer_pair.right)] # Trim
            return
        else: # Doesn't match perfectly, check if there's few enough mismatches
            left_mismatch = compare_seqs(primer_pair.left, seq_left)
            if left_mismatch > max_mismatch:
                continue
            right_mismatch = compare_seqs(primer_pair.right, seq_right)
            if right_mismatch > max_mismatch:
                continue

            # It's good enough, take it
            seq.locus = locus # Sort
            seq.bases = seq.bases[len(primer_pair.left):-len(primer_pair.right)] # Trim
            return

def deoligo_seqs(seqs, oligos, bdiffs, ldiffs, pdiffs):
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
    return SortedReads(counts)

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

