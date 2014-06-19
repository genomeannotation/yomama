import sys
from collections import namedtuple
from src.sequtil import reverse_complement
from src.sequence import add_sample_name_from_header

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
            left_mismatch = 0
            for i in range(0, len(primer_pair.left)):
                if primer_pair.left[i] != seq_left[i]:
                    left_mismatch += 1
                    if left_mismatch > max_mismatch:
                        break
            if left_mismatch > max_mismatch:
                continue
            right_mismatch = 0
            for i in range(0, len(primer_pair.right)):
                if primer_pair.right[i] != seq_right[i]:
                    right_mismatch += 1
                    if right_mismatch > max_mismatch:
                        break
            if right_mismatch > max_mismatch:
                continue

            # It's good enough, take it
            seq.locus = locus # Sort
            seq.bases = seq.bases[len(primer_pair.left):-len(primer_pair.right)] # Trim
            return

def deoligo_seqs(seqs, oligos):
    counts = {}
    for seq in seqs:
        # Get sample name for each sequence
        add_sample_name_from_header(seq)
        
        # Skip if no sample name found
        if not seq.sample:
            sys.stderr.write("Failed to find sample for sequence:\t" +
                             seq.header + "\t" + seq.bases + "\n")
            continue

        # Deprimer each sequence
        sort_seq(oligos, seq)

        # Skip if deprimering didn't work
        if not seq.locus:
            sys.stderr.write("Failed to deprimer sequence:\t" + seq.header +
                             "\t" + seq.bases + "\n")
            continue

        # Build dictionary of counts of unique reads for each locus/sample
        update_counts_dict(counts, seq)
    return counts

def update_counts_dict(counts_dict, seq):
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
        sample_dict[seq.bases] = 0
    sample_dict[seq.bases] += 1

