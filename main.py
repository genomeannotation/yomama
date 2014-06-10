#!/usr/bin/env python

import sys
from src.fastq import read_fastq
from src.sequence import add_sample_name_from_header
from src.oligos import read_oligos, sort_seq

def main():
    # Read fastq file
    fastq = open("foo.fastq", "r")
    seqs = read_fastq(fastq)
    if not seqs:
        sys.stderr.write("Oh snap, failed to read fastq.\n")
        exit()

    # Read oligos file
    with open("foo.oligos", "r") as oligos:
        oligos = read_oligos(oligos)
    if not oligos:
        sys.stderr.write("Aww naw, couldn't read oligos.\n")
        exit()

    # Create dictionary to hold counts
    counts = {}

    for seq in seqs:
        # Get sample name for each sequence
        add_sample_name_from_header(seq)

        # Deprimer each sequence
        sort_seq(oligos, seq)

        # Build dictionary of counts of unique reads for each locus/sample
        update_counts_dict(counts, seq)

    # Write contents of dictionary to stdout
    # TODO
    print(counts)

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


###################

if __name__ == "__main__":
    main()
