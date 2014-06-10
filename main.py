#!/usr/bin/env python

import sys
from src.fastq import read_fastq
from src.sequence import get_sample_name

def main():
    # Read fastq file
    with open("foo.fastq", "r") as fastq:
        seqs = read_fastq(fastq)
    if not seqs:
        sys.stderr.write("Oh snap, failed to read fastq.\n")
        exit()

    # Create dictionary to hold counts
    counts_dict = {}

    for seq in seqs:
        # Get sample name for each sequence
        add_sample_name_from_header(seq)

        # Deprimer each sequence
        # TODO 

        # Build dictionary of counts of unique reads for each locus/sample
        update_counts_dict(counts_dict, seq)


###################

if __name__ == "__main__":
    main()
