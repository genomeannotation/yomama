#!/usr/bin/env python

import sys
from src.fastq import read_fastq
from src.oligos import read_oligos, deoligo_seqs
from src.consensus import call_consensus_for_yohan # ;-)
from src.reports import write_top_n_counts, write_read_counts

def main(args):
    if len(args) == 1:
        # Print usage
        print("Usage: yomama.py <mode> [path] [primer_diffs]")
        exit()

    mode = args[1]

    if mode == "consensus":
        BARCODE_DIFFS = 1
        LINKER_DIFFS = 3
        PRIMER_DIFFS = 3

        # Read oligos file
        with open("foo.oligos", "r") as oligos:
            oligos = read_oligos(oligos)
        if not oligos:
            sys.stderr.write("Aww naw, couldn't read oligos.\n")
            exit()

        # Read fastq file and deoligo seqs
        with open("foo.fastq", "r") as fastq:
            print("Reading fastq...")
            seqs = read_fastq(fastq)
            if not seqs:
                sys.stderr.write("Oh snap, failed to read fastq.\n")
                exit()
            print("Sorting reads...")
            sorted_reads = deoligo_seqs(seqs, oligos, BARCODE_DIFFS, LINKER_DIFFS, PRIMER_DIFFS)

        if not sorted_reads:
            sys.stderr.write("Failed to deoligo seqs, I'm out.\n")

        # Now call consensus
        print("Calling consensus...")
        with open("homozygous.fasta", "w") as ones,\
             open("heterozygous.fasta", "w") as twos:
            call_consensus_for_yohan(sorted_reads, ones, twos)

        # Generate a reports
        print("Generating reports...")
        with open("top_3_counts.tsv", "w") as report:
            write_top_n_counts(sorted_reads, 3, report)
        with open("reads_matrix.tsv", "w") as matrix:
            write_read_counts(sorted_reads.reads, matrix)
    elif mode == "heteroplasmy":
        print("Not implemented yet")
        exit()
    else:
        print("Invalid mode '"+mode+"', use 'consensus' or 'heteroplasmy'")



###################

if __name__ == "__main__":
    main(sys.argv)
