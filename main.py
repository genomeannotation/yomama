#!/usr/bin/env python

import sys
from src.fastq import read_fastq
from src.oligos import read_oligos, deoligo_seqs
from src.consensus import call_consensus_for_yohan # ;-)

def main():
    PRIMER_DIFFS = 2

    # Read oligos file
    with open("foo.oligos", "r") as oligos:
        oligos = read_oligos(oligos)
    if not oligos:
        sys.stderr.write("Aww naw, couldn't read oligos.\n")
        exit()

    # Read fastq file and deoligo seqs
    with open("foo.fastq", "r") as fastq:
        seqs = read_fastq(fastq)
        if not seqs:
            sys.stderr.write("Oh snap, failed to read fastq.\n")
            exit()
        counts = deoligo_seqs(seqs, oligos, PRIMER_DIFFS)

    if not counts:
        sys.stderr.write("Failed to deoligo seqs, I'm out.\n")

    # Now call consensus
    with open("homozygous.fasta", "w") as ones,\
         open("heterozygous.fasta", "w") as twos:
        call_consensus_for_yohan(counts, ones, twos)



###################

if __name__ == "__main__":
    main()
