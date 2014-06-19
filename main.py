#!/usr/bin/env python

import sys
from src.fastq import read_fastq
from src.oligos import read_oligos, deoligo_seqs
from src.consensus import call_consensus_for_yohan # ;-)

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
    counts = deoligo_seqs(seqs, oligos)

    # Now call consensus
    call_consensus_for_yohan(counts)



###################

if __name__ == "__main__":
    main()
