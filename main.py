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

    # Create dictionary to hold counts, lists to hold sample and locus names
    counts = {}
    samples = []
    loci = []

    for seq in seqs:
        # Get sample name for each sequence
        add_sample_name_from_header(seq)
        
        # Skip if no sample name found
        if not seq.sample:
            continue

        # Add sample name to list
        if seq.sample not in samples:
            samples.append(seq.sample)

        # Deprimer each sequence
        sort_seq(oligos, seq)

        # Skip if deprimering didn't work
        if not seq.locus:
            continue

        # Add locus name to list
        if seq.locus not in loci:
            loci.append(seq.locus)

        # Build dictionary of counts of unique reads for each locus/sample
        update_counts_dict(counts, seq)

    fastq.close()
    # Write contents of dictionary to stdout
    #print_summary(counts)
    samples = sorted(samples)
    loci = sorted(loci)
    print_read_counts(counts, loci, samples)

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

def print_summary(counts_dict):
    for locus, locus_dict in counts_dict.items():
        for sample, sample_dict in locus_dict.items():
            for seq, count in sample_dict.items():
                print(locus+"\t"+sample+"\t"+seq+"\t"+str(count))

def print_read_counts(counts_dict, loci, samples):
    # Print header
    sys.stdout.write("\t")
    for sample in samples[:-1]:
        sys.stdout.write(sample + "\t")
    sys.stdout.write(samples[-1] + "\n")

    # For each locus, write read counts for each sample
    for locus in loci:
        if locus not in counts_dict:
            write_zero_row(locus, samples)
        else:
            write_counts_per_sample(locus, counts_dict[locus], samples)

def write_zero_row(locus, samples):
    sys.stdout.write(locus + "\t")
    for sample in samples[:-1]:
        sys.stdout.write("0\t")
    sys.stdout.write("0\n")

def write_counts_per_sample(locus, locus_dict, samples):
    sys.stdout.write(locus + "\t")
    for sample in samples[:-1]:
        if sample not in locus_dict:
            sys.stdout.write("0\t")
        else:
            sys.stdout.write(total_reads(locus_dict[sample]) + "\t")
    # Take care of last sample and newline
    if samples[-1] not in locus_dict:
        sys.stdout.write("0\n")
    else:
        sys.stdout.write(total_reads(locus_dict[samples[-1]]) + "\n")

def total_reads(sample_dict):
    total = 0
    for count in sample_dict.values():
        total += count
    return str(count)


###################

if __name__ == "__main__":
    main()
