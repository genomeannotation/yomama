#!/usr/bin/env python

import sys
from src.fastq import read_fastq
from src.sequence import add_sample_name_from_header
from src.oligos import read_oligos, sort_seq
from src.consensus import call_consensus

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
    total_reads = 0
    deprimered_reads = 0

    for seq in seqs:
        total_reads += 1

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

        deprimered_reads += 1

        # Add locus name to list
        if seq.locus not in loci:
            loci.append(seq.locus)

        # Build dictionary of counts of unique reads for each locus/sample
        update_counts_dict(counts, seq)

    fastq.close()
    #print_summary(counts)
    samples = sorted(samples)
    loci = sorted(loci)
    #print_read_counts(counts, loci, samples)
    sys.stdout.write("Total reads: " + str(total_reads) + "\n")
    sys.stdout.write("Deprimered reads: " + str(deprimered_reads) + "\n")
    call_consensus_sequences(counts, 30, 0.25, dry_run=True)

def call_consensus_sequences(counts_dict, min_count, min_percentage, dry_run=False):
    # dry_run generates a summary of read counts at each locus/sample;
    # non-dry_run writes fasta files of consensus calls
    if dry_run:
        sys.stdout.write("locus\tsample\ttotal_reads\ttop_3_seq_counts\n")
    else:
        twofers = open("two_consensus_seqs.fasta", "w")
        onefers = open("one_consensus_seq.fasta", "w")
    for locus, locus_dict in counts_dict.items():
        for sample, sample_dict in locus_dict.items():
            if dry_run:
                sys.stdout.write(locus + "\t" + sample + "\t")
                seq_counts = sample_dict.values()
                if len(seq_counts) < 3:
                    sys.stdout.write("(fewer than 3 uniq seqs; counts are: "+
                                    str(seq_counts) + ")\n")
                    continue
                total_reads = sum(seq_counts)
                top_3_counts = sorted(seq_counts)[::-1][:3]
                sys.stdout.write(str(total_reads) + "\t")
                for count in top_3_counts[:-1]:
                    sys.stdout.write(str(count) + ",")
                sys.stdout.write(str(top_3_counts[-1]) + "\n")
            else:
                consensus = call_consensus(sample_dict,
                                min_count, min_percentage)
                if consensus:
                    seqs = consensus.keys()
                    if len(seqs) == 1:
                        seq = seqs[0]
                        percent = str(consensus[seq])
                        # write it twice
                        onefers.write(consensus_fasta(locus, sample, seq, percent))
                        onefers.write(consensus_fasta(locus, sample, seq, percent))
                    elif len(seqs) == 2:
                        for seq in seqs:
                            percent = str(consensus[seq])
                            twofers.write(consensus_fasta(locus, sample, seq, percent))
    if not dry_run:
        onefers.close()
        twofers.close()

def consensus_fasta(locus, sample, seq, percent):
    result = ">" + locus + "_" + sample + "_" + percent + "\n"
    result += seq + "\n"
    return result

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
    for sample in samples:
        sys.stdout.write(sample + "\t")
    sys.stdout.write("total\n")

    # For each locus, write read counts for each sample
    for locus in loci:
        if locus not in counts_dict:
            write_zero_row(locus, samples)
        else:
            write_counts_per_sample(locus, counts_dict[locus], samples)

def write_zero_row(locus, samples):
    sys.stdout.write(locus + "\t")
    for sample in samples:
        sys.stdout.write("0\t")
    sys.stdout.write("0\n")

def write_counts_per_sample(locus, locus_dict, samples):
    sys.stdout.write(locus + "\t")
    total = 0
    for sample in samples:
        if sample not in locus_dict:
            sys.stdout.write("0\t")
        else:
            read_count = total_reads(locus_dict[sample])
            total += read_count
            sys.stdout.write(str(read_count) + "\t")
    sys.stdout.write(str(total) + "\n")

def total_reads(sample_dict):
    total = 0
    for count in sample_dict.values():
        total += count
    return total


###################

if __name__ == "__main__":
    main()
