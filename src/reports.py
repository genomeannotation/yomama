def write_top_n_counts(sorted_reads, num_counts, outfile):
    """Generates a table giving total reads and top N read counts.

    Format of output is locus  sample  total_reads  top_n_counts

    Args:
        sorted_reads: a SortedReads object
        num_counts: the number of sequence counts to include in the output
    """
    outfile.write("locus\tsample\ttotal_reads\ttop_" +
                        str(num_counts) + "_seq_counts_proportions\n")
    for read in sorted_reads.reads_by_locus_sample():
        locus = read[0]
        sample = read[1]
        seqs_counts = read[2]
        # seqs_counts is a list of (seq, [qual], count) tuples
        outfile.write(locus + "\t" + sample + "\t")
        seq_counts = [s[2] for s in seqs_counts]
        if len(seq_counts) < num_counts:
            outfile.write("(fewer than " + str(num_counts) +
                            " uniq seqs; counts are: "+
                            str(seq_counts) + ")\n")
            continue
        total_reads = sum(seq_counts)
        top_n_counts = sorted(seq_counts)[::-1][:3]
        outfile.write(str(total_reads) + "\t")
        proportions = []
        for count in top_n_counts:
            proportions.append(str(float(count) / total_reads))
        outfile.write(",".join(proportions) + "\n")


def print_summary(counts_dict):
    """Prints a table of the counts of each seq at each locus/sample point."""
    for locus, locus_dict in counts_dict.items():
        for sample, sample_dict in locus_dict.items():
            for seq, count in sample_dict.items():
                print(locus+"\t"+sample+"\t"+seq+"\t"+str(count))


def print_read_counts(counts_dict):
    """Writes a matrix of total read counts for each locus/sample."""
    # Generate list of loci, samples
    loci = sorted(counts_dict.keys())
    samples = []
    for locus in loci:
        for sample in counts_dict[locus].keys():
            if sample not in samples:
                samples.append(sample)
    samples = sorted(samples)

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
