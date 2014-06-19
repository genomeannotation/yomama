def top_n_counts(counts_dict, num_counts):
    """Generates a table giving total reads and top N read counts.

    Format of output is locus  sample  total_reads  top_n_counts

    Args:
        counts_dict: a dictionary that maps loci to locus dicts,
                     which map samples to sample dicts,
                     which map sequences to counts :)
        num_counts: the number of sequence counts to include in the output
    """
    sys.stdout.write("locus\tsample\ttotal_reads\ttop_" +
                        str(num_counts) + "_seq_counts\n")
    for locus, locus_dict in counts_dict.items():
        for sample, sample_dict in locus_dict.items():
            sys.stdout.write(locus + "\t" + sample + "\t")
            seq_counts = sample_dict.values()
            if len(seq_counts) < num_counts:
                sys.stdout.write("(fewer than " + str(num_counts) +
                                " uniq seqs; counts are: "+
                                str(seq_counts) + ")\n")
                continue
            total_reads = sum(seq_counts)
            top_n_counts = sorted(seq_counts)[::-1][:3]
            sys.stdout.write(str(total_reads) + "\t")
            for count in top_n_counts[:-1]:
                sys.stdout.write(str(count) + ",")
            sys.stdout.write(str(top_n_counts[-1]) + "\n")


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
