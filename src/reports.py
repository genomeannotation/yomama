import sys

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


def write_most_abundant_fasta(sorted_reads, outfile):
    """Writes fasta file of most abundant reads at each locus
    """
    print("Writing consensus")
    for read in sorted_reads.reads_by_locus_sample():
        print("enter")
        locus = read[0]
        sample = read[1]
        seqs_counts = read[2]
        # seqs_counts is a list of (seq, [qual], count) tuples
        print(seq_counts)
        if len(seq_counts) == 0:
            continue
        print("yayaya")
        # Sort by count and take last (largest count) seq
        consensus = sorted(seq_counts, key=lambda x: x[2])[-1]
        # Write fasta
        outfile.write(">"+locus + " " + sample + " count:"+consensus[2]+"\n"+consensus[0]+"\n")


def write_summary(counts_dict, io_buffer=sys.stdout):
    """Prints a table of the counts of each seq at each locus/sample point."""
    for locus, locus_dict in counts_dict.items():
        for sample, sample_dict in locus_dict.items():
            for seq, count in sample_dict.items():
                io_buffer.write(locus+"\t"+sample+"\t"+seq+"\t"+str(count)+"\n")


def write_read_counts(counts_dict, io_buffer=sys.stdout):
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
    io_buffer.write("\t")
    for sample in samples:
        io_buffer.write(sample + "\t")
    io_buffer.write("total\n")

    # For each locus, write read counts for each sample
    for locus in loci:
        if locus not in counts_dict:
            write_zero_row(locus, samples, io_buffer)
        else:
            write_counts_per_sample(locus, counts_dict[locus], samples, io_buffer)


def write_zero_row(locus, samples, io_buffer):
    io_buffer.write(locus + "\t")
    for sample in samples:
        io_buffer.write("0\t")
    io_buffer.write("0\n")


def write_counts_per_sample(locus, locus_dict, samples, io_buffer):
    io_buffer.write(locus + "\t")
    total = 0
    for sample in samples:
        if sample not in locus_dict:
            io_buffer.write("0\t")
        else:
            read_count = total_reads(locus_dict[sample])
            total += read_count
            io_buffer.write(str(read_count) + "\t")
    io_buffer.write(str(total) + "\n")


def total_reads(sample_dict):
    total = 0
    for info in sample_dict.values():
        total += info[1]
    return total
