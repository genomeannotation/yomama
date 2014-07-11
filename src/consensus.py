import sys
import operator

def call_consensus_sequences(counts_dict, min_count, min_percentage,
                             one_seq_file, two_seqs_file):
    # writes fasta files of consensus calls
    onefers = one_seq_file
    twofers = two_seqs_file
    for locus, locus_dict in counts_dict.items():
        for sample, sample_dict in locus_dict.items():
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

def call_consensus_for_yohan(sorted_reads, one_seq_file, two_seqs_file):
    onefers = one_seq_file
    twofers = two_seqs_file
    for read in sorted_reads.reads_by_locus_sample():
        locus = read[0]
        sample = read[1]
        seqs_counts = read[2]
        # TODO not sample dict...
        consensus = yohan_consensus(sample_dict)
        if consensus:
            if len(consensus) == 1:
                seq = consensus[0][0]
                count = consensus[0][1]
                # write it twice
                onefers.write(consensus_fasta(locus, sample, seq, count))
                onefers.write(consensus_fasta(locus, sample, seq, count))
            elif len(consensus) == 2:
                for entry in consensus:
                    seq = entry[0]
                    count = entry[1]
                    twofers.write(consensus_fasta(locus, sample, seq, count))

def consensus_fasta(locus, sample, seq, count):
    result = ">" + locus + "_" + sample + "_" + str(count) + "\n"
    result += seq + "\n"
    return result


def call_consensus(seqs_dict, min_count=0, min_percent=0.0):
    """Calls consensus seq(s) and returns in fasta format.

    Args:
        seqs_dict: a dictionary mapping sequences to their counts
        min_count: an integer representing the minimum number of copies
                   of a sequence necessary for it to be considered
        min_percent: a number between 0.0 and 1.0 representing the minimum
                     percentage of reads that should contain the seq

    Applies filters to sequences; if one or two remain after filters,
    they are returned. If zero or three or more remain, writes an error
    to stderr and returns None.
    """
    total_count = sum(seqs_dict.values())
    passing_seqs = {}
    for seq, count in seqs_dict.items():
        seq_percent = float(count)/total_count
        if count >= min_count and seq_percent >= min_percent:
            passing_seqs[seq] = seq_percent
    if len(passing_seqs) == 0:
        sys.stderr.write("call_consensus found no seqs passing filters...\n")
        return None
    elif len(passing_seqs) >= 3:
        sys.stderr.write("call_consensus found > 2 consensus seqs\n")
        return None
    else:
        return passing_seqs

def yohan_consensus(seqs_dict):
    MIN_COUNT_HIGH = 50
    MIN_COUNT_LOW = 10
    sorted_seqs = sorted(seqs_dict.items(), key=operator.itemgetter(1), reverse=True)
    if not sorted_seqs:
        return []
    elif len(sorted_seqs) == 1:
        seq = sorted_seqs[0]
        if seq[1] > MIN_COUNT_HIGH:
            return [seq]
        else:
            return []
    elif len(sorted_seqs) == 2:
        first_seq = sorted_seqs[0]
        if first_seq[1] > MIN_COUNT_HIGH:
            second_seq = sorted_seqs[1]
            if second_seq[1] > 0.5 * first_seq[1]:
                return [first_seq, second_seq]
            else:
                return [first_seq]
        else:
            return []
    else:
        first_seq = sorted_seqs[0]
        second_seq = sorted_seqs[1]
        third_seq = sorted_seqs[2]
        first_seq_count = first_seq[1]
        second_seq_count = second_seq[1]
        third_seq_count = third_seq[1]
        if first_seq_count < MIN_COUNT_LOW:
            # fail
            return []

        homozygote = True
        if second_seq_count > 0.5 * first_seq_count:
            homozygote = False

        if first_seq_count < MIN_COUNT_HIGH:
            valid_seq = False
            # noise test
            if homozygote:
                if second_seq_count <= 0.2 * first_seq_count:
                    valid_seq = True
            else:
                if third_seq_count <= 0.2 * first_seq_count:
                    valid_seq = True
            if not valid_seq:
                return []
        else: # first_seq_count >= MIN_COUNT_HIGH:
            if homozygote:
                return [first_seq]
            else:
                return [first_seq, second_seq]

