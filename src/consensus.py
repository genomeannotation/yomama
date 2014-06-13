import sys

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
