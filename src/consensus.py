def call_consensus(seqs_dict, min_count=0, min_percent=0.0):
    """Takes a dictionary of sequences to counts, returns consensus seq(s).

    Applies filters to sequences; if one or two remain after filters,
    they are returned. If zero or three or more remain, writes an error
    to stderr and returns None.
    """
    total_count = sum(seqs_dict.values())
    passing_seqs = []
    for seq, count in seqs_dict.items():
        if count >= min_count and count/total_count >= min_percent:
            passing_seqs.append(seq)
    if len(passing_seqs) == 0 or len(passing_seqs) >= 3:
        return None
    else:
        return passing_seqs
