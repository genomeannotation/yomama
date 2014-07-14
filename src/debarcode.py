def read_samples(io_buffer):
    oligos = {} # Map of locus names to primer pairs
    for line in io_buffer:
        if not line:
            continue
        columns = line.strip("\t\n ").split("\t")
        if len(columns) != 4:
            continue
        oligos[columns[3]] = PrimerPair(columns[1], columns[2])
    return oligos
