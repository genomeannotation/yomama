from collections import namedtuple

PrimerPair = namedtuple("Oligo", "left right")

def read_oligos(io_buffer):
    oligos = {} # Map of sample names to primer pairs
    for line in io_buffer:
        if not line:
            continue
        columns = line.strip("\t\n ").split("\t")
        if len(columns) != 4:
            continue
        oligos[columns[3]] = PrimerPair(columns[1], columns[2])
    return oligos
