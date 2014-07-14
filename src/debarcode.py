from src.sequtil import compare_seqs

def read_samples(io_buffer):
    samples = {}
    for line in io_buffer:
        columns = line.strip().split("\t")
        if len(columns) != 3:
            continue
        samples[columns[2]] = columns[1]
    return samples

def debarcode_seqs(seqs, samples, mismatch_limit):
    for seq in seqs:
        for sample, barcode in samples.items():
            mismatches = compare_seqs(barcode, seq.bases[:len(barcode)])
            if mismatches <= mismatch_limit:
                seq.sample = sample
                yield seq
