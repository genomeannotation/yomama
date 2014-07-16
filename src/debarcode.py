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
            seq_barcode = seq.bases[:len(barcode)]
            if barcode == seq_barcode:
                seq.bases = seq.bases[len(barcode):] # Trim
                seq.sample = sample
                break
            else:
                mismatches = compare_seqs(barcode, seq_barcode)
                if mismatches <= mismatch_limit:
                    seq.bases = seq.bases[len(barcode):] # Trim
                    seq.sample = sample
                    break
        yield seq
    print("Debarcoded")
