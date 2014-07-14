def read_samples(io_buffer):
    samples = {}
    for line in io_buffer:
        columns = line.strip().split("\t")
        if len(columns) != 3:
            continue
        samples[columns[2]] = columns[1]
    return samples
