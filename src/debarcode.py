def read_samples(io_buffer):
    samples = {}
    for line in io_buffer:
        columns = line.strip().split("\t")
        samples[columns[0]] = (columns[2], columns[1])
    return samples
