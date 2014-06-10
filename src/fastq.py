from src.sequence import Sequence

def read_fastq(io_buffer):
    """Returns a generator of Sequence objects"""
    scores_are_next = False
    header = ''
    bases = ''
    scores = []
    for line in io_buffer:
        if line[0] == '@':
            if len(header) > 0:
                # Yield the previous Sequence
                yield Sequence(header, bases, scores)
            header = line[1:].strip()  # Get the next header
        elif line[0] == '+':
            scores_are_next = True
        else:
            if scores_are_next:
                scores = translate_scores(line.strip())
                scores_are_next = False
            else:
                bases += line.strip()
    # Add the last sequence
    yield Sequence(header, bases, scores)

def translate_scores(scorestring):
    """Translates a string of phred scores to a list of integers"""
    # TODO
    return []
