from .CommonBlocks import CommonBlocks
from .biotools import sequence_to_record


def _turn_sequence_into_record_if_necessary(sequence, record_id="id"):
    if hasattr(sequence, "seq"):
        return sequence
    else:
        return sequence_to_record(sequence, record_id=record_id)


def sequences_are_circularly_equal(sequences):
    """Return whether all the sequences represent the same circular sequence.

    This means that the sequences are differently rotated versions of a same
    circular sequence, and for any pair s1, s2 in the sequences, there is an
    index i such that s1 = s2[i:] + s2[:i].

    The ``sequences`` parameter should be a list of "ATGC" strings or
    SeqRecords.
    """
    sequences = [
        _turn_sequence_into_record_if_necessary(seq, record_id="REC_%d" % i)
        for i, seq in enumerate(sequences)
    ]
    if len(sequences) < 2:
        raise ValueError("Provide at least 2 sequences")
    elif len(sequences) > 2:
        first_equal = sequences_are_circularly_equal(sequences[:2])
        return first_equal and sequences_are_circularly_equal(sequences[1:])
    s1, s2 = sequences[:2]
    if s1.id == s2.id:
        s1.id += "_b"
    if len(s1) != len(s2):
        return False
    blocks = CommonBlocks.from_sequences(sequences=[s1, s2], min_block_size=2)
    if len(blocks.common_blocks) > 2:
        return False
    potential_pivot_indices = [
        index
        for data in blocks.common_blocks.values()
        for (origin, (start, end, _)) in data["locations"]
        for index in [start, end]
    ]
    s1, s2 = str(s1.seq), str(s2.seq)
    for index in potential_pivot_indices:
        new_s1 = s1[index:] + s1[:index]
        if new_s1 == s2:
            return True
    return False
