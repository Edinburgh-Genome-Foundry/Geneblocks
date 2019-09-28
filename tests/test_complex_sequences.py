from geneblocks import DiffBlocks, CommonBlocks, random_dna_sequence
import geneblocks.sequence_modification_utils as smu
import matplotlib.pyplot as plt


def test_complex_sequences():
    seq1 = random_dna_sequence(50000, seed=123)
    seq1 = smu.copy(seq1, 25000, 30000, 50000)

    seq2 = seq1
    seq2 = smu.insert(seq2, 39000, random_dna_sequence(100))
    seq2 = smu.insert(seq2, 38000, random_dna_sequence(100))
    seq2 = smu.reverse(seq2, 30000, 35000)
    seq2 = smu.swap(seq2, (30000, 35000), (45000, 480000))
    seq2 = smu.delete(seq2, 20000, 2000)
    seq2 = smu.insert(seq2, 10000, random_dna_sequence(2000))
    seq2 = smu.insert(seq2, 0, 1000 * "A")
    diff_blocks = DiffBlocks.from_sequences(seq1, seq2).merged()

    b = diff_blocks.blocks
    assert len(b) == 15
    assert b[0].operation == "insert"
    start, end, _ = b[0].s2_location.to_tuple()
    assert end - start == 1000
    assert b[1].operation == "equal"
    assert b[2].operation == "insert"
    start, end, _ = b[2].s2_location.to_tuple()
    assert end - start == 2000
    assert sorted([b[6].operation, b[7].operation]) == ["change", "transpose"]
    assert sorted([b[-1].operation, b[-2].operation]) == ["change", "reverse"]

    s1, s2 = diff_blocks.reconstruct_sequences_from_blocks(diff_blocks.blocks)
    assert s1 == seq1
    assert s2 == seq2
