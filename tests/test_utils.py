from geneblocks.utils import sequences_are_circularly_equal
from Bio import SeqIO
import os

this_directory = os.path.dirname(os.path.realpath(__file__))


def test_sequences_are_circularly_equal():

    block_1 = "ATGTGCACACGCACCGTGTGTGCACACACGTGTGCACACACGTGCACACGGTGT"
    block_2 = "ACACACATATACGCGTGCGTGCAAAACACATTTTACACGGCACGTGCA"
    block_3 = "ACCCACACTTTGTGTCGCGCACACGTGTG"

    # Three rotated sequences
    seq_1 = block_1 + block_2 + block_3
    seq_2 = block_2 + block_3 + block_1
    seq_3 = block_3 + block_1 + block_2

    # Sequences not equivalent to the previous
    seq_4 = block_2 + block_3 + block_1 + "A"
    seq_5 = block_1 + block_2 + block_2 + block_3

    assert sequences_are_circularly_equal([seq_1, seq_2])
    assert sequences_are_circularly_equal([seq_1, seq_2, seq_3])
    assert not sequences_are_circularly_equal([seq_1, seq_2, seq_3, seq_4])
    assert not sequences_are_circularly_equal([seq_1, seq_5])


def test_long_equivalents():
    """Test with two 12kb sequences with only a 10bp shift"""
    path = os.path.join(this_directory, "sequences", "long_equivalents.fa")
    seq_a, seq_b = SeqIO.parse(path, "fasta")
    assert sequences_are_circularly_equal([seq_a, seq_b])
