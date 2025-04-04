from geneblocks.DiffBlocks import DiffBlock, DiffBlocks
from geneblocks.biotools import random_dna_sequence


def diff_string(seq1, seq2, contract_under=0):
    return str(DiffBlocks.from_sequences(seq1, seq2))


flanks_dict = {"short": 10 * "A", "long": random_dna_sequence(10000)}


def flanked(seq, flanks="short"):
    flank = flanks_dict[flanks]
    return flank + seq + flank


def test_1():
    assert (
        diff_string(flanked("T", "short"), flanked("C", "short"))
        == "replace 10-11|10-11"
    )


def test_2():
    assert diff_string(flanked("T", "long"), flanked("C", "long")) == ", ".join(
        [
            "equal 0-10000|0-10000",
            "replace 10000-10001|10000-10001",
            "equal 10001-20001|10001-20001",
        ]
    )


def test_3():
    assert (
        diff_string("AATAATAAT", "AAAAAAAAA")
        == "replace 2-3|2-3, replace 5-6|5-6, replace 8-9|8-9"
    )


def test_4():
    assert diff_string(
        flanked("AATAATAAT", "long"), flanked("AAAAAAAAA", "long")
    ) == ", ".join(
        [
            "equal 0-10002|0-10002",
            "replace 10002-10003|10002-10003",
            "replace 10005-10006|10005-10006",
            "replace 10008-10009|10008-10009",
            "equal 10009-20009|10009-20009",
        ]
    )


def test_5():
    assert (
        diff_string("AAAAATTTTGGAAA", "AAATTTTGGAAAAA")
        == "delete 3-5|3-3, insert 11-11|9-11"
    )


def test_6():
    assert diff_string(
        flanked("AAAAATTTTGGAAA", "long"), flanked("AAATTTTGGAAAAA", "long")
    ) == ", ".join(
        [
            "equal 0-10003|0-10003",
            "delete 10003-10005|10003-10003",
            "insert 10011-10011|10009-10011",
            "equal 10011-20014|10011-20014",
        ]
    )
