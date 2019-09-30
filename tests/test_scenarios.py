import os
import matplotlib
from geneblocks import CommonBlocks, DiffBlocks, load_record
from geneblocks.biotools import reverse_complement, random_dna_sequence
from geneblocks.DiffBlocks import DiffBlock
from geneblocks.Location import Location

matplotlib.use("Agg")


def test_CommonBlocks_basics(tmpdir):
    C1, A, B, C2, C3, D, E, F = [
        random_dna_sequence(100 * L) for L in range(1, 9)
    ]

    sequences = {
        "a": C1 + A + C2,
        "b": B + D + C2 + C3,
        "c": E + C1 + C2,
        "d": C2 + C1 + F + C3,
        "e": C3 + reverse_complement(C2 + C1),
    }

    common_blocks = CommonBlocks.from_sequences(sequences)
    assert len(common_blocks.common_blocks) == 3

    axes = common_blocks.plot_common_blocks()
    fig_path = os.path.join(str(tmpdir), "basic_example.png")
    axes[0].figure.savefig(fig_path, bbox_inches="tight")

    # GET ALL COMMON BLOCKS AS BIOPYTHON RECORDS
    _ = common_blocks.common_blocks_records()
    _ = common_blocks.unique_blocks_records()

    # WRITE ALL COMMON BLOCKS INTO A CSV SPREADSHEET
    csv_path = os.path.join(str(tmpdir), "basic_example.csv")
    common_blocks.common_blocks_to_csv(target_file=csv_path)


def test_DiffBlocks_basics(tmpdir):
    seq_1 = load_record(os.path.join("tests", "sequences", "sequence1.gb"))
    seq_2 = load_record(os.path.join("tests", "sequences", "sequence2.gb"))

    diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2).merged()
    # next line is just to cover separate_axes=false
    diff_blocks.plot(figure_width=8, separate_axes=False)
    ax1, __ax2 = diff_blocks.plot(figure_width=8)
    fig_path = os.path.join(str(tmpdir), "diff_blocks.png")
    ax1.figure.savefig(fig_path, bbox_inches="tight")
    assert list(map(str, diff_blocks.blocks)) == [
        "insert 0-0|0-120",
        "equal 0-1000|120-1120",
        "replace 1000-1004|1120-1124",
        "equal 1004-1503|1124-1623",
        "insert 1503-1503|1623-1723",
        "equal 1503-2304|1723-2524",
        "delete 2304-2404|2524-2524",
        "equal 2404-3404|2524-3524",
    ]

def test_features_transfer():
    seq_folder = os.path.join("tests", "sequences", "features_transfer")
    insert = load_record(os.path.join(seq_folder, "insert.gb"), name="insert")
    plasmid = load_record(
        os.path.join(seq_folder, "plasmid_to_annotate.gb"), name="plasmid"
    )
    blocks = CommonBlocks.from_sequences([insert, plasmid])
    records = blocks.copy_features_between_common_blocks(inplace=False)
    assert len(records["plasmid"].features) == 6
    assert len(plasmid.features) == 2
    blocks.copy_features_between_common_blocks(inplace=True)
    assert len(plasmid.features) == 6


def test_good_management_of_homologies():
    """This checks for a former obscure bug where a sequence with 2 homologies
    in seq2 corresponding to a single sequence in s1 used to cause an index
    error due to the "*" added by the algorithm to the end of homologies."""
    b1 = random_dna_sequence(4000, seed=123)
    b2 = random_dna_sequence(4000, seed=234)
    b3 = random_dna_sequence(4000, seed=345)
    seq1 = b1 + "A" + "T" + b2 + b3
    seq2 = "T" + b1 + "T" + b3 + b2 + b1 + b1

    blocks = DiffBlocks.from_sequences(seq1, seq2).merged()
    assert len(blocks.blocks) == 10
