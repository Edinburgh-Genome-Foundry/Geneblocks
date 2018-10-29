import os
import matplotlib
matplotlib.use('Agg')

from geneblocks import CommonBlocks, DiffBlocks, load_record
from geneblocks.biotools import reverse_complement, random_dna_sequence

def test_CommonBlocks_basics(tmpdir):
    C1, A, B, C2, C3, D, E, F = [
        random_dna_sequence(100 * L) for L in range(1, 9)]

    sequences = {
        "a": C1 + A + C2,
        "b": B + D + C2 + C3,
        "c": E + C1 + C2,
        "d": C2 + C1 + F + C3,
        "e": C3 + reverse_complement(C2 + C1)
    }

    common_blocks = CommonBlocks(sequences)
    assert len(common_blocks.common_blocks) == 3

    axes = common_blocks.plot_common_blocks()
    fig_path = os.path.join(str(tmpdir), "basic_example.png")
    axes[0].figure.savefig(fig_path, bbox_inches="tight")

    # GET ALL COMMON BLOCKS AS BIOPYTHON RECORDS
    blocks_records = common_blocks.common_blocks_records()
    antiblock_records = common_blocks.unique_blocks_records()

    # WRITE ALL COMMON BLOCKS INTO A CSV SPREADSHEET
    csv_path = os.path.join(str(tmpdir), "basic_example.csv")
    common_blocks.common_blocks_to_csv(target_file=csv_path)


def test_DiffBlocks_basics(tmpdir):
    seq_1 = load_record(os.path.join("tests", "sequences", "sequence1.gb"))
    seq_2 = load_record(os.path.join("tests", "sequences", "sequence2.gb"))

    diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2)
    # next line is just to cover separate_axes=false
    diff_blocks.plot(figure_width=8, separate_axes=False)
    ax1, ax2 = diff_blocks.plot(figure_width=8)
    fig_path = os.path.join(str(tmpdir), "diff_blocks.png")
    ax1.figure.savefig(fig_path, bbox_inches='tight')
    assert list(map(str, diff_blocks.blocks)) == [
        "insert 0-0|0-120",
        "equal 0-1000|120-1120",
        "replace 1000-1004|1120-1124",
        "equal 1004-1503|1124-1623",
        "insert 1503-1503|1623-1723",
        "equal 1503-2304|1723-2524",
        "delete 2304-2404|2524-2524",
        "equal 2404-3404|2524-3524"
    ]

def test_features_transfer():
    seq_folder = os.path.join("tests", "sequences", 'features_transfer')
    insert = load_record(os.path.join(seq_folder, 'insert.gb'), name='insert')
    plasmid = load_record(os.path.join(seq_folder, 'plasmid_to_annotate.gb'),
                          name='plasmid')
    blocks = CommonBlocks([insert, plasmid])
    records = blocks.copy_features_between_common_blocks(inplace=False)
    assert (len(records['plasmid'].features) == 6)
    assert (len(plasmid.features) == 2)
    blocks.copy_features_between_common_blocks(inplace=True)
    assert (len(plasmid.features) == 6)

