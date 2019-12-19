import os
from geneblocks import CommonBlocks, load_record
from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt

# LOAD THE TWO RECORDS

part_path = os.path.join("sequences", "features_transfer", "part.gb")
part = load_record(part_path, name="part")
plasmid_path = os.path.join(
    "sequences", "features_transfer", "plasmid_to_annotate.gb"
)
plasmid = load_record(plasmid_path, name="plasmid")


# TRANSFER THE FEATURES

blocks = CommonBlocks.from_sequences([part, plasmid])
new_records = blocks.copy_features_between_common_blocks(inplace=False)
annotated_plasmid = new_records["plasmid"]  # record with all features


# PLOT ALL RECORDS


fig, (ax0, ax1, ax2) = plt.subplots(3, 1, sharex=True, figsize=(7, 5))
ax0.set_title("Record 1: Annotated part", loc="left")
ax1.set_title("Record 2: Partially annotated plasmid", loc="left")
ax2.set_title("Record 2 (after) with all annotations", loc="left")


class SpecialBiopythonTranslator(BiopythonTranslator):
    def compute_feature_color(self, f):
        original_record = "".join(f.qualifiers.get("original_record", ""))
        return "#ffafaf" if (original_record == "part") else "#afafff"


tr = SpecialBiopythonTranslator()
tr.translate_record(part).plot(ax0, with_ruler=False, x_lim=(0, len(plasmid)))
tr.translate_record(plasmid).plot(ax1, with_ruler=False)
tr.translate_record(annotated_plasmid).plot(ax2, with_ruler=False)

fig.tight_layout()
fig.subplots_adjust(hspace=0.6)

fig.savefig("features_transfer.png", dpi=125)
