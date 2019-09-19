from copy import deepcopy
from difflib import SequenceMatcher
from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
import networkx as nx

from .Location import Location
from .biotools import sequence_to_record, sequences_differences
from .CommonBlocks import CommonBlocks


class DiffBlock:
    def __init__(self, operation, s1_location, s2_location):
        self.operation = operation
        self.s1_location = s1_location
        self.s2_location = s2_location

    def to_feature(self):
        s1_length, s2_length = len(self.s1_location), len(self.s2_location)
        max_length = max([s1_length, s2_length])

        if self.operation == "insert":
            if max_length <= 7:
                label = "+%s" % self.s2_location.extract_sequence()
            else:
                label = "+ %d nuc." % s2_length
        elif self.operation == "delete":
            if max_length <= 7:
                label = "-%s" % self.s1_location.extract_sequence()
            else:
                label = "- %d nuc." % s1_length
        elif self.operation == "replace":
            if max([s1_length, s2_length]) <= 6:
                label = "%s ➤ %s" % (
                    self.s1_location.extract_sequence(),
                    self.s2_location.extract_sequence(),
                )
            else:
                if s1_length == s2_length:
                    sub_s1 = self.s1_location.extract_sequence()
                    sub_s2 = self.s2_location.extract_sequence()
                    diffs = sequences_differences(sub_s1, sub_s2)
                    label = "%d/%dn changed" % (diffs, s1_length)
                else:
                    label = "%sn ➤ %sn change" % (s1_length, s2_length)
        elif self.operation == "reverse":
            label = "was reversed at %d-%d" % (
                self.s1_location.start,
                self.s1_location.end,
            )
        elif self.operation == "transpose":
            label = "was at %d-%d" % (
                self.s1_location.start,
                self.s1_location.end,
            )
        elif self.operation == "equal":
            label = "Equal"

        return self.s2_location.to_biopython_feature(
            feature_type="diff_" + self.operation, label=label
        )

    def __str__(self):
        return "%s %s|%s" % (
            self.operation,
            self.s1_location,
            self.s2_location,
        )

    def __repr__(self):
        return str(self)

    def trim_replace_block(self):
        s1 = self.s1_location.extract_sequence()
        s2 = self.s2_location.extract_sequence()
        start, end = 0, 0
        for start in range(1, min(len(s1), len(s2))):
            if s1[:start] != s2[:start]:
                start = start - 1
                break
        for end in range(1, min(len(s1), len(s2))):
            if s1[-end:] != s2[-end:]:
                end = end - 1
                break
        if start > 0:
            left_block = DiffBlock(
                operation="equal",
                s1_location=Location(
                    self.s1_location.start,
                    self.s1_location.start + start,
                    sequence=self.s1_location.sequence,
                ),
                s2_location=Location(
                    self.s2_location.start,
                    self.s2_location.start + start,
                    sequence=self.s2_location.sequence,
                ),
            )
            self.s1_location.start += start
            self.s2_location.start += start
        else:
            left_block = None
        if end > 0:
            right_block = DiffBlock(
                operation="equal",
                s1_location=Location(
                    self.s1_location.end - end,
                    self.s1_location.end,
                    sequence=self.s1_location.sequence,
                ),
                s2_location=Location(
                    self.s2_location.end - end,
                    self.s2_location.end,
                    sequence=self.s2_location.sequence,
                ),
            )
            self.s1_location.end -= end
            self.s2_location.end -= end
        else:
            right_block = None
        return left_block, right_block


class DiffRecordTranslator(BiopythonTranslator):

    ignored_features_types = ("diff_equal",)

    @staticmethod
    def compute_feature_color(f):
        return dict(
            diff_delete="red",
            diff_insert="green",
            diff_replace="orange",
            diff_reverse="white",
            diff_transpose="white",
        ).get(f.type, "white")

    @staticmethod
    def compute_feature_box_linewidth(f):
        return 1 if f.type.startswith("diff_") else 0

    @staticmethod
    def compute_feature_fontdict(f):
        return {"fontsize": 12 if f.type.startswith("diff_") else 8}


class DiffBlocks:
    """Class to generate and represent DiffBlocks.

    Usage:

    >>> DiffBlocks.from_sequences(s1, s2)
    """

    def __init__(self, s1, s2, blocks):
        self.s1 = s1
        self.s2 = s2
        self.blocks = blocks

    @staticmethod
    def from_sequences(
        s1, s2, contract_under=3, blast_over=5000, use_junk_over=100
    ):

        # Simple case to eliminate the trivial case of equality
        seq_s1 = str(s1.seq) if hasattr(s1, "seq") else str(s1)
        seq_s2 = str(s2.seq) if hasattr(s2, "seq") else str(s2)
        if seq_s1.upper() == seq_s2.upper():
            return DiffBlocks(s1, s2, [])

        if (blast_over is not None) and (len(s1) + len(s2)) > blast_over:
            sequences = {"s1": s1, "s2": s2}
            common_blocks = CommonBlocks(
                sequences,
                ignore_self_homologies=True,
                block_selection="larger_first",
            ).common_blocks
            blocks_in_seqs, remarks = get_optimal_common_blocks(common_blocks)
            diffblocks = []
            for b1, b2 in zip(blocks_in_seqs["s1"], blocks_in_seqs["s2"]):
                diffblocks.append(
                    DiffBlock(
                        "equal",
                        s1_location=Location(*b1[:2], sequence=s1),
                        s2_location=Location(*b2[:2], sequence=s2),
                    )
                )
            for seq, blocks in blocks_in_seqs.items():
                blocks_in_seqs[seq] = (
                    [(0, 0, "START")]
                    + blocks_in_seqs[seq]
                    + [(len(sequences[seq]), len(sequences[seq]), "END")]
                )

            for i in range(len(blocks_in_seqs["s2"]) - 1):
                _, end1, _ = blocks_in_seqs["s1"][i]
                next_start1, _, _ = blocks_in_seqs["s1"][i + 1]
                _, end2, _ = blocks_in_seqs["s2"][i]
                next_start2, _, _ = blocks_in_seqs["s2"][i + 1]
                if next_start2 < end2:
                    subdiffblocks = [
                        DiffBlock(
                            "delete",
                            s1_location=Location(
                                end1, next_start1, sequence=s1
                            ),
                            s2_location=Location(
                                next_start2, next_start2, sequence=s2
                            ),
                        )
                    ]
                else:
                    subsequence_1 = s1[end1:next_start1]
                    subsequence_2 = s2[end2:next_start2]
                    subdiffblocks = DiffBlocks.from_sequences(
                        subsequence_1,
                        subsequence_2,
                        contract_under=contract_under,
                        blast_over=None,
                    ).blocks
                    for block in subdiffblocks:
                        block.s1_location.start += end1
                        block.s1_location.end += end1
                        block.s1_location.sequence = s1
                        block.s2_location.start += end2
                        block.s2_location.end += end2
                        block.s2_location.sequence = s2
                diffblocks += subdiffblocks

            sorted_blocks = sorted(
                diffblocks + remarks, key=lambda b: b.s2_location.to_tuple()
            )
            return DiffBlocks(s1, s2, sorted_blocks)

        upper_s1, upper_s2 = s1.upper(), s2.upper()
        matcher = SequenceMatcher(a=upper_s1, b=upper_s2, autojunk=False)
        ratio = matcher.quick_ratio()
        if (max(len(s1), len(s2)) * (1.0 - ratio)) > use_junk_over:
            matcher = SequenceMatcher(a=upper_s1, b=upper_s2, autojunk=True)
        blocks = [
            DiffBlock(
                operation,
                Location(s1s, s1e, sequence=s1),
                Location(s2s, s2e, sequence=s2),
            )
            for operation, s1s, s1e, s2s, s2e in matcher.get_opcodes()
        ]
        diffblocks = DiffBlocks(s1, s2, blocks)
        if contract_under:
            diffblocks.contract_subblocks_as_replace(
                max_replace_size=contract_under
            )
        diffblocks.trim_all_replace_blocks()
        return diffblocks

    def trim_all_replace_blocks(self):
        for block in list(self.blocks):
            if block.operation == "replace":
                left, right = block.trim_replace_block()
                if left is not None:
                    self.blocks.append(left)
                if right is not None:
                    self.blocks.append(right)
        self.sort_blocks()
        self.merge_successive_equal_blocks()

    def sort_blocks(self):
        self.blocks = sorted(
            self.blocks, key=lambda b: b.s2_location.to_tuple()
        )

    def merge_successive_equal_blocks(self):
        def contract_two_blocks(blocks):
            blocks = list(blocks)
            for i in range(len(blocks) - 1):
                b1, b2 = blocks[i], blocks[i + 1]
                if (b1.operation == "equal") and (b2.operation == "equal"):
                    b1.s1_location.end = b2.s1_location.end
                    b1.s2_location.end = b2.s2_location.end
                    blocks.remove(b2)
                    return blocks
            return blocks

        self.sort_blocks()
        blocks = self.blocks
        while True:
            new_blocks = contract_two_blocks(blocks)
            if len(new_blocks) == len(blocks):
                break
            blocks = new_blocks
        self.blocks = blocks

    def diffs_as_features(self):
        return [block.to_feature() for block in self.blocks]

    def plot(self, translator_class="default", separate_axes=True, **plot_kw):

        if translator_class == "default":
            translator_class = DiffRecordTranslator
        translator = translator_class()
        record = deepcopy(self.s2)
        if not hasattr(record, "features"):
            record = sequence_to_record(record)

        diff_features = self.diffs_as_features()

        if separate_axes:
            gr_record = translator.translate_record(record)
            record.features = diff_features
            gr_diffrecord = translator.translate_record(record)
            width = plot_kw.get("figure_width", 8)
            if "axes" in plot_kw:
                ax1, ax2 = plot_kw["axes"]
                fig = ax1.figure
            else:
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(width, 6))
            _, stats1 = gr_record.plot(ax=ax1, **plot_kw)
            _, stats2 = gr_diffrecord.plot(ax=ax2, with_ruler=False, **plot_kw)
            max_level_1 = max(
                [0] + [v["annotation_y"] for v in stats1[1].values()]
            )
            max_features_1 = max(
                [0] + [v["feature_y"] for v in stats1[1].values()]
            )
            max_level_2 = max(
                [0] + [v["annotation_y"] for v in stats2[1].values()]
            )
            n_levels = max_level_1 + max_level_2
            fig.set_size_inches((width, 0.6 * n_levels))
            ax2.set_ylim(bottom=-0.5)
            ax2.invert_yaxis()
            for f in gr_diffrecord.features:
                ax1.fill_between(
                    [f.start, f.end],
                    y1=max_features_1 + 1,
                    y2=-1,
                    facecolor=f.color,
                    alpha=0.07,
                    zorder=1000,
                )
            return (ax1, ax2)

        else:
            record.features += diff_features
            gr_record = translator.translate_record(record)
            ax, _ = gr_record.plot(**plot_kw)
            return ax

    @staticmethod
    def reconstruct_sequences_from_blocks(blocks):
        s1, s2 = "", ""
        blocks = sorted(blocks, key=lambda b: b.s2_location.to_tuple())
        for block in blocks:
            if block.operation in ("equal", "replace", "delete"):
                s1 = s1 + block.s1_location.extract_sequence()
            if block.operation in ("equal", "replace", "insert"):
                s2 = s2 + block.s2_location.extract_sequence()
        return s1, s2

    def contract_subblocks_as_replace(self, max_replace_size=4):
        def contract_one_block_tuple(blocks, tuple_size):
            for i in range(len(blocks) - tuple_size):
                subblocks = blocks[i : i + tuple_size]
                if subblocks[-1].operation == "equal":
                    subblocks = subblocks[:-1]
                if subblocks[0].operation == "equal":
                    subblocks = subblocks[1:]
                if len(subblocks) < 2:
                    continue
                s1_location = Location(
                    min([b.s1_location.start for b in subblocks]),
                    max([b.s1_location.end for b in subblocks]),
                    sequence=self.s1,
                )
                s2_location = Location(
                    min([b.s2_location.start for b in subblocks]),
                    max([b.s2_location.end for b in subblocks]),
                    sequence=self.s2,
                )
                if max(len(s1_location), len(s2_location)) <= max_replace_size:
                    new_block = DiffBlock("replace", s1_location, s2_location)
                    # left_block, right_block = new_block.trim_replace_block()
                    return blocks[:i] + [new_block] + blocks[i + tuple_size :]
            return blocks

        blocks = self.blocks
        for tuple_size in (4, 3, 2):
            while True:
                new_blocks = contract_one_block_tuple(blocks, tuple_size)
                if len(new_blocks) == len(blocks):
                    break
                blocks = new_blocks

        self.blocks = blocks


def get_optimal_common_blocks(common_blocks):
    print("there")
    common_blocks = deepcopy(common_blocks)
    remarks = []

    # Make so that there is never an anitsense block in s1 and a + block in s2.
    # If it is so, flip the block in s2. It will become antisense and be later
    # removed
    for block_name, data in common_blocks.items():
        locations = data["locations"]
        s1_strands = [strand for (s, (_, _, strand)) in locations if s == "s1"]
        if 1 not in s1_strands:
            for i, location in enumerate(locations):
                seq, (start, end, strand) = location
                if seq == "s2":
                    locations[i] = (seq, (start, end, -strand))

    # Remove every antisense blocks now. For the ones in s2, log this
    # with a remark.
    for block_name, data in common_blocks.items():
        locations = data["locations"]
        for i, location in enumerate(locations):
            (seq, (start, end, strand)) = location
            if seq == "s2" and (strand == -1):
                locations.remove(location)
                _, (start1, end1, strand1) = locations[0]
                remarks.append(
                    DiffBlock(
                        "reverse",
                        s1_location=Location(start1, end1, strand1),
                        s2_location=Location(start, end),
                    )
                )

    # We start the structure that will be returned in the end
    blocks_in_seqs = {
        seq: sorted(
            [
                (start, end, bname)
                for bname, data in common_blocks.items()
                for (s, (start, end, strand)) in data["locations"]
                if s == seq
            ]
        )
        for seq in ("s1", "s2")
    }

    # Identify blocks appearing only in one of the two sequences
    blocks_in_s1 = set(b[-1] for b in blocks_in_seqs["s1"])
    blocks_in_s2 = set(b[-1] for b in blocks_in_seqs["s2"])
    uniblocks = (blocks_in_s1.union(blocks_in_s2)).difference(
        blocks_in_s1.intersection(blocks_in_s2)
    )

    # Remove blocks appearing only in one of the two sequences
    # as they are useless for sequences comparison
    # this should be very rare but you never know.
    for block_list in blocks_in_seqs.values():
        for b in block_list:
            if b[-1] in uniblocks:
                block_list.remove(b)

        for b1, b2 in zip(block_list, block_list[1:]):
            start1, end1, __name1 = b1
            start2, end2, __name2 = b2
            if end2 <= end1:
                block_list.remove(b2)

    # If a block appears several time in a sequence (self-homology)
    # give unique names to each occurence: block_1, block_1*, etc.
    blocks_in_seqs_dicts = dict(s1={}, s2={})
    for seq, blocks_list in list(blocks_in_seqs.items()):
        seen_blocks = set()
        for i, (start, end, block_name) in enumerate(blocks_list):
            while block_name in seen_blocks:
                block_name = block_name + "*"
            blocks_list[i] = start, end, block_name
            blocks_in_seqs_dicts[seq][block_name] = dict(
                rank=i, location=(start, end)
            )
            seen_blocks.add(block_name)

    # Find and retain the largest sequence of blocks which is in the right order
    # in both sequences. We will remove every other blocks

    if len(blocks_in_seqs_dicts["s2"]) < 2:
        retained_blocks = list(blocks_in_seqs_dicts["s2"])
    else:
        s1_dict = blocks_in_seqs_dicts["s1"]
        graph = nx.DiGraph(
            [
                (b1, b2)
                for b1, data1 in blocks_in_seqs_dicts["s2"].items()
                for b2, data2 in blocks_in_seqs_dicts["s2"].items()
                if (b2 in s1_dict)
                and (b1 in s1_dict)
                and (s1_dict[b2]["rank"] > s1_dict[b1]["rank"])
                and (data2["rank"] > data1["rank"])
            ]
        )
        retained_blocks = nx.dag_longest_path(graph)

    # remove any "misplaced" block that is not in the retained list.
    # log a remark for the ones in s2.
    for seq in ("s1", "s2"):
        blocks_list = blocks_in_seqs[seq]
        for block in list(blocks_list):  # copy cause we will remove elements
            start, end, block_name = block
            if block_name not in retained_blocks:
                blocks_list.remove(block)
                if seq == "s2":
                    new_block_name = block_name.strip("*")
                    s1_blocks = blocks_in_seqs_dicts["s1"]
                    s1_loc = s1_blocks[new_block_name]["location"]
                    start1, end1 = s1_loc
                    if (
                        len(remarks)
                        and (start == remarks[-1].s2_location.end)
                        and (start1 == remarks[-1].s1_location.end)
                    ):
                        remarks[-1].s1_location.end = end1
                        remarks[-1].s2_location.end = end
                    else:
                        remarks.append(
                            DiffBlock(
                                "transpose",
                                s1_location=Location(start1, end1),
                                s2_location=Location(start, end),
                            )
                        )

    # Reduce blocks when there is overlap
    blocks_to_reduce = {}
    for seq in ("s1", "s2"):
        blocks_list = blocks_in_seqs[seq]
        for b1, b2 in zip(blocks_list, blocks_list[1:]):
            start1, end1, block_name1 = b1
            start2, end2, __block_name2 = b2
            diff = end1 - start2
            if diff > 0:
                if block_name1 not in blocks_to_reduce:
                    blocks_to_reduce[block_name1] = 0
                blocks_to_reduce[block_name1] = max(
                    blocks_to_reduce[block_name1], diff
                )

    for seq in ("s1", "s2"):
        blocks_list = blocks_in_seqs[seq]
        for i, (start, end, block_name) in enumerate(blocks_list):
            if block_name in blocks_to_reduce:
                diff = blocks_to_reduce[block_name]
                blocks_list[i] = (start, end - diff, block_name)

    return blocks_in_seqs, remarks
