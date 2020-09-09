from copy import deepcopy

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from ..Location import Location
from ..biotools import sequence_to_record
from ..CommonBlocks import CommonBlocks

from .DiffBlock import DiffBlock
from .DiffRecordTranslator import DiffRecordTranslator
from .diffblocks_tools import (
    compute_levenshtein_blocks,
    get_optimal_common_blocks,
    merge_blocs_by_location,
    merge_successive_blocks,
    compute_sorted_blocks,
)


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
    def from_sequences(s1, s2, blast_over=500, max_complexity=1e8):
        """Create DiffBlocks by comparing two sequences.

        Parameters
        ----------

        s1, s2
          Two sequences, either "ATGC..." strings or Biopython records.

        blast_over
          A blast will be triggered to accelerate homology finding if
          len(s1) + len(s2) > blast_over.

        max_complexity
          If len(s1) * len(s2) is over max_complexity, no analysis is done and
          s1 is just labeled as a "change" of s2 (useful internally during the
          recursions of this method).
        """
        # Note: the sequences will always be upperized before they are
        # compared. however we also need to keep the initial sequences to
        # create the final blocks (possibly with upper/lowercase nucleotides)

        # If the sequences are records,  convert to string

        seq_s1 = str(s1.seq) if hasattr(s1, "seq") else str(s1)
        seq_s2 = str(s2.seq) if hasattr(s2, "seq") else str(s2)

        # Simple case to eliminate the trivial case of equality

        if seq_s1.upper() == seq_s2.upper():
            return DiffBlocks(s1, s2, [])

        # If the sequences are too big for straight-on Levenshtein, first
        # find the large sub-blocks that are identical, and the ones that
        # differ.

        if (blast_over is not None) and (len(s1) + len(s2)) > blast_over:

            diffblocks = []

            # Use CommonBlocks to find all big common blocks

            sequences = {"s1": s1, "s2": s2}
            common_blocks = CommonBlocks.from_sequences(
                sequences,
                min_block_size=100,
                include_self_homologies=False,
                block_selection_method="larger_first",
            ).common_blocks
            blocks_in_seqs, remarks = get_optimal_common_blocks(common_blocks)

            # First, each common block is added as an "equal" diffblock

            for b1, b2 in zip(blocks_in_seqs["s1"], blocks_in_seqs["s2"]):
                diffblocks.append(
                    DiffBlock(
                        "equal",
                        s1_location=Location(*b1[:2], sequence=s1),
                        s2_location=Location(*b2[:2], sequence=s2),
                    )
                )

            # for sequence in s1, s2, complete the sequence's list of blocks
            # with a (0, 0, "START") on the left, (L, L, "END") on the right.
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
                            s1_location=Location(end1, next_start1, sequence=s1),
                            s2_location=Location(next_start2, next_start2, sequence=s2),
                        )
                    ]
                else:
                    subsequence_1 = s1[end1:next_start1]
                    subsequence_2 = s2[end2:next_start2]
                    subdiffblocks = DiffBlocks.from_sequences(
                        subsequence_1,
                        subsequence_2,
                        blast_over=None,
                        max_complexity=max_complexity,
                    )
                    for block in subdiffblocks.blocks:
                        block.s1_location.start += end1
                        block.s1_location.end += end1
                        block.s1_location.sequence = s1
                        block.s2_location.start += end2
                        block.s2_location.end += end2
                        block.s2_location.sequence = s2
                diffblocks += subdiffblocks.blocks
            diffblocks = [
                b for b in diffblocks if len(b.s1_location) or len(b.s2_location)
            ]
            sorted_blocks = compute_sorted_blocks(diffblocks + remarks)
            return DiffBlocks(s1, s2, sorted_blocks)

        s1_std = str(s1.seq if hasattr(s1, "seq") else s1).upper()
        s2_std = str(s2.seq if hasattr(s2, "seq") else s2).upper()
        levenshtein_blocks = compute_levenshtein_blocks(
            s1_std, s2_std, max_complexity=max_complexity
        )
        blocks = [
            DiffBlock(
                operation,
                Location(s1s, s1e, sequence=s1),
                Location(s2s, s2e, sequence=s2),
            )
            for operation, (s1s, s1e), (s2s, s2e) in levenshtein_blocks
        ]
        return DiffBlocks(s1, s2, blocks)

    def merged(
        self, blocks_per_span=(3, 600), change_gap=100, replace_gap=10, reference="s2",
    ):
        blocks = [b for b in self.blocks if b.operation not in ["reverse", "transpose"]]
        remarks = [b for b in self.blocks if b.operation in ["reverse", "transpose"]]
        if blocks_per_span is not None:
            max_blocks, span = blocks_per_span
            blocks = merge_blocs_by_location(
                blocks=blocks,
                max_blocks=max_blocks,
                max_span=span,
                reference=reference,
            )
        if change_gap is not None:
            blocks = merge_successive_blocks(
                blocks=blocks,
                change_gap=change_gap,
                replace_gap=replace_gap,
                reference="s2",
            )
        blocks = compute_sorted_blocks(blocks + remarks)
        return DiffBlocks(s1=self.s1, s2=self.s2, blocks=blocks)

    def sort_blocks(self):
        self.blocks = compute_sorted_blocks(self.blocks)

    def diffs_as_features(self, sequence="s2"):
        return [block.to_feature(sequence=sequence) for block in self.blocks]

    def plot(
        self, translator_class="default", separate_axes=True, sequence="s2", **plot_kw
    ):

        if translator_class == "default":
            translator_class = DiffRecordTranslator
        translator = translator_class()
        record = deepcopy(self.s2 if sequence == "s2" else self.s1)
        if not hasattr(record, "features"):
            record = sequence_to_record(record)

        diff_features = self.diffs_as_features(sequence=sequence)

        if separate_axes:
            gr_record = translator.translate_record(record)
            record.features = diff_features
            gr_diffrecord = DiffRecordTranslator().translate_record(record)
            width = plot_kw.get("figure_width", 8)
            if "axes" in plot_kw:
                ax1, ax2 = plot_kw["axes"]
                fig = ax1.figure
            else:
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(width, 6))
            plot_kw["annotate_inline"] = plot_kw.get("annotate_inline", True)
            _, stats1 = gr_record.plot(ax=ax1, **plot_kw)
            _, stats2 = gr_diffrecord.plot(ax=ax2, with_ruler=False, **plot_kw)
            max_features_1 = gr_record.feature_level_height * max(
                [0] + [v for v in stats1[0].values()]
            )
            max_level_1 = max(
                [max_features_1] + [v["annotation_y"] for v in stats1[1].values()]
            )
            max_level_2 = max([1] + [v["annotation_y"] for v in stats2[1].values()]) + 2
            max_level_1 = int(np.round(max_level_1))
            max_level_2 = int(np.round(max_level_2))
            # print (stats2)
            n_levels = max_level_1 + max_level_2
            if max_level_1 and max_level_2:
                plt.close(fig)

                ##
                easing = 3
                gs = gridspec.GridSpec(n_levels + 2 * easing, 1)
                fig = plt.figure(figsize=(width, 1 + 0.5 * n_levels), facecolor="w")
                ax1 = fig.add_subplot(gs[: max_level_1 + easing])
                ax2 = fig.add_subplot(gs[max_level_1 + easing :])
                _, stats1 = gr_record.plot(ax=ax1, **plot_kw)
                _, stats2 = gr_diffrecord.plot(ax=ax2, with_ruler=False, **plot_kw)

            # fig.set_size_inches((width, 3 + 0.4 * n_levels))
            ax2.set_ylim(bottom=-2)
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
            if block.operation in ("equal", "replace", "change", "delete"):
                s1 = s1 + block.s1_location.extract_sequence()
            if block.operation in ("equal", "replace", "change", "insert"):
                s2 = s2 + block.s2_location.extract_sequence()
        return s1, s2

    def __str__(self):
        return ", ".join([str(b) for b in self.blocks])
