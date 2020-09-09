"""Defines central class BlockFinder."""
import itertools
from collections import OrderedDict
from copy import deepcopy


from .CommonBlocksRecordTranslator import CommonBlocksRecordTranslator
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from ..biotools import annotate_record

from .commonblocks_tools import (
    format_sequences_as_dicts,
    select_common_blocks,
    find_homologies_between_sequences,
)

# TODO: Simplify the code by using the new Location class in Location.py


class CommonBlocks:
    """Class to represent a set of common blocks from different sequences.

    Create with ``CommonBlocks.from_sequences``:

    >>> common_blocks = CommonBlocks.from_sequences({'s1': 'ATGC...'})

    Parameters
    ----------

    common_blocks
      A dictionary of the sequences to compare, of the form
      {sequence_name: ATGC_sequence_string} or a list of records, all with
      different IDs.

    records
      A dictionary of the Biopython records of the sequences
      {record_id: record}.
    """

    def __init__(self, common_blocks, records):
        """Initialize, compute best blocks."""
        self.common_blocks = common_blocks
        self.records = records

    @staticmethod
    def from_sequences(
        sequences,
        block_selection_method="most_coverage_first",
        include_self_homologies=True,
        min_block_size=80,
        max_block_size=None,
    ):
        sequences_dict, records_dict = format_sequences_as_dicts(sequences)
        homologies_dict = find_homologies_between_sequences(
            sequences_dict,
            min_size=min_block_size,
            max_size=max_block_size,
            include_self_homologies=include_self_homologies,
        )
        common_blocks = select_common_blocks(
            homologies_dict,
            sequences_dict,
            min_size=min_block_size,
            method=block_selection_method,
        )
        return CommonBlocks(common_blocks=common_blocks, records=records_dict)

    def compute_unique_blocks(self):
        """Return a dictionary listing unique blocks by sequence.

        The unique blocks are the blocks between the selected common blocks.

        The result is of the form {seq: [(start, end), (start2, end2), ...]}
        """

        unique_blocks = OrderedDict()
        for seqname, rec in self.sequences_with_annotated_blocks().items():
            blocks_locations = (
                [(0, 0)]
                + sorted(
                    [
                        (f.location.start, f.location.end)
                        for f in rec.features
                        if f.qualifiers.get("is_block", False)
                    ]
                )
                + [(len(rec), len(rec))]
            )
            unique_blocks[seqname] = [
                (end1, start2)
                for (_, end1), (start2, _) in zip(
                    blocks_locations, blocks_locations[1:]
                )
                if (start2 - end1) > 1
            ]
        return unique_blocks

    def common_blocks_to_csv(self, target_file=None):
        """Write the common blocks into a CSV file.

        If a target CSV file is provided the result is written to that file.
        Otherwise the result is returned as a string.

        The columns of the CSV file are "block", "size", "locations", and
        sequence.
        """
        csv_content = "\n".join(
            ["block;size;locations;sequence"]
            + [
                ";".join(
                    [
                        block_name,
                        str(len(data["sequence"])),
                        " ".join(
                            [
                                "%s(%d, %d, %d)" % (cst, start, end, strand)
                                for (cst, (start, end, strand)) in data["locations"]
                            ]
                        ),
                        data["sequence"],
                    ]
                )
                for block_name, data in self.common_blocks.items()
            ]
        )
        if target_file:
            with open(target_file, "w+") as f:
                f.write(csv_content)
        else:
            return csv_content

    def common_blocks_records(self):
        """Return all common blocks as a list of Biopython records.
        """
        if self.records is None:
            raise ValueError("")
        records = []
        for block_name, data in self.common_blocks.items():
            cst, (start, end, strand) = data["locations"][0]
            record = self.records[cst][start:end]
            if strand == -1:
                record = record.reverse_complement()
            record.id = record.name = block_name
            records.append(record)
        return records

    def unique_blocks_records(self, target_file=None):
        """Return all unique blocks as a list of Biopython records."""
        if self.records is None:
            raise ValueError("")
        records = []
        for seqname, locations in self.compute_unique_blocks().items():
            for i, (start, end) in enumerate(locations):
                record = self.records[seqname][start:end]
                record.id = "%s_%03d" % (seqname, i)
                records.append(record)
        return records

    def sequences_with_annotated_blocks(self, colors="auto"):
        """Return a list of Biopython records representing the sequences
        with annotations indicating the common blocks.

        Parameter ``colors`` is either a list of colors or "auto" for the
        default.
        """
        records = deepcopy(self.records)
        if colors == "auto":
            colors = itertools.cycle([cm.Paired(0.21 * i % 1.0) for i in range(30)])
        blocks_and_colors = zip(self.common_blocks.items(), colors)
        for (name, data), color in blocks_and_colors:
            for (seqname, location) in data["locations"]:
                annotate_record(
                    records[seqname],
                    location,
                    feature_type="misc_feature",
                    is_block=True,
                    label=name,
                    color=color,
                )
        return records

    def plot_common_blocks(
        self, colors="auto", axes=None, figure_width=10, ax_height=2
    ):
        """Plot the common blocks found on vertically stacked axes.

        The axes on which the plots are drawn are returned at the end.

        Parameters
        ----------

        colors
          Either a list of colors to use for blocks or "auto" for the default.

        axes
          A list of matplotlib axes on which to plot, or None for new axes.

        figure_width
          Width of the final figure in inches.

        ax_eight
          Height of each plot.
        """

        translator = CommonBlocksRecordTranslator()
        records = self.sequences_with_annotated_blocks(colors=colors)
        if axes is None:
            fig, axes = plt.subplots(
                len(self.records),
                1,
                facecolor="white",
                sharex=True,
                figsize=(figure_width, ax_height * len(self.records)),
            )
        else:
            fig = axes[0].figure
        for (ax, (seqname, record)) in zip(axes, records.items()):
            gr_record = translator.translate_record(record)
            gr_record.plot(
                ax,
                x_lim=(0, max([len(rec) for rec in self.records.values()])),
                with_ruler=(ax == axes[-1]),
            )
            ax.set_ylim(top=ax.get_ylim()[1])
            ax.set_title(seqname, loc="left", fontdict=dict(weight="bold"))
        # fig.tight_layout()
        return axes

    def copy_features_between_common_blocks(self, inplace=False):
        def extract_subrecord(record, location):
            start, end, strand = location
            record = record[start:end]
            if strand == -1:
                record = record.reverse_complement()
            return record

        def extract_features(record, offset, reverse=False):
            if reverse:
                record = record.reverse_complement()
            new_features = [deepcopy(f) for f in record.features]
            for f in new_features:
                f.qualifiers["original_record"] = record.id
            for f in new_features:
                f.location += offset
            return new_features

        if inplace:
            records = self.records
        else:
            records = deepcopy(self.records)
        for data in self.common_blocks.values():
            locations = data["locations"]
            subrecords = {
                rec_id: extract_subrecord(records[rec_id], location)
                for rec_id, location in data["locations"]
            }
            for l1, l2 in itertools.combinations(locations, 2):
                for ((id1, loc1), (id2, __loc2)) in ((l1, l2), (l2, l1)):
                    start1, __end1, strand1 = loc1
                    # start2, end2, strand2 = loc2
                    records[id1].features += extract_features(
                        subrecords[id2], offset=start1, reverse=(strand1 == -1)
                    )
        return records
