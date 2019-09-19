"""Defines central class BlockFinder."""
import subprocess
import itertools
from collections import defaultdict, OrderedDict
import tempfile
import numpy as np
import re
from copy import deepcopy


from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from .biotools import reverse_complement, sequence_to_record, annotate_record

# TODO: Simplify the code by using the new Location class in Location.py


class CommonBlocks:
    """CommonBlocks.

    Parameters
    ----------

    sequences
      A dictionnary of the sequences to compare, of the form
      {sequence_name: ATGC_sequence_string} or a list of records, all with
      different IDs.

    blocks_selection
      Either 'most_coverage_first' or 'larger_first'. If 'larger_first', the
      largest blocks (in number of nucleotides). If 'most_coverage_first',
      we pick first the blocks whose size times number of occurences is larger

    min_block_size
      Blocks with a size below this threshold will be ignored and the blocks
      search stops when all blocks are below this threshold.
    """

    def __init__(
        self,
        sequences,
        block_selection="most_coverage_first",
        ignore_self_homologies=False,
        min_block_size=80,
        max_block_size=None,
    ):
        """Initialize, compute best blocks."""
        self.block_selection = block_selection
        self.min_block_size = min_block_size
        self.max_block_size = max_block_size
        self.ignore_self_homologies = ignore_self_homologies
        if isinstance(sequences, (list, tuple)):
            if hasattr(sequences[0], "seq"):
                self.records = OrderedDict(
                    [(record.id, record) for record in sequences]
                )
                sequences = OrderedDict(
                    [
                        (record.id, str(record.seq).upper())
                        for record in sequences
                    ]
                )
            else:
                sequences = OrderedDict(sequences)
                if isinstance(list(sequences.values())[0], str):
                    self.records = OrderedDict(
                        [
                            (name, sequence_to_record(seq, name=name))
                            for name, seq in sequences.items()
                        ]
                    )
                else:
                    self.records = sequences
        elif hasattr(list(sequences.values())[0], "seq"):
            self.records = sequences
            sequences = OrderedDict(
                [
                    (record_id, str(record.seq).upper())
                    for record_id, record in sequences.items()
                ]
            )
        else:
            self.records = OrderedDict(
                [
                    (name, sequence_to_record(seq, name=name))
                    for name, seq in sequences.items()
                ]
            )
        if isinstance(sequences, dict):
            sequences = OrderedDict(sorted(sequences.items()))
        self.sequences = sequences
        self.max_sequence_length = max(
            [len(seq) for seq in self.sequences.values()]
        )

        self._find_intermatches()
        self._find_common_blocks()
        self._find_unique_blocks()

    def _find_intermatches(self):
        temp_fasta_path = tempfile.mktemp(".fa")

        with open(temp_fasta_path, "w+") as f:
            f.write(
                "\n\n".join(
                    [
                        "> %s\n%s" % (name, seq.upper())
                        for name, seq in self.sequences.items()
                    ]
                )
            )
        proc = subprocess.Popen(
            [
                "blastn",
                "-query",
                temp_fasta_path,
                "-subject",
                temp_fasta_path,
                "-perc_identity",
                "100",
                "-dust",
                "no",
                "-evalue",
                "100000000000",
                "-culling_limit",
                "10",
                "-ungapped",
                "-outfmt",
                "6 qseqid qstart qend sseqid sstart send",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        result, __blast_err = proc.communicate()
        parsing = [
            line.split("\t") for line in result.decode("utf-8").splitlines()
        ]
        self.intermatches = {
            name: defaultdict(lambda *a: [])
            for name, seq in self.sequences.items()
        }
        for query, qstart, qend, subject, sstart, send in parsing:
            if self.ignore_self_homologies and (
                (query == subject) and qstart != sstart
            ):
                continue
            qstart, qend = int(qstart) - 1, int(qend)
            sstart, send = int(sstart) - 1, int(send)
            if qend - qstart < self.min_block_size:
                continue
            max_size = self.max_block_size
            if (max_size is not None) and (qend - qstart > max_size):
                continue
            location = (subject, sstart, send)
            self.intermatches[query][(qstart, qend)].append(location)

    def _compute_intersections(self, matches_dict):
        intersections = {}
        if len(matches_dict) == 1:
            intersection = list(matches_dict.keys())[0]
            intersections[intersection] = 1
        matches_list = sorted(matches_dict.keys())
        for i, match1 in enumerate(matches_list):
            for match2 in matches_list[i + 1 :]:
                intersection = start, end = (
                    match2[0],
                    min(match1[1], match2[1]),
                )
                if end < start:
                    break
                elif (end - start > self.min_block_size) and (
                    intersection not in intersections
                ):
                    intersections[intersection] = len(
                        [
                            matching
                            for (match_start, match_end) in matches_list
                            for matching in matches_dict[
                                (match_start, match_end)
                            ]
                            if match_start <= start <= end <= match_end
                        ]
                    )
        return intersections

    def _best_intersection(self, intersections):
        def intersection_score(intersection):
            start, end = intersection
            if self.block_selection == "most_coverage_first":
                f = intersections[intersection]
            else:
                f = 1
            return f * (end - start)

        return max(
            [(0, (None, None))]
            + [
                (intersection_score(intersection), intersection)
                for intersection in intersections
            ]
        )

    @staticmethod
    def _subtract_from_segment(segment, subtracted):
        """Return the difference between ``segment`` and ``subtracted``."""
        seg_start, seg_end = segment
        sub_start, sub_end = subtracted
        result = []
        if sub_start > seg_start:
            result.append((seg_start, min(sub_start, seg_end)))
        if sub_end < seg_end:
            result.append((max(seg_start, sub_end), seg_end))
        return sorted(list(set(result)))

    def _find_common_blocks(self):
        """Find the largest common blocks, iteratively."""
        common_blocks = []
        best_score = 1
        all_intersections = {
            seqname: self._compute_intersections(self.intermatches[seqname])
            for seqname in self.sequences
        }
        while True:
            global_best_intersection = max(
                [
                    (
                        self._best_intersection(all_intersections[seqname]),
                        seqname,
                    )
                    for seqname in self.sequences
                ]
            )
            (best_score, (start, end)), seqname = global_best_intersection
            if best_score == 0:
                break
            best_sequence = self.sequences[seqname][start:end]
            best_sequence_rev = reverse_complement(best_sequence)
            locations = []
            for seqname, sequence in self.sequences.items():
                for bestseq, strand in [
                    (best_sequence, 1),
                    (best_sequence_rev, -1),
                ]:
                    for match in re.finditer(bestseq, sequence):
                        start, end = match.start(), match.end()
                        locations.append((seqname, (start, end, strand)))
                        match_segment = tuple(sorted([start, end]))
                        intersections = list(all_intersections[seqname].keys())
                        for intersection in intersections:
                            score = all_intersections[seqname].pop(
                                intersection
                            )
                            for diff in self._subtract_from_segment(
                                intersection, match_segment
                            ):
                                diff_start, diff_end = diff
                                if diff_end - diff_start > self.min_block_size:
                                    all_intersections[seqname][diff] = score
            common_blocks.append((best_sequence, locations))

        # trick to remove self-homologous sequences
        common_blocks = [
            (seq, locations)
            for (seq, locations) in common_blocks
            if len(locations) >= 2
        ]

        self.common_blocks = OrderedDict()
        if len(common_blocks) > 0:
            number_size = int(np.log10(len(common_blocks))) + 1
            for i, (sequence, locations) in enumerate(common_blocks):
                block_name = "block_%s" % (str(i + 1).zfill(number_size))
                self.common_blocks[block_name] = {
                    "sequence": sequence,
                    "locations": locations,
                }

    def _find_unique_blocks(self):
        self.unique_blocks = OrderedDict()
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
            self.unique_blocks[seqname] = [
                (end1, start2)
                for (_, end1), (start2, _) in zip(
                    blocks_locations, blocks_locations[1:]
                )
                if (start2 - end1) > 1
            ]

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
                                for (cst, (start, end, strand)) in data[
                                    "locations"
                                ]
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
        for seqname, locations in self.unique_blocks.items():
            for i, (start, end) in enumerate(locations):
                record = self.records[seqname][start:end]
                record.id = "%s_%03d" % (seqname, i)
                records.append(record)
        return records

    def sequences_with_annotated_blocks(self, colors="auto"):
        """Return a list of Biopython records representing the sequences
        with annotations indicating the common blocks

        Parameter ``colors`` is either a list of colors or "auto" for the
        default.
        """
        records = deepcopy(self.records)
        if colors == "auto":
            colors = itertools.cycle(
                [cm.Paired(0.21 * i % 1.0) for i in range(30)]
            )
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
          A list of matplotlib axes on which to plot, or None for new axes

        figure_width
          Width of the final figure in inches

        ax_eight
          Height of each plot
        """

        translator = BiopythonTranslator()
        records = self.sequences_with_annotated_blocks(colors=colors)
        if axes is None:
            fig, axes = plt.subplots(
                len(self.sequences),
                1,
                facecolor="white",
                sharex=True,
                figsize=(figure_width, ax_height * len(self.sequences)),
            )
        else:
            fig = axes[0].figure
        for (ax, (seqname, record)) in zip(axes, records.items()):
            gr_record = translator.translate_record(record)
            gr_record.plot(
                ax,
                x_lim=(0, self.max_sequence_length),
                with_ruler=(ax == axes[-1]),
            )
            ax.set_title(seqname, loc="left", fontdict=dict(weight="bold"))
        fig.tight_layout()
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
