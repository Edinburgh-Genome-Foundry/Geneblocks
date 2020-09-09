"""Algorithmic methods for the selection of common blocks in DiffBlocks

- select_common_blocks
  -
  - segments_difference
"""

import re
import tempfile
import subprocess
from collections import defaultdict, OrderedDict

import numpy as np
from ..biotools import reverse_complement, sequence_to_record


def format_sequences_as_dicts(sequences):
    """Standardize different formats into a single one.


    The ``sequences`` can be either:
    - A list [('sequence_id', 'ATGC...'), ('sequence_2', ...)]
    - A list of Biopython records (all with different IDs)
    - A dict {'sequence_id': "ATGC..."}
    - A dict {'sequence_id': biopython_record}

    The output is a tuple (sequences_dict, records_dict), where
    - sequences_dict is of the form {'sequence_id': 'ATGC...'}
    - sequences_dict is of the form {'sequence_id': 'ATGC...'}
    """
    if isinstance(sequences, (list, tuple)):
        if hasattr(sequences[0], "seq"):
            # SEQUENCES = LIST OF RECORDS
            records_dict = OrderedDict([(record.id, record) for record in sequences])
            sequences_dict = OrderedDict(
                [(record.id, str(record.seq).upper()) for record in sequences]
            )
        else:
            # SEQUENCES = LIST OF ATGC STRINGS
            sequences_dict = OrderedDict(sequences)
            if isinstance(list(sequences_dict.values())[0], str):
                records_dict = OrderedDict(
                    [
                        (name, sequence_to_record(seq, name=name))
                        for name, seq in sequences_dict.items()
                    ]
                )
            else:
                records_dict = sequences
    elif hasattr(list(sequences.values())[0], "seq"):
        # SEQUENCES = DICT {SEQ_ID: RECORD}
        records_dict = OrderedDict(sorted(sequences.items()))
        sequences_dict = OrderedDict(
            [
                (record_id, str(record.seq).upper())
                for record_id, record in sequences.items()
            ]
        )
    else:
        # SEQUENCES = DICT {SEQ_ID: ATGC}
        sequences_dict = OrderedDict(sorted(sequences.items()))
        records_dict = OrderedDict(
            [
                (name, sequence_to_record(seq, name=name))
                for name, seq in sequences.items()
            ]
        )
    return sequences_dict, records_dict


def segments_difference(segment, subtracted):
    """Return the difference between segment (start, end) and subtracted.

    The result is a list containing either zero, one, or two segments of the
    form (start, end).

    Examples
    --------

    >>> segment=(10, 100), subtracted=(0, 85) => [(85, 100)]
    >>> segment=(10, 100), subtracted=(40, 125) => [(10, 40)]
    >>> segment=(10, 100), subtracted=(30, 55) => [(10, 30), (55, 100)]
    >>> segment=(10, 100), subtracted=(0, 150) => []
    """
    seg_start, seg_end = segment
    sub_start, sub_end = subtracted
    result = []
    if sub_start > seg_start:
        result.append((seg_start, min(sub_start, seg_end)))
    if sub_end < seg_end:
        result.append((max(seg_start, sub_end), seg_end))
    return sorted(list(set(result)))


def find_homologies_between_sequences(
    sequences, min_size=0, max_size=None, include_self_homologies=True
):
    """Return a dict listing the locations of all homologies between sequences.

    The result is a dict of the form below, where the sequence identifiers
    are used as keys.

    >>> {
    >>>     'seq_1': {
    >>>        (start1, end1): [('seq2_5', _start, _end), ('seq_3', )...]
    >>>        (start2, end2): ...
    >>>     }
    >>> }

    Parameters
    ----------

    sequences
      A dict {'sequence_id': 'ATTGTGCAG...'}.

    min_size, max_size
      Minimum and maximum size outside which homologies will be ignored.

    include_self_homologies
      If False, self-homologies will be removed from the list.
    """

    # BLAST THE SEQUENCES USING NCBI-BLAST

    temp_fasta_path = tempfile.mktemp(".fa")

    with open(temp_fasta_path, "w+") as f:
        f.write(
            "\n\n".join(
                ["> %s\n%s" % (name, seq.upper()) for name, seq in sequences.items()]
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
            "1000000000000000",
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

    # PARSE THE RESULT FROM BLAST

    parsing = [line.split("\t") for line in result.decode("utf-8").splitlines()]
    homologies = {name: defaultdict(lambda *a: []) for name, seq in sequences.items()}

    # FILTER THE RESULTS (MIN_SIZE, MAX_SIZE, SELF-HOMOLOGIES)

    for query, qstart, qend, subject, sstart, send in parsing:
        is_self_homology = (query == subject) and (qstart != sstart)
        if is_self_homology and (not include_self_homologies):
            continue
        qstart, qend = int(qstart) - 1, int(qend)
        sstart, send = int(sstart) - 1, int(send)
        if qend - qstart < min_size:
            continue
        if (max_size is not None) and (qend - qstart > max_size):
            continue
        location = (subject, sstart, send)
        homologies[query][(qstart, qend)].append(location)
    return homologies


def count_homologies(matches, min_size):
    """Return a dict {(start, end): number_of_homologies_count}.
    """
    homologies_counts = {}
    if len(matches) == 1:
        segment = list(matches.keys())[0]
        homologies_counts[segment] = 1
    matches_list = sorted(matches.keys())
    for i, match1 in enumerate(matches_list):
        for match2 in matches_list[i + 1 :]:
            segment = start, end = (match2[0], min(match1[1], match2[1]))
            if end < start:
                # The segment is empty, match1 and match2 as disjunct.
                break
            elif (end - start > min_size) and (segment not in homologies_counts):
                homologies_counts[segment] = len(
                    [
                        matching
                        for (match_start, match_end) in matches_list
                        for matching in matches[(match_start, match_end)]
                        if match_start <= start <= end <= match_end
                    ]
                )
    return homologies_counts


def segment_with_most_homologies(homologies_counts, method="most_coverage_first"):
    """Select the "best" segment, that should be selected next as a common
    block."""

    def segment_score(segment):
        if method == "most_coverage_first":
            factor = homologies_counts[segment]
        else:
            factor = 1
        start, end = segment
        return factor * (end - start)

    return max(
        [(0, (None, None))]
        + [(segment_score(segment), segment) for segment in homologies_counts]
    )


def select_common_blocks(
    homologies, sequences, min_size=0, method="most_coverage_first"
):
    """Select a collection of the largest common blocks, iteratively."""
    common_blocks = []
    homologies_counts = {
        seqname: count_homologies(matches=homologies[seqname], min_size=min_size)
        for seqname in sequences
    }

    # ITERATIVELY SELECT A COMMON BLOCK AND REMOVE THAT BLOCK FROM THE
    # homologies IN VARIOUS SEQUENCES, UNTIL THERE IS NO HOMOLOGY

    while True:

        # FIND THE HOMOLOGY WITH THE BEST OVERALL SCORE ACROSS ALL SEQS

        (best_score, (start, end)), seqname = max(
            [
                (
                    segment_with_most_homologies(
                        homologies_counts[seqname], method=method
                    ),
                    seqname,
                )
                for seqname in sequences
            ]
        )

        # IF NO HOMOLOGY WAS FOUND AT ALL, STOP

        if best_score == 0:
            break

        # FIND WHERE THE SELECTED SUBSEQUENCE APPEARS IN OTHER SEQUENCES.
        # AT EACH LOCATION, "EXTRUDE" THE SUBSEQUENCE FROM THE CURRENT
        # LOCATIONS IN homologies_counts

        best_subsequence = sequences[seqname][start:end]
        locations = []
        for seqname, sequence in sequences.items():
            seq_n_intersections = homologies_counts[seqname]

            # we look for both the subsequence and its reverse complement:
            for strand in [1, -1]:
                if strand == 1:
                    matches = re.finditer(best_subsequence, sequence)
                else:
                    matches = re.finditer(
                        reverse_complement(best_subsequence), sequence
                    )
                for match in matches:

                    # add the location to the list for this subsequence...
                    start, end = match.start(), match.end()
                    locations.append((seqname, (start, end, strand)))

                    # ...then subtract the location from the sequence's
                    # homologies list
                    match_as_segment = tuple(sorted([start, end]))
                    for intersection in list(seq_n_intersections.keys()):
                        score = seq_n_intersections.pop(intersection)
                        for diff in segments_difference(intersection, match_as_segment):
                            diff_start, diff_end = diff
                            if diff_end - diff_start > min_size:
                                seq_n_intersections[diff] = score
        common_blocks.append((best_subsequence, locations))

    # REMOVE SELF-HOMOLOGOUS SEQUENCES

    common_blocks = [
        (seq, locations) for (seq, locations) in common_blocks if len(locations) >= 2
    ]

    # CREATE THE FINAL COMMON_BLOCKS_DICT

    common_blocks_dict = OrderedDict()
    if len(common_blocks) > 0:
        number_size = int(np.log10(len(common_blocks))) + 1
        for i, (sequence, locations) in enumerate(common_blocks):
            block_name = "block_%s" % (str(i + 1).zfill(number_size))
            common_blocks_dict[block_name] = {
                "sequence": sequence,
                "locations": locations,
            }

    return common_blocks_dict
