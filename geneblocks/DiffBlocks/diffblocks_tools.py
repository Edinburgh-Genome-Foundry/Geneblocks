from copy import deepcopy

import networkx as nx
import Levenshtein

from ..Location import Location

from .DiffBlock import DiffBlock


def compute_levenshtein_blocks(seq1, seq2, max_complexity=1e8):
    """Compute the Levenshtein blocks of insertion, deletion, replacement.
    """
    # TODO: better method for dealing with long sequences?
    l1, l2 = len(seq1), len(seq2)
    if l1 * l2 > max_complexity:
        return [("change", (0, l1), (0, l2))]

    def block_format(op, s1, e1, s2, e2):
        if op == "delete":
            return (op, (s1, e1 + 1), (s2, e2))
        if op == "insert":
            return (op, (s1, e1), (s2, e2 + 1))
        else:
            return (op, (s1, e1 + 1), (s2, e2 + 1))

    edits = Levenshtein.editops(seq1, seq2)
    if len(edits) == 0:
        return []
    bop, s1, s2 = edits[0]
    e1, e2 = s1, s2
    blocks = []
    for (op, _e1, _e2) in edits[1:]:
        continuity = any(
            [
                all([op == "delete", _e1 == e1 + 1, e2 == _e2]),
                all([op == "insert", _e1 == e1, _e2 == e2 + 1]),
                all([op == "replace", _e1 == e1 + 1, _e2 == e2 + 1]),
            ]
        )
        if op == bop and continuity:
            e1, e2 = _e1, _e2
        else:
            blocks.append(block_format(bop, s1, e1, s2, e2))
            bop, s1, s2 = op, _e1, _e2
            e1, e2 = s1, s2
    blocks.append(block_format(bop, s1, e1, s2, e2))
    return blocks


def merge_subblocks(subblocks):
    s1_location = Location(
        min([b.s1_location.start for b in subblocks]),
        max([b.s1_location.end for b in subblocks]),
        sequence=subblocks[0].s1_location.sequence,
    )
    s2_location = Location(
        min([b.s2_location.start for b in subblocks]),
        max([b.s2_location.end for b in subblocks]),
        sequence=subblocks[0].s2_location.sequence,
    )
    if len(s1_location) == len(s2_location):
        return DiffBlock("replace", s1_location, s2_location)
    else:
        return DiffBlock("change", s1_location, s2_location)


def merge_blocs_by_location(blocks, max_blocks, max_span, reference="s2"):
    while 1:
        for i in range(len(blocks) - max_blocks):
            subblocks = blocks[i : i + max_blocks]
            b1, b2 = subblocks[0], subblocks[-1]
            if reference == "s2":
                span = b2.s2_location.end - b1.s2_location.start
            else:
                span = b2.s1_location.end - b1.s1_location.start
            if span < max_span:
                new_block = merge_subblocks(subblocks)
                blocks = blocks[:i] + [new_block] + blocks[i + max_blocks :]
                break
        else:
            break
    return blocks


def merge_successive_blocks(blocks, change_gap=10, replace_gap=5, reference="s2"):
    while 1:
        for i in range(len(blocks) - 1):
            b1, b2 = blocks[i], blocks[i + 1]
            operations = (b1.operation, b2.operation)
            if reference == "s2":
                gap = b2.s2_location.start - b1.s2_location.end
            else:
                gap = b2.s1_location.start - b1.s1_location.end
            if "equal" in operations:
                continue
            if "change" in operations and (gap < change_gap):
                new_block = merge_subblocks([b1, b2])
                blocks = blocks[:i] + [new_block] + blocks[i + 2 :]
                break
            if operations == ("replace", "replace") and gap < replace_gap:
                new_block = merge_subblocks([b1, b2])
                blocks = blocks[:i] + [new_block] + blocks[i + 2 :]
                break
        else:
            break
    return blocks


def compute_sorted_blocks(blocks, reference="s2"):
    def sort_key(block):
        if reference == "s2":
            return block.s2_location.to_tuple()
        else:
            return block.s1_location.to_tuple()

    return sorted(blocks, key=sort_key)


def get_optimal_common_blocks(common_blocks):
    common_blocks = deepcopy(common_blocks)
    remarks = []

    # Make so that there is never an antisense block in s1 and a + block in s2.
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
            blocks_in_seqs_dicts[seq][block_name] = dict(rank=i, location=(start, end))
            seen_blocks.add(block_name)

    # Find and retain the largest sequence of blocks which is in the right
    # order in both sequences. We will remove every other block.
    # We do that by looking for the longest path in a graph

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
                blocks_to_reduce[block_name1] = max(blocks_to_reduce[block_name1], diff)

    for seq in ("s1", "s2"):
        blocks_list = blocks_in_seqs[seq]
        for i, (start, end, block_name) in enumerate(blocks_list):
            if block_name in blocks_to_reduce:
                diff = blocks_to_reduce[block_name]
                blocks_list[i] = (start, end - diff, block_name)

    return blocks_in_seqs, remarks
