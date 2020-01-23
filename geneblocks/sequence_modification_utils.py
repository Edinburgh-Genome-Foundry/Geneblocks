"""These methods are only useful to build examples and tests for Geneblocks."""

from .biotools import reverse_complement


def change(seq, start, end, change):
    """Return the sequence with ``seq[start:end]`` replaced by ``change``"""
    return seq[:start] + change + seq[end:]


def insert(seq, pos, inserted):
    """Return the sequence with ``inserted`` inserted, starting at index 'pos'
    """
    return seq[:pos] + inserted + seq[pos:]


def delete(seq, pos, deletions):
    """Return the sequence with a number of deletions from position pos."""
    return seq[:pos] + seq[pos + deletions :]


def reverse(seq, start, end):
    """Return the sequence with segment seq[start:end] reverse-complemented."""
    return seq[:start] + reverse_complement(seq[start:end]) + seq[end:]


def move(seq, start, end, diff):
    """Move a subsequence by "diff" nucleotides the left or the right."""
    sub = seq[start:end]
    if diff > 0:
        return seq[:start] + seq[end : end + diff] + sub + seq[end + diff :]
    else:
        return (
            seq[: start + diff] + sub + seq[start + diff : start] + seq[end:]
        )


def swap(seq, pos1, pos2):
    """Return a new sequence with segments at position pos1 and pos2 swapped.

    pos1, pos2 are both of the form (start1, end1), (start2, end2)
    """
    (start1, end1), (start2, end2) = sorted([pos1, pos2])
    return (
        seq[:start1]
        + seq[start2:end2]
        + seq[end1:start2]
        + seq[start1:end1]
        + seq[end2:]
    )


def copy(seq, start, end, new_start):
    """Return the sequence with segment [start, end] also copied elsewhere,
    starting in new_start."""
    return insert(seq, new_start, seq[start:end])
