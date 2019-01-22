from .biotools import reverse_complement
def change(seq, start, end, change):
    return seq[:start] + change + seq[end:]

def insert(seq, pos, inserted):
    return seq[:pos] + inserted + seq[pos:]

def delete(seq, pos, deletions):
    return seq[:pos] + seq[pos + deletions:]

def reverse(seq, start, end):
    return seq[:start] + reverse_complement(seq[start:end]) + seq[end:]

def move(seq, start, end, diff):
    sub = seq[start: end]
    if diff > 0:
        return seq[:start] + seq[end:end + diff] + sub + seq[end + diff:]
    else:
        return seq[:start + diff] + sub + seq[start + diff: start] + seq[end:]

def swap(seq, pos1, pos2):
    (start1, end1), (start2, end2) = sorted([pos1, pos2])
    return (seq[:start1] + seq[start2:end2] + seq[end1: start2] +
            seq[start1:end1] + seq[end2:])

def copy(seq, start, end, new_start):
    return insert(seq, new_start, seq[start:end])