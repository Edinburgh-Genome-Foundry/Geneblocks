from .biotools import reverse_complement
from Bio.SeqFeature import SeqFeature, FeatureLocation


class Location:
    def __init__(
        self, start, end, strand=None, sequence=None, sequence_id=None
    ):

        self.start = start
        self.end = end
        self.strand = strand
        self.sequence = sequence
        self.sequence_id = sequence_id

    def extract_sequence(self, sequence=None):
        """Return the subsequence read at the given location.

        If sequence is None, ``self.sequence`` is used.
        """
        if sequence is None:
            sequence = self.sequence
        if hasattr(sequence, "seq"):
            sequence = str(sequence.seq)
        result = sequence[self.start : self.end]
        if self.strand == -1:
            return reverse_complement(result)
        else:
            return result

    def __repr__(self):
        """Represent"""
        result = "%d-%d" % (self.start, self.end)
        if self.strand is not None:
            result += {1: "(+)", -1: "(-)", 0: ""}[self.strand]
        if self.sequence_id is not None:
            result = self.sequence_id + "|" + result
        return result

    def __len__(self):
        """Size of the location"""
        return abs(self.end - self.start)

    def to_tuple(self):
        return self.start, self.end, self.strand

    def to_biopython_location(self):
        """Return a Biopython FeatureLocation equivalent to the location."""
        start, end, strand = [
            None if e is None else int(e)
            for e in [self.start, self.end, self.strand]
        ]
        return FeatureLocation(start, end, strand)

    def to_biopython_feature(self, feature_type="misc_feature", **qualifiers):
        """Return a Biopython SeqFeature with same location and custom
        qualifiers."""
        return SeqFeature(
            self.to_biopython_location(),
            type=feature_type,
            qualifiers=qualifiers,
        )
