from ..Location import Location
from ..biotools import sequences_differences


class DiffBlock:
    """Class to represent a segment that differs between sequences.

    Parameters
    ----------
    operation
      One of "insert", "delete", "replace", or "equal".

    s1_location
      The Location(start, end) of the region in sequence sequence s1.

    s2_location
      The Location(start, end) of the region in sequence sequence s2.
    """

    def __init__(self, operation, s1_location, s2_location):
        self.operation = operation
        self.s1_location = s1_location
        self.s2_location = s2_location

    def to_feature(self, sequence="s2"):
        s1_length, s2_length = len(self.s1_location), len(self.s2_location)
        max_length = max([s1_length, s2_length])

        if sequence == "s1":
            return DiffBlock(
                self.operation, self.s2_location, self.s1_location
            ).to_feature()

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
                sub_s1 = self.s1_location.extract_sequence()
                sub_s2 = self.s2_location.extract_sequence()
                diffs = sequences_differences(sub_s1, sub_s2)
                label = "%d mutations" % diffs

        elif self.operation == "change":
            if max([s1_length, s2_length]) <= 6:
                label = "%s ➤ %s" % (
                    self.s1_location.extract_sequence(),
                    self.s2_location.extract_sequence(),
                )
            else:
                label = "%sn ➤ %sn change" % (s1_length, s2_length)
        elif self.operation == "reverse":
            label = "was reversed at %d-%d" % (
                self.s1_location.start,
                self.s1_location.end,
            )
        elif self.operation == "transpose":
            label = "was at %d-%d" % (self.s1_location.start, self.s1_location.end,)
        elif self.operation == "equal":
            label = "Equal"

        return self.s2_location.to_biopython_feature(
            feature_type="diff_" + self.operation, label=label
        )

    def __str__(self):
        """Represent a diffbloc, for instance: "insert 10-30|105-105" """
        return "%s %s|%s" % (self.operation, self.s1_location, self.s2_location,)

    def __repr__(self):
        return str(self)
