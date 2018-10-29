from copy import deepcopy
from difflib import SequenceMatcher

from .Location import Location

from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from .biotools import sequence_to_record


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
                label = "%s ➤ %s" % (self.s1_location.extract_sequence(),
                                     self.s2_location.extract_sequence())
            else:
                if s1_length == s2_length:
                    label = "%dn changed" % s1_length
                else:
                    label = "%sn ➤ %sn change" % (s1_length, s2_length)
        elif self.operation == "equal":
            label = "Equal"

        return self.s2_location.to_biopython_feature(
            feature_type="diff_" + self.operation, label=label)

    def __str__(self):
        return ("%s %s|%s" % (self.operation,
                              self.s1_location,
                              self.s2_location))

    def __repr__(self):
        return str(self)


class DiffRecordTranslator(BiopythonTranslator):

    ignored_features_types = ('diff_equal',)

    @staticmethod
    def compute_feature_color(f):
        return dict(
            diff_delete='red',
            diff_insert='green',
            diff_replace='orange'
        ).get(f.type, 'white')

    @staticmethod
    def compute_feature_box_linewidth(f):
        return 1 if f.type.startswith('diff_') else 0

    @staticmethod
    def compute_feature_fontdict(f):
        return {'fontsize': 12 if f.type.startswith('diff_') else 8}


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
    def from_sequences(s1, s2):
        matcher = SequenceMatcher(a=s1.upper(), b=s2.upper(), autojunk=False)
        blocks = [
            DiffBlock(operation,
                      Location(s1s, s1e, sequence=s1),
                      Location(s2s, s2e, sequence=s2))
            for operation, s1s, s1e, s2s, s2e in matcher.get_opcodes()
        ]
        return DiffBlocks(s1, s2, blocks)

    def diffs_as_features(self):
        return [block.to_feature() for block in self.blocks]

    def plot(self, translator_class='default', separate_axes=True,  **plot_kw):

        if translator_class == 'default':
            translator_class = DiffRecordTranslator
        translator = translator_class()

        record = deepcopy(self.s2)
        if not hasattr(record, 'features'):
            record = sequence_to_record(record)

        diff_features = self.diffs_as_features()

        if separate_axes and len(record.features):
            gr_record = translator.translate_record(record)
            record.features = diff_features
            gr_diffrecord = translator.translate_record(record)
            width = plot_kw.get('figure_width', 8)
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(width, 6))
            _, stats1 = gr_record.plot(ax=ax1, **plot_kw)
            _, stats2 = gr_diffrecord.plot(ax=ax2, with_ruler=False, **plot_kw)
            max_level_1 = max([v['annotation_y'] for v in stats1[1].values()])
            max_features_1 = max([v['feature_y'] for v in stats1[1].values()])
            max_level_2 = max([v['annotation_y'] for v in stats2[1].values()])
            n_levels = max_level_1 + max_level_2
            fig.set_size_inches((width, 0.6*n_levels))
            ax2.set_ylim(ymin=-0.5)
            ax2.invert_yaxis()
            for f in gr_diffrecord.features:
                ax1.fill_between([f.start, f.end], y1=max_features_1 + 1,
                                  y2=-1,
                                 facecolor=f.color, alpha=0.07, zorder=1000)
            return (ax1, ax2)

        else:
            record.features += diff_features
            gr_record = translator.translate_record(record)
            ax, _ = gr_record.plot(**plot_kw)
            return ax
