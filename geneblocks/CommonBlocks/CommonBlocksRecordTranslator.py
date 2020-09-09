from dna_features_viewer import BiopythonTranslator


class CommonBlocksRecordTranslator(BiopythonTranslator):

    ignored_features_types = ("diff_equal",)
    default_box_color = None

    def compute_feature_color(self, f):
        if f.qualifiers.get("is_block", False):
            return BiopythonTranslator.compute_feature_color(self, f)
        else:
            return "white"

    @staticmethod
    def compute_feature_box_linewidth(f):
        return 1 if f.qualifiers.get("is_block", False) else 0

    @staticmethod
    def compute_feature_fontdict(f):
        return {"fontsize": 12 if f.qualifiers.get("is_block", False) else 9}
