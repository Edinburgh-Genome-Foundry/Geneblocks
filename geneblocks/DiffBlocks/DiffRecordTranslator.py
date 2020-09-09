from dna_features_viewer import BiopythonTranslator


class DiffRecordTranslator(BiopythonTranslator):

    ignored_features_types = ("diff_equal",)
    default_box_color = None

    @staticmethod
    def compute_feature_color(f):
        return dict(
            diff_delete="#E76F51",  # RED
            diff_insert="#2A9D8F",  # GREEN
            diff_replace="#E9C46A",  # YELLOW
            diff_change="#F4A261",  # ORANGE
            diff_reverse="white",
            diff_transpose="white",
        ).get(f.type, "white")

    @staticmethod
    def compute_feature_box_linewidth(f):
        return 1 if f.type.startswith("diff_") else 0

    @staticmethod
    def compute_feature_fontdict(f):
        return {"fontsize": 12 if f.type.startswith("diff_") else 9}
