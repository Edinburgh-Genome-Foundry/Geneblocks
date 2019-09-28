from dna_features_viewer import BiopythonTranslator

class DiffRecordTranslator(BiopythonTranslator):

    ignored_features_types = ("diff_equal",)

    @staticmethod
    def compute_feature_color(f):
        return dict(
            diff_delete="red",
            diff_insert="green",
            diff_replace="orange",
            diff_reverse="white",
            diff_transpose="white",
        ).get(f.type, "white")

    @staticmethod
    def compute_feature_box_linewidth(f):
        return 1 if f.type.startswith("diff_") else 0

    @staticmethod
    def compute_feature_fontdict(f):
        return {"fontsize": 12 if f.type.startswith("diff_") else 7}