import subprocess
import itertools
from collections import defaultdict, OrderedDict
import tempfile
import numpy as np
import re
from copy import deepcopy

try:
    from dna_features_viewer import BiopythonTranslator
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    PLOTS_AVAILABLE = True
except:
    PLOTS_AVAILABLE = False

from .biotools import reverse_complement, sequence_to_record, annotate_record

class BlocksFinder:

    def __init__(self, sequences, factor="coverage", min_block_size=80):
        self.factor = factor
        self.min_block_size = min_block_size
        if isinstance(sequences, (list, tuple)):
            if hasattr(sequences[0], 'seq'):
                sequences = OrderedDict([
                    (record.id, str(record.seq).upper())
                    for record in sequences
                ])

            sequences = OrderedDict(sequences)
        if isinstance(sequences, dict):
            sequences = OrderedDict(sorted(sequences.items()))
        self.sequences = sequences
        self.max_sequence_length = max([
            len(seq)
            for seq in self.sequences.values()
        ])

        self.find_intermatches()
        self.compute_match_profiles()
        self.find_common_blocks()

    def find_intermatches(self):
        temp_fasta_path = tempfile.mktemp(".fa")

        with open(temp_fasta_path, "w+") as f:
            f.write("\n\n".join([
                "> %s\n%s" % (name, seq)
                for name, seq in self.sequences.items()
            ]))
        proc = subprocess.Popen([
            "blastn",
            "-query", temp_fasta_path,
            "-subject", temp_fasta_path,
            "-perc_identity", '100',
            "-ungapped",
            "-outfmt", '6 qseqid qstart qend sseqid sstart send'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        result, blast_err = proc.communicate()
        parsing = [
            line.split('\t')
            for line in result.decode('utf-8').splitlines()
        ]
        self.intermatches = {
            name: defaultdict(lambda *a: [])
            for name, seq in self.sequences.items()
        }
        for query, qstart, qend, subject, sstart, send in parsing:
            qstart, qend = int(qstart), int(qend)
            sstart, send = int(sstart), int(send)
            if qend - qstart < self.min_block_size:
                continue
            location = (subject, sstart, send)
            self.intermatches[query][(qstart, qend)].append(location)

    def compute_match_profiles(self):
        self.match_profiles = OrderedDict([
            (name, np.zeros(len(sequence)))
            for name, sequence in self.sequences.items()
        ])
        for seqname, matches in self.intermatches.items():
            for (start, end), matchlist in matches.items():
                self.match_profiles[seqname][start+1:end-1] += len(matchlist)
        self.max_frequency = max([
            profile.max()
            for profile in self.match_profiles.values()
        ])
        self.original_match_profiles = deepcopy(self.match_profiles)

    def plot_match_profiles(self, axes=None, color='#8eb7f9', original=True):
        match_profiles = (self.original_match_profiles if original else
                          self.match_profiles)
        if axes is None:
            fig, axes = plt.subplots(len(self.sequences), 1,
                                     figsize=(10, 2*len(self.sequences)),
                                     facecolor="white", sharex=True)
        for ax, (seq, match_profile) in zip(axes, match_profiles.items()):
            ax.fill_between(range(len(match_profile)), match_profile,
                            facecolor=color)
            ax.set_ylabel(seq)
        return axes


    def find_most_common_block_in_profile(self, profile):
        def block_score(block):
            start, end = block
            block_length = end - start
            if block_length < self.min_block_size:
                return 0
            frequency = profile[int((start+end)/2)]
            if frequency == 1:
                return 0
            factor = block_length if self.factor == "coverage" else 1
            return factor * frequency

        changes = np.diff(profile).nonzero()[0]
        return max([
            (block_score(block), block)
            for block in zip(changes, changes[1:])
        ])

    def find_common_blocks(self):
        self.common_blocks = []
        best_score = 1
        while True:
            global_best_block = max([
                (self.find_most_common_block_in_profile(profile), seqname)
                for seqname, profile in self.match_profiles.items()
            ])
            (best_score, (start, end)), seqname = global_best_block
            if best_score == 0:
                break
            # print ("global_best_block", global_best_block)
            best_sequence = self.sequences[seqname][start: end]
            best_sequence_rev = reverse_complement(best_sequence)
            locations = []
            for seqname, sequence in self.sequences.items():
                for bestseq, strand in [(best_sequence, 1),
                                        (best_sequence_rev, -1)]:
                    for match in re.finditer(bestseq, sequence):
                        start, end = match.start(), match.end()
                        self.match_profiles[seqname][start: end] = 0
                        locations.append((seqname, (start, end, strand)))
            self.common_blocks.append(locations)

    def common_blocks_to_records(self, colors="auto"):
        records = OrderedDict([
            (seqname, sequence_to_record(seq))
            for seqname, seq in self.sequences.items()
        ])
        if colors == "auto":
            colors = itertools.cycle([cm.Paired(0.21 * i % 1.0)
                                      for i in range(30)])
        for i, (blocks, color) in enumerate(zip(self.common_blocks, colors)):
            for (seqname, location) in blocks:
                annotate_record(records[seqname], location,
                                feature_type='misc_feature',
                                label="block_%d" % (i + 1), color=color)
        return records

    def plot_common_blocks(self, colors="auto", axes=None, figure_width=10,
                           ax_height=2):
        if not PLOTS_AVAILABLE:
            raise ImportError("Plotting requires Matplotlib and "
                              "DNA Features Viewer installed. See docs.")

        translator = BiopythonTranslator()
        records = self.common_blocks_to_records(colors=colors)
        if axes is None:
            fig, axes = plt.subplots(
                len(self.sequences), 1, facecolor="white", sharex=True,
                figsize=(figure_width, ax_height * len(self.sequences)),
            )
        for (ax, (seqname, record)) in zip(axes, records.items()):
            gr_record = translator.translate_record(record)
            gr_record.plot(ax, x_lim=(0, self.max_sequence_length),
                           with_ruler=(ax == axes[-1]))
            ax.set_title(seqname, loc="left", fontdict=dict(weight="bold"))
        fig.tight_layout()
        return axes
