from geneblocks import DiffBlocks, CommonBlocks, random_dna_sequence
import geneblocks.sequence_modification_utils as smu
import matplotlib.pyplot as plt
import numpy

numpy.random.seed(1)  # ensures the sequences will be the same at each run

# GENERATE 2 "SISTER" SEQUENCES FOR THE EXAMPLE
seq1 = random_dna_sequence(50000)
seq1 = smu.copy(seq1, 25000, 30000, 50000)

seq2 = seq1
seq2 = smu.insert(seq2, 39000, random_dna_sequence(100))
seq2 = smu.insert(seq2, 38000, random_dna_sequence(100))
seq2 = smu.reverse(seq2, 30000, 35000)
seq2 = smu.swap(seq2, (30000, 35000), (45000, 480000))
seq2 = smu.delete(seq2, 20000, 2000)
seq2 = smu.insert(seq2, 10000, random_dna_sequence(2000))
seq2 = smu.insert(seq2, 0, 1000*"A")

# FIND COMMON BLOCKS AND DIFFS
common_blocks = CommonBlocks.from_sequences({'seq1': seq1, 'seq2': seq2})
diff_blocks = DiffBlocks.from_sequences(seq1, seq2).merged()

# PLOT EVERYTHING
fig, axes = plt.subplots(3, 1, figsize=(16, 8))
common_blocks.plot_common_blocks(axes=axes[:-1])
diff_blocks.plot(ax=axes[-1], separate_axes=False)
axes[-1].set_xlabel("Changes in seq2 vs. seq1")
fig.savefig("complex_sequences.png", bbox_inches='tight')