import os
from geneblocks import DiffBlocks, load_record

seq_1 = load_record(os.path.join("sequences", "sequence1.gb"))
seq_2 = load_record(os.path.join("sequences", "sequence2.gb"))

diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2)
ax1, ax2 = diff_blocks.plot(figure_width=8)
ax1.figure.savefig("diff_blocks.png", bbox_inches='tight')
