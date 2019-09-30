# Code organisation

This document walks you trough the Geneblocks code. Please request changes if anything is unclear.

- **CommonBlocks/** implements the CommonBlocks class. The code is separated in two files:
  - **commonblocks_tools.py** implements the (rather long) algorithms to extract common blocks from a set of sequences
  - **CommonBlocks.py** implements the CommonBlocks class, which stores a set of common blocks and provides methods to display them or export them as Genbank.

- **DiffBlocks/** implements the DiffBlocks class. The code is separated in several files:
  - **DiffBlock.py** implements a class to store a unit diff block (sequence location and edit operation), as well as methods to pretty-print the block or export it as a Genbank feature.
  - **diffblocks_tools.py** implements the (rather long) algorithms to identify different diffblocks using the Levenshtein algorithm, and merge various diffblocks together.
  - **DiffBlocks.py** implements the CommonBlocks class, which stores a set of common blocks and provides methods to display them or export them as Genbank.

- **biotools.py** implements methods related to sequence and Genbank record analysis and manipulation (reverse_complemtn, BLAST)

- **sequence_modification_utils.py** implements a collection opf sequence manipulation routines (to insert, delete, replace fragments), only to create examples and tests.

