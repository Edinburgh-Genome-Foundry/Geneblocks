.. raw:: html

    <p align="center">
    <img alt="DNA Chisel Logo" title="DNA Chisel" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Geneblocks/master/docs/_static/images/title.png" width="450">
    <br /><br />
    </p>

GeneBlocks
===========

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/Geneblocks.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/Geneblocks
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/Geneblocks/badge.svg
  :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/Geneblocks

GeneBlocks is a Python library to compare DNA sequences. It can be used to:

- Find common blocks in a group of DNA sequences, to factorize them (e.g. only analyze or synthetize each common block once)
- Highlight differences between sequences (insertions, deletions, mutations)
- Transfer Genbank features from one record to another sharing similar subsequences.

At the Edinburgh Genome Foundry, we use GeneBlocks to optimize sequence assembly, explore sets of non-annotated sequences, or visualize the differences
between different versions of a sequence, and re-annotate records coming from third parties such as DNA manufacturers.

Live demos
----------

`Find common sequence blocks <https://cuba.genomefoundry.org/find-common-blocks>`_

`Plot sequence diffs <https://cuba.genomefoundry.org/compare-two-sequences>`_

`Transfer Genbank features between records <https://cuba.genomefoundry.org/transfer-features>`_

Installation
-------------

The CommonBlocks feature requires NCBI BLAST+. On Ubuntu, install it with

.. code:: shell

    (sudo) apt-get install ncbi-blast+


You can install GeneBlocks through PIP

.. code:: shell

    (sudo) pip install geneblocks

Alternatively, you can unzip the sources in a folder and type

.. code:: shell

    (sudo) python setup.py install


Usage
------


Finding common blocks in a set of sequences:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from geneblocks import CommonBlocks

    # Input sequences are in a dictionnary as follows:
    sequences = {'seq1': 'ATTTGCGT...', 'seq2': 'ATGCCCGCACG...', ...}

    common_blocks = CommonBlocks.from_sequences(sequences)

    # PLOT THE BLOCKS
    axes = common_blocks.plot_common_blocks()
    axes[0].figure.savefig("basic_example.png", bbox_inches="tight")

    # GET ALL COMMON BLOCKS AS BIOPYTHON RECORDS
    blocks_records = common_blocks.common_blocks_records()

    # WRITE ALL COMMON BLOCKS INTO A CSV SPREADSHEET
    common_blocks.common_blocks_to_csv(target_file="blocks.csv")

Result:

.. raw:: html
   
   <img src='https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/examples/common_blocks.png'
    width='600px'/>

Transfering features between genbank records:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this snippet we assume that we have two genbank records:

- A record of an annotated part, containing an expression module.
- A record of a plasmid which contains the part but the part was not properly annotated

We will use Geneblocks to automatically detect where the part is located in
the plasmid and automatically copy the features from the part record to the
plasmid record.

.. code:: python

    from geneblocks import CommonBlocks, load_record
    part = load_record('part.gb', name='insert')
    plasmid = load_record('part.gb', name='plasmid')
    blocks = CommonBlocks.from_sequences([part, plasmid])
    new_records = blocks.copy_features_between_common_blocks(inplace=False)
    annotated_plasmid = new_records['plasmid'] # Biopython record


The resulting annotated plasmids has annotations from both the original plasmid and the annotated part:

.. raw:: html
   
   <img src='https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/examples/features_transfer.png'
    width='600px'/>

Highlighting the differences between two sequences:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    seq_1 = load_record("sequence1.gb")
    seq_2 = load_record("sequence2.gb")

    diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2)
    ax1, ax2 = diff_blocks.merged().plot(figure_width=8)
    ax1.figure.savefig("diff_blocks.png")

Result:

.. raw:: html
   
     <img src='https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/examples/diff_blocks.png'
      width='700px'/>

Example on more complex sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the example below we build two ~50kb "sister" sequences with many insertions,
deletions, transpositions between them, and even a self-homology. Then we ask
Geneblocks to compare the two sequences, in term of common blocks, and in terms
of modifications brought to the second sequence:

.. code:: python

    from geneblocks import DiffBlocks, CommonBlocks, random_dna_sequence
    import geneblocks.sequence_modification_utils as smu
    import matplotlib.pyplot as plt

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
    fig, axes = plt.subplots(3, 1, figsize=(15, 8))
    common_blocks.plot_common_blocks(axes=axes[:-1])
    diff_blocks.plot(ax=axes[-1], separate_axes=False)
    axes[-1].set_xlabel("Changes in seq2 vs. seq1")
    fig.savefig("complex_sequences.png", bbox_inches='tight')

.. raw:: html
   
     <img src='https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/examples/complex_sequences.png'
      width='700px'/>

Licence
--------

Geneblocks is an open-source software originally written at the `Edinburgh Genome Foundry
<http://www.genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/Geneblocks>`_ under the MIT licence (copyright Edinburgh Genome Foundry).
Everyone is welcome to contribute !
