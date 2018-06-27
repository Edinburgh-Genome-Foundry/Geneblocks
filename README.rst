GeneBlocks
=============
.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/Geneblocks
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/Geneblocks/badge.svg
  :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/Geneblocks

GeneBlocks is a Python library to compare different DNA sequences. It can either find common blocks in a group of DNA sequences, or highlight the differences between two sequences.

At the Edinburgh Genome Foundry, we use GeneBlocks to optimize sequence assembly, explore sets of non-annotated sequences, or visualize the differences
between different versions of a sequence.

Usage
------


Finding common blocks in a set of sequences:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from geneblocks import CommonBlocks

    # Input sequences are in a dictionnary as follows:
    sequences = {'seq1': 'ATTTGCGT...', 'seq2': 'ATGCCCGCACG...', ...}

    common_blocks = CommonBlocks(sequences)

    # PLOT THE BLOCKS
    axes = common_blocks.plot_common_blocks()
    axes[0].figure.savefig("basic_example.png", bbox_inches="tight")

    # GET ALL COMMON BLOCKS AS BIOPYTHON RECORDS
    blocks_records = common_blocks.common_blocks_records()

    # WRITE ALL COMMON BLOCKS INTO A CSV SPREADSHEET
    common_blocks.common_blocks_to_csv(self, target_file=None)

Result:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/examples/common_blocks.png
   :alt: [illustration]
   :align: center
   :width: 700px


Highlighting the differences between two sequences:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    seq_1 = load_record("sequence1.gb")
    seq_2 = load_record("sequence2.gb")

    diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2)
    ax1, ax2 = diff_blocks.plot(figure_width=8)
    ax1.figure.savefig("diff_blocks.png")

Result:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/examples/diff_blocks.png
   :alt: [illustration]
   :align: center
   :width: 700px


Installation
-------------

The CommonBlocks feature requires NCBI BLAST+. On Ubuntu, install it with

.. code:: shell
    (sudo) apt-get install ncbi-blast+


You can install GeneBlocks through PIP

.. code:: shell

    sudo pip install geneblocks

Alternatively, you can unzip the sources in a folder and type

.. code:: shell

    sudo python setup.py install


Licence
--------

Geneblocks is an open-source software originally written at the `Edinburgh Genome Foundry
<http://www.genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/Geneblocks>`_ under the MIT licence (copyright Edinburgh Genome Foundry).
Everyone is welcome to contribute !
