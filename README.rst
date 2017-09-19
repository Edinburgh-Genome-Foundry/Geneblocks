GeneBlocks
=============

GeneBlocks is a Python library to find common blocks in several DNA fragments.

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/illustration.jpeg
   :alt: [illustration]
   :align: center
   :width: 500px

At the Edinburgh Genome Foundry, we use GeneBlocks to:

- Gain insight on the structure non-annotated sequences that we get from different sources.
- Speed up sequence screening by avoiding to BLAST several times the same sequence
- Plan the assembly of several sequences at once.

Usage
------

.. code:: python

    from geneblocks import BlocksFinder

    # Input sequences are in a dictionnary as follows:
    sequences = {'seq1': 'ATTTGCGT...', 'seq2': 'ATGCCCGCACG...', etc.}

    blocks_finder = BlocksFinder(sequences)
    axes = blocks_finder.plot_common_blocks()
    axes[0].figure.savefig("basic_example.png", bbox_inches="tight")

Installation
-------------

Geneblocks requires NCBI BLAST+. On Ubuntu, install it with
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
