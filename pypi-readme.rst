GeneBlocks
=============

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel.svg?branch=master
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

Infos
-----

**PIP installation:**

.. code:: bash

  pip install geneblocks

**Github Page**

`<https://github.com/Edinburgh-Genome-Foundry/geneblocks>`_

**Live demos**


`Find common sequence blocks <https://cuba.genomefoundry.org/find-common-blocks>`_

`Plot sequence diffs <https://cuba.genomefoundry.org/compare-two-sequences>`_

`Transfer Genbank features between records <https://cuba.genomefoundry.org/transfer-features>`_

**License:** MIT, Copyright Edinburgh Genome Foundry

More biology software
-----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

Geneblocks is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.