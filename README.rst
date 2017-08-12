GeneBlocks
=============
.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/GeneBlocks.svg?branch=master
  :target: https://travis-ci.org/Edinburgh-Genome-Foundry/GeneBlocks
  :alt: Travis CI build status

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GeneBlocks/master/docs/_static/images/title.png
   :alt: [logo]
   :align: center
   :width: 500px



GeneBlocks (complete documentation `here <https://edinburgh-genome-foundry.github.io/GeneBlocks/>`_)
is a Python library to find common blocks in several DNA fragments.

At the Edinburgh Genome Foundry, we use GeneBlocks to:

- Gain insight on the structure non-annotated sequences that we get from different sources.
- Speed up sequence screening by avoiding to BLAST several times the same sequence
- Plan the assembly of several sequences at once.

Usage
------


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
<http://www.genomefoundry.io>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaCauldron>`_ under the MIT licence (copyright Edinburgh Genome Foundry).
Everyone is welcome to contribute !
