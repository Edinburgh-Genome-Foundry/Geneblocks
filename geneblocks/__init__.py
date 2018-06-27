""" geneblocks/__init__.py """

# __all__ = []

from .CommonBlocks import CommonBlocks
from .DiffBlocks import DiffBlocks, DiffRecordTranslator
from .biotools import load_record

from .version import __version__
