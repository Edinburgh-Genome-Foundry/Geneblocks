import tempfile
import subprocess

import numpy as np

try:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False


complements_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}


def random_dna_sequence(length, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    proba
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


def load_record(filename, linear=True, name="id", upperize=True):
    formt = "genbank" if filename.endswith(("gb", "gbk")) else "fasta"
    record = SeqIO.read(filename, formt)
    if upperize:
        record.seq = record.seq.upper()
    record.linear = linear
    if name != "id":
        record.id = name
    record.name = record.id.replace(" ", "_")[:20]

    return record


def complement(sequence):
    return "".join(complements_dict[c] for c in sequence)


def reverse_complement(sequence):
    return complement(sequence)[::-1]


def sequence_to_record(sequence, record_id=None, name="unnamed", features=()):
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("Creating records requires Biopython installed.")
    if has_dna_alphabet:  # Biopython <1.78
        sequence = Seq(sequence, alphabet=DNAAlphabet())
    else:
        sequence = Seq(sequence)

    seqrecord = SeqRecord(sequence, name=name, id=record_id, features=list(features),)
    seqrecord.annotations["molecule_type"] = "DNA"

    return seqrecord


def annotate_record(
    seqrecord, location="full", feature_type="feature", margin=0, **qualifiers
):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """

    if not BIOPYTHON_AVAILABLE:
        raise ImportError("Creating records requires Biopython installed.")

    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type,
        )
    )


def sequences_differences_array(seq1, seq2):
    """Return an array [0, 0, 1, 0, ...] with 1s for sequence differences.

    seq1, seq2 should both be ATGC strings.
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            "Only use on same-size sequences (%d, %d)" % (len(seq1), len(seq2))
        )
    arr1 = np.fromstring(seq1, dtype="uint8")
    arr2 = np.fromstring(seq2, dtype="uint8")
    return arr1 != arr2


def sequences_differences(seq1, seq2):
    """Return the number of nucleotides that differ in the two sequences.

    seq1, seq2 should be strings of DNA sequences e.g. "ATGCTGTGC"
    """
    return sequences_differences_array(seq1, seq2).sum()
