from enum import Enum
from pathlib import Path
from typing import Union

def normalize_whitespace(string: str) -> str:
    return ' '.join(string.split())

def identity(x):
    return x

def find_base_index_path(alignment_index: Union[Path, str]) -> Path:
    """
    When BWA indexes a reference genome, it produces five files:
      reference.fasta -> [
        reference.fasta.amb
        reference.fasta.ann
        reference.fasta.bwt
        reference.fasta.pac
        reference.fasta.sa
      ]

    When aligning reads with BWA, the "base" name of the reference
    (reference.fasta above) must be given as the index parameter --
    BWA appends the other file extensions automatically. The genome
    itself, in FASTA format, is not required for alignment, just
    BWA's index files.

    If given a directory, this function constructs the "base" path
    to the reference genome, suitable for passing as a command-line
    argument to BWA. The resulting `Path` may or may not exist on
    disk; we don't care and we don't check for this.

    If `alignment_index` does not refer to a directory, we assume
    it's already the base name of the index, and return unchanged.

    :param alignment_index: BWA alignment index, either the directory
      or the actual "base path" given to BWA
    :return: "Base" path to alignment reference. Note that this path
      doesn't actually (need to) exist on disk
    """
    # Could use any of the files expected to be part of the BWA index,
    # and could check for the presence of all of them, but may as well
    # just use one
    alignment_index = Path(alignment_index)
    if alignment_index.is_dir():
        bwt_files = list(alignment_index.glob('*.bwt'))
        assert len(bwt_files) == 1
        return alignment_index / bwt_files[0].stem
    else:
        # Don't care whether it actually exists -- assume it's usable
        # as input to BWA
        return alignment_index
