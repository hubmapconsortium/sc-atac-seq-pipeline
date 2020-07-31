import bz2
from dataclasses import dataclass
from enum import Enum
import gzip
import lzma
from os import PathLike
from pathlib import Path
from typing import Iterable, Tuple, Union

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

class FileType(Enum):
    def __new__(cls, filetype, open_function):
        obj = object.__new__(cls)
        obj._value_ = filetype
        obj.open_function = open_function
        return obj

    GZ = ("gz", gzip.open)
    BZ2 = ("bz2", bz2.open)
    XZ = ("xz", lzma.open)
    TEXT = ("txt", open)

def get_file_type_by_extension(file_path: Path) -> FileType:
    suffix = file_path.suffix.lstrip(".")
    try:
        return FileType(suffix)
    except ValueError:
        # No special suffix, assume text
        return FileType.TEXT

def smart_open(file_path: PathLike, mode="rt", *args, **kwargs):
    file_type = get_file_type_by_extension(Path(file_path))
    return file_type.open_function(file_path, mode, *args, **kwargs)

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
    'fq',
    'fq.gz',
]

@dataclass
class Read:
    read_id: str
    seq: str
    unused: str
    qual: str

    def serialize(self):
        unused = self.unused or '+'
        return '\n'.join([self.read_id, self.seq, unused, self.qual])

def fastq_reader(fastq_file: Path) -> Iterable[Read]:
    with smart_open(fastq_file) as f:
        while True:
            lines = [f.readline().strip() for _ in range(4)]
            if not all(lines):
                return
            yield Read(*lines)
