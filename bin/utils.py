from enum import Enum


class Assay(Enum):
    def __new__(cls, key: str, fastq_count: int):
        obj = object.__new__(cls)
        obj._value_ = key
        obj.fastq_count = fastq_count
        return obj

    def __str__(self):
        return self.value

    SNARESEQ = "snareseq", 3
    SCISEQ = "sciseq", 2
    SNSEQ = "snseq", 3
    MULTIOME_10X = "multiome_10x", 3
