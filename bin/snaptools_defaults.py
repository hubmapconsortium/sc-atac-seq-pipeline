from pathlib import Path

SUPPLEMENTARY_DATA_PATH = Path('/opt/supplementary-data')

# Some of these are only used by SnapATAC.
# TODO: refactor so the same data is usable both in Python and R

DEFAULT_ALIGNMENT_INDEX = SUPPLEMENTARY_DATA_PATH / 'bwa-index'
DEFAULT_GENE_ANNOTATION = SUPPLEMENTARY_DATA_PATH / 'gencode.v32.annotation.bed'
DEFAULT_SIZE_INDEX = SUPPLEMENTARY_DATA_PATH / 'grch38.fasta.fai'
DEFAULT_ENCODE_BLACKLIST = SUPPLEMENTARY_DATA_PATH / 'hg38.blacklist.bed'
DEFAULT_PROMOTERS = SUPPLEMENTARY_DATA_PATH / 'hg38.promoters.bed'
