#!/usr/bin/env python3
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from subprocess import run
from typing import Callable, Dict, List, Tuple, Union

import snaptools_defaults
from utils import identity

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

# 3-tuples:
#  [0] command-line argument to snaptools
#  [1] default path
#  [2] function used to adjust the path (both default path and user-provided)
#      (if no adjustment required, use 'identity')

CommandData = Tuple[
    str,
    Path,
    Callable[
        [Path],
        Path,
    ],
]

SNAPTOOLS_COMMAND_DEFAULTS: Dict[str, List[CommandData]] = {
    'align_paired_end': [
        ('--input-reference', snaptools_defaults.DEFAULT_ALIGNMENT_INDEX, find_base_index_path),
    ],
    'snap-pre': [
        ('--genome-size', snaptools_defaults.DEFAULT_SIZE_INDEX, identity),
    ],
}

SNAPTOOLS_ALIGN_COMMAND = [
    'snaptools',
    'align-paired-end',
    '--input-reference',
    '{base_index_path}',
]

def run_with_defaults(snaptools_command: str, other_args: List[Union[Path, str]]):
    args = other_args.copy()

    # For each command, check known options/arguments with default values.
    # If the option is present, transform with the provided function (e.g.
    # finding the "base" path of a genome index). If the option is not present,
    # append with the default value
    for option, default, func in SNAPTOOLS_COMMAND_DEFAULTS[snaptools_command]:
        if snaptools_command in args:
            existing_option_index = args.index(option)
            value_index = existing_option_index + 1
            args[value_index] = func(args[value_index])
        else:
            args.append(option)
            args.append(func(default))

    command = ['snaptools']
    command.extend(fspath(arg) for arg in args)

    print('Running:', ' '.join(command))
    run(command, check=True)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('snaptools_command')
    p.add_argument(
        'other_arg',
        nargs='*',
        help='Arguments to the chosen snaptools command',
    )
    args = p.parse_args()

    run_with_defaults(args.snaptools_command, args.other_arg)
