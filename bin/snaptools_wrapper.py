#!/usr/bin/env python3
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from subprocess import run
from typing import Callable, Dict, List, Tuple, Union

import snaptools_defaults
from utils import find_base_index_path, identity

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
    "mem": [
        ("--alignment-index", snaptools_defaults.DEFAULT_ALIGNMENT_INDEX, find_base_index_path),
    ]
    #    "align-paired-end": [
    #        ("--input-reference", snaptools_defaults.DEFAULT_ALIGNMENT_INDEX, find_base_index_path),
    #    ],
    #    "snap-pre": [
    #        ("--genome-size", snaptools_defaults.DEFAULT_SIZE_INDEX, identity),
    #    ],
}


def run_with_defaults(snaptools_command: str, other_args: List[Union[Path, str]]):
    args = other_args.copy()

    # For each command, check known options/arguments with default values.
    # If the option is present, transform with the provided function (e.g.
    # finding the "base" path of a genome index). If the option is not present,
    # append with the default value
    for option, default, func in SNAPTOOLS_COMMAND_DEFAULTS[snaptools_command]:
        print("option {} default {} func {}".format(option, default, func))
        print("args {}".format(args))
        if option in args:
            existing_option_index = args.index(option)
            print("existing option index: {}".format(existing_option_index))
            value_index = existing_option_index + 1
            # args[value_index] = func(args[value_index])
            # insert the reference genome index path just after 'bwa mem'
            # and don't use the option  --input-reference
            existing_option_value = args[value_index]
            args.remove(option)
            args.remove(existing_option_value)
            print("args before insertion: {}".format(args))
            args.insert(0, func(existing_option_value))
        else:
            # insert the reference genome index path just after 'bwa mem'
            args.insert(0, (func(default)))

    command = ["bwa", snaptools_command]
    # command = ["snaptools", snaptools_command]
    command.extend(fspath(arg) for arg in args)

    print("Running:", " ".join(command))
    run(command, check=True)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("snaptools_command")
    p.add_argument(
        "other_arg",
        nargs="*",
        help="Arguments to the chosen snaptools command",
    )
    args = p.parse_args()
    print("args {}".format(args))
    print("snaptools command: {}".format(args.snaptools_command))
    print("other args: {}".format(args.other_arg))

    run_with_defaults(args.snaptools_command, args.other_arg)
