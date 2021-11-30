#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_create_cell_by_bin_matrix
label: snaptools create cell by bin matrix
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.4

  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.snap_file)
        writable: true

inputs:
  snap_file:
    type: File
    inputBinding:
      position: 1
      prefix: --snap-file
    doc: The SNAP file in which to place the cell by bin matrix.

  bin_size_list:
    type: int[]?
    inputBinding:
      prefix: --bin-size-list
      position: 2
    default: [50, 1000, 5000, 10000]
    doc: A list of bin size(s) to create in the cell-by-bin count matrix.

  tmp_folder:
    type: string?
    inputBinding:
      position: 11
      prefix: --tmp-folder
    default: "/tmp"
    doc: Directory to store temporary files.

  verbose:
    type: string?
    inputBinding:
      position: 13
      prefix: --verbose
    default: "TRUE"
    doc: A boolean tag; if true output the progress.

outputs:
  snap_file_w_cell_by_bin:
    type: File
    outputBinding:
      glob: $(inputs.snap_file.basename)

baseCommand: [snaptools, snap-add-bmat]
