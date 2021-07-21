#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_add_pmat_tool
label: snapanalysis add pmat
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.3

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
    doc: The SNAP file to be processed.

  peak_file:
    type: File
    inputBinding:
      position: 2
      prefix: --peak-file
    doc: The BED file that contains peak information.

  buffer_size:
    type: int?
    inputBinding:
      position: 3
      prefix: --buffer-size
    doc: Max number of barcodes be stored in memory (default 1000).

  tmp_folder:
    type: string?
    inputBinding:
      position: 4
      prefix: --tmp-folder
    default: "/tmp"
    doc: Directory to store temporary files.

  verbose:
    type: string?
    inputBinding:
      position: 5
      prefix: --verbose
    default: "TRUE"
    doc: A boolean tag; if true output the progress.

outputs:
  snap_file_w_peaks:
    type: File
    outputBinding:
      glob: $(inputs.snap_file.basename)

baseCommand: [snaptools, snap-add-pmat]
