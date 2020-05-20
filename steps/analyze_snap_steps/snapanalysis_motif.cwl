#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_motif
label: snap analysis motif
cwlVersion: v1.0

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2019-11-15"
s:license: https://spdx.org/licenses/Apache-2.0

s:keywords: edam:topic_0091 , edam:topic_0622
s:programmingLanguage: Python

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/

$schemas:
  - https://schema.org/docs/schema_org_rdfa.html
  - http://edamontology.org/EDAM_1.18.owl


dct:creator:
  '@id':  https://orcid.org/0000-0001-5173-4627
  foaf:name: Walter Shands
  foaf:mbox: jshands@ucsc.edu

requirements:
  DockerRequirement:
    dockerPull: "hubmap/sc-atac-seq-grch38"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000
  # Set the environment variable for the location of
  # the R environment file. In this file will be set
  # the location of the R temporary working directory
  # to a directory on the host and not in
  # the Docker container. If this is not done R
  # can run out of memory due to memory limitations
  # in the container
  # https://stackoverflow.com/questions/27302983/how-to-change-directory-for-temporary-files-problems-with-huge-temporary-raste
  # https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html
  # https://www.commonwl.org/user_guide/12-env/index.html
  #EnvVarRequirement:
  #  envDef:
  #    #R_ENVIRON  : $(inputs.r_environ_path)
  #    TEMP  : $(inputs.tmpdir)
  #    TMP  : $(inputs.tmpdir)
  #    TMPDIR  : $(inputs.tmpdir)

inputs:
  snap_rds:
    type: File
    inputBinding:
      position: 1
      prefix: --snap_rds
    doc: The snap object RDS file.

  snap_file:
    type: File
    inputBinding:
      position: 2
      prefix: --snap_file
    doc: The snap file.

  #tmpdir:
  #  type: string?
  #  inputBinding:
  #    position: 7
  #    prefix: --tmpdir
  #  doc: Path for the temporary directory


outputs:
  CSV_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.csv"

  RData_file:
    type: File
    outputBinding:
      glob: "chromvar_data.RData"

baseCommand: [Rscript, /opt/snapMotifAnalysis.R]
