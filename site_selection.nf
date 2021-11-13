#!/usr/bin/env nextflow

// Copyright (C) 2021 Sur Herrera Paredes

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.midas_dir = ''
params.meta = 'map.txt'
params.outdir = 'output'

midas_dir = file(params.midas_dir)

Channel.fromPath("$params.midas_dir/*", type: 'dir')
  .map{ dirname -> tuple(dirname.name, file(dirname)) }
  .into{MIDAS1; MIDAS2}

process FIT {
  label 'r'
  tag "$spec"
  publishDir "$params.outdir/FIT", mode: 'rellink'

  input:
  tuple spec, file(midas_dir) from MIDAS1
  file meta from file(params.meta)

  output:
  tuple spec, file("${spec}.tsv") optional true

  """
  Rscript $workflow.projectDir/fit_test.r \
    $midas_dir \
    $meta \
    --output ${spec}.tsv
  """

}

process S_COEF {
  label 'r'
  tag "$spec"
  publishDir "$params.outdir/s_coef", mode: 'rellink'

  input:
  tuple spec, file(midas_dir) from MIDAS2
  file meta from file(params.meta)

  output:
  tuple spec, file("${spec}.tsv") optional true

  """
  Rscript $workflow.projectDir/s_coef.r \
    $midas_dir \
    $meta \
    --output ${spec}.tsv
  """

}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '4h'
  memory = '2G'
  withLabel: 'r'{
    module = 'R/4.1.0'
    // module = "R/4.0.2:v8/8.4.371.22" // Make sure you have ~/.R/Makevars with CXX14
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
