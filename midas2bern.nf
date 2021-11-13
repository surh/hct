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

// Params
params.midas_dir = ''
params.map = '' // TO DO: We need one map per species in the end
params.outdir = 'output'
params.depth_thres = 5
params.npop_thres = 3
params.max_sites = 0

// Process params
midas_dir = file(params.midas_dir)
mapfile = file(params.map)

MIDAS = Channel.fromPath("$midas_dir/*", type: 'dir', maxDepth: 1)
  .map{ midas_dir -> tuple(midas_dir.name,
    file(midas_dir)) }

process midas2bern {
  tag "$spec"
  label 'r'
  publishDir "$params.outdir/sites", mode: 'rellink',
    pattern: "output/sites.tsv", saveAs: {"${spec}.tsv"}
  publishDir "$params.outdir/pops", mode: 'rellink',
    pattern: "output/pops.tsv", saveAs: {"${spec}.tsv"}

  input:
  tuple spec, file(midas_dir) from MIDAS
  file mapfile from mapfile
  val depth_thres from params.depth_thres
  val npop_thres from params.npop_thres
  val max_sites from params.max_sites

  output:
  file "output/sites.tsv"
  file "output/pops.tsv"

  """
  Rscript ${workflow.projectDir}/midas2bern.r \
    $midas_dir \
    $mapfile \
    --depth_thres $depth_thres \
    --n_thres $npop_thres \
    --outdir output \
    --max_sites $max_sites
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '24h'
  memory = '10G'
  withLabel: 'r'{
    module = 'R/3.6.1'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
