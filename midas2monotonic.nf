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

// Parameters
params.midas_dir = ''
params.map = ''
params.outdir = 'output'
params.brms_cpus = 4
params.min_sites = 5
params.iter = 3000
params.warmup = 2000

// Process params
midas_dir = file(params.midas_dir)
mapfile = file(params.map)

MIDAS = Channel.fromPath("$midas_dir/*", type: 'dir', maxDepth: 1)
  .map{ midas_dir -> tuple(midas_dir.name,
    file(midas_dir)) }

process preprocess{
  tag "$spec"
  label 'r'
  publishDir "$params.outdir/site_data", mode: 'rellink',
    pattern: "site_data.tsv.gz", saveAs: {"${spec}.tsv"}

  input:
  tuple spec, file(midas_dir) from MIDAS
  file mapfile from mapfile

  output:
  tuple spec, file("site_data.tsv.gz") optional true into DAT

  """
  Rscript ${workflow.projectDir}/midasmerge2sitefreqs.r \
    $midas_dir \
    $mapfile \
    --output site_data.tsv.gz
  """
}

process model{
  tag "$spec"
  label 'r'
  label 'long'
  publishDir "$params.outdir/model", mode: 'rellink',
    pattern: "output/", saveAs: {"$spec"}
  cpus params.brms_cpus
  // afterScript 'rm -r Rtmp*'

  input:
  tuple spec, file(site_data) from DAT
  val cpus from params.brms_cpus
  val min_sites from params.min_sites
  val iter from params.iter
  val warmup from params.warmup

  output:
  tuple spec, file("output") optional true

  """
  TMPDIR=${workflow.workDir}
  Rscript ${workflow.projectDir}/model_monotonic_gene_maf.r \
    $site_data \
    --outdir output \
    --min_sites $min_sites \
    --iter $iter \
    --warmup $warmup \
    --cores $cpus
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
  withLabel: 'long'{
    time = '300h'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
