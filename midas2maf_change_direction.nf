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
params.map = ''
params.outdir = 'output'
params.prop_thres = 0.8
params.recode = "no"

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
    pattern: "site_data.tsv.gz", saveAs: {"${spec}.tsv.gz"}

  input:
  tuple spec, file(midas_dir) from MIDAS
  file mapfile from mapfile
  val recode from params.recode

  output:
  tuple spec, file("site_data.tsv.gz") optional true into DAT

  """
  Rscript ${workflow.projectDir}/midasmerge2sitefreqs.r \
    $midas_dir \
    $mapfile \
    --recode $recode \
    --output site_data.tsv.gz
  """
}

process dists{
  tag "$spec"
  label 'r'
  publishDir "$params.outdir/pts_dist", mode: 'rellink',
    pattern: "output/pts_dist.tsv.gz", saveAs: {"${spec}.tsv.gz"}
  publishDir "$params.outdir/sites_dist", mode: 'rellink',
    pattern: "output/sites_dist.tsv.gz", saveAs: {"${spec}.tsv.gz"}
  publishDir "$params.outdir/figs/depthchange_vs_binompval/", mode: 'rellink',
    pattern: "output/output/depthchange_vs_binompval.jpeg",
    saveAs: {"${spec}.jpeg"}
  publishDir "$params.outdir/figs/hexbin_binomial_inc_vs_dec/", mode: 'rellink',
    pattern: "output/output/hexbin_binomial_inc_vs_dec.jpeg",
    saveAs: {"${spec}.jpeg"}
  publishDir "$params.outdir/figs/hexbin_days_vs_change/", mode: 'rellink',
    pattern: "output/output/hexbin_days_vs_change.jpeg",
    saveAs: {"${spec}.jpeg"}

  input:
  tuple spec, file(data) from DAT
  val prop_thres from params.prop_thres

  output:
  file "output/depthchange_vs_binompval.jpeg"
  file "output/hexbin_binomial_inc_vs_dec.jpeg"
  file "output/hexbin_days_vs_change.jpeg"
  file "output/pts_dist.tsv.gz"
  file "output/sites_dist.tsv.gz"

  """
  Rscript $workflow.projectDir/site_and_pt_maf_direction.r \
    $data \
    --outdir output \
    --prop_thres $prop_thres
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
