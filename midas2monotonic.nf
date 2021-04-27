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

// Paramters
params.midas_dir = ''
params.map = ''
params.outdir = 'output'
params.brms_cpus = 4

// Process params
midas_dir = file(params.midas_dir)
mapfile = file(params.map)

MIDAS = Channel.fromPath("$midas_dir/*", type: 'dir', maxdepth: 1)
  .map{ midas_dir -> tuple(midas_dir.name,
    file(midas_dir)) }

process preprocess{
  tag "$spec"
  label 'r'
  publishDir "$params.outdir/site_data", mode: 'rellink',
    pattern: "site_date.tsv", saveAs: {"${spec}.tsv"}

  input:
  tuple spec, file(midas_dir) from MIDAS
  file mapfile from mapfile

  output:
  tuple spec, file("site_data.tsv") into DAT

  """
  Rscript ${workflow.projectDir}/midasmerge2sitefreqs.r \
    $midas_dir \
    $mapfile \
    --output site_data.tsv
  """
}
