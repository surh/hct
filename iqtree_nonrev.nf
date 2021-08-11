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

params.indir = "input/"


indir = file(params.indir)

INPUTS = Channel.fromPath("$indir/MGY*{[.]partitions}")
// INPUTS = Channel.fromFilePairs("$indir/MGY*{.fasta, .partitions}")

INPUTS.view()

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '4h'
  memory = '2G'
  withLabel: 'iqtree2'{
    time = '24h'
    cpus = 8
    module = 'iq-tree/2'
  }

}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
