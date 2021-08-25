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
params.outdir = "output/"
params.cpus = 8

indir = file(params.indir)

// Doesn't like dots in pattern
INPUTS = Channel.fromFilePairs("$indir/*.{fasta,partitions}", flat: true)
// INPUTS.view()


process best_scheme {
  label 'iqtree2'
  label 'big'
  tag "$spec"
  publishDir "$params.outdir/find_best", mode: 'rellink',
    pattern: "$spec"
  cpus params.cpus

  input:
  tuple spec, file(aln), file(part) from INPUTS
  val cpus from params.cpus

  output:
  tuple spec, file(aln), file("${spec}/find_best.best_scheme.nex") into BEST
  file "$spec"

  """
  # Convert partitions file to nexus format
  echo -e "#nexus\nbegin sets;" > ${spec}.nex
  cat $part | awk '{print " charset " \$2 "dna = " \$4 ";"}' >> ${spec}.nex
  echo "end;" >> ${spec}.nex

  # Prepare output directory
  mkdir $spec

  # Run IQ-TREE to find best partition scheme
  iqtree2 \
    -s $aln \
    -p ${spec}.nex \
    -m MF+MERGE \
    -T $cpus \
    --prefix $spec/find_best
  """
}

process nonrev {
  label 'iqtree2'
  tag "$spec"
  publishDir "$params.outdir/nonrev", mode: 'rellink'
  cpus params.cpus

  input:
  tuple spec, file(aln), file(part) from BEST
  val cpus from params.cpus

  output:
  file "$spec"

  """
  # Create directory for output
  mkdir $spec

  # Run
  iqtree2 \
    -s $aln \
    -p $part \
    --model-joint 12.12 \
    -B 1000 \
    -T $cpus \
    --prefix $spec/nonrev_dna
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
  withLabel: 'iqtree2'{
    time = '36h'
    module = 'iq-tree/2'
  }
  withLabel: 'big'{
    memory = '72G'
  }

}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
