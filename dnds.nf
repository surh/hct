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

// Command line params
params.gff = ''
params.enogg = ''
params.fna = ''
params.midas = ''
params.outdir = 'output'

// Process params
gff_dir = file(params.gff)
enogg_dir = file(params.enogg)
fna_dir = file(params.fna)
midas_dir = file(params.midas)

GFF = Channel.fromPath("$gff_dir/*.gff")
  .map{ infile -> tuple(infile.name.replaceAll(/\.gff$/, ''),
    file(infile))}
ENOGG = Channel.fromPath("$enogg_dir/*.tsv")
  .map{ infile -> tuple(infile.name.replaceAll(/\.tsv$/, ''),
    file(infile))}
FNA = Channel.fromPath("$fna_dir/*.fna")
  .map{ infile -> tuple(infile.name.replaceAll(/\.fna$/, ''),
    file(infile))}
MIDAS = Channel.fromPath("$enogg_dir/*", type: 'dir', maxdepth: 1)
  .map{ infile -> tuple(infile.name,
    file(infile))}

BUILDIN = MIDAS.join(ENOGG).join(FNA).join(GFF)


process buildref{
  label 'r'
  tag "$spec"
  publishDir "$params.outdir/buildref/", mode: 'rellink'
    pattern: "^$spec$"

  input:
  tuple spec, file("midas_dir"), file("eggnog.tsv"),
      file("genome.fna"), file("genes.gff") from BUILDIN

  output:
  file "$spec"
  tuple spec, file("$spec/cdsfile.tsv")
  tuple spec, file("$spec/reference.rda")

  """
  Rscript ${workflow.projectDir}/dndscv_buildref.r \
    genes.gff \
    eggnog.tsv \
    genome.fna \
    midas_dir/ \
    --outdir $spec \
    --enogg_format uhgg
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
