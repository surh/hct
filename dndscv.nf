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
params.pdir = ''
params.outdir = 'output'
params.maf_thres = 0.8
params.max_coding_muts_per_sample = 100
params.max_muts_per_gene_per_sample = 10
genetic_code = 1

// Process params
gff_dir = file(params.gff)
enogg_dir = file(params.enogg)
fna_dir = file(params.fna)
midas_dir = file(params.midas)
pdir_dir = file(params.pdir)

GFF = Channel.fromPath("$gff_dir/*")
  .map{ infile -> tuple(infile.name.replaceAll(/\.gff$/, ''),
    file(infile))}
ENOGG = Channel.fromPath("$enogg_dir/*")
  .map{ infile -> tuple(infile.name.replaceAll(/\.tsv$/, ''),
    file(infile))}
PDIR = Channel.fromPath("$pdir_dir/*")
  .map{ infile -> tuple(infile.name.replaceAll(/\.tsv\.gz$/, ''),
    file(infile))}
FNA = Channel.fromPath("$fna_dir/*")
  .map{ infile -> tuple(infile.name.replaceAll(/\.fna$/, ''),
    file(infile))}
Channel.fromPath("$midas_dir/*", type: 'dir', maxDepth: 1)
  .map{ infile -> tuple(infile.name,
    file(infile)) }
  .into{MIDAS1; MIDAS2}
INFO = Channel.fromPath("$midas_dir/*", type: 'dir', maxDepth: 1)
  .map{ infile -> tuple(infile.name,
    file("$infile/snps_info.txt")) }

BUILDIN = MIDAS1.join(ENOGG).join(FNA).join(GFF).view()

process buildref{
  label 'r'
  tag "$spec"
  publishDir "$params.outdir/buildref/", mode: 'rellink',
    pattern: "$spec"

  input:
  tuple spec, file("midas_dir"), file("eggnog.tsv"),
      file("genome.fna"), file("genes.gff") from BUILDIN

  output:
  file "$spec"
  tuple spec, file("$spec/cdsfile.tsv") into CDS
  tuple spec, file("$spec/reference.rda") into REF

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

mode = 'all_dummy'
if (mode == 'all_dummy' || mode == "dummy_singletons"){
  curr_max_coding_muts_per_sample = "Inf"
  curr_max_muts_per_gene_per_sample = "Inf"
}

process dndscv{
  label 'r'
  tag "$mode $spec"
  publishDir "$params.outdir/dndscv/", mode: 'rellink',
    pattern: "$spec"

  input:
  val mode from mode
  val maf_thres from params.maf_thres
  val max_coding_muts_per_sample from curr_max_coding_muts_per_sample
  val max_muts_per_gene_per_sample from curr_max_muts_per_gene_per_sample
  genetic_code = 1
  tuple spec, file('midas_dir'), file('reference.rda') from MIDAS2.join(REF)

  output:
  tuple spec, mode, file("$spec/dnds_cv.tsv") into DNDSCV
  file "$spec"

  """
  Rscript ${workflow.projectDir}/dndscv_run.r \
    midas_dir/ \
    reference.rda \
    --mode $mode \
    --maf_thres \
    --max_coding_muts_per_sample $max_coding_muts_per_sample \
    --max_muts_per_gene_per_sample $max_muts_per_gene_per_sample \
    --genetic_code \
    --outdir $spec
  """
}


proces plot{
  label 'r'
  tag "$mode $spec"
  publishDir "$params.outdir/plot/$mode", mode: 'rellink'

  input:
  tuple spec, mode, file('dnds_csv.tsv'),
    file('cdsfile.tsv'),
    file('snps_info.txt'),
    file('p_directional.tsv.gz') from DNDSCV.join(CDS).join(INFO).join(PDIR)

  output:
  file "$spec"

  """
  Rscript ${workflow.projectDir}/dndscv_vs_pdirectional.r \
    dnds_csv.tsv \
    cdsfile.tsv \
    snps_info.txt \
    p_directional.tsv.gz \
    --outdir $spec
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
