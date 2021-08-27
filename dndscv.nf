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
params.genetic_code = 1
params.mut_modes = ["dummy_singletons", "all_dummy", "random_top_thres",
  "all_thres", "all", "true_singletons", "thres_singletons"]

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

BUILDIN = MIDAS1.join(ENOGG).join(FNA).join(GFF)

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

DNDSIN = MIDAS2.join(REF).combine(params.mut_modes)
  .map{ elems ->

    spec = elems[0]
    midas_dir = elems[1]
    ref_rda = elems[2]
    mut_mode = elems[3]

    if (mut_mode == 'all_dummy' || mut_mode == "dummy_singletons"){
      max_coding_muts_per_sample = "Inf"
      max_muts_per_gene_per_sample = "Inf"
    }else{
      max_coding_muts_per_sample = params.max_coding_muts_per_sample
      max_muts_per_gene_per_sample = params.max_muts_per_gene_per_sample
    }

    tuple(spec, file(midas_dir), file(ref_rda), mut_mode,
      max_coding_muts_per_sample, max_muts_per_gene_per_sample) }

process dndscv{
  label 'dndscv'
  label 'r'
  tag "$mut_mode $spec"
  publishDir "$params.outdir/dndscv/$mut_mode", mode: 'rellink',
    pattern: "$spec"

  input:
  tuple spec, file('midas_dir'), file('reference.rda'),
    mut_mode,
    max_coding_muts_per_sample,
    max_muts_per_gene_per_sample from DNDSIN
  val maf_thres from params.maf_thres
  val genetic_code from params.genetic_code

  output:
  tuple spec, mut_mode, file("$spec/dnds_cv.tsv") optional true into DNDSCV
  file "$spec" optional true

  """
  Rscript ${workflow.projectDir}/dndscv_run.r \
    midas_dir/ \
    reference.rda \
    --mode $mut_mode \
    --maf_thres $maf_thres \
    --max_coding_muts_per_sample $max_coding_muts_per_sample \
    --max_muts_per_gene_per_sample $max_muts_per_gene_per_sample \
    --genetic_code $genetic_code \
    --outdir $spec
  """
}

process compare_pdir{
  label 'r'
  tag "$mut_mode $spec"
  publishDir "$params.outdir/plot/$mut_mode", mode: 'rellink'

  input:
  tuple spec, mut_mode, file('dnds_csv.tsv'),
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
  time = '20m'
  memory = '2G'
  withLabel: 'r'{
    module = 'R/4.1.0'
    // module = "R/4.0.2:v8/8.4.371.22" // Make sure you have ~/.R/Makevars with CXX14
  }
  withLabel: 'dndscv'{
    errorStrategy = { task.attempt < 3 ? 'retry' : 'finish' }
    time = { 60.m + (task.attempt - 1) * 120.m }
    memory = { 2.GB + (task.attempt - 1) * 4.GB }
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
