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

// A problem that is not easy to solve is how to deal with not all
// inputs present for process

// Parametrs
params.bern = ''
params.midas_merge = ''
params.outdir = 'output'
params.manhattan = false
params.contig_sizes = ''
params.locus_test = false
params.ns_test = false
params.annot = false
params.min_size = 5
params.OR_trans = false
params.bool_labs = false // Must be quoted and space-separated, e.g "a b"

// Process Params
bern = file(params.bern)
midas_merge = file(params.midas_merge)

// Prepare command
opt_pars = ''

// Getting speclist
Channel.fromPath("$bern/*")
 .map{bernfile -> bernfile.name[ 0..<bernfile.name.indexOf('.') ]}
 .into{SPECSCTGS; SPECSANNOT}

//  Reading basic files
BERN = Channel.fromPath("$bern/*")
  .map{bernfile -> tuple(bernfile.name[ 0..<bernfile.name.indexOf('.') ],
    file(bernfile))}
// bernfile.name[ 0..<bernfile.name.indexOf('.') ]
// BERN.subscribe{println it}

INFO = Channel.fromPath("$midas_merge/*", type: 'dir', maxDepth: 1)
  .map{midas_dir -> tuple(midas_dir.name,
    file("$midas_dir/snps_info.txt"))}
// INFO.subscribe{println it}

if(params.manhattan){
  contig_sizes = file(params.contig_sizes)
  // println("Reading contig sizes files")
  CTGS = Channel.fromPath("$contig_sizes/*")
    .map{ctgfile -> tuple(ctgfile.name[ 0..<ctgfile.name.indexOf('.') ],
      file(ctgfile))}
  opt_pars = opt_pars + ' --manhattan --contig_sizes contig_sizes'
  // CTGS.subscribe{println it}
}else{
  CTGS = SPECSCTGS
    .map{ spec -> tuple(spec, "dummy_ctg")}
}

if(params.locus_test){
  opt_pars = opt_pars + ' --locus_test'
}

if(params.ns_test){
  opt_pars = opt_pars + ' --ns_test'
}

if(params.annot){
  annot = file(params.annot)
  ANNOT = Channel.fromPath("$annot/*")
    .map{ annotfile -> tuple(annotfile.name[ 0..<annotfile.name.indexOf('.') ],
      file(annotfile))}
  // ANNOT.subscribe{println it}
  opt_pars = opt_pars + ' --annot gene_annotations'
}else{
  ANNOT = SPECSANNOT
    .map{spec -> tuple(spec, 'dummy_annot')}
}

if(params.OR_trans){
  opt_pars = opt_pars + ' --OR_trans'
}

if(params.bool_labs){
  opt_pars = opt_pars + " --bool_labs $params.bool_labs"
}

// println(opt_pars)
ALLIN = BERN.join(INFO).join(CTGS).join(ANNOT)

process fgsea {
  label 'r'
  tag "$spec"
  publishDir params.outdir, mode: 'rellink',
    pattern: "output/bern_manhattan.png",
    saveAs: {"${spec}.manhattan.png"}
  publishDir params.outdir, mode: 'rellink',
    pattern: "output/locus_test.tsv",
    saveAs: {"${spec}.locus.tsv"}
  publishDir params.outdir, mode: 'rellink',
    pattern: "output/ns_test.tsv",
    saveAs: {"${spec}.ns.tsv"}
  publishDir params.outdir, mode: 'rellink',
    pattern: "output/annot_test.tsv",
    saveAs: {"${spec}.annot.tsv"}

  input:
  tuple spec, file("bern_file"), file("snps_info"), file("contig_sizes"), file("gene_annotations") from ALLIN
  val min_size from params.min_size
  val opt_pars

  output:
  file "output/bern_manhattan.png" optional true
  file "output/locus_test.tsv" optional true
  file "output/ns_test.tsv" optional true
  file "output/annot_test.tsv" optional true

  """
  Rscript $workflow.projectDir/fgsea_bern.r bern_file snps_info \
    --min_size $min_size \
    --outdir output \
    $opt_pars
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
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
