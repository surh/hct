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

// Process Params
bern = file(params.bern)
midas_merge = file(params.midas_merge)
contig_sizes = file(params.contig_sizes)

// Prepare command
opt_pars = ''

// Getting speclist
Channel.fromPath("$bern/*")
 .map{bernfile -> bernfile.name[ 0..<bernfile.name.indexOf('.') ]}
 .into{SPECSINFO; SPECSCTGS; SPECSANNOT}

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
  contig_sizes = file(contig_sizes)
  println("Reading contig sizes files")
  CTGS = Channel.fromPath("$contig_sizes/*")
    .map{ctgfile -> tuple(ctgfile.name[ 0..<ctgfile.name.indexOf('.') ],
      file(ctgfile))}
  opt_pars = opt_pars + ' --manhattan --contig_sizes contig_sizes.txt'
  // CTGS.subscribe{println it}
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
  opt_pars = opt_pars + ' --annot gene_annotations.txt'
}

if(params.OR_trans){
  opt_pars = opt_pars + ' --OR_trans'
}

println(opt_pars)
