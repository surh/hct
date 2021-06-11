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

// Parametrs
params.bern = ''
params.midas_merge = ''
params.outdir = 'output'
params.manhattan = false
params.contig_sizes = ''
params.locus_test = false
params.ns_test = false


// Process Params
bern = file(params.bern)
midas_merge = file(params.midas_merge)
println("$bern")

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
      file(ctgsize))}
}
