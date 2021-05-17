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
params.input = ''
params.outdir = 'output'
params.q_thres = 0.1
params.min_patients = 5


// Process params
indir = file(params.input)

SITES = Channel.fromPath("$indir/*")
  .map{ sites_file -> tuple(sites_file.name.replaceAll(/\.tsv\.gz$/, ''),
    file(sites_file)) }
