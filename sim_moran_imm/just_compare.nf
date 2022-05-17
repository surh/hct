#!/usr/bin/env nextflow

// Copyright (C) 2022 Sur Herrera Paredes

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

// parameters
params.sims_dir = ''
params.scoef_dir = ''
params.FIT_dir = ''
params.pdir_dir = ''
params.outdir = 'output'

// process params
sims_dir = file(params.sims_dir)
scoef_dir = file(params.scoef_dir)
FIT_dir = file(params.FIT_dir)
pdir_dir = file(params.pdir_dir)


SIMS = Channel.fromPath("$sims_dir/*", type: 'dir', maxDepth: 1)
  .map{simdir -> tuple(simdir.name,
    file("$simdir/maf_changes.tsv.gz"),
    file("$simdir/snps_info.txt.gz"))}

SCOEFS = Channel.fromPath("$scoef_dir/*", type: 'file')
  .map{resfile -> tuple(resfile.name.replaceAll(/\.tsv/. ''),
    file(resfile))}

FITS = Channel.fromPath("$FIT_dir/*", type: 'file')
  .map{resfile -> tuple(resfile.name.replaceAll(/\.tsv$/. ''),
    file(resfile))}

PDIRS = Channel.fromPath("$pdir_dir/*", type: 'file')
  .map{resfile -> tuple(resfile.name.replaceAll(/\.tsv.gz$/. ''),
    file(resfile))}


process compare{
  tag "$sim_id"
  label 'r'
  publishDir "$params.outdir/comp_htmls/", mode: 'rellink',
    pattern: "compare_tests_sim_moran_imm.html",
    saveAs:{"${sim_id}.html"}
  publishDir "$params.outdir/comparisons/", mode: 'rellink',
    pattern: "output",
    saveAs:{"$sim_id"}

  input:
  tuple val(sim_id), file(maf), file(info), file(s_coef), file(fit),
    file(pdir) from SIMS.join(SCOEFS).join(FITS).join(PDIRS)

  output:
  file "output/"
  file "compare_tests_sim_moran_imm.html"

  """
  Rscript $workflow.projectDir/render_compare_tests_sim_moran_imm.r \
    --s_coef $s_coef \
    --FIT $fit \
    --pdir $pdir\
    --maf $maf \
    --info $info \
    --outdir output
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
  withLabel: 'long'{
    time = '48h'
  }
  withLabel: 'bern'{
    time = '48h'
    memory = '8G'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
