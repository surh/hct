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

// Parameters
params.sims_params = ""
params.outdir = "output"


SIMS = Channel
  .fromPath(params.sims_params)
  .splitCsv(header:true, sep:"\t")
  .map{ row -> tuple(row.sim_id,
    row.f,
    row.g,
    row.x_0,
    row.nsites,
    row.npops,
    row.p_imm,
    row.prop_selected,
    row.popsize,
    row.T,
    row.seed)}

process sim_moran_imm {
  tag "$sim_id"
  publishDir "$params.outdir/", mode: 'rellink',
    saveAs:{"$sim_id/"}
  label 'r'

  input:
  tuple val(sim_id), f, g, x_0, nsites, npops, p_imm, prop_selected,
    popsize, T, seed from SIMS

  output:
  file "output/"

  """
  ${workflow.projectDir}/sim_moran_imm.r \
    --nsites $nsites \
    --npops $npops \
    --pop_size $popsize \
    --p_imm $p_imm \
    --x_0 $x_0 \
    --f $f \
    --g $g \
    --T $T \
    --prop_selected $prop_selected \
    --seed $seed \
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
    // module = "R/4.2.0" // Make sure you have ~/.R/Makevars with CXX14
  }
  withLabel: 'py3'{
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/fraserconda'
    # conda = '/home/groups/hbfraser/modules/packages/conda/4.6.14/envs/fraserconda'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
