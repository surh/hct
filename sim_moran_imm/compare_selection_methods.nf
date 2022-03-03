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

// Parameters
params.simdirs = 'sims/'
params.sims_params = ''
params.outdir = 'output/'
// FIT parameters

params.fit_n_timepoints = 3

// Bern parameters

// Load simulation parameters
SIMPARS = Channel
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

// Get list of simulation directories
SIMDIRS = Channel.fromPath("$params.simdirs/*", type:'dir', maxDepth: 1)
  .map{simdir -> tuple(simdir.name, file(simdir))}

// Splitting for each selection method
SIMDIRS.join(SIMPARS).into{SIMS1; SIMS2; SIMS3}

process s_coef{
  label 'r'
  tag "$sim_id"
  publishDir "$params.outdir/s_coef/", mode: 'rellink',
    saveAs: {"${sim_id}.tsv"}

  input:
  tuple val(sim_id), file(simdir), f, g, x_0, nsites,
    npops, p_imm, prop_selected, popsize, T, seed from SIMS1

  output:
  tuple sim_id, file("s_coef.tsv") into SCOEFS

  """
  $workflow.projectDir/s_coef_sim_moran_imm.r \
    $simdir \
    --N $popsize \
    --T $T \
    --x_0 $x_0 \
    --output s_coef.tsv
  """
}


process FIT{
  label 'r'
  tag "$sim_id"
  publishDir "$params.outdir/FIT/", mode: 'rellink',
    saveAs: {"${sim_id}.tsv"}

  input:
  tuple val(sim_id), file(simdir), f, g, x_0, nsites,
    npops, p_imm, prop_selected, popsize, T, seed from SIMS1
  val n_timepoints from params.fit_n_timepoints

  output:
  tuple sim_id, file("FIT.tsv") into FITS

  """
  $workflow.projectDir/s_coef_sim_moran_imm.r \
    $simdir \
    --n_timepoints $n_timepoints \
    --output FIT.tsv
  """
}
