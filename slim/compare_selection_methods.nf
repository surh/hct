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
params.q_thres = 0.1
params.min_patients = 5
params.iter = 3000
params.warmup = 2000
params.chains = 4
params.vp = 5
params.vq = 5

// Load simulation parameters
Channel
  .fromPath(params.sims_params)
  .splitCsv(header:true, sep:"\t")
  .map{ row -> tuple(row.sim_id,
    row.run_id_short,
    row.Ne,
    row.Mu,
    row.Rho,
    row.genome_size,
    row.gcBurning,
    row.tractlen,
    row.n_generations,
    row.scoef,
    row.prop_selection,
    row.n_pops,
    row.sample_size,
    row.seed_x1,
    row.seed_x2,
    row.seed_x3,
    row.seed_sim
    )}
  .into{SIMPARS_SCOEF; SIMPARS}

// Get list of simulation directories
Channel.fromPath("$params.simdirs/*", type:'dir', maxDepth: 1)
  .map{simdir -> tuple(simdir.name, file(simdir))}
  .into{SIMDIRS_TEMP1; SIMDIRS_TEMP2; SIMDIRS_FIT}


// Getting additional files needed for comparisons
INFOS = SIMDIRS_TEMP1
  .map{sim_id, simdir -> tuple(sim_id, file("$simdir/snps_info.txt"))}
SIMDIRS_TEMP2
  .map{sim_id, simdir -> tuple(sim_id, file("$simdir/maf_changes.tsv"))}
  .into(MAFS_SCOEF; MAFS2)


//
// // Splitting for each selection method
// SIMDIRS.join(SIMPARS).into{SIMS1; SIMS2; SIMS3}


process s_coef{
  label 'r'
  tag "$sim_id"
  publishDir "$params.outdir/s_coef/", mode: 'rellink',
    saveAs: {"${sim_id}.tsv"}

  input:
  tuple sim_id, file(maf_changes),
    run_id_short, Ne, Mu, Rho, genome_size,
    gcBurnin, tractlen, n_generations, scoef,
    prop_selection, n_pops, sample_size, seed_x1,
    seed_x2, seed_x3,
    seed_sim from MAFS_SCOEF.join(SIMPARS_SCOEF)

  output:
  tuple sim_id, file("s_coef.tsv") into SCOEFS

  """
  Rscript $workflow.projectDir/s_coef.r \
    $maf_changes \
    --time $n_generations \
    --output s_coef.tsv
  """
}

process FIT{
  label 'r'
  label 'long'
  tag "$sim_id"
  publishDir "$params.outdir/FIT/", mode: 'rellink',
    saveAs: {"${sim_id}.tsv"}

  input:
  tuple sim_id, file(simdir) from SIMDIRS_FIT
  val n_timepoints from params.fit_n_timepoints

  output:
  tuple sim_id, file("FIT.tsv") into FITS

  """
  Rscript $workflow.projectDir/FIT.r \
    $simdir \
    --n_timepoints $n_timepoints \
    --output FIT.tsv
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
