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
Channel.fromPath("$params.simdirs/*", type:'dir', maxDepth: 1)
  .map{simdir -> tuple(simdir.name, file(simdir))}
  .into{SIMDIRS; SIMDIRS_TEMP}

// Splitting for each selection method
SIMDIRS.join(SIMPARS).into{SIMS1; SIMS2; SIMS3}

INFOS = SIMDIRS_TEMP
  .map{sim_id, simdir -> tuple(sim_id, file("$simdir/snps_info.txt"))}

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
  Rscript $workflow.projectDir/s_coef_sim_moran_imm.r \
    $simdir \
    --N $popsize \
    --T $T \
    --x_0 $x_0 \
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
  tuple val(sim_id), file(simdir), f, g, x_0, nsites,
    npops, p_imm, prop_selected, popsize, T, seed from SIMS2
  val n_timepoints from params.fit_n_timepoints

  output:
  tuple sim_id, file("FIT.tsv") into FITS

  """
  Rscript $workflow.projectDir/FIT_sim_moran_imm.r \
    $simdir \
    --n_timepoints $n_timepoints \
    --output FIT.tsv
  """
}

process bern_mix{
  cpus params.chains
  tag "$sim_id"
  label 'r'
  label 'bern'
  publishDir "$params.outdir/stan_models", mode: 'rellink',
    pattern: "output/m1.stan.rdat", saveAs: {"${sim_id}.stan.rdat"}
  publishDir "$params.outdir/p_directional", mode: 'rellink',
    pattern: "output/p_directional.tsv.gz", saveAs: {"${sim_id}.tsv.gz"}
  publishDir "$params.outdir/model_summaries", mode: 'rellink',
    pattern: "output/model_summaries.tsv.gz", saveAs: {"${sim_id}.tsv.gz"}
  publishDir "$params.outdir/CHECK_RHAT", mode: 'rellink',
    pattern: "CHECK_RHAT", saveAs: {"$sim_id"}

  input:
  tuple val(sim_id), file(simdir), f, g, x_0, nsites,
    npops, p_imm, prop_selected, popsize, T, seed from SIMS3
  val q_thres from params.q_thres
  val min_patients from params.min_patients
  val iter from params.iter
  val warmup from params.warmup
  val chains from params.chains
  val vp from params.vp
  val vq from params.vq

  output:
  file "output/m1.stan.rdat" optional true
  tuple sim_id, "output/p_directional.tsv.gz" optional true into PDIRS
  file  "output/model_summaries.tsv.gz" optional true
  file "CHECK_RHAT" optional true

  """
  Rscript $workflow.projectDir/../bern_mix.r \
    $simdir/sites.tsv \
    --q_thres $q_thres \
    --min_patients $min_patients \
    --outdir output \
    --iter $iter \
    --warmup $warmup \
    --chains $chains \
    --vp $vp \
    --vq $vq
  """
}

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
  tuple val(sim_id), file(s_coef), file(fit),
    file(pdir), file(info) from SCOEFS.join(FITS).join(PDIRS).join(INFOS)

  output:
  file "output/"
  file "compare_tests_sim_moran_imm.html"

  """
  Rscript $workflow.projectDir/render_compare_tests_sim_moran_imm.r \
    --s_coef $s_coef \
    --FIT $fit \
    --pdir $pdir\
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
    // module = "R/4.0.2:v8/8.4.371.22:pandoc/2.7.3" // Make sure you have ~/.R/Makevars with CXX14
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
