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
params.sims_params = ""
params.outdir = "output"

sims_params = file(params.sims_params)

Channel
    .fromPath(sims_params)
    .splitCsv(header:true, sep: "\t")
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
      row.seed_sim)}
    .into{SIMS1; SIMS2}

process standing_variation {
  tag "$sim_id"
  label 'slim'
  // publishDir "$params.outdir/standing_variation", mode: 'rellink'

  input:
  tuple sim_id,
    run_id_short,
    Ne,
    Mu,
    Rho,
    genome_size,
    gcBurnin,
    tractlen,
    n_generations,
    scoef,
    prop_selection,
    n_pops,
    sample_size,
    seed_x1,
    seed_x2,
    seed_x3,
    seed_sim from SIMS1

  output:
  tuple sim_id, file("${sim_id}.ms") into STD_VAR

  """
  slim \
    -t \
    -m \
    -define Ne=$Ne \
    -define Mu=$Mu \
    -define Rho=$Rho \
    -define genomeSize=$genome_size \
    -define gcBurnin=$gcBurnin \
    -define tractlen=$tractlen \
    -define "runId='$sim_id'" \
    -define "runIdShort='$run_id_short'" \
    -define N_generations=$n_generations \
    -define x1=$seed_x1 \
    -define x2=$seed_x2 \
    -define x3=$seed_x3 \
    $workflow.projectDir/generate_standing_variation.slim
  """
}

process process_standing {
  tag "$sim_id"
  label 'r'
  publishDir "$params.outdir/standing_variation", mode: 'rellink',
    saveAs: {"$sim_id"}

  input:
  tuple sim_id,
    file("standing_variation/standing_variation.ms"),
    run_id_short,
    Ne,
    Mu,
    Rho,
    genome_size,
    gcBurnin,
    tractlen,
    n_generations,
    scoef,
    prop_selection,
    n_pops,
    sample_size,
    seed_x1,
    seed_x2,
    seed_x3,
    seed_sim from STD_VAR.join(SIMS2)

  output:
  tuple sim_id, file("standing_variation/"),
    run_id_short, Ne, Mu, Rho, genome_size,
    gcBurnin, tractlen, n_generations, scoef,
    n_pops, sample_size, seed_sim into STDVAR2SLIM

  """
  Rscript $workflow.projectDir/process_standing_variation.r \
    standing_variation/standing_variation.ms \
    $genome_size \
    --prop_selected $prop_selection \
    --min_maf 0.05 \
    --outdir standing_variation

  """
}



process sim_pops {
  tag "$sim_id"
  label 'slim'
  publishDir "$params.outdir/sims", mode: 'rellink'

  input:
  tuple sim_id, file("$sim_id"),
    run_id_short, Ne, Mu, Rho, genome_size,
    gcBurnin, tractlen, n_generations, scoef,
    n_pops, sample_size, seed_sim from STDVAR2SLIM

  output:
  tuple sim_id, file("$sim_id") into RES

  """
  $workflow.projectDir/sim_pops.py \
    --info_file $sim_id/snps_info.txt \
    --standing_variation $sim_id/standing_variation.ms \
    --sim_id $sim_id \
    --slim_script $workflow.projectDir/single_pop.slim \
    --n_pops $n_pops \
    --Ne $Ne \
    --Mu $Mu \
    --Rho $Rho \
    --genome_size $genome_size \
    --gcBurnin $gcBurnin \
    --tractlen $tractlen \
    --run_id_short $run_id_short \
    --n_generations $n_generations \
    --sample_size $sample_size \
    --scoef $scoef \
    --prop_selection 0.0 \
    --sim_seed $seed_sim \
    --print_period 10 \
    --outdir $sim_id
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
  withLabel: 'slim'{
    conda = '/opt/modules/pkgs/conda/4/envs/bactslim'
    module = 'slim/3.7.1'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
