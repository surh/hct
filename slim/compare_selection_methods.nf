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
    row.run_id_short,
    row.Ne,
    row.Mu,
    row.Rho,
    row.genome_size,
    row.gcBurning,
    row.tractlen,
    row.n_generations
    )}
