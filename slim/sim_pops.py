#!/usr/bin/env python

# Copyright (C) 2022 Sur Herrera Paredes
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import random
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Run multiple replicates of Bacterial SLiMulations")

    # Define required arguments
    required.add_argument("--standing_variation",
                          help=("ms file with standing variation"),
                          required=True,
                          type=str)
    required.add_argument("--info_file",
                          help=("tsv file with header & sites for selection "
                                "in the 5-th column."),
                          required=True,
                          type=str)
    required.add_argument("--sim_id",
                          help="Sim ID",
                          type=str,
                          required=True)
    required.add_argument("--slim_script",
                          help="SLiM script",
                          type=str,
                          required=True)

    # Define other arguments
    parser.add_argument("--n_pops",
                        help="Number of replicate populations",
                        type=int,
                        default=1)
    parser.add_argument("--Ne",
                        help=("Effective population size"),
                        type=int,
                        default=1000)
    parser.add_argument("--Mu",
                        help="Mutation rate",
                        type=float,
                        default=1e-7)
    parser.add_argument("--Rho",
                        help="Recombination rate",
                        type=float,
                        default=1e-7)
    parser.add_argument("--genome_size",
                        help="Genome size in bp",
                        type=int,
                        default=3e6)
    parser.add_argument("--gcBurnin",
                        help="gcBurnin",
                        type=float,
                        default=0)
    parser.add_argument("--tractlen",
                        help="Average HGT tract length in bp.",
                        type=int,
                        default=1e5)
    parser.add_argument("--run_id_short",
                        help="Short run_id",
                        type=str)
    parser.add_argument("--n_generations",
                        help="Number of generations.",
                        type=int,
                        default=50)
    parser.add_argument("--sample_size",
                        help="Sample size for each time point.",
                        type=int,
                        default=20)
    parser.add_argument("--scoef",
                        help="Selection advantage.",
                        type=float,
                        default=0.01)
    parser.add_argument("--prop_selection",
                        help="Proportion of de novo sites under selection.",
                        type=float,
                        default=0.00)
    parser.add_argument("--sim_seed",
                        help="Seed for simulation",
                        type=int)
    parser.add_argument("--print_period",
                        help="How often to print",
                        type=int,
                        default=10)
    parser.add_argument("--outdir",
                        help="output directory",
                        type=str,
                        default="output")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.sim_seed is None:
        print("Creating automatic seed")
        args.sim_seed = random.randrange(1, 10000000)

    if args.run_id_short is None:
        print("Setting run_id_short to sim_id")
        args.run_id_short = args.sim_id

    return args


if __name__ == '__main__':
    args = process_arguments()
    print(args)

    random.seed(args.sim_seed)
    Seeds = []
    for pop_i in range(1, args.n_pops + 1):

        # Seed for specific pop simulation
        seed = random.randrange(1, 10000000)
        Seeds.append(seed)

        # outdir for specific pop
        outdir = os.path.join(args.outdir, "pop_" + str(pop_i))
        os.mkdir(outdir)

        cmd = ['slim',
               "-t",  # print SLiM's total execution time (in user clock time)
               "-m",  # print SLiM's peak memory usage
               "-define", """"runId='{}'" """.format(args.sim_id),
               "-define", """"runIdShort='{}'" """.format(args.run_id_short),
               "-define", """"outdir='{}'" """.format(outdir),
               "-define", """"info_file='{}'" """.format(args.info_file),

               "-define", """"standing_variation='{}'" """
               .format(args.standing_variation),

               "-define", "Ne={}".format(args.Ne),
               "-define", "Mu={}".format(args.Mu),
               "-define", "Rho={}".format(args.Rho),
               "-define", "genomeSize={}".format(args.genome_size),
               "-define", "gcBurnin={}".format(args.gcBurnin),
               "-define", "tractlen={}".format(args.tractlen),
               "-define", "N_generations={}".format(args.n_generations),
               "-define", "sampleSize={}".format(args.sample_size),
               "-define", "scoef1={}".format(args.scoef),
               "-define", "prop_selection={}".format(args.prop_selection),
               "-define", "seed={}".format(seed),
               "-define", "print_period={}".format(args.print_period),

               args.slim_script]
        cmd = " ".join(cmd)
        print("Running\n>" + cmd)
        os.system(cmd)

    with open(args.outdir + "pop_seeds.txt", 'w') as oh:
        pop_i = 0
        for seed in Seeds:
            pop_i = pop_i + 1
            line = ["pop_" + str(pop_i), str(seed)]
            line = "\t".join(line)
            oh.write(line + "\n")
    oh.close()
