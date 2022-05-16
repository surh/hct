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
    required.add_argument("--run_id",
                          help="Run ID",
                          type=str,
                          required=True)

    # Define other arguments
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
    parser.add_argument("--gcBurning",
                        help="gcBurnin",
                        type=float,
                        default=0)
    parser.add_argument("--tractlen",
                        help="Average HGT tract length in bp.",
                        type=int,
                        default=1e5)
    parser.add_argument("--rund_id_short",
                        help="Short run_id",
                        type=str,
                        default='short_id')
    parser.add_argument("--n_generations",
                        help="Number of generations.",
                        type=int,
                        default=50)
    parser.add_argument("--sample_size",
                        help="Sample size for each time point.",
                        type=int,
                        default=20)
    parser.add_argument("--scoef",
                        help="Selection coefficient.",
                        type=float,
                        default=0.01)





    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args
