#!/usr/bin/env Rscript

# (C) Copyright 2022 Sur Herrera Paredes
# This file is part of This program.
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
# along with This program.  If not, see <https://www.gnu.org/licenses/>.

# library(knitr)
library(tidyverse)

cmd_args <- commandArgs(trailingOnly = FALSE)
ii <- cmd_args %>% str_detect("^--file=")
if(sum(ii) == 0){
  stop("ERROR: file not found")
}else if(sum(ii) > 1){
  stop("ERROR: more than one file parameters")
}
r_script <- cmd_args[ ii ] %>% str_remove("^--file=")
script_dir <- dirname(r_script)
script_to_render <- file.path(script_dir, "compare_tests_sim_moran_imm.r")

rmarkdown::render(input = script_to_render,
                  output_format = "html_document",
                  output_dir = getwd(),
                  intermediates_dir = getwd())
