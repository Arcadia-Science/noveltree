#!/usr/bin/env Rscript

# Wrapper to generate html from an R markdown file.

library(optparse)
library(rmarkdown)

option_list <- list(
  make_option(c("--input_dir"), type = "character", default = NULL,
              help = "Path to the top-level directory of the results to summarize"),
  make_option(c("--html_output_file"), type = "character", default = NULL,
              help = "Path to the file to write the HTML summary to"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Path to the directory to write additional summary files to"),
  make_option(c("--rmd_file"), type = "character", default = NULL,
              help = "Path to Rmd file"),
  make_option(c("--outgroups"), type = "character", default = NULL,
              help = "Optional: outgroups used to root Asteroid tree. Defaults to 'none'")
)

args <- parse_args(OptionParser(option_list = option_list))

rmarkdown::render(input = args$rmd_file,
                  output_format = "html_document",
                  output_file = args$html_output_file,
                  params = list(input_dir = args$input_dir,
                                output_dir = args$output_dir,
                                outgroups = args$outgroups))
