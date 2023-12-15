#!/usr/bin/env Rscript

library("argparse")
library("styler")

parser <- ArgumentParser(description = "Call styler")
parser$add_argument("dir",
  metavar = "DIR", type = "character",
  help = "File to parse"
)
parser$add_argument("--dry",
  choices = c("off", "on", "fail"), default = "on"
)
args <- parser$parse_args()
result <- style_dir(args$dir, dry = args$dry)
