#!/usr/bin/env Rscript

####################################
# This script is intended to be run from the command line.
# Specify the input csv and two output csvs.
# Example:
# ./targeted_learning.R input.csv outfile1.csv outfile2.csv
# outfile1.csv will contain the average treatment effect on the treated.
# outifle2.csv will contain unit-level treatment effect estimates.

####################################
# Based on code by:
# Susan Gruber and Mark van der Laan
# from the 2016 Atlantic Causal Inference competition.
####################################

conf = list(
  # Set to T for extra output during execution.
  verbose = T,

  # Set auto-install to T for code to install any missing packages.
  auto_install = F,

  # Set to T to use simple SL libraries for testing purposes.
  debug = T,

  # Use up to this many cores if available.
  max_cores = 4
)

args <- commandArgs(TRUE)

if (length(args) != 3) {
  cat("Detected", length(args), "arguments but expected 3.\n")
  stop("usage: targeted_learning.R inFile outFile1 outFile2")
}

source("lib/function_library.R")

load_all_packages(auto_install = conf$auto_install, verbose = conf$verbose)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = conf$verbose)

##############################
# Parse command-line arguments.

inFile <- args[[1]]

# outfile1 is for the overall ATT inference.
outFile1 <- args[[2]]

# outfile2 should be for individual treatment effects.
outFile2 <- args[[2]]

if (!file.exists(inFile)) stop("Cannot find '", inFile, "'")

d <- read.csv(inFile)
colnames(d)[1:2] <- c("A", "Y")

if (!require(SuperLearner)) stop("Cannot find package SuperLearner")
if (!require(gam)) stop("Cannot find package gam")
if (!require(glmnet)) stop("Cannot find package glmnet")
if (!require(earth)) stop("Cannot find package earth")
#if (!require(BayesTree)) stop("Cannot find package BayesTree")
#if (!require(dbarts)) stop("Cannot find package dbarts")

# Setup parallelization? Use up to 4 cores.
num_cores = RhpcBLASctl::get_num_cores()

if (conf$verbose) {
  cat("Cores detected:", num_cores, "\n")
}

use_cores = min(num_cores, conf$max_cores)
options("mc.cores" = use_cores)

if (conf$verbose) {
  # Check how many parallel workers we are using:
  cat("Cores used:", getOption("mc.cores"), "\n")
}

# This is the Q and g library.
if (conf$debug) {
  SL.library = c("SL.glm", "SL.mean")
} else {
  SL.library <- list(c("SL.glm", "All",  "prescreen.nosq"),
                     # Not working, can we fix it?
                     # c("SL.gam", "All", "prescreen.nosq"),
                     #c("sg.gbm.2500", "prescreen.nocat"),
                     "SL.xgboost",
                     "SL.randomForest",
                     "SL.glmnet",
                     c("SL.earth", "prescreen.nosq"),
                     c("SL.bartMachine", "prescreen.nocat"),
                     "SL.mean")
}

# Just use the same library for g and Q.
g.SL.library = SL.library

#if (conf$debug) {
  #g.SL.library = c("SL.glm", "SL.mean")
#} else {
  #g.SL.library <- list(c("SL.glm", "All", "prescreen.nosq"),
                       #c("sg.gbm.2500", "prescreen.nocat"),
                       #"SL.glmnet",
                       # Not working:
                       #c("SL.gam", "All", "prescreen.nosq"),
  #                     c("SL.earth", "prescreen.nosq"))#,
                      # c("SL.bartMachine", "prescreen.nocat"))
#}

set.seed(10, "L'Ecuyer-CMRG")

results = estimate_att(A = d$A,
                       Y = d$Y,
                       W = d[, -(1:2)],
                       verbose = conf$verbose,
                       SL.library = SL.library,
                       g.SL.library = g.SL.library)

cat("\nResults:\n")
print(results)

df_result <- data.frame(est = results$est,
                     ci_lower = results$ci_lower,
                     ci_upper = results$ci_upper)

# Save point estimates and CI to outFile1
write.csv(df_result, file = outFile1, row.names = FALSE)

# Save unit-level estimates to outFile2.
write.csv(results$unit_estimates, file = outFile2, row.names = FALSE)

