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
  debug = F,

  # Set to F to disable parallelism.
  parallel = F,

  # Use up to this many cores if available.
  max_cores = 4,
  
  # Maximum amount of memory to allow rJava heap to use for bartMachine.
  java_mem = "4g",
  
  # Number of rows to trigger larger java memory allocation.
  large_data_rows = 1000,
  
  # Higher java memory setting when dataset has more rows than
  # conf$large_data_rows . This is intended to avoid out-of-memory errors.
  java_large_mem = "16g"
)

args <- commandArgs(TRUE)

if (length(args) != 3) {
  cat("Detected", length(args), "arguments but expected 3.\n")
  stop("usage: targeted_learning.R inFile outFile1 outFile2")
}

source("lib/function_library.R")

##############################
# Parse command-line arguments.

inFile <- args[[1]]

# outfile1 is for the overall ATT inference.
outFile1 <- args[[2]]

# outfile2 should be for individual treatment effects.
outFile2 <- args[[3]]

if (!file.exists(inFile)) stop("Cannot find '", inFile, "'")

d <- read.csv(inFile)

# Check for larger dataset and allocate more RAM to Java for bartMachine.
if (nrow(d) > conf$large_data_rows) {
  
  # Switch to large memory heap size for rJava.
  conf$java_mem = conf$java_large_mem
  
  if (conf$verbose) {
    cat("Found larger data file",
        paste0("(", prettyNum(nrow(d), big.mark = ","), " rows)."),
        "Allocating", conf$java_mem, "ram to java for bartMachine.\n")
  }
}

#####################################

load_all_packages(auto_install = conf$auto_install,
                  verbose = conf$verbose,
                  java_mem = conf$java_mem)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = conf$verbose)



# The first column is treatment indicator A and the second column is the
# continuous outcome.
colnames(d)[1:2] <- c("A", "Y")

if (!require(SuperLearner)) stop("Cannot find package SuperLearner")
#if (!require(gam)) stop("Cannot find package gam")
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

if (conf$parallel) {
  options("mc.cores" = use_cores)
}

if (conf$verbose) {
  # Check how many parallel workers we are using:
  cat("Cores used:", getOption("mc.cores"), "\n")
}

if (conf$verbose) {
  cat("Setting bartMachine cores to", use_cores, "\n")
}

# Speed up bartMachine by enabling multicore execution.
bartMachine::set_bart_machine_num_cores(use_cores)

# This is the Q and g library.
if (conf$debug) {
  q_lib = c("SL.glm", "SL.mean")
  g_lib = q_lib
} else {
  q_lib = c(list(# speedglm doesn't work :/ just use plain ol' glm.
               c("SL.glm", "All", "screen.corRank4", "screen.corRank8", "prescreen.nosq")#,
               #c("SL.mgcv", "All", "prescreen.nosq"),
               #c("sg.gbm.2500", "prescreen.nocat"),
               #"SL.xgboost",
               #"SL.xgboost_threads_4"
               # Effect modification learners can't be used with g, only Q.
            ),
            # create.Learner() grids.
            sl_glmnet_em15$names,
            sl_xgb$names,
            # Temporarily turn off SVM due to errors Vince is getting.
            #sl_ksvm$names, 
            list(
               #"SL.randomForest_fast",
               "SL.ranger_fast",
               c("SL.glmnet_fast", "All", "screen.corRank4", "screen.corRank8"),
               c("SL.nnet", "All", "screen.corRank4", "screen.corRank8"),
               c("SL.earth", "prescreen.nosq"),
               # Works only if parallel = F. Do not use with mcSuperlearner!
               "SL.bartMachine",
               "SL.mean"))
  
  # Need a separate g lib that does not include effect modification learners.
  g_lib = c(list(c("SL.glm", "All", "screen.corRank4", "screen.corRank8", "prescreen.nosq"),
               #c("SL.mgcv", "All", "prescreen.nosq"),
               #c("sg.gbm.2500", "prescreen.nocat"),
               #"SL.xgboost",
               #"SL.xgboost_threads_4",
               "SL.ranger_fast"#,
              ), # create.Learner() grids.
            sl_xgb$names,
            # Temporarily turn off SVM due to errors Vince is getting.
            #sl_ksvm$names, 
            list(
               c("SL.glmnet_fast", "All", "screen.corRank4", "screen.corRank8"),
               c("SL.nnet", "All", "screen.corRank4", "screen.corRank8"),
               c("SL.earth", "prescreen.nosq"),
               # Works only if parallel = F. Do not use with mcSuperlearner!
               "SL.bartMachine",
               "SL.mean"))
}

# Just use the same library for g and Q.
# g.SL.library = SL.library

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
                       SL.library = q_lib,
                       #parallel = conf$parallel,
                       parallel = F,
                       g.SL.library = g_lib)


# Extract unit-level estimates before printing.
unit_estimates = results$unit_estimates
results$unit_estimates = NULL

cat("\nResults:\n")
print(results)

cat("\nSummary of unit-level estimates:\n")
print(summary(unit_estimates))

df_result <- data.frame(est = results$est,
                     ci_lower = results$ci_lower,
                     ci_upper = results$ci_upper)

# Save point estimates and CI to outFile1
write.csv(df_result, file = outFile1, row.names = FALSE)

# Save unit-level estimates to outFile2.
write.csv(unit_estimates, file = outFile2, row.names = FALSE)

