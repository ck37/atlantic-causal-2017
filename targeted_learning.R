#!/usr/bin/env Rscript
####
####################################
# Based on code by:
# Susan Gruber and Mark van der Laan
# from the 2016 Atlantic Causal Inference competition.
####################################

args <- commandArgs(TRUE)
if (length(args) != 3)
  stop("usage: targeted_learning.R inFile outFile1 outFile2")

source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F,
                  verbose = T)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = T)

##############################
# Parse command-line arguments.

inFile <- args[[1]]

# outfile1 is for the overall ATT inference.
outFile1 <- args[[2]]

# outfile2 should be for individual treatment effects.
outFile2 <- args[[2]]

if (!file.exists(inFile)) stop("cannot find '", inFile, "'")

d <- read.csv(inFile)
colnames(d)[1:2] <- c("A", "Y")

if (!require(SuperLearner)) stop("Cannot find package SuperLearner")
if (!require(gam)) stop("Cannot find package gam")
if (!require(glmnet)) stop("Cannot find package glmnet")
if (!require(earth)) stop("Cannot find package earth")
#if (!require(BayesTree)) stop("Cannot find package BayesTree")
#if (!require(dbarts)) stop("Cannot find package dbarts")


SL.library <- list(c("SL.glm", "All",  "prescreen.nosq"),
                   c("SL.gam", "All", "prescreen.nosq"),
                   "SL.glmnet",
                   c("sg.gbm.2500", "prescreen.nocat"),
                   c("SL.earth", "prescreen.nosq"),
                   c("SL.bartMachine", "prescreen.nocat"))

g.SL.library <- list(c("SL.glm", "All", "prescreen.nosq"),
                     c("sg.gbm.2500", "prescreen.nocat"),
                     c("SL.gam", "All", "prescreen.nosq"),
                     c("SL.earth", "prescreen.nosq"),
                     c("SL.bartMachine", "prescreen.nocat"))

set.seed(10, "L'Ecuyer-CMRG")

results = estimate_att(A = d$A,
                      Y = d$Y,
                      W = d[, -(1:2)],
                      SL.library = SL.library,
                      g.SL.library = g.SL.library)


df_result <- data.frame(est = results$est,
                     ci_lower = results$ci_lower,
                     ci_upper = results$ci_upper)

write.csv(df_result, file = outFile1, row.names = FALSE)


