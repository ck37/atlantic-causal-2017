#!/usr/bin/env Rscript

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
outFile <- args[[2]]
# TODO: outfile2 should be for individual treatment effects.


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
                   c("SL.bart", "prescreen.nocat"))

g.SL.library <- list(c("SL.glm", "All", "prescreen.nosq"),
                     c("sg.gbm.2500", "prescreen.nocat"),
                     c("SL.gam", "All", "prescreen.nosq"),
                     c("SL.earth", "prescreen.nosq"),
                     c("SL.bart", "prescreen.nocat"))

gbounds <- c(0.01, 0.99)
depsilon <-  0.001
V <- 5
alpha<- c(.0005, .9995)
# convert factors to dummies
n <- nrow(d)
W <- model.matrix(as.formula(Y ~ .), data = d[,-1])[,-1]
nonbinary <- which(colMeans(W) > 1)  # this isn't general, but true for these data
Wcat <- matrix(as.integer(W[,nonbinary] < rep(colMeans(W[,nonbinary]), each = n)), nrow = n, byrow = FALSE)
colnames(Wcat) <- paste0(colnames(W[,nonbinary]), "cat")
W <-cbind(W, x_3aug = as.integer(W[,"x_3"] > 0), x_4aug = as.integer(W[,"x_4"] > 0),
 				Wcat
 				 )

set.seed(10)

a <- min(d$Y)
b <- max(d$Y)
Ystar <- (d$Y - a)/(b-a)
 keep <- which(prescreen.uni(d$Y, d$A, W, alpha = .05))
 keep.nonbinary <- nonbinary[nonbinary %in% keep]
 if(length(keep.nonbinary) > 0){
 	keep.sq <- keep.nonbinary[prescreen.uni(d$Y, d$A, W[,keep.nonbinary, drop = FALSE]^2, min = 0)]
 	 if (sum(keep.sq) > 0) {
 	 	Wsq <- W[,keep.sq, drop = FALSE]^2
 	 	colnames(Wsq) <- paste0(colnames(Wsq), "sq")
 	 }
 }
 X <- cbind(W[,keep],Wsq)
 n.columns <- ncol(X)
 g.SL <- try(SuperLearner(Y=d$A, X = data.frame(X),
 				SL.library = g.SL.library, family = "binomial", cvControl = list(V = V)) )
if(class(g.SL) == "try-error") {
	 	g1W <- predict(glm(d$A ~ X, family = "binomial"), type = "response")
} else {
		g1W <- g.SL$SL.predict
}
g1W <- .bound(g1W, gbounds)
A0 <- d$A == 0
m.SL.A0 <- SuperLearner(Y=d$Y[A0], X = as.data.frame(X[A0,]), newX = as.data.frame(X), SL.library = SL.library, family = "gaussian", cvControl = list(V = V))
m.SL.A1 <- SuperLearner(Y=d$Y[!A0], X = as.data.frame(X[!A0,]), newX = as.data.frame(X), SL.library = SL.library, family = "gaussian", cvControl = list(V = V))

Q.unbd <- cbind(QAW = d$A * m.SL.A1$SL.predict + (1-d$A) * m.SL.A0$SL.predict,
							Q0W = m.SL.A0$SL.predict, Q1W = m.SL.A1$SL.predict)
colnames(Q.unbd) <- c("QAW", "Q0W", "Q1W")
Q <- .bound((.bound(Q.unbd, c(a,b) )- a) / (b-a), alpha)

results.oneStep <- one_step_att(Y = Ystar, A = d$A, Q = Q, g1W = g1W, depsilon = depsilon, max_iter = max(1000, 2/depsilon), gbounds = gbounds, Qbounds = alpha)

est <- results.oneStep[1]*(b-a)
se <- sqrt(results.oneStep[2])*(b-a)

result <- data.frame(est = est, ci_lower = est - 1.96*se, ci_upper = est + 1.96*se)
write.csv(result, file = outFile, row.names = FALSE)


