#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 2)
  stop("usage: Gruber_vDLv2 inFile outFile")

inFile <- args[[1]]
outFile <- args[[2]]

# Susan Gruber and Mark van der Laan
# sgruber@hsph.harvard.edu, laan@berkeley.edu
# June 5, 2016 - Add Bart to SL library for g
# May 11, 2016
# Entry for ACIC 2016 Competition, Option1

#This program has been tested using older versions of some R packages
# We've included these in the submission zip file.
# SuperLearner 2.0-19
# nnls 1.4
# gam 1.09
# glmnet 1.9-5
# gbm 2.1
# earth 3.2-7
# plotmo 1.3-3
# plotrix 3.5-3
# Matrix (>= 1.0-6) 1.1-2
 # survival 2.37-7
# lattice  0.20-27


if (!file.exists(inFile)) stop("cannot find '", inFile, "'")
d <- read.csv(inFile)
colnames(d)[1:2] <- c("A", "Y")

if (!require(SuperLearner)) stop("Cannot find package SuperLearner")
if (!require(gam)) stop("Cannot find package gam")
if (!require(glmnet)) stop("Cannot find package glmnet")
if (!require(earth)) stop("Cannot find package earth")
#if (!require(BayesTree)) stop("Cannot find package BayesTree")
if (!require(dbarts)) stop("Cannot find package dbarts")

#-----------------------------------One-Step TMLE for ATT parameter  ----------------------------------------
oneStepATT <- function(Y, A, Q, g1W, depsilon, max_iter, gbounds, Qbounds){
n <- length(Y)
q <- mean(A)
calcLoss <- function(Q, g1W){
		-mean(Y * log(Q[,"QAW"]) + (1-Y) * log(1 - Q[,"QAW"]) + A * log(g1W) + (1-A) * log(1 - g1W))
}
psi.prev <- psi  <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q) 
H1.AW =  (A -(1-A) * g1W / (1-g1W))/q
IC.prev <- IC.cur <- H1.AW* (Y-Q[, "QAW"]) + A*(Q[,"Q1W"]-Q[,"Q0W"] - psi)/q 
deriv <-  mean(H1.AW* (Y-Q[, "QAW"]) + A*(Q[,"Q1W"]-Q[,"Q0W"] - psi)/q )
if (deriv > 0) { depsilon <- -depsilon}
loss.prev <- Inf
 	loss.cur <-  calcLoss(Q, g1W)
 	if(is.nan(loss.cur) | is.na(loss.cur) | is.infinite(loss.cur)) { 
 		loss.cur <- Inf
 		loss.prev <- 0
 	}
 	iter <-  0
 	while (loss.prev > loss.cur & iter < max_iter){
	IC.prev <-  IC.cur		
	Q.prev <- Q
	g1W.prev <- g1W	
	g1W <- .bound(plogis(qlogis(g1W.prev) - depsilon  * (Q.prev[,"Q1W"] - Q.prev[,"Q0W"] - psi.prev)/q), gbounds) 
 		H1 <- cbind(HAW = A/q - (1-A) * g1W.prev / (q * (1-g1W.prev)),
				  H0W =   - g1W.prev/(q * (1-g1W.prev)),
				  H1W =   1/q) 
 		Q  <- .bound(plogis(qlogis(Q.prev) - depsilon * H1), Qbounds)
 		psi.prev <- psi
 		psi <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q) 
 		loss.prev <- loss.cur
 		loss.cur <- calcLoss(Q, g1W)
 		IC.cur <- ((A - (1-A) * g1W / (1-g1W)) * (Y-Q[, "QAW"]) + A *(Q[,"Q1W"]-Q[,"Q0W"] - psi))/q
 		if(is.nan(loss.cur) | is.infinite(loss.cur) | is.na(loss.cur)) {loss.cur <- Inf}
 		iter <- iter + 1
 		#print(psi.prev)
 	}
 	 	return(c(psi = psi.prev, var.psi = var(IC.prev)/n,  conv = loss.prev < loss.cur))
 }

#------------- bound function --------------
.bound <- function(x, bds){
x[x > max(bds)] <- max(bds)
x[x < min(bds)] <- min(bds)
x
}

#------- SL wrappers --------
SL.bart <- function(Y, X, newX, family, printevery = 100000, keeptrainfits = FALSE, ...) {
    SuperLearner:::.SL.require("dbarts")	
    fit.bart <- bart(x.train = X, y.train = Y, x.test = newX, printevery = printevery, keeptrainfits = keeptrainfits, 
    		printcutoffs = 0, verbose = FALSE)
    if (family == "gaussian"){
    		pred <- fit.bart$yhat.test.mean
    } else {
   		 pred <- pnorm(colMeans(fit.bart$yhat.test))
   	}
     fit <- list(object = fit.bart)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.dbarts")
    return(out)
}

sg.gbm.2500 <- function (Y, X, newX, family, obsWeights, gbm.trees = 2500, 
    interaction.depth = 2, ...) {
    SuperLearner:::.SL.require("gbm")
    gbm.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+")))
    if (family$family == "gaussian") {
        fit.gbm <- gbm(formula = gbm.model, data = X, distribution = "gaussian", 
            n.trees = gbm.trees, interaction.depth = interaction.depth, 
            cv.folds = 5, keep.data = TRUE, weights = obsWeights, 
            verbose = FALSE)
    }
    if (family$family == "binomial") {
        fit.gbm <- gbm(formula = gbm.model, data = X, distribution = "bernoulli", 
            n.trees = gbm.trees, interaction.depth = interaction.depth, 
            cv.folds = 5, keep.data = TRUE, verbose = FALSE, class.stratify.cv = TRUE,
            weights = obsWeights)
    }
    best.iter <- gbm.perf(fit.gbm, method = "cv", plot.it = FALSE)
    pred <- predict(fit.gbm, newdata = newX, best.iter, type = "response")
    fit <- list(object = fit.gbm, n.trees = best.iter)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.gbm")
    return(out)
}

prescreen.nocat <- function(Y, X, ...){
 	whichVariable <- rep(TRUE, NCOL(X))
 	 omit <- grep("cat", colnames(X))
if(length(omit) > 0){
		whichVariable[omit] <- FALSE
}
return(whichVariable)
}


prescreen.nosq <- function(Y, X, ...){
 	whichVariable <- rep(TRUE, NCOL(X))
 	 omit <- grep("sq", colnames(X))
if(length(omit) > 0){
		whichVariable[omit] <- FALSE
}
return(whichVariable)
}

SL.library <- list(c("SL.glm", "All",  "prescreen.nosq"), c("SL.gam", "All", "prescreen.nosq"),"SL.glmnet", c("sg.gbm.2500", "prescreen.nocat"), c("SL.earth", "prescreen.nosq"), c("SL.bart", "prescreen.nocat"))

g.SL.library <- list(c("SL.glm", "All", "prescreen.nosq"), c("sg.gbm.2500", "prescreen.nocat"), c("SL.gam", "All", "prescreen.nosq"), 		c("SL.earth", "prescreen.nosq"), c("SL.bart", "prescreen.nocat"))

# keep covariates with univariate associations
prescreen.uni <- function(Y, A, X, alpha = .05, min = 5, ...){
pvalues <- rep(NA, ncol(X))
for (i in 1:ncol(X)){
 		m <- lm(Y~ A+ X[,i])
 		p <- try(summary(m)$coef[3,4], silent = TRUE)
 		if (class(p) == "try-error") {
 			pvalues[i] <- 1
 		} else {
 			pvalues[i] <- p
 		} 			
 	}
 	keep <- pvalues <= alpha
 	if(sum(keep) < min){
 		keep[order(pvalues)[1:min]] <- TRUE
 	}
 	return(keep)
 }

 
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

results.oneStep <- oneStepATT(Y = Ystar, A = d$A, Q = Q, g1W = g1W, depsilon = depsilon, max_iter = max(1000, 2/depsilon), gbounds = gbounds, Qbounds = alpha)
	
est <- results.oneStep[1]*(b-a)
se <- sqrt(results.oneStep[2])*(b-a)
result <- data.frame(est = est, ci_lower = est - 1.96*se, ci_upper = est + 1.96*se)
write.csv(result, file = outFile, row.names = FALSE)


