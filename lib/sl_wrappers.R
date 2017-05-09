#------- SL wrappers --------

#Don't put both mgcv and gam in the library!
SL.mgcv<-function(Y, X, newX, family, deg.gam = 2, cts.num = 4, ...) 
{
  SuperLearner:::.SL.require("mgcv")
  if ("gam" %in% loadedNamespaces()) 
    warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", 
                                                    colnames(X[, cts.x, drop = FALSE]), ", k=", deg.gam, 
                                                    ")", sep = ""), collapse = "+"), "+", paste(colnames(X[, 
                                                                                                           !cts.x, drop = FALSE]), collapse = "+")))
  }
  else {
    gam.model <- as.formula(paste("Y~", paste(paste("s(", 
                                                    colnames(X[, cts.x, drop = FALSE]), ", k=", deg.gam, 
                                                    ")", sep = ""), collapse = "+")))
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("Y~", paste(colnames(X), 
                                              collapse = "+"), sep = ""))
  }
  fit.gam <- mgcv::gam(gam.model, data = X, family = family)
  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.mgcv")
  return(out)
}


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
  whichVariable <- apply(X, 2, FUN = function(x) {
    if (var(x)==0|sum(x==0)>=(length(x)-1)|sum(x==1)>=(length(x)-1)) return(FALSE) else return(TRUE)
  })
  omit <- grep("sq", colnames(X))
  if(length(omit) > 0){
    whichVariable[omit] <- FALSE
  }
  return(whichVariable)
}


# keep covariates with univariate associations
prescreen.uni <- function(Y, A, X, alpha = .05, min = 5, ...){
  pvalues <- rep(NA, ncol(X))
  for (i in 1:ncol(X)){
    x=X[,i]
    if (var(x)==0|sum(x==0)>=(length(x)-1)|sum(x==1)>=(length(x)-1)) pvalues[i]=1 else {
      m <- lm(Y~ A+ X[,i])
      p <- try(summary(m)$coef[3,4], silent = TRUE)
      if (class(p) == "try-error") {
        pvalues[i] <- 1
      } else {
        pvalues[i] <- p
      }
    }
  }
  keep <- pvalues <= alpha
  if(sum(keep) < min){
    keep[order(pvalues)[1:min]] <- TRUE
  }
  return(keep)
}

screen.corRank4 = function(...) {
  screen.corRank(..., rank = 4)
}

screen.corRank8 = function(...) {
  screen.corRank(..., rank = 8)
}
