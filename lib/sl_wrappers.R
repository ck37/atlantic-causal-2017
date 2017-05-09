#------- SL wrappers --------

# Don't put both mgcv and gam in the library!
SL.mgcv<-function(Y, X, newX, family, deg.gam = 2, cts.num = 4, ...) {
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


# Not being used, should be moved into archive.
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

# Not being used, should be moved into archive.
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


# screens out vars with pvals below .15 inc ems 
SL.glmnet_em15 = function (Y, X, newX, family, obsWeights, id, nfolds = 10,
                           nlambda = 100, useMin = TRUE, alpha, ...)
{
  require("glmnet")
  
  # create the formula of all ems
  form_em = paste0("~A*(",paste(colnames(X),collapse="+"),")")
  form_em = formula(form_em)
  X <- model.matrix(form_em, as.data.frame(X))
  newX <- model.matrix(form_em, as.data.frame(newX))
  X = X[,-1]
  newX = newX[,-1]
  
  # simplify colnames except A to insure no problems with names
  colnames(X)[colnames(X)!="A"]=vapply(1:(ncol(X)-1),FUN=function(x) paste0("X",x),FUN.VALUE="cc")
  colnames(newX)=colnames(X)
  
  # cut out those with 0 var, keep min number of vars with p-val below .25, keep A
  cutoff = .15
  min = 6
  pvalues = vapply(1:ncol(X),FUN = function(x){
    x=3
    V = X[,x]
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=0 else {
      m <- glm(Y~ V,family='gaussian')
      p <- try(summary(m)$coef[2,4], silent = TRUE)
      if (class(p) == "try-error") p=1}
    return(p)},FUN.VALUE = 1)
  
  keep <- pvalues <= cutoff
  if(sum(keep) < min){
    keep[order(pvalues)[1:min]] <- TRUE}
  
  X = X[,keep]
  newX = newX[,keep]
  
  # proceed as in normal glmnet wrapper
  fitCV <- glmnet::cv.glmnet(x = X, y = Y,
                             lambda = NULL, type.measure = "deviance", nfolds = nfolds,
                             family = 'gaussian', alpha = alpha, nlambda = nlambda)
  pred <- predict(fitCV$glmnet.fit, newx = newX,
                  s = ifelse(useMin,fitCV$lambda.min,fitCV$lambda.1se),
                  type = "response")
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

# function to create glmnets with various alpha values (ridge-lasso balance)
create.SL.glmnet_em15 <- function(alpha) {
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet_em15_', alpha[mm], '<- function(..., alpha = ', alpha[mm],
                            ') SL.glmnet_em15(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

# TODO: convert to create.Learner()
create.SL.glmnet_em15(c(0,.5,1))

environment(SL.glmnet_em15_0) <- asNamespace("SuperLearner")
environment(SL.glmnet_em15_0.5) <- asNamespace("SuperLearner")
environment(SL.glmnet_em15_1) <- asNamespace("SuperLearner")
sl_glmnet_em15 = list(names = c("SL.glmnet_em15_0", "SL.glmnet_em15_0.5", "SL.glmnet_em15_1"))

SL.glm_em05 = function (Y, X, newX, family, obsWeights, ...)
{
  
  form_em = paste0("A*(",paste(colnames(X),"",collapse="+"),")")
  
  # create sq terms on the cont vars
  form_em = formula(paste0("~",form_em))
  
  X <- model.matrix(form_em, as.data.frame(X))
  newX <- model.matrix(form_em, as.data.frame(newX))
  X = X[,-1]
  newX = newX[,-1]
  
  # simplify colnames except A to insure no problems with names
  colnames(X)[colnames(X)!="A"]=vapply(1:(ncol(X)-1),FUN=function(x) paste0("X",x),FUN.VALUE="cc")
  colnames(newX)=colnames(X)
  
  # cut out those with 0 var, keep min number of vars (inc A) with p-val below .25, keep A
  alpha = .05
  min = 6
  pvalues = vapply(1:ncol(X),FUN = function(x){
    V = X[,x]
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=0 else {
      m <- glm(Y~ V,family=family)
      p <- try(summary(m)$coef[2,4], silent = TRUE)
      if (class(p) == "try-error") p=1}
    return(p)},FUN.VALUE = 1)
  
  keep <- pvalues <= alpha
  if(sum(keep) < min){
    keep[order(pvalues)[1:min]] <- TRUE}
  
  X = as.data.frame(X[,keep])
  newX = as.data.frame(newX[,keep])
  
  fit.glm <- glm(Y ~ ., data = as.data.frame(X), family = family, weights = obsWeights)
  pred <- predict(fit.glm, newdata = as.data.frame(newX), type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}
environment(SL.glm_em05) <- asNamespace("SuperLearner")

SL.glm_em20 = function (Y, X, newX, family, obsWeights, ...)
{
  
  form_em = paste0("A*(",paste(colnames(X),"",collapse="+"),")")
  
  # create sq terms on the cont vars
  form_em = formula(paste0("~",form_em))
  
  X <- model.matrix(form_em, as.data.frame(X))
  newX <- model.matrix(form_em, as.data.frame(newX))
  X = X[,-1]
  newX = newX[,-1]
  
  # simplify colnames except A to insure no problems with names
  colnames(X)[colnames(X)!="A"]=vapply(1:(ncol(X)-1),FUN=function(x) paste0("X",x),FUN.VALUE="cc")
  colnames(newX)=colnames(X)
  
  # cut out those with 0 var, keep min number of vars (inc A) with p-val below .25, keep A
  alpha = .20
  min = 6
  pvalues = vapply(1:ncol(X),FUN = function(x){
    V = X[,x]
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=0 else {
      m <- glm(Y~ V,family=family)
      p <- try(summary(m)$coef[2,4], silent = TRUE)
      if (class(p) == "try-error") p=1}
    return(p)},FUN.VALUE = 1)
  
  keep <- pvalues <= alpha
  if(sum(keep) < min){
    keep[order(pvalues)[1:min]] <- TRUE}
  
  X = as.data.frame(X[,keep])
  newX = as.data.frame(newX[,keep])
  
  fit.glm <- glm(Y ~ ., data = as.data.frame(X), family = family, weights = obsWeights)
  pred <- predict(fit.glm, newdata = as.data.frame(newX), type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}
environment(SL.glm_em20) <- asNamespace("SuperLearner")

