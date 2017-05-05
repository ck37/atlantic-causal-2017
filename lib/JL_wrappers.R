


SL.glmnet_em = function (Y, X, newX, family, obsWeights, id, nfolds = 10,
                        nlambda = 100, useMin = TRUE, alpha, ...)
{
  require("glmnet")

  # create the formula of all ems
  form_em = paste0("A*(",paste(colnames(X),"",collapse="+"),")")

  # capture continuous vars--
  whichCont = apply(X,2,FUN = function(col) length(unique(col))>2)
  contVars = colnames(X)[whichCont]

  # create sq terms on the cont vars
  formCont = paste(paste0(paste0(paste0("I(",contVars),"^2"),")"),"",collapse="+")
  formSq = formula(paste0("~",form_em,"+",formCont))

  X <- model.matrix(formSq, as.data.frame(X))
  newX <- model.matrix(formSq, as.data.frame(newX))
  X = X[,-1]
  newX = newX[,-1]

  # simplify colnames except A to insure no problems with names
  colnames(X)[colnames(X)!="A"]=vapply(1:(ncol(X)-1),FUN=function(x) paste0("X",x),FUN.VALUE="cc")
  colnames(newX)=colnames(X)

  # get rid of obs with either no variance or only 1 obs that is different from 0 or 1
  listp <- apply(X, 2, FUN = function(x) {
    if (var(x)==0|sum(x==0)>=(length(x)-1)|sum(x==1)>=(length(x)-1)) return(0) else return(1)
  })
  keep = (listp==1)

  X = X[,keep]
  newX = newX[,keep]

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

environment(SL.glmnet_em) <- asNamespace("SuperLearner")

create.SL.glmnet_em <- function(alpha) {
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet_em', alpha[mm], '<- function(..., alpha = ', alpha[mm],
                            ') SL.glmnet_em(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
    }
  invisible(TRUE)
}

create.SL.glmnet_em(seq(0,1,.2))
environment(SL.glmnet_em0) <- asNamespace("SuperLearner")
environment(SL.glmnet_em0.2) <- asNamespace("SuperLearner")
environment(SL.glmnet_em0.4) <- asNamespace("SuperLearner")
environment(SL.glmnet_em0.6) <- asNamespace("SuperLearner")
environment(SL.glmnet_em0.8) <- asNamespace("SuperLearner")
environment(SL.glmnet_em1) <- asNamespace("SuperLearner")

# screens out vars with pvals below .25 inc ems and sq terms
SL.glmnet_em25 = function (Y, X, newX, family, obsWeights, id, nfolds = 10,
                         nlambda = 100, useMin = TRUE, alpha, ...)
{
  require("glmnet")

  # create the formula of all ems
  form_em = paste0("A*(",paste(colnames(X),"",collapse="+"),")")

  # capture continuous vars--
  whichCont = apply(X,2,FUN = function(col) length(unique(col))>2)
  contVars = colnames(X)[whichCont]

  # create sq terms on the cont vars
  formCont = paste(paste0(paste0(paste0("I(",contVars),"^2"),")"),"",collapse="+")
  formSq = formula(paste0("~",form_em,"+",formCont))

  X <- model.matrix(formSq, as.data.frame(X))
  newX <- model.matrix(formSq, as.data.frame(newX))
  X = X[,-1]
  newX = newX[,-1]

  # simplify colnames except A to insure no problems with names
  colnames(X)[colnames(X)!="A"]=vapply(1:(ncol(X)-1),FUN=function(x) paste0("X",x),FUN.VALUE="cc")
  colnames(newX)=colnames(X)

  # cut out those with 0 var, keep min number of vars with p-val below .25, keep A
  cutoff = .25
  min = 6
  pvalues = vapply(1:ncol(X),FUN = function(x){
    V = X[,x]
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=1 else {
      m <- glm(Y~ V,family='binomial')
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
create.SL.glmnet_em25 <- function(alpha) {
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet_em25', alpha[mm], '<- function(..., alpha = ', alpha[mm],
                            ') SL.glmnet_em25(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

create.SL.glmnet_em25(seq(0,1,.2))

environment(SL.glmnet_em250) <- asNamespace("SuperLearner")
environment(SL.glmnet_em250.2) <- asNamespace("SuperLearner")
environment(SL.glmnet_em250.4) <- asNamespace("SuperLearner")
environment(SL.glmnet_em250.6) <- asNamespace("SuperLearner")
environment(SL.glmnet_em250.8) <- asNamespace("SuperLearner")
environment(SL.glmnet_em251) <- asNamespace("SuperLearner")

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
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=1 else {
      m <- glm(Y~ V,family='binomial')
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

SL.bayesglm_em05 = function (Y, X, newX, family, obsWeights, ...)
{
  require("arm")
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
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=1 else {
      m <- glm(Y~ V,family='binomial')
      p <- try(summary(m)$coef[2,4], silent = TRUE)
      if (class(p) == "try-error") p=1}
    return(p)},FUN.VALUE = 1)

  keep <- pvalues <= alpha
  if(sum(keep) < min){
    keep[order(pvalues)[1:min]] <- TRUE}

  X = as.data.frame(X[,keep])
  newX = as.data.frame(newX[,keep])

  fit.glm <- arm::bayesglm(Y ~ ., data = X, family = family,
                           weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.bayesglm")
  return(out)
}
environment(SL.bayesglm_em05) <- asNamespace("SuperLearner")


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
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=1 else {
      m <- glm(Y~ V,family='binomial')
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

SL.bayesglm_em20 = function (Y, X, newX, family, obsWeights, ...)
{
  require("arm")
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
    if ((var(V) <= 0)|(colnames(X)[x]=="A")) p=1 else {
      m <- glm(Y~ V,family='binomial')
      p <- try(summary(m)$coef[2,4], silent = TRUE)
      if (class(p) == "try-error") p=1}
    return(p)},FUN.VALUE = 1)

  keep <- pvalues <= alpha
  if(sum(keep) < min){
    keep[order(pvalues)[1:min]] <- TRUE}

  X = as.data.frame(X[,keep])
  newX = as.data.frame(newX[,keep])

  fit.glm <- arm::bayesglm(Y ~ ., data = X, family = family,
                           weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.bayesglm")
  return(out)
}
environment(SL.bayesglm_em20) <- asNamespace("SuperLearner")

SL.library = list("SL.bayesglm_em05","SL.bayesglm_em20","SL.glm_em05","SL.glm_em20",
                  "SL.glmnet_em0","SL.glmnet_em0.2","SL.glmnet_em0.4","SL.glmnet_em0.6",
                  "SL.glmnet_em0.8","SL.glmnet_em1","SL.glmnet_em250","SL.glmnet_em250.2",
                  "SL.glmnet_em250.4","SL.glmnet_em250.6","SL.glmnet_em250.8","SL.glmnet_em251")


# For a quick test use the code below which uses a model matrix as susan did
if (F) {
  # Run manually if desired.
  X = read.csv("inbound/pre_data/X_subset_y.csv",header=FALSE)
  W = read.csv("inbound/pre_data/X_subset_z.csv",header=FALSE)
  A = read.csv("inbound/pre_data/4.100.z.csv")[,1]
  Y = read.csv("inbound/pre_data/4.101.y.csv")[,1]
  Y = (Y-min(Y))/(max(Y)-min(Y))
  X$A = A
  X = as.data.frame(model.matrix(~ ., data = X))
  X = X[,-1]
  X1 = X
  X1$A = 1
  X0 = X
  X0$A = 0
  newX = rbind(X,X1,X0)

  testSL = SuperLearner(Y,X,newX=newX,family=binomial(),SL.library=SL.library)
  testSL$library.predict[1:10,]
}
