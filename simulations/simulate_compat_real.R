source("lib/weighted_update.R")
source("lib/sl_wrappers.R")
source("lib/bound.R")
source("lib/estimate_att.R")
library(ggplot2)
library(parallel)
library(SuperLearner)
library(glmnet)

# bring in covariate data from last year
X_f = read.csv("inbound/pre_data/X_subset_y.csv",header=FALSE)

# to be used with various functional forms
g0=function(v,coeff) plogis(v %*% coeff)
Q0=function(v,coeff) v %*% coeff

# function to create functional forms to simulate, using data from last year
# currently works for formg="linear" and formQ="linear"
create_siminfo = function(numvarsg,numvarsQ, formg, formQ){
  formg=formQ="linear"
  numvarsg=10
  numvarsQ=5
  # select the covariates at random
  Wz_names = sample(colnames(X_f),numvarsg)
  Wy_names = sample(colnames(X_f),numvarsQ)
  
  # create the linear form if linear main terms if linear is specified
  if (formg=="linear"){
  form_z = paste0("~",paste(Wz_names,collapse = "+"))
  form_z = formula(form_z)
  }
  
  if (formQ=="linear"){
    form_y = paste0("~",paste(c("A",Wy_names),collapse = "+"))
    form_y = formula(form_y)
  }

  # create treatment design
  Xz = model.matrix(form_z,X_f[,Wz_names])
  
  # generate A, making sure to have enough A=1 rep
  A=1
  while (mean(A)>=.66|mean(A)<=.33){
    coeff_z = runif(length(colnames(Xz)),-5/numvarsg,5/numvarsg)/apply(Xz,2,FUN = function(x) {
      ifelse(var(x)==0,1,sd(x))})
    A = rbinom(250,1,g0(Xz,coeff_z))
  }
  
  # generate y design and true mean outcomes as well as Y
  Xy = cbind(A,X_f[,Wy_names])
  Xy = model.matrix(form_y,Xy)
  
  coeff_y = runif(length(colnames(Xy)),-3,3)/apply(Xy,2,FUN = function(x) {
    ifelse(var(x)==0,1,sd(x))})
  return(list(Xz=Xz,Xy=Xy,coeff_z=coeff_z,coeff_y=coeff_y))
}

# Y
# data.frame(QAW,Y)
# max(QAW)-min(QAW)
# mean(Q1W-Q0W)

# function just to give estimates for now
sim_ATT = function(useSL = F,
                   # Bounds used for Qbar when rescaled to [0, 1].
                   alpha = c(.0005, .9995),
                   # Bounds used for g to bound away from 0, 1.
                   gbounds = c(0.01, 0.99),
                   # SL library for outcome regression (Qbar).
                   SL.library = c("SL.mean", "SL.glmnet",
                                  "SL.glm"#,
                                  #"SL.bartMachine"#,
                                  # "SL.ranger"#,
                                  #"SL.xgboost",
                                  #"SL.speedglm"
                   ),
                   verbose = T) {
  
  
  # pull the sim info generated before--this should not change within this sim function
  Xz = siminfo$Xz
  Xy = siminfo$Xy
  coeff_z = siminfo$coeff_z
  coeff_y = siminfo$coeff_y

  # generate A, making sure to have enough A=1 rep
  A=99
  while (mean(A)>=.66|mean(A)<=.33){
    A = rbinom(250,1,g0(Xz,coeff_z))
  }
  A
  mean(A)
  # generate y design and true mean outcomes as well as Y
  Xy1 = Xy0 = Xy
  Xy1[,"A"] = 1
  Xy0[,"A"] = 0
  
  QAWtrue = Q0(Xy,coeff_y)
  Q1Wtrue = Q0(Xy1,coeff_y)
  Q0Wtrue = Q0(Xy0,coeff_y)
  Y = rnorm(250,QAW,sd(QAW)/3)
  
  
  # Target parameter: sample average treatment effect on treated units (SATT).
  Psi_0 = sum((A == 1) * (Q1Wtrue - Q0Wtrue)) / sum(A == 1)
  
  sim_results = estimate_att(A = A, 
                             Y = Y, 
                             W = X_f,
                             SL.library = SL.library,
                             g.SL.library = SL.library,
                             useSL = useSL,
                             verbose = verbose)
  
  # Indicator for our CI containing the true parameter.
  covered = (Psi_0 >= sim_results$ci_lower) && (Psi_0 <= sim_results$ci_upper)
  
  # Compile results into a named vector.
  sim_output = c(Psi_0 = Psi_0,
              Psi = sim_results$est,
              covered = covered,
              conf_interval = c(sim_results$ci_lower,sim_results$ci_upper),
              standard_error = sim_results$se)
  
  return(sim_output)
}

# TODO: this is missing a call to create_siminfo().
sim_ATT()
# Set multicore-compatible seed.
# set.seed(1, "L'Ecuyer-CMRG")

# Run B sims of n=1000 observations, then compile results.
B = 10

# Takes only a few seconds without using SL.
# With SL (useSL = T) will take 10 minutes or more.
res = mclapply(1:B, FUN = function(x) sim_ATT(useSL = F), mc.cores = 2)

if (F) {
  # Run non-parallel version manually if extra output is useful.
  res = lapply(1:B, FUN = function(x) sim_ATT(useSL = F))
}

if (F) {
  # Run this if desired for debugging.
  debugonce(sim_ATT)
  
  # Run version with SL fit for Q.
  # This takes a very long time to run 1000 times, like 30 minutes to 1+ hrs
  # depending on the SL library.
  res = mclapply(1:B, FUN = function(x) sim_ATT(n, useSL = T), mc.cores = 4)
  
  # Single-run version for testing:
  res = lapply(1:50, FUN = function(x) sim_ATT(n, useSL = T))
}

res = t(sapply(res, FUN = function(x) x))

# Review coverage. We want this to be close to 95%.
mean(res[, 3])

# Check bias in our estimates. We want this to be close to 0.
mean(res[, 2]) - mean(res[, 1])
hist((res[,2]-res[,1]))
