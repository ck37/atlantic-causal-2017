source("lib/weighted_update.R")
source("lib/sl_wrappers.R")
source("lib/bound.R")
source("lib/estimate_att.R")
source("lib/estimate_att_jl.R")
library(ggplot2)
library(parallel)
library(SuperLearner)
library(glmnet)

# bring in covariate data from last year
X_f = read.csv("inbound/pre_data/X_subset_y.csv",header=FALSE)

# to be used with various functional forms
g0=function(v,coeff) plogis(v %*% coeff)
Q0=function(v,coeff) v %*% coeff

whichCont = which(apply(X_f,2,
                        FUN = function(x) length(unique(x))>5&!(is.factor(x))))
whichBin = which(!(1:58 %in% whichCont))

# function to create functional forms to simulate, using data from last year
# allowed forms are "linear","squared" and "trans"
create_siminfo = function(numvarsg=5,numvarsQ=5, formg="linear", formQ="linear"){

  # create linear main terms form if linear is specified
  if (formg=="linear"){
    Wz_names = sample(colnames(X_f),numvarsg)
    form_z = paste0("~",paste(Wz_names,collapse = "+"))
    form_z = formula(form_z)
  }
  
  if (formQ=="linear"){
    Wy_names = sample(colnames(X_f),numvarsQ)
    form_y = paste0("~",paste(c("A",Wy_names),collapse = "+"))
    form_y = formula(form_y)
  }
  
  if (formg=="squared"){
    contCols = colnames(X_f)[sample(whichCont,floor(numvarsg/2))]
    others = colnames(X_f)[sample(whichBin,ceiling(numvarsg/2))]
    formCont = paste(paste0(paste0("I(",contCols),"^2)"),collapse="+")
    formOthers = paste(others,collapse="+")
    form_z = paste0("~",formOthers,"+",formCont)
    form_z = formula(form_z)
    Wz_names=c(contCols,others)
  }
  
  if (formQ=="squared"){
    contCols = colnames(X_f)[sample(whichCont,floor(numvarsQ/2))]
    others = colnames(X_f)[sample(whichBin,ceiling(numvarsQ/2))]
    formCont = paste(paste0(paste0("I(",contCols),"^2)"),collapse="+")
    formOthers = paste0("A*(",paste(others,collapse="+"),")")
    form_y = paste0("~",formOthers,"+",formCont)
    form_y = formula(form_y)
    Wy_names=c(contCols,others)
  }
  
  if (formg=="trans"){
    contCols = colnames(X_f)[sample(whichCont,floor(numvarsg/2))]
    others = colnames(X_f)[sample(whichBin,ceiling(numvarsg/2))]
    formCont = paste(paste0(paste0("I(cos(",contCols),"))"),collapse="+")
    formOthers = paste(others,collapse="+")
    form_z = paste0("~",formOthers,"+",formCont)
    form_z = formula(form_z)
    Wz_names=c(contCols,others)
  }
  
  if (formQ=="trans"){
    contCols = colnames(X_f)[sample(whichCont,floor(numvarsQ/2))]
    others = colnames(X_f)[sample(whichBin,ceiling(numvarsQ/2))]
    formCont = paste(paste0(paste0("I(sin(",contCols),"))"),collapse="+")
    formOthers = paste0("A*(",paste(others,collapse="+"),")")
    form_y = paste0("~",formOthers,"+",formCont)
    form_y = formula(form_y)
    Wy_names=c(contCols,others)
  }

  # the design matrix for treatment mech
  Xz = model.matrix(form_z,X_f[,Wz_names])
  # generate A, making sure to have enough A=1 rep. Keep the betas that work
  A=1
  while (mean(A)>=.66|mean(A)<=.33){
    coeff_z = runif(length(colnames(Xz)),-5/ncol(Xz),5/ncol(Xz))/apply(Xz,2,FUN = function(x) {
      ifelse(var(x)==0,1,max(abs(x)))})
    A = rbinom(250,1,g0(Xz,coeff_z))
  }

  return(list(Xz=Xz,coeff_z=coeff_z,form_y=form_y,Wy_names=Wy_names))
}

# plug siminfo into sim function
sim_ATT_jl = function(siminfo, useSL = F,
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
  coeff_z = siminfo$coeff_z
  Wy_names = siminfo$Wy_names
  form_y = siminfo$form_y
  # generate A, making sure to have enough A=1 rep btwn 33 and 66 percent
  A=99
  while (mean(A)>=.66|mean(A)<=.33){
    A = rbinom(250,1,g0(Xz,coeff_z))
  }
  
  # generate y design and true mean outcomes as well as Y
  Xy = cbind.data.frame(A,X_f[,Wy_names])
  Xy1 = Xy0 = Xy
  Xy[,"A"] = A
  Xy1[,"A"] = 1
  Xy0[,"A"] = 0
  
  # generate designs for outcome under A=1,A=0 and sample
  Xy = model.matrix(form_y,Xy)
  Xy1 = model.matrix(form_y,Xy1)
  Xy0 = model.matrix(form_y,Xy0)
  
  #generate betas for outcome that are tame and then generate true means
  coeff_y = runif(length(colnames(Xy)),-3,3)/apply(Xy,2,FUN = function(x) {
    ifelse(var(x)==0,1,max(abs(x)))})
  QAWtrue = Q0(Xy,coeff_y)
  Q1Wtrue = Q0(Xy1,coeff_y)
  Q0Wtrue = Q0(Xy0,coeff_y)
  # add noise to the true mean
  Y = rnorm(250,QAWtrue,sd(QAWtrue)/3)
  Qtrue = data.frame(Q0 = Q0Wtrue,Q1 = Q1Wtrue)
  
  # Target parameter: sample average treatment effect on treated units (SATT).
  Psi_0 = sum((A == 1) * (Q1Wtrue - Q0Wtrue)) / sum(A == 1)

  sim_results = estimate_att(A = A, 
                             Y = Y, 
                             W = X_f,
                             SL.library = SL.library,
                             g.SL.library = SL.library,
                             useSL = useSL,
                             verbose = verbose,
                             pooled_outcome = T,
                             prescreen=c(.10,10))

  # Indicator for our CI containing the true parameter.
  covered = (Psi_0 >= sim_results$ci_lower) && (Psi_0 <= sim_results$ci_upper)
  
  # Proportion of unit-level CIs containing the true unit-level effect.
  # Coverage pretty bad right now, ~36%. We shouldn't expect to be able to make this very reliable, 
  # but conservative would be better than anti-conservative.
  covered_units = mean((sim_results$unit_estimates$ci_lower <= Q1Wtrue - Q0Wtrue) &
                  (sim_results$unit_estimates$ci_upper >= Q1Wtrue - Q0Wtrue))
  
  # Compile results into a named vector.
  sim_output = c(Psi_0 = Psi_0,
              Psi = sim_results$est,
              covered = covered,
              conf_interval = c(sim_results$ci_lower,sim_results$ci_upper),
              standard_error = sim_results$se,
              covered_units = covered_units)
  
  return(sim_output)
}

# NOTE: We do not want to create siminfo within the sim
SL.library=list(c("SL.glm","prescreen.nosq","All"), c("SL.glmnet","prescreen.nosq"),
                c("SL.nnet","prescreen.nosq","All")) 
siminfo = create_siminfo(numvarsg=5,numvarsQ=9,"linear","trans") 
siminfo

test = sim_ATT_jl(siminfo,useSL=T,SL.library=SL.library)
test

# Set multicore-compatible seed.
set.seed(1, "L'Ecuyer-CMRG")

############
############
### TEST HERE--RUN IT! 
############
############

# NOTE: siminfo takes arguments:
# numvarsg, numvarsQ, formg and formQ
# formg and formQ can each be one of "trans","linear","squared"
# trans offers transcendental forms, linear is main terms, squared
# is interactions and squares of some continuous variables

# Run B sims on the actual covs n = 250 from last year's data

SL.library=list(c("SL.glm","prescreen.nosq","All"), c("SL.glmnet","prescreen.nosq"),
                c("SL.nnet","prescreen.nosq","All")) 

siminfo = create_siminfo(numvarsg=5,numvarsQ=9,formg="linear",formQ="trans")
B = 100

res = mclapply(1:B,
               FUN = function(x) sim_ATT_jl(siminfo,useSL=T,SL.library=SL.library),
               mc.cores = 2)

res = t(sapply(res, FUN = function(x) x))

# Review coverage. We want this to be close to 95%.
mean(res[, 3])

