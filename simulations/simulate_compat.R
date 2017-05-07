source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F, verbose = T)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = T)

# generate conditional means

Q0=function(A,W1,W2,W3,W4) return(A+2*A*W4+3*W1+1*W2^2+.5*W3*W4+.25*W4)
Q0=function(A,W1,W2,W3,W4) return(15*A*(W4<-1)+1*A*(W3^2<1)+.63*W1+1*W2^2+A*.5*cos(W3)+
                                    .1*W4*W1+A*(W2>1))
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+5*W2^2*W4+.08*W3+5*abs(W4)-1)}

# well-spec
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+1*W2+.08*W3-.3*W4-1)}
Q0=function(A,W1,W2,W3,W4) return(A-2*W3+3*W1+1*W2+.25*W4)

# random draw of sample size n--continuous unscaled Y
gendata = function(n) {

  U1 = runif(n,0,1)
  W1 = -1*(U1<=.5) + (U1>.5)
  W2 = rnorm(n)
  W3 = rnorm(n)
  W4 = rnorm(n)
  A = rbinom(n,1,g0(W1,W2,W3,W4))
  Y = rnorm(n,Q0(A,W1,W2,W3,W4),2)
  mean(Q0(A,W1,W2,W3,W4))
  mean(Y)
  Q0Wtrue = Q0(A = rep(0, n),W1,W2,W3,W4)
  Q1Wtrue = Q0(A = rep(1, n),W1,W2,W3,W4)
  data.frame(A, W1, W2, W3, W4, Y, Q0Wtrue, Q1Wtrue)
}

# function just to give estimates for now
sim_ATT = function(n,
                   # If F use glm(), if T use an SL library for Q.
                   useSL = F,
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
                   )) {

  # Draw sample from data-generating process.
  data = gendata(n)

  # Pull true potential outcomes out of the dataframe.
  Q0Wtrue = data$Q0Wtrue
  Q1Wtrue = data$Q1Wtrue

  # Now remove these true values from the dataframe so they aren't used
  # in the glm() when we do ~ .
  data$Q0Wtrue = data$Q1Wtrue = NULL

  # Target parameter: sample average treatment effect on treated units (SATT).
  Psi_0 = sum((data$A == 1) * (Q1Wtrue - Q0Wtrue)) / sum(data$A == 1)

  sim_results = estimate_att(A = data$A,
                             Y = data$Y,
                             W = as.data.frame(cbind(data$W1,data$W2,data$W3,data$W4)),
                             SL.library = SL.library,
                             g.SL.library = SL.library,
                             useSL = useSL)

  # Indicator for our CI containing the true parameter.
  covered = (Psi_0 >= sim_results$ci_lower) && (Psi_0 <= sim_results$ci_upper)
  
  # Evaluate unit-level effect estimates
  pehe          = sd(sim_results$unit_estimates$est - Q1Wtrue + Q0Wtrue)
  covered_units = mean((sim_results$unit_estimates$ci_lower <= Q1Wtrue - Q0Wtrue) & 
                         (sim_results$unit_estimates$ci_upper >= Q1Wtrue - Q0Wtrue))
  
  # Compile results into a named vector.
  sim_output = c(Psi_0 = Psi_0,
              Psi = sim_results$est,
              covered = covered,
              conf_interval = c(sim_results$ci_lower,sim_results$ci_upper),
              standard_error = sim_results$se,
              pehe = pehe,
              covered_units = covered_units)

  return(sim_output)
}

# Set multicore-compatible seed.
# set.seed(1, "L'Ecuyer-CMRG")

# Run B sims of n=1000 observations, then compile results.
B = 1000
n = 1000

# Takes only a few seconds without using SL.
# With SL (useSL = T) will take 10 minutes or more.
res = mclapply(1:100, FUN = function(x) sim_ATT(n, useSL = T,SL.library="SL.glm"), mc.cores = 4)

if (F) {
  # Run non-parallel version manually if extra output is useful.
  res = lapply(1:B, FUN = function(x) sim_ATT(n))
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
# Currently getting 94.5% so we're looking pretty good!
# But 93.6 - 93.3% with SuperLearner :/
mean(res[, 3])

# Check bias in our estimates. We want this to be close to 0.
mean(res[, 2]) - mean(res[, 1])
hist((res[,2]-res[,1]))
