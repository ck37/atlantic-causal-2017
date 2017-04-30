source("fluctuate.R")
library(ggplot2)
library(parallel)

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
  W1= -1*(U1<=.5)+(U1>.5)
  W2=rnorm(n)
  W3=rnorm(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rnorm(n,Q0(A,W1,W2,W3,W4),2)
  mean(Q0(A,W1,W2,W3,W4))
  mean(Y)
  Q0Wtrue = Q0(A=rep(0,n),W1,W2,W3,W4)
  Q1Wtrue = Q0(A=rep(1,n),W1,W2,W3,W4)
  data.frame(A,W1,W2,W3,W4,Y,Q0Wtrue,Q1Wtrue)
}

# function just to give estimates for now
sim_ATT = function(n) {

  # Draw sample.
  data = gendata(n)

  # Truth
  Q0Wtrue = data$Q0Wtrue
  Q1Wtrue = data$Q1Wtrue

  # Target parameter: sample average treatment effect on treated units (SATT).
  Psi_0 = sum((data$A == 1) * (Q1Wtrue - Q0Wtrue)) / sum(data$A == 1)

  # Max and min for scaling Y.
  a = min(data$Y)
  b = max(data$Y)

  # Scale Y to [0, 1] interval.
  data$Y = (data$Y - a) / (b - a)

  data$Q0Wtrue = data$Q1Wtrue = NULL

  # covariates including A
  X = data
  X$Y = NULL

  # Estimate smoothed outcome under observed treatment status.
  QAWfit = suppressWarnings(glm(Y~.,data=data,family='binomial'))

  # Estimate propensity score.
  gfit = glm(A~.,data=X, family='binomial')

  if (F) {
    # Run this manually to debug into the update() call.
    debug(update)
  }

  # Predict Q0W and g
  data0 = data
  data0$A = 0

  # Predict smoothed potential outcome (Q0_bar) under A = control.
  Q0W = suppressWarnings(predict(QAWfit,newdata=data0,type='response'))

  # Predict propensity score.
  g = predict(gfit, type = 'response')

  # Compile columns needed for fluctuation step.
  initdata = data.frame(A = data$A, Y = data$Y, Q0W = Q0W, g = g)

  # Fluctuate by logistic regression.
  update_results = suppressWarnings(update(initdata))
  Q0Wstar = update_results$Q0Wstar

  Psi = with(initdata, sum((A == 1) * (Y - Q0Wstar)) / sum(A))

  Dstar = with(data, (A == 0) * g / (mean(A) * (1 - g)) * (Y - Q0Wstar))

  # Calculate rescaled standard error with degrees of freedom correction.
  std_err = (b - a) * sd(Dstar) * sqrt(n - 1) / n

  # Return parameter estimate to original scal.
  Psi = (b - a) * Psi

  # Calculate confidence interval.
  CI = Psi + c(-1, 1) * 1.96 * std_err

  # Indicator for our CI containing the true parameter.
  covered = (Psi_0 >= CI[1]) && (Psi_0 <= CI[2])

  # Compile results into a named vector.
  results = c(Psi_0 = Psi_0,
              Psi = Psi,
              covered = covered,
              standard_error = std_err)

  # PsiLS=(b-a)*PsiLS
  return(results)
}

# Set multicore-compatible seed.
set.seed(1, "L'Ecuyer-CMRG")

# Run B sims of n=1000 observations, then compile results.
B = 1000
n = 1000

res = mclapply(1:B, FUN = function(x) sim_ATT(n), mc.cores = 4)
if (F) {
  # Run non-parallel version manually if extra output is useful.
  res = lapply(1:B, FUN = function(x) sim_ATT(n))
}
res = t(sapply(res, FUN = function(x) x))

# Review coverage. We want this to be close to 95%.
# Currently getting 86% so we're undercovering.
mean(res[, 3])

# Check bias in our estimates. We want this to be close to 0.
mean(res[, 2]) - mean(res[, 1])

# Outdated - could be updated.
#colnames(res) = c("true","tmle")

# we can see if there's a difference, should be very little
df=data.frame(type=c(rep("true",B),rep("tmle",B)),est=c(res[,1],res[,2]))
gghist = ggplot(df,aes(x=est,fill=type))+geom_density(alpha=.3)
gghist

hist((res[,2]-res[,1]),breaks=100)
