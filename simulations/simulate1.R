source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F, verbose = T)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = T)

source("simulations/fluctuate1.R")

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

  # Max and min for scaling Y.
  a = min(data$Y)
  b = max(data$Y)

  # Scale Y to [0, 1] interval. Call this rescaled Y "Ystar" per Susan Gruber
  # code, and keep the unscaled Y to allow either to be modeled.
  Ystar = (data$Y - a) / (b - a)

  # covariates including A
  X = data

  # Remove Y so it isn't included in any regressions.
  X$Y = NULL

  # Estimate smoothed outcome under observed treatment status.
  if (!useSL) {
    # Use Ystar with family = binomial and Y with family = gaussian.
    # Feel free to switch back to the Ystar version.
    #QAWfit = suppressWarnings(glm(Ystar ~ ., data = data, family = 'binomial'))
    QAWfit = suppressWarnings(glm(Y ~ ., data = data, family = 'gaussian'))
  } else {
    QAWfit = SuperLearner(Y = data$Y,
                          X = X,
                          SL.library = SL.library,
                          family = 'gaussian')
  }

  # Estimate propensity score.
  if (!useSL) {
    gfit = glm(A ~ ., data = X, family = 'binomial')
  } else {
    gfit = SuperLearner(X = subset(X, select = -A), Y = X$A, family = 'binomial',
                        SL.library = SL.library)
  }

  if (F) {
    # Run this manually to debug into the update() call.
    debug(update)
  }

  # Predict Q0W and g with all units set to A = control.
  data0 = X
  data0$A = 0

  # Predict smoothed potential outcome (Q0_bar) under A = control.
  if (!useSL) {
    Q0W = suppressWarnings(predict(QAWfit, newdata = data0, type = 'response'))
  } else {
    cat("QAWfit from SL:\n")
    print(QAWfit)
    Q0W = predict(QAWfit, newdata = data0, onlySL = T)$pred

    num_nas = sum(is.na(Q0W))
    if (num_nas > 0) {
      cat("Found", num_nas, "missing values in Q0W after SL.\n")
    }
  }

  # Predict Q1W with all units set to A = treated.
  data1 = X
  data1$A = 1

  if (!useSL) {
    Q1W = suppressWarnings(predict(QAWfit, newdata = data1, type = 'response'))
  } else {
    cat("QAWfit from SL:\n")
    print(QAWfit)
    Q1W = predict(QAWfit, newdata = data1, onlySL = T)$pred

    num_nas = sum(is.na(Q1W))
    if (num_nas > 0) {
      cat("Found", num_nas, "missing values in Q1W after SL.\n")
    }
  }

  # NOTE: change this bounding below if Q0W is modeled on Ystar rather than Y.

  # Make sure to bound prediction, as it can go out of bounds esp.
  # with linear terms.
  # Here we bound to the original scale.
  Q0W_orig_scale = .bound(Q0W, c(a, b))
  Q1W_orig_scale = .bound(Q1W, c(a, b))

  # Rescale to [0, 1] and bound away from 0, 1 by alpha.
  Q0W = .bound((Q0W_orig_scale - a) / (b - a), alpha)
  Q1W = .bound((Q1W_orig_scale - a) / (b - a), alpha)

  # TODO: confirm Q0W is within [0, 1].

  # Predict propensity score.
  if (!useSL) {
    g = predict(gfit, type = 'response')
  } else {
    g = gfit$SL.predict
  }

  # Bound g away from 0, 1.
  g = .bound(g, gbounds)

  # Compile columns needed for fluctuation step. Here we use Ystar (rescaled).
  initdata = data.frame(A = data$A, Y = Ystar, Q0W = Q0W, Q1W = Q1W, g = g)

  # Fluctuate by logistic regression.
  update_results = suppressWarnings(update(initdata))
  Q0Wstar = update_results$Q0Wstar
  Q1Wstar = update_results$Q1Wstar
  QAWstar = with(initdata, ifelse(A == 1, Q1Wstar, Q0Wstar))

  Psi = with(initdata, sum((A == 1) * (Q1Wstar - Q0Wstar)) / sum(A))

  Dstar = with(data, ((A == 1) - (A == 0) * g / (1 - g)) / mean(A) * (Ystar - QAWstar))

  num_nas = sum(is.na(Dstar))
  if (num_nas > 0) {
    cat("Dstar num NAs is:", num_nas, "\n")
  }

  # Calculate rescaled standard error and cancel sd()'s degrees of freedom
  # correction.
  std_err = (b - a) * sd(Dstar) * sqrt(n - 1) / n

  if (is.na(std_err)) {
    cat("Std_err is na!\n")
  }

  # Return parameter estimate to original scale.
  Psi = (b - a) * Psi

  # Calculate confidence interval.
  CI = Psi + c(-1, 1) * qnorm(.975) * std_err

  # Indicator for our CI containing the true parameter.
  covered = (Psi_0 >= CI[1]) && (Psi_0 <= CI[2])

  # Compile results into a named vector.
  results = c(Psi_0 = Psi_0,
              Psi = Psi,
              covered = covered,
              conf_interval = CI,
              standard_error = std_err)

  # PsiLS=(b-a)*PsiLS
  return(results)
}

# Set multicore-compatible seed.
# set.seed(1, "L'Ecuyer-CMRG")

# Run B sims of n=1000 observations, then compile results.
B = 1000
n = 1000

# Takes only a few seconds.
res = mclapply(1:B, FUN = function(x) sim_ATT(n, useSL = FALSE), mc.cores = 4)

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

hist((res[,2]-res[,1]))
