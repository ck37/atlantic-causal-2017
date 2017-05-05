source("lib/function_library.R")

# Set auto-install to T for code to install any missing packages.
load_all_packages(auto_install = F, verbose = T)

# Load all .R files in the lib directory.
ck37r::load_all_code("lib", verbose = T)

source("simulations/fluctuate.R")

# generate conditional means
Q0=function(A,W1,W2,W3,W4) return(A+2*A*W4+3*W1+1*W2^2+.5*W3*W4+.25*W4)
Q0=function(A,W1,W2,W3,W4) return(15*A*(W4<-1)+1*A*(W3^2<1)+.63*W1+1*W2^2+A*.5*cos(W3)+.1*W4*W1+A*(W2>1))
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+5*W2^2*W4+.08*W3+5*abs(W4)-1)}
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+5*W2 + W4+.08*W3 -1)}

# random draw of sample size n--continuous unscaled Y
gendata = function(n) {
  U1 = runif(n,0,1)
  W1= -1*(U1<=.5)+(U1>.5)
  W2=rnorm(n)
  W3=rnorm(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rnorm(n,Q0(A,W1,W2,W3,W4),2)
  Q0Wtrue = Q0(A=rep(0,n),W1,W2,W3,W4)
  Q1Wtrue = Q0(A=rep(1,n),W1,W2,W3,W4)
  data.frame(A,W1,W2,W3,W4,Y,Q0Wtrue,Q1Wtrue)
}

# function just to give estimates for now
sim_fluctuate = function(n) {

  # draw sample
  data = gendata(n)

  # Truth
  Q0Wtrue = data$Q0Wtrue
  Q1Wtrue = data$Q1Wtrue

  # Target parameter: sample average treatment on treated units (SATT).
  Psi_0 = sum((data$A == 1) * (Q1Wtrue - Q0Wtrue)) / sum(data$A == 1)

  # Max and min for scaling Y.
  a = min(data$Y)
  b = max(data$Y)

  # Scale Y to [0, 1] interval.
  data$Y = (data$Y - a) / (b - a)

  data$Q0Wtrue = NULL

  # Covariates including A.
  X = data
  X$Y = NULL

  # Estimate smoothed outcome under observed treatment status.
  QAWfit = suppressWarnings(glm(Y ~ ., data = data, family = 'binomial'))

  # Estimate propensity score.
  gfit = glm(A ~ ., data = X, family = 'binomial')

  # Predict Q0W and g
  X0 = data
  X0$A = 0

  # Predict smoothed potential outcome (Q_bar) under A = control.
  Q0W = suppressWarnings(predict(QAWfit, newdata = X0, type = 'response'))

  # Predict propensity score.
  g_hat = predict(gfit, type = 'response')

  # Compile data needed for fluctuation.
  initdata = data.frame(Q0W = Q0W, Y = data$Y, g = g_hat, A = data$A)

  # Fluctuate by logistic regression
  update_results = suppressWarnings(update(initdata))
  Q0Wstar = update_results$Q0Wstar

  Psi = with(initdata, 1 / sum(A) * sum((A == 1) * (Y - Q0Wstar)))

  # fluctuate by least squares regression--either fluc is fine for simulations
  # Q0WstarLS = updateLS(initdata)
  # PsiLS = with(initdata,1/sum(A)*sum((A==1)*(Y-Q0WstarLS)))

  # Return parameter estimate to original scale.
  Psi = (b - a) * Psi
  # PsiLS=(b-a)*PsiLS

  Psi_0 <= Psi
  return(c(Psi_0 = Psi_0, Psi = Psi))
}

# Set multicore-compatible seed.
set.seed(101, "L'Ecuyer-CMRG")

# run B sims of n=1000 compile results
B = 200
n = 1000
res = mclapply(1:B, FUN = function(x) sim_fluctuate(n), mc.cores = 4)

if (F) {
  # Non-parallel version for running manually if desired.
  # This may show additional output compared to the parallel version.
  res = lapply(1:B, FUN = function(x) sim_fluctuate(1000))
}

res = t(sapply(res, FUN = function(x) x))
colnames(res) = c("true", "logistic")

# we can see if there's a difference, should be very little
df = data.frame(type = c(rep("true", B), rep("logistic", B)),
              est = c(res[,1], res[,2]))
bins = 50
e = (max(df$est) - min(df$est)) / bins

gghist = ggplot(df, aes(x = est, fill = type)) +
  geom_density(alpha = .3) + theme_minimal()
gghist
