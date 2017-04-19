W = rnorm(1000)
A = rbinom(1000,1,plogis(W+.4))
Y = rnorm(1000,.1*A+.2*W+1,2)
Y = (Y-min(Y))/(max(Y)-min(Y))
Q0Wfit = glm(Y~A+W,family='binomial')
gfit = glm(A~W,family='binomial')
Q0W = predict(Q0Wfit,type='response')
g = predict(gfit, type='response')

initdata = data.frame(A=A,Y=Y,Q0W=Q0W,g=g)
# takes initdata, dataframe or list containing elements, Q0W,Y,g,A 
update <- function(initdata) {
  H = with(initdata, (A==0)*g/(mean(A)*(1-g)))
  fit = glm(Y ~ offset(qlogis(Q0W)),
               weights=H,family="binomial")
  Q0Wstar = with(initdata, qlogis(plogis(Q0W-fit$coef)))
  return(Q0Wstar)
}

# basic truncate function
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
  pmin(pmax(x, lower), upper)
}


