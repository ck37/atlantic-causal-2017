Q0=function(A,W1,W2,W3,W4) return(A+2*A*W4+3*W1+1*W2^2+.5*W3*W4+.25*W4)
Q0=function(A,W1,W2,W3,W4) return(A+2*A*W4^2+.63*W1+1*W2^2+A*.5*cos(W3)+.25*W4)
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+.1*W2+.08*W3+1*W4-1)}

gendata=function(n){
  U1 = runif(n,0,1)
  W1= -1*(U1<=.5)+(U1>.5)
  W2=rnorm(n)
  W3=rnorm(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rnorm(n,Q0(A,W1,W2,W3,W4),2)
  Q0Wtrue = Q0(A=rep(0,n),W1,W2,W3,W4)
  data.frame(A,W1,W2,W3,W4,Y,Q0Wtrue)
}

n=1000
data = gendata(n)

# Truth
Q0Wtrue = data$Q0Wtrue 
Psi_0 = sum((data$A==1)*(data$Y-Q0Wtrue))/sum(data$A==1)
Psi_0

a = min(data$Y)
b = max(data$Y)
data$Y = (data$Y-a)/(b-a)
data$Q0Wtrue = NULL

X= data
X$Y = NULL

QAWfit = glm(Y~.,data=data,family='binomial')
gfit = glm(A~.,data=X, family='binomial')

# set up to predict Q0W
X0 = data
X0$A = 0

# predict Q0W and g
Q0W = predict(QAWfit,newdata=X0,type='response')
g = predict(gfit, type='response')

# put in data.frame or list
initdata = data.frame(Q0W=Q0W,Y=data$Y,g=g,A=data$A)

# takes initdata, dataframe or list containing elements, Q0W,Y,g,A 
update <- function(initdata) {
  H = with(initdata, (A==0)*g/(mean(A)*(1-g)))
  # fit a glm with the weight
  fit = glm(Y ~ offset(qlogis(Q0W)),data=initdata,
               weights=H,family="binomial")

  # update
  Q0Wstar = with(initdata, plogis(qlogis(Q0W)+fit$coef))
  return(Q0Wstar)
}


# takes initdata, dataframe or list containing elements, Q0W,Y,g,A 
updateLS <- function(initdata) {
  H = with(initdata, (A==0)*g/(mean(A)*(1-g)))
  # fit a glm with the weight
  fit = glm(Y ~ offset(Q0W),data=initdata,
            weights=H,family="gaussian")
  
  # update yo momma
  Q0Wstar = with(initdata, truncate(Q0W+fit$coef,lower=1e-4))
  return(Q0Wstar)
}

Q0Wstar = update(initdata)
Psi = with(initdata,1/sum(A)*sum((A==1)*(Y-Q0Wstar)))

Q0WstarLS = updateLS(initdata)
PsiLS = with(initdata,1/sum(A)*sum((A==1)*(Y-Q0WstarLS)))

Psi=(b-a)*Psi
PsiLS=(b-a)*PsiLS
c(Psi,PsiLS)

hist(update(initdata)-updateLS(initdata))
# basic truncate function
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
  pmin(pmax(x, lower), upper)
}


