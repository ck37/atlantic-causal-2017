source("fluctuate.R")
library(ggplot2)
library(parallel)

# generate conditional means
Q0=function(A,W1,W2,W3,W4) return(A+2*A*W4+3*W1+1*W2^2+.5*W3*W4+.25*W4)
Q0=function(A,W1,W2,W3,W4) return(15*A*(W4<-1)+1*A*(W3^2<1)+.63*W1+1*W2^2+A*.5*cos(W3)+.1*W4*W1+A*(W2>1))
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+5*W2^2*W4+.08*W3+5*abs(W4)-1)}

# random draw of sample size n--continuous unscaled Y
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

# function just to give estimates for now
sim_ATT = function(n){
  
  # draw sample
  data = gendata(n)
  
  # Truth
  Q0Wtrue = data$Q0Wtrue 
  Psi_0 = sum((data$A==1)*(data$Y-Q0Wtrue))/sum(data$A==1)
  
  # max and mins for scaling, adjust Y
  a = min(data$Y)
  b = max(data$Y)
  data$Y = (data$Y-a)/(b-a)
  data$Q0Wtrue = NULL
  
  # covariates including A
  X= data
  X$Y = NULL
  
  # fits
  QAWfit = suppressWarnings(glm(Y~.,data=data,family='binomial'))
  gfit = glm(A~.,data=X, family='binomial')
  
  # predict Q0W and g
  X0 = data
  X0$A = 0
  Q0W = suppressWarnings(predict(QAWfit,newdata=X0,type='response'))
  g = predict(gfit, type='response')
  
  # put in data.frame or list
  initdata = data.frame(Q0W=Q0W,Y=data$Y,g=g,A=data$A)
  
  # fluctuate by logistic regression
  Q0Wstar = suppressWarnings(update(initdata))
  Psi = with(initdata,1/sum(A)*sum((A==1)*(Y-Q0Wstar)))
  
  # fluctuate by least squares regression--either fluc is fine for simulations
  # Q0WstarLS = updateLS(initdata)
  # PsiLS = with(initdata,1/sum(A)*sum((A==1)*(Y-Q0WstarLS)))
  
  # scale the outcome back 
  Psi=(b-a)*Psi
  # PsiLS=(b-a)*PsiLS
  return(c(Psi_0=Psi_0,Psi=Psi))
}

# run B sims of n=1000 compile results
B=200
n=1000
res=mclapply(1:B,FUN = function(x) sim_ATT(1000),mc.cores=4)
res=t(sapply(res,FUN=function(x) x))
colnames(res)=c("true","tmle")

# we can see if there's a difference, should be very little 
df=data.frame(type=c(rep("true",B),rep("tmle",B)),est=c(res[,1],res[,2])) 
gghist = ggplot(df,aes(x=est,fill=type))+geom_density(alpha=.3) 
gghist
