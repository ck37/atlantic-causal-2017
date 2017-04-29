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
gendata=function(n){
  
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
n=1000
# function just to give estimates for now
sim_ATT = function(n){
  # draw sample
  data = gendata(n)
  
  # Truth
  Q0Wtrue = data$Q0Wtrue 
  Q1Wtrue = data$Q1Wtrue
  Psi_0 = sum((data$A==1)*(Q1Wtrue-Q0Wtrue))/sum(data$A==1)
  
  # max and mins for scaling, adjust Y
  a = min(data$Y)
  b = max(data$Y)
  data$Y = (data$Y-a)/(b-a)
  data$Q0Wtrue = data$Q1Wtrue = NULL
  
  # covariates including A
  X= data
  X$Y = NULL
  
  # fits
  QAWfit = suppressWarnings(glm(Y~.,data=data,family='binomial'))
  gfit = glm(A~.,data=X, family='binomial')
  # debug(update)
  # predict Q0W and g
  data0 = data
  data0$A = 0
  Q0W = suppressWarnings(predict(QAWfit,newdata=data0,type='response'))
  g = predict(gfit, type='response')
  
  # fluctuate by logistic regression
  initdata = data.frame(A=data$A, Y=data$Y, Q0W=Q0W, g=g)
  Q0Wstar = suppressWarnings(update(initdata))
  Psi = with(initdata,sum((A==1)*(Y-Q0Wstar))/sum(A))
  
  Dstar = with(data, (A==0)*g/(mean(A)*(1-g))*(Y-Q0Wstar))
  SE = (b-a)*sd(Dstar)*sqrt(n-1)/n
  
  # scale the outcome back 
  Psi=(b-a)*Psi
  CI = c(Psi-1.96*SE,Psi+1.96*SE)
  cover = (Psi_0>=CI[1])&(Psi_0<=CI[2])
  # PsiLS=(b-a)*PsiLS
  return(c(Psi_0=Psi_0,Psi=Psi,cover))
}

# run B sims of n=1000 compile results
B=1000
n=1000
res=mclapply(1:B,FUN = function(x) sim_ATT(1000),mc.cores=4)
res=t(sapply(res,FUN=function(x) x))

mean(res[,3])
mean(res[,2])-mean(res[,1])
colnames(res)=c("true","tmle")

# we can see if there's a difference, should be very little 
df=data.frame(type=c(rep("true",B),rep("tmle",B)),est=c(res[,1],res[,2])) 
gghist = ggplot(df,aes(x=est,fill=type))+geom_density(alpha=.3) 
gghist

hist((res[,2]-res[,1]),breaks=100)
