# Generate data
Q0=function(A,W1,W2,W3,W4) return(A+2*A*W4+3*W1+1*W2^2+.5*W3*W4+.25*W4)
Q0=function(A,W1,W2,W3,W4) return(A+2*A*W4+.63*W1+1*W2^3+A*.5*cos(W3)+.25*W4)
g0=function(W1,W2,W3,W4) {plogis(-.28*W1+1*W2+.08*W3-.12*W4-1)}

gendata=function(n){
  U1 = runif(n,0,1)
  W1= -1*(U1<=.5)+(U1>.5)
  W2=rnorm(n)
  W3=rnorm(n)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rnorm(n,Q0(A,W1,W2,W3,W4),2)
  data.frame(A,W1,W2,W3,W4,Y)
}
big = gendata(1000000)
mean((big$Y-min(big$Y))/(max(big$Y)-min(big$Y)))

# define custom logistic reg objective (log lik loss)
logregobj <- function(preds, dtrain) { 
  labels <- getinfo(dtrain, "label") 
  
  # interestingly, the model constraint by the logistic function not present here
  # if I don't transform and then try log-lik loss it is less stable due to the
  # xgboost starting point, perhaps, gives similar results to logistic
  
  # grad <- -labels/preds-(labels-1)/(1-preds)
  # hess <- labels/preds^2-(labels-1)/(1-preds)^2
  
  preds <- 1/(1+exp(-preds))
  grad <- preds - labels
  hess <- preds * (1 - preds)
  return(list(grad = grad, hess = hess)) }

# eval error can be any good error function to compare models
# here I did MSE of course
evalerror <- function(preds, dtrain) { 
  labels <- getinfo(dtrain, "label") 
  err <- sqrt(mean((preds-labels)^2))
  return(list(metric = "MSE", value = err)) }

# set up their matrix object, notice scaled continuous Y
m=as.matrix(gendata(1000))
a = min(m[,"Y"])
b = max(m[,"Y"])
dtest <- xgb.DMatrix(m[,1:5], label = (m[,"Y"]-a)/(b-a))
# hyperparams, also note two examples below give same answer so not necessary to program it
# can use built-in function "binary:logistic" to override default "reg:linear"
param_logistic <- list(max_depth = 2, eta = .01, silent = 1, objective = "binary:logistic")
param_logistic1 <- list(max_depth = 2, eta = .01, silent = 1, objective = logregobj)
param_leastsqLinear <- list(max_depth = 2, eta = .01, silent = 1, objective = "reg:linear",booster="gblinear")
param_leastsq<- list(max_depth = 2, eta = .01, silent = 1, objective = "reg:linear")

bst_logistic <- xgb.train(param_logistic, data=dtest, nrounds = 1000, maximize  =FALSE,watchlist=list())
bst_logistic1 <- xgb.train(param_logistic1, data=dtest, nrounds = 1000, maximize  =FALSE,watchlist=list())
bst_leastsqLinear <- xgb.train(param_leastsqLinear, data=dtest, nrounds = 1000, maximize  =FALSE,watchlist=list())
bst_leastsq <- xgb.train(param_leastsq, data=dtest, nrounds = 1000, maximize  =FALSE,watchlist=list())

# compare yhat with Y
yhat_logistic= predict(bst_logistic, dtest)
yhat_logistic1= plogis(predict(bst_logistic1, dtest))
# yhat_logistic1= predict(bst_logistic1, dtest)
yhat_leastsqLinear= predict(bst_leastsqLinear, dtest)
yhat_leastsq= predict(bst_leastsq, dtest)

# results are very similar here but somehow the handmade logistic regression is slightly different

# built-in logistic
mean(yhat_logistic)

# hard coded
mean(yhat_logistic1)

mean(yhat_leastsqLinear)
mean(yhat_leastsq)

# very close to each other--hard coded vs built-in logistic
hist(yhat_logistic-yhat_logistic1,breaks=100)

hist(yhat_logistic-yhat_leastsq,breaks=100)
hist(yhat_logistic1-yhat_leastsq,breaks=100)
hist(yhat_leastsqLinear-yhat_leastsq,breaks=100)

# anything with trees respects param bounds but linear regression does not
test_maxmin = data.frame(sampleY=c(a,b),logistic=c(min((b-a)*yhat_logistic)+a,(b-a)*max(yhat_logistic)+a),
                         logistic1=c(min((b-a)*yhat_logistic1)+a,(b-a)*max(yhat_logistic1)+a),
                         leastsq=c((b-a)*min(yhat_leastsq)+a,(b-a)*max(yhat_leastsq)+a),
                         leastsqLinear=c(min((b-a)*yhat_leastsqLinear)+a,(b-a)*max(yhat_leastsqLinear)+a))
rownames(test_maxmin)=c("min","max")
t(test_maxmin)
