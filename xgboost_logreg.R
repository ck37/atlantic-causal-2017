# Generate data
W = rnorm(1000)
A = rbinom(1000,1,plogis(W+.4))
Y = rnorm(1000,.1*A+W^2+10*(W>1),2)
Y = (Y-min(Y))/(max(Y)-min(Y))
# Y = rbinom(1000,1,plogis(.1*A+W^2+10*(W>1)))
m = matrix(c(W,A,Y), ncol=3)
colnames(m)=c("W","A","Y")

# define custom logistic reg objective (log lik loss)
logregobj <- function(preds, dtrain) { 
  labels <- getinfo(dtrain, "label") 
  # preds <- 1/(1+exp(-preds))
  grad <- preds-labels
  hess <- preds * (1 - preds)
  return(list(grad = grad, hess = hess)) }

# eval error can be any good error function to compare models
# here I did MSE of course
evalerror <- function(preds, dtrain) { 
  labels <- getinfo(dtrain, "label") 
  err <- sqrt(mean((preds-labels)^2))
  return(list(metric = "MSE", value = err)) }

# set up their matrix object
dtest <- xgb.DMatrix(m[,1:2], label = m[,"Y"])
# hyperparams, also note two examples below give same answer so not necessary to program it
# can use built-in function "binary:logistic" to override default "reg:linear"
param_logistic <- list(max_depth = 2, eta = .01, silent = 1, objective = "binary:logistic")
param_logistic1 <- list(max_depth = 2, eta = .01, silent = 1, objective = logregobj,eval_error="evalerror")
param_leastsq <- list(max_depth = 2, eta = .01, silent = 1, objective = "reg:linear")

bst_logistic <- xgb.train(param_logistic, data=dtest, nrounds = 1000, maximize  =FALSE,watchlist=list())
bst_logistic1 <- xgb.train(param_logistic1, data=dtest, nrounds = 1000, maximize  =FALSE,watchlist=list())
bst_leastsq <- xgb.train(param_leastsq, data=dtest, nrounds = 1000, maximize  =FALSE,watchlist=list())

# compare yhat with Y
yhat_logistic= predict(bst_logistic, dtest, type='response')
yhat_logistic1= predict(bst_logistic1, dtest,type='response')
yhat_leastsq= predict(bst_leastsq, dtest, type='response')

# results are very similar here but somehow the handmade logistic regression is slightly different
mean(yhat_logistic)
mean(yhat_logistic1)
mean(yhat_leastsq)

hist(yhat_logistic-yhat_logistic1,breaks=100)
hist(yhat_logistic-yhat_leastsq,breaks=100)
hist(yhat_logistic1-yhat_leastsq,breaks=100)

# check if leastsq stays within 0 and 1--of course logistic regressions will do so
test_maxmin = data.frame(true_Yminmax=c(min(Y),max(Y)),leastsq_Yminmax=c(min(yhat_leastsq),max(yhat_leastsq)))
test_maxmin
