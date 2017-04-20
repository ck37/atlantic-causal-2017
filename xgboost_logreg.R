# Generate data
W = rnorm(1000)
A = rbinom(1000,1,plogis(W+.4))
Y = rnorm(1000,.1*A+.2*W+1,2)
Y = (Y-min(Y))/(max(Y)-min(Y))
m = matrix(c(W,A,Y), ncol=3)
colnames(m)=c("W","A","Y")

# define custom logistic reg objective(log lik loss)
logregobj <- function(preds, dtrain) { 
  labels <- getinfo(dtrain, "label") 
  grad <- preds-labels
  hess <- preds * (1 - preds)
  return(list(grad = grad, hess = hess)) }

# eval error can be any good error function to compare models
# here I did MSE of course
evalerror <- function(preds, dtrain) { 
  labels <- getinfo(dtrain, "label") 
  err <- sqrt(mean((preds-labels)^2))
  return(list(metric = "MSE", value = err)) }
xgboost(m,label = Y,nrounds=10,params = list(objective="logregobj"))

# set up their matrix object
dtest <- xgb.DMatrix(m[,1:2], label = m[,"Y"])
# hyperparams
param <- list(max_depth = 2, eta = .01, silent = 1)
bst <- xgb.train(param, data=dtest, nrounds = 10000, logregobj, evalerror, 
                 maximize  =FALSE,watchlist=list())

# compare yhat with Y
data.frame(yhat= predict(bst, dtest, type='response'), Y=Y)



library(xgboost)
data(agaricus.train, package='xgboost')
data(agaricus.test, package='xgboost')
train <- agaricus.train
test <- agaricus.test
bst <- xgboost(data = train$data, label = train$label, max_depth = 2, eta = 1,
               nrounds = 2, objective = "binary:logistic")