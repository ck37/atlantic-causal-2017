
# Multithreaded version of XGBoost when using sequential SuperLearner.
SL.xgboost_threads_4 = function(...) SL.xgboost(..., nthread = 4)

# Faster glmnet.
SL.glmnet_fast = function(...) SL.glmnet(..., nlambda = 20, nfolds = 5)

# Faster randomForest.
SL.randomForest_fast = function(...) SL.randomForest(..., ntree = 200, verbose = F)

# Faster ranger (itself a faster version of RF).
SL.ranger_fast = function(...) SL.ranger(..., num.trees = 200, num.threads = 4)

# Restrict to top 4 variables based on univariate correlation.
# TODO: keep treatment indicator in if we're running outcome regression.
screen.corRank4 = function(...) {
  screen.corRank(..., rank = 4)
}

# Restrict to top 8 variables based on univariate correlation.
# TODO: keep treatment indicator in if we're running outcome regression.
screen.corRank8 = function(...) {
  screen.corRank(..., rank = 8)
}