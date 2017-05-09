
# Multithreaded version of XGBoost when using sequential SuperLearner.
SL.xgboost_threads_4 = function(...) SL.xgboost(..., nthread = 4)

# Faster glmnet.
SL.glmnet_fast = function(...) SL.glmnet(..., nlambda = 20, nfolds = 5)

# Faster randomForest.
SL.randomForest_fast = function(...) SL.randomForest(..., ntree = 200, verbose = F)

# Faster ranger (itself a faster version of RF).
SL.ranger_fast = function(...) SL.ranger(..., num.trees = 200, num.threads = 4)


screen.corRank4 = function(...) {
  screen.corRank(..., rank = 4)
}

screen.corRank8 = function(...) {
  screen.corRank(..., rank = 8)
}