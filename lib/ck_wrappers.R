
# Multiple versions of XGBoost if we can afford the extra computation.
# Keep the grid pretty small: 6 learners.
sl_xgb = create.Learner("SL.xgboost", detailed_names = T,
                        params = list(nthread = 4, ntrees = 1000),
                        tune = list(max_depth = c(2, 4),
                                    shrinkage = c(0.05, 0.1, 0.2)))
length(sl_xgb$names)

# Small SVM grid for cost parameter and kernel.
sl_ksvm = create.Learner("SL.ksvm", detailed_names = T,
                         tune = list(
                            # Regularization parameter, could be 2^-5 to 2^15.
                            C = 2^c(-3, 0, 4),
                            kernel = c("rbfdot", "laplacedot"))) 
length(sl_ksvm$names)
                             

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