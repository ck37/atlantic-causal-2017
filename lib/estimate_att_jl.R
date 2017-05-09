#' TODO: document function.
estimate_att_jl =
  function(# Binary treatment indicator.
           A,
           # Outcome
           Y,
           # Adjustment covariates.
           W,
           # Only gaussian (continuous outcomes) is supported currently.
           family = "gaussian",
           # SL library for outcome regression.
           SL.library,
           # SL library for propensity score estimation.
           g.SL.library,
           # Bounding of estimated propensity score.
           gbounds = c(0.01, 0.99),
           # Cross-validation folds used by SuperLearner.
           V = 10,
           # If T estimate a single outcome regression. If F
           # estimate treatment and control outcomes separately.
           pooled_outcome = T,
           # Set to F to disable parallelism.
           parallel = T,
           # Bounds used for Y when rescaled to [0, 1].
           alpha = c(.0005, .9995),
           verbose = F,
           # Added optional prescreening: (0,0) if you don't want it
           # First term is p-value cut-off, second is minimum # of terms.
           prescreen = c(0.05, 10),
           # If T use SL estimation, otherwise use GLM.
           useSL = T
           ) {

  time_start = proc.time()

  n <- nrow(W)

  if (parallel) {
    if (verbose) cat("estimate_att: Using mcSuperLearner\n")
    # This assumes we're on OSX or Linux; won't work in Windows.
    sl_fn = mcSuperLearner
  } else {
    if (verbose) cat("estimate_att: Using SuperLearner\n")
    sl_fn = SuperLearner
  }

  if (verbose) {
    cat("estimate_att: Preprocessing data.\n")
  }

  # Convert factors to dummies
  W <- model.matrix(~ ., data = data.frame(W))
 
  # Remove intercept that model.matrix() added.
  W = W[, -1]
  num_cols = ncol(W) 
  # Remove constant columns from W.
  nonconstant_columns = which(apply(W, MARGIN = 2, var) != 0)
  # Make sure we retain matrix format.
  W = W[, nonconstant_columns, drop = F]
  if (verbose) {
    cat("Removed", num_cols-length(nonconstant_columns), "constant columns from W.\n")
  }
  
  # Remove linearly correlated columns from W before running gee.
  # Use caret to identify collinearity.
  linear_combos = caret::findLinearCombos(cbind(Y, A, W))
  
  remove_columns = linear_combos$remove
  
  if (length(linear_combos$remove) > 0) {
    # Remove Y and A from columns if they are included.
    # linear_combos$remove = setdiff(linear_combos$remove, c("Y", "A"))
    
    if (verbose) {
      cat("Removing", length(linear_combos$remove), "W vars due to collinearity:\n")
      cat(paste0(colnames(W)[linear_combos$remove - 2], collapse = ", "), "\n")
      cat("Indices:", paste(linear_combos$remove, collapse = ", "), "\n") 
    }
    
    # Make sure we don't switch to a vector if only 1 column remains.
    W = W[, !colnames(W) %in% colnames(W)[linear_combos$remove - 2], drop = F]
    
    # Check if Y, A, W is full rank.
    # Check for full-rank covariate matrix so that glm() can run without error.
    data_mat = cbind(Y, A, W)
    
    # Compute covariance matrix.
    cov_mat = cov(data_mat)
    
    # Compute QR decomp of covariance matrix.
    qr_cov = qr(cov_mat)
    
    # These need to be equal for the covariance matrix to be full rank.
    if (ncol(cov_mat) != qr_cov$rank && verbose) {
      cat("Warning: covariance of (Y, A, W) matrix is not full rank.\n")
      cat("Covariance columns:", ncol(cov_mat), "QR rank:", qr_cov$rank, "\n")
    }
  }

  # TODO: remove collinear columns from W?

  #  Identify range of outcome variable.
  a <- min(Y)
  b <- max(Y)

  # Rescale Y to [0, 1]
  Ystar <- (Y - a) / (b - a)

  if (!useSL) {
    # Estimate outcome regression
    QAWfit <- suppressWarnings(glm(Y ~ W + A, family = 'gaussian'))

    # Estimate propensity score.
    gfit <- glm(A ~ W, family = 'binomial')

    # Predict smoothed potential outcome (Q0_bar) under A = control.
    new0        <- data.frame(W,0)
    names(new0) <- c(colnames(W),"A")
    Q0W         <- suppressWarnings(predict(QAWfit, newdata = new0, type = 'response'))

    # Predict Q1W with all units set to A = treated.
    new1   <- new0
    new1$A <- 1
    Q1W    <- suppressWarnings(predict(QAWfit, newdata = new1, type = 'response'))

    # Predict propensity score.
    g1W <- predict(gfit, type = 'response')
    g1W <- .bound(g1W, gbounds)

    Q.unbd <- cbind(QAW = ifelse(A == 1, Q1W, Q0W),
                    Q0W = Q0W,
                    Q1W = Q1W)
  } else {

  # Prescreening

  #Make prescreening optional
  if (sum(prescreen)>0) {
    if (verbose) cat("Keep covariates with univariate associations. \n")

    # Identify non-binary variables.
    nonbinary <- apply(W,2,function(x) { length(unique(x))>2 })
    
    #Keep for all variables

    keep <- which(prescreen.uni(Y, A, W, alpha = prescreen[1], min=prescreen[2]))
    
    keep.nonbinary<-data.frame(t(subset(t(nonbinary),select=keep)))
    names(keep.nonbinary)<-"val"
    
    keep.nonbin_sub<-subset(keep.nonbinary, val=="TRUE")
    keep.nonbinary<-names(data.frame(W)) %in% row.names(keep.nonbin_sub)
    
    # Initialize to NULL so that cbind() will still work.
    Wsq = NULL

    if (length(which(keep.nonbinary)) > 0) {
      
      new<-cbind.data.frame(names(data.frame(W[, keep.nonbinary])),prescreen.uni(Y, A, W[, keep.nonbinary]^2, alpha=prescreen[1], min = 0))
      names(new)<-c("name","val")
      
      keep.nonbin_sub<-subset(new, val=="TRUE")
      keep.sq<-names(data.frame(W)) %in% keep.nonbin_sub$name
      
      if (sum(keep.sq) > 0) {
        Wsq <- W[, keep.sq, drop = FALSE]^2
        colnames(Wsq) <- paste0(colnames(Wsq), "sq")
      }
    }

    # Add squared terms to X.
    X <- cbind(W[, keep], Wsq)
    n.columns <- ncol(X)

  } else {
    if (verbose) cat("Keep all covariates. \n")

    Wsq <- W^2
    colnames(Wsq) <- paste0(colnames(Wsq), "sq")

    # Add squared terms to X.
    X<-cbind(W, Wsq)
    n.columns <- ncol(X)
  }

  if (verbose) {
    cat("\nestimate_att: Estimating g.\n")
  }

  g.SL <- try(sl_fn(Y = A,
                    X = data.frame(X),
                    SL.library = g.SL.library,
                    verbose = verbose,
                    family = "binomial",
                    # Stratify the CV folds to maximize power.
                    cvControl = list(V = V, stratifyCV = T)))

  if (class(g.SL) == "try-error") {
    if (verbose) {
      cat("SL failed for g, trying simple SL.\n")
    }
      g.SL2 <- try(sl_fn(Y = A,
                        X = data.frame(X),
                        SL.library = c("SL.mean", "SL.glmnet"),
                        verbose = verbose,
                        family = "binomial",
                        # Stratify the CV folds to maximize power.
                        cvControl = list(V = V, stratifyCV = T)))
      
      if(class(g.SL2) == "try-error") {
        if (verbose) {
          cat("Both simple and complex SL failed for g, falling back to glm().\n")
        }
          g1W <- predict(glm(A ~ X, family = "binomial"), type = "response")
      }
        } else {
          cat("Propensity score SuperLearner results:\n")
          print(g.SL)

          g1W <- g.SL$SL.predict
        }

  # Bound g away from 0, 1.
  g1W <- .bound(g1W, gbounds)

  # Create indicator for the control group.
  A0 <- A == 0

  if (!pooled_outcome) {
    # Stratified outcome regression version.

    if (verbose) {
      cat("\nestimate_att: Estimating control regression.\n")
    }

    # Outcome regression for control units.
    # Note that we are modeling Y (original scale) rather than Ystar (scaled).
    m.SL.A0 <- try(sl_fn(Y = Y[A0],
                         X = as.data.frame(X[A0, ]),
                         # Predicted potential outcome for all observations.
                         newX = as.data.frame(X),
                         SL.library = SL.library,
                         family = family,
                         verbose = verbose,
                         cvControl = list(V = V)))

    if (class(m.SL.A0) != "try-error") {
      cat("Outcome regression for controls:\n")
      print(m.SL.A0)
    } else {
      cat("Outcome regression for controls failed! Using simple SL instead.")
      
      m.SL.A0 <- try(sl_fn(Y = Y[A0],
                           X = as.data.frame(X[A0, ]),
                           # Predicted potential outcome for all observations.
                           newX = as.data.frame(X),
                           SL.library = c("SL.mean", "SL.glmnet"),
                           family = family,
                           verbose = verbose,
                           cvControl = list(V = V)))
      
      cat("Outcome regression for controls:\n")
      print(m.SL.A0)
    }

    if (verbose) {
      cat("\nestimate_att: Estimating treated regression.\n")
    }

    # Outcome regression for treated units.
    # Note that we are modeling Y (original scale) rather than Ystar (scaled)
    m.SL.A1 <- sl_fn(Y = Y[!A0],
                     X = as.data.frame(X[!A0,]),
                     # Predicted potential outcome for all observations.
                     newX = as.data.frame(X),
                     SL.library = SL.library,
                     family = family,
                     verbose = verbose,
                     cvControl = list(V = V))

    if (class(m.SL.A1) != "try-error") {
      cat("Outcome regression for treatment:\n")
      print(m.SL.A1)
    } else {
      cat("Outcome regression for treatment failed! Using simple SL instead.")
      
      m.SL.A1 <- sl_fn(Y = Y[!A0],
                       X = as.data.frame(X[!A0,]),
                       # Predicted potential outcome for all observations.
                       newX = as.data.frame(X),
                       SL.library = c("SL.mean", "SL.glmnet"),
                       family = family,
                       verbose = verbose,
                       cvControl = list(V = V))
      
      cat("Outcome regression for treatment:\n")
      print(m.SL.A1)      
      
    }

    Q0W = m.SL.A0$SL.predict
    Q1W = m.SL.A1$SL.predict
    QAW = ifelse(A == 1, Q1W, Q0W)

    Q.unbd = cbind(QAW = QAW, Q1W = Q1W, Q0W = Q0W)

  } else {
    # Pooled outcome regression version.
    if (verbose) {
      cat("\nestimate_att: Estimating outcome regression.\n")
    }

    # Create dataframe used for prediction.
    # We want to get Q0W and Q1W.
    # We can then create QAW based on those and speed up
    # prediction just a tad.
    pred_X = rbind(#data.frame(A = A, X),
                   data.frame(A = 0, X),
                   data.frame(A = 1, X))

    # Confirm that we have n * 2 rows, otherwise fail out.
    stopifnot(nrow(pred_X) == n * 2)

    # Outcome regression for all units.
    # Note that we are modeling Y (original scale) rather than Ystar (scaled).
    sl_outcome = try(sl_fn(Y = Y,
                           X = data.frame(A = A, X),
                           # Predicted potential outcome for all observations.
                           newX = pred_X,
                           SL.library = SL.library,
                           family = family,
                           verbose = verbose,
                           cvControl = list(V = V)))

    if (class(sl_outcome) != "try-error") {
      cat("Outcome regression:\n")
      print(sl_outcome)
    } else {
      cat("Outcome regression failed! Using simple SL instead.")
      
      sl_outcome = try(sl_fn(Y = Y,
                             X = data.frame(A = A, X),
                             # Predicted potential outcome for all observations.
                             newX = pred_X,
                             SL.library = c("SL.mean", "SL.glmnet"),
                             family = family,
                             verbose = verbose,
                             cvControl = list(V = V)))
      
      cat("Outcome regression:\n")
      print(sl_outcome) 
    }

    # Order of predictions is Q0W then Q1W.
    Q0W = sl_outcome$SL.predict[1:n]
    Q1W = sl_outcome$SL.predict[n + 1:n]
    QAW = ifelse(A == 1, Q1W, Q0W)

    Q.unbd <- cbind(QAW = QAW, Q0W = Q0W, Q1W = Q1W)

  } # Done with pooled regression version.
  } # Done with glm vs SL option.

  colnames(Q.unbd) <- c("QAW", "Q0W", "Q1W")

  # Bound Q on the original scale.
  Q_orig_scale = .bound(Q.unbd, c(a, b))

  # Then convert to [0, 1] scale and bound within alpha away from 0, 1.
  Q <- .bound((Q_orig_scale - a) / (b - a), alpha)

  if (verbose) {
    cat("\nestimate_att: fluctuating via update().\n")
  }

  # Fluctuation via weighted logistic regression.
  initdata = data.frame(A = A, Y = Ystar, Q = Q, g = g1W)
  results.wtupdate <- suppressWarnings(update(initdata))

  # Convert estimate back to original scale.
  est <- results.wtupdate[1] * (b - a)

  # Convert se back to original scale.
  se <- sqrt(results.wtupdate[2]) * (b - a)

  ##############
  # Difference in our unit-level estimated potential outcomes, or \hat{tau_i}.
  # These are already unscaled so we don't need to multiply by (b - a)
  unit_est = Q1W - Q0W
  
  
  # Remove near-zero variance columns.
  # Not working :/ "The following pre-processing methods were eliminated: 'nzv'"
  if (F) {
    preproc = caret::preProcess(data.frame(W), method = "nzv")
    W = predict(preproc, W)
  }
  
  # Sandwich variance estimate for linear model with all effect modifiers
  glm_fit = suppressWarnings(glm(Y ~ A + X, family = 'gaussian'))
  Sigma   = vcovHC(glm_fit, type = "HC1")
  
  # Generate matrix of differences in model matrices for treated and controls
  dat0     = dat1 = data.frame(A,X)
  dat0$A   = 0
  dat1$A   = 1
  new0     = model.matrix(~ A + X, data = dat0)
  new1     = model.matrix(~ A + X, data = dat1)
  
  # If vcovHC returns NaNs, run a more stringent screener on X
  if(sum(is.nan(Sigma)) > 0){
    keep        = which(prescreen.uni(Y, A, X, alpha = .05, min=1))
    Xscr        = X[, keep]
    glm_fit_scr = suppressWarnings(glm(Y ~ A + Xscr, family = 'gaussian'))
    Sigma       = vcovHC(glm_fit, type = "HC1")
    dat0     = dat1 = data.frame(A,Xscr)
    dat0$A   = 0
    dat1$A   = 1
    new0     = model.matrix(~ A + Xscr, data = dat0)
    new1     = model.matrix(~ A + Xscr, data = dat1)
    
    # If still getting NaNs, return NULL for CI
    if(sum(is.nan(Sigma)) > 0){
      cat("Unit-level effect inference failed.\n")
      ci_lower = NULL
      ci_upper = NULL
    }
  }
  
  new_diff = new1 - new0
  
  # Create vector of variances for unit-level effect estimates
  unit_var = apply(new_diff, 1, function(x){t(x) %*% Sigma %*% x})
  
  ci_lower = unit_est - qnorm(.975) * sqrt(unit_var)
  ci_upper = unit_est + qnorm(.975) * sqrt(unit_var)

  # Compile unit-level effects, preferable with a ci_lower and ci_upper.
  unit_estimates = data.frame(est = unit_est,
                              ci_lower = ci_lower,
                              ci_upper = ci_upper)

  time_end = proc.time()

  results = list(est = est,
                 se = se,
                 b = b,
                 a = a,
                 ci_lower = est - qnorm(.975) * se,
                 unit_estimates = unit_estimates,
                 ci_upper = est + qnorm(.975) * se,
                 time = time_end - time_start)

  return(results)
}
