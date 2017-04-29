#' TODO: document function.
estimate_att = function(A,
                        Y,
                        W,
                        # Only gaussian is supported currently.
                        family = "gaussian",
                        SL.library,
                        g.SL.library,
                        gbounds = c(0.01, 0.99),
                        depsilon =  0.001,
                        V = 5,
                        # Set to F to disable parallelism.
                        parallel = T,
                        alpha = c(.0005, .9995),
                        verbose = F
                        ) {

  time_start = proc.time()

  n <- nrow(W)

  if (parallel) {
    if (verbose) cat("estimate_att: Using mcSuperLearner\n")
    sl_fn = mcSuperLearner
  } else {
    if (verbose) cat("estimate_att: Using SuperLearner\n")
    sl_fn = SuperLearner
  }

  if (verbose) {
    cat("estimate_att: Preprocessing data.\n")
  }

  # Convert factors to dummies
  W <- model.matrix(~ ., data = W)

  # Remove intercept that model.matrix() added.
  W = W[, -1]

  # SG: This isn't general, but true for these data.
  nonbinary <- which(colMeans(W) > 1)

  Wcat <- matrix(as.integer(W[,nonbinary] < rep(colMeans(W[,nonbinary]), each = n)), nrow = n, byrow = FALSE)
  colnames(Wcat) <- paste0(colnames(W[,nonbinary]), "cat")

  W <- cbind(W, x_3aug = as.integer(W[,"x_3"] > 0), x_4aug = as.integer(W[,"x_4"] > 0),
            Wcat
  )

  #  Identify range of outcome variable.
  a <- min(Y)
  b <- max(Y)

  # Rescale Y to [0, 1]
  Ystar <- (Y - a) / (b - a)

  # Prescreening
  keep <- which(prescreen.uni(Y, A, W, alpha = .05))
  keep.nonbinary <- nonbinary[nonbinary %in% keep]

  # Initialize to NULL so that cbind() will still work.
  Wsq = NULL

  if (length(keep.nonbinary) > 0) {
    keep.sq <- keep.nonbinary[prescreen.uni(Y, A, W[, keep.nonbinary, drop = FALSE]^2, min = 0)]
    if (sum(keep.sq) > 0) {
      Wsq <- W[, keep.sq, drop = FALSE]^2
      colnames(Wsq) <- paste0(colnames(Wsq), "sq")
    }
  }

  # Add squared terms to X.
  X <- cbind(W[, keep], Wsq)
  n.columns <- ncol(X)

  if (verbose) {
    cat("\nestimate_att: Estimating g.\n")
  }

  g.SL <- try(sl_fn(Y = A,
                    X = data.frame(X),
                    SL.library = g.SL.library,
                    verbose = verbose,
                    family = "binomial",
                    cvControl = list(V = V)))

  if (class(g.SL) == "try-error") {
    # TODO: fall back to a simpler SL algorithm before using glm().
    if (verbose) {
      cat("SL failed for g, falling back to glm().\n")
    }
    g1W <- predict(glm(A ~ X, family = "binomial"), type = "response")
  } else {
    cat("Propensity score SuperLearner results:\n")
    print(g.SL)

    g1W <- g.SL$SL.predict
  }

  g1W <- .bound(g1W, gbounds)

  # Create indicator for the control group.
  A0 <- A == 0

  if (verbose) {
    cat("\nestimate_att: Estimating control regression.\n")
  }

  # TODO: convert these two to a pooled regression.

  # Outcome regression for control units.
  # TODO: check for errors and fall back to simpler library like we do for g.
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
    cat("Outcome regression for controls failed!")
  }

  if (verbose) {
    cat("\nestimate_att: Estimating treated regression.\n")
  }

  # Outcome regression for treated units.
  # TODO: check for errors and fall back to simpler library like we do for g.
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
    cat("Outcome regression for treatment failed!")
  }

  Q.unbd <- cbind(QAW = A * m.SL.A1$SL.predict + (1 - A) * m.SL.A0$SL.predict,
                  Q0W = m.SL.A0$SL.predict,
                  Q1W = m.SL.A1$SL.predict)
  colnames(Q.unbd) <- c("QAW", "Q0W", "Q1W")

  Q <- .bound((.bound(Q.unbd, c(a, b)) - a) / (b - a), alpha)

  if (verbose) {
    cat("\nestimate_att: fluctuating via one_step_att().\n")
  }

  # Fluctuation via the one-step algorithm.
  # TODO: convert this to an interated TMLE fluctuation.
  results.oneStep <- one_step_att(Y = Ystar,
                                  A = A,
                                  Q = Q,
                                  g1W = g1W,
                                  depsilon = depsilon,
                                  max_iter = max(1000, 2 / depsilon),
                                  verbose = verbose,
                                  gbounds = gbounds,
                                  Qbounds = alpha)

  # Convert estimate back to original scale.
  est <- results.oneStep[1] * (b - a)

  # Convert se back to original scale.
  se <- sqrt(results.oneStep[2]) * (b - a)

  ##############
  # Unit-level estimate uses observed outcome and the complementary modeled outcome.

  # NOTE: if we fit a pooled outcome regression for both treatment & control,
  # we should update this to not use separate models for A1 and A0.

  # Best Y1 is observed Y if A = 1 and SL.predict otherwise.
  best_y1 = A * Y + (1 - A) * m.SL.A1$SL.predict

  # Best Y0 is observed Y if A = 0 and SL.predict otherwise.
  best_y0 = (1 - A) * Y + A * m.SL.A0$SL.predict

  # Difference is our unit-level estimated treatment effect.
  # These are already unscaled so we don't need to multiply by (b - a)
  unit_est = best_y1 - best_y0

  # Compile unit-level effects, preferable with a ci_lower and ci_upper.
  unit_estimates = data.frame(est = unit_est,
                              ci_lower = NA,
                              ci_upper = NA)

  time_end = proc.time()

  results = list(est = est,
                 se = se,
                 b = b,
                 a = a,
                 ci_lower = est - 1.96 * se,
                 unit_estimates = unit_estimates,
                 ci_upper = est + 1.96 * se,
                 time = time_end - time_start)

  return(results)
}
