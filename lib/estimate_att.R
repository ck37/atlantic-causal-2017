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
                        alpha = c(.0005, .9995),
                        verbose = F
                        ) {

  time_start = proc.time()

  n <- nrow(W)

  # Convert factors to dummies
  W <- model.matrix(~ ., data = W)

  # Remove intercept that model.matrix() added.
  W = W[, -1]

  # This isn't general, but true for these data.
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

  keep <- which(prescreen.uni(Y, A, W, alpha = .05))
  keep.nonbinary <- nonbinary[nonbinary %in% keep]

  if (length(keep.nonbinary) > 0) {
    keep.sq <- keep.nonbinary[prescreen.uni(Y, A, W[,keep.nonbinary, drop = FALSE]^2, min = 0)]
    if (sum(keep.sq) > 0) {
      Wsq <- W[, keep.sq, drop = FALSE]^2
      colnames(Wsq) <- paste0(colnames(Wsq), "sq")
    }
  }
  X <- cbind(W[,keep],Wsq)
  n.columns <- ncol(X)

  g.SL <- try(SuperLearner(Y = A,
                           X = data.frame(X),
                           SL.library = g.SL.library,
                           family = "binomial",
                           cvControl = list(V = V)))

  if (class(g.SL) == "try-error") {
    g1W <- predict(glm(A ~ X, family = "binomial"), type = "response")
  } else {
    g1W <- g.SL$SL.predict
  }

  g1W <- .bound(g1W, gbounds)

  A0 <- A == 0

  m.SL.A0 <- SuperLearner(Y = Y[A0],
                          X = as.data.frame(X[A0,]),
                          newX = as.data.frame(X),
                          SL.library = SL.library,
                          family = family,
                          cvControl = list(V = V))

  m.SL.A1 <- SuperLearner(Y = Y[!A0],
                          X = as.data.frame(X[!A0,]),
                          newX = as.data.frame(X),
                          SL.library = SL.library,
                          family = family,
                          cvControl = list(V = V))

  Q.unbd <- cbind(QAW = A * m.SL.A1$SL.predict + (1 - A) * m.SL.A0$SL.predict,
                  Q0W = m.SL.A0$SL.predict, Q1W = m.SL.A1$SL.predict)
  colnames(Q.unbd) <- c("QAW", "Q0W", "Q1W")

  Q <- .bound((.bound(Q.unbd, c(a,b)) - a) / (b - a), alpha)

  results.oneStep <- one_step_att(Y = Ystar,
                                  A = A,
                                  Q = Q,
                                  g1W = g1W,
                                  depsilon = depsilon,
                                  max_iter = max(1000, 2/depsilon),
                                  gbounds = gbounds,
                                  Qbounds = alpha)

  # Convert estimate back to original scale.
  est <- results.oneStep[1] * (b - a)

  # Convert se back to original scale.
  se <- sqrt(results.oneStep[2]) * (b - a)


  time_end = proc.time()

  results = list(est = est,
                 se = se,
                 b = b,
                 a = a,
                 ci_lower = est - 1.96 * se,
                 ci_upper = est + 1.96 * se,
                 time = time_end - time_start)
  return(results)
}
