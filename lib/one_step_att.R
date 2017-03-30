#-----------------------------------One-Step TMLE for ATT parameter  ----------------------------------------
# TODO: document arguments.
one_step_att <- function(Y, A, Q,
                         g1W,
                         depsilon,
                         max_iter,
                         gbounds,
                         Qbounds,
                         verbose = F) {
  n <- length(Y)
  q <- mean(A)

  calcLoss <- function(Q, g1W) {
    -mean(Y * log(Q[, "QAW"]) + (1 - Y) * log(1 - Q[, "QAW"]) +
            A * log(g1W) + (1 - A) * log(1 - g1W))
  }

  psi.prev <- psi  <- mean((Q[, "Q1W"] - Q[, "Q0W"]) * g1W / q)
  H1.AW =  (A - (1 - A) * g1W / (1 - g1W)) / q
  IC.prev <- IC.cur <- H1.AW * (Y - Q[, "QAW"]) + A*(Q[,"Q1W"] - Q[,"Q0W"] - psi) / q

  deriv <-  mean(H1.AW * (Y - Q[, "QAW"]) + A * (Q[,"Q1W"] - Q[,"Q0W"] - psi) / q)

  if (deriv > 0) {
    depsilon <- -depsilon
  }

  loss.prev <- Inf
  loss.cur <-  calcLoss(Q, g1W)

  if (is.nan(loss.cur) | is.na(loss.cur) | is.infinite(loss.cur)) {
    loss.cur <- Inf
    loss.prev <- 0
  }

  iter <-  0
  while (loss.prev > loss.cur & iter < max_iter) {
    IC.prev <-  IC.cur
    Q.prev <- Q

    g1W.prev <- g1W
    g1W <- .bound(plogis(qlogis(g1W.prev) - depsilon  * (Q.prev[,"Q1W"] - Q.prev[,"Q0W"] - psi.prev) / q), gbounds)

    H1 <- cbind(HAW = A/q - (1 - A) * g1W.prev / (q * (1 - g1W.prev)),
                H0W = -g1W.prev / (q * (1 - g1W.prev)),
                H1W = 1 / q)

    Q  <- .bound(plogis(qlogis(Q.prev) - depsilon * H1), Qbounds)

    psi.prev <- psi
    psi <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W / q)

    loss.prev <- loss.cur
    loss.cur <- calcLoss(Q, g1W)

    IC.cur <- ((A - (1 - A) * g1W / (1 - g1W)) * (Y - Q[, "QAW"]) + A * (Q[,"Q1W"] - Q[,"Q0W"] - psi)) / q

    if (is.nan(loss.cur) | is.infinite(loss.cur) | is.na(loss.cur)) {
      loss.cur <- Inf
    }

    iter <- iter + 1

    if (verbose) print(psi.prev)
  }

  results = c(psi = psi.prev,
              var.psi = var(IC.prev) / n,
              conv = loss.prev < loss.cur)

  return(results)
}
