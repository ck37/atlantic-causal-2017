# basic truncate function
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
  pmin(pmax(x, lower), upper)
}


# Takes initdata, dataframe or list containing elements: Q0W, Q1W, Y, g, A
update <- function(initdata) {

  # Create clever covariate, which is 0 for treated units.
  HAW_1 = with(initdata, (- (A==0)*g / (1 - g)) / mean(A))
  HAW_2 = with(initdata, (A==1) / mean(A))
  # HAW = with(initdata, ifelse(A == 1, H1W, H0W))
  
  QAW = with(initdata, ifelse(A == 1, Q1W, Q0W))

  # Fit a glm with H as clever covariate.
  fit = glm(Y ~ offset(qlogis(QAW)) + HAW_1 +HAW_2- 1,
            data = initdata,
            family = "binomial")

  # Extract epsilon from the MLE.
  epsilon_hat = fit$coefficients

  if (epsilon_hat == 0) {
    warning("Epsilon hat estimated to be zero.")
  }

  if (is.na(epsilon_hat)) {
    warning("Epsilon hat is NA.")
  }

  H0W = with(initdata, (- g / (1 - g)) / mean(A))
  H1W = with(initdata, 1 / mean(A))
  # Update estimated Q(0, W) and Q(A, W)
  Q0Wstar = with(initdata, plogis(qlogis(Q0W) + epsilon_hat[1] * H0W))
  Q1Wstar = with(initdata, plogis(qlogis(Q1W) + epsilon_hat[2] * H1W))

  num_nas0 = sum(is.na(Q0Wstar))
  num_nas1 = sum(is.na(Q1Wstar))
  if (num_nas0 + num_nas1 > 0) {
    cat("Updated Q0Wstar num NAs is:", num_nas0, "\n","Updated Q1Wstar num NAs is:", num_nas1, "\n")
    print(summary(initdata$Q0W))
    print(summary(initdata$Q1W))
    print(summary(qlogis(initdata$Q0W)))
    print(summary(qlogis(initdata$Q1W)))
    cat("Epsilon hat:", epsilon_hat, "\n")
  }

  # Calculate percentage of treated units that changed after fluctuation update.
  pct_changed = mean(initdata$Q0W[initdata$A == 1] != Q0Wstar[initdata$A == 1] |
                       initdata$Q1W[initdata$A == 1] != Q1Wstar[initdata$A == 1])

  # If the treated units didn't change there is some sort of issue.
  if (pct_changed == 0) {
    warning("Fluctuation did not change any potential outcomes for treated units.")
  }

  # Compile results.
  results = list(Q0Wstar = Q0Wstar,
                 Q1Wstar = Q1Wstar,
                 pct_changed = pct_changed,
                 epsilon_hat = epsilon_hat)

  return(results)
}


# takes initdata, dataframe or list containing elements, Q0W,Y,g,A
updateLS <- function(initdata) {
  H = with(initdata, (A==0)*g/(mean(A)*(1-g)))
  # fit a glm with the weight
  fit = glm(Y ~ offset(Q0W),data=initdata,
            weights=H,family="gaussian")

  eps = fit$coef

  # update and truncate to prevent overstepping outcome bounds
  Q0Wstar = with(initdata, truncate(Q0W + eps,lower=1e-4))

  # Compile results.
  results = list(Q0Wstar = Q0Wstar,
                 epsilon = eps)
  return(results)
}

