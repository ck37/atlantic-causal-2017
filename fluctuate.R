# basic truncate function
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
  pmin(pmax(x, lower), upper)
}


# Takes initdata, dataframe or list containing elements: Q0W, Y, g, A
update <- function(initdata) {

  # Create clever covariate, which is 0 for treated units.
  H = with(initdata, (A == 0) * g / (mean(A) * (1 - g)))

  # Fit a glm putting the clever covariate in the weight.
  fit = glm(Y ~ offset(qlogis(Q0W)),
            weights = H, data = initdata,
            family = "binomial")

  # Extract epsilon from the MLE.
  epsilon_hat = fit$coefficients[1]

  if (epsilon_hat == 0) {
    warning("Epsilon hat estimated to be zero.")
  }

  # Update estimated Q(0, W) - we only care about the treated units though.
  Q0Wstar = with(initdata, plogis(qlogis(Q0W) + epsilon_hat))

  # Calculate percentage of treated units that changed after fluctuation update.
  pct_changed = mean(initdata$Q0W[initdata$A == 1] != Q0Wstar[initdata$A == 1])

  # If the treated units didn't change there is some sort of issue.
  if (pct_changed == 0) {
    warning("Fluctuation did not change any potential outcomes for treated units.")
  }

  # Compile results.
  results = list(Q0Wstar = Q0Wstar,
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

