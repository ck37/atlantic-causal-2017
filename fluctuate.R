# basic truncate function
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
  pmin(pmax(x, lower), upper)
}


# takes initdata, dataframe or list containing elements, Q0W,Y,g,A
update <- function(initdata) {
  H = with(initdata, (A==0)*g/(mean(A)*(1-g)))

  # fit a glm with the weight
  fit = glm(Y ~ offset(qlogis(Q0W)),
            weights = H, data = initdata,
            family = "binomial")

  eps = fit$coefficients[1]

  # update
  Q0Wstar = with(initdata, plogis(qlogis(Q0W) + eps))

  # Compile results.
  results = list(Q0Wstar = Q0Wstar,
                 epsilon = eps)
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

