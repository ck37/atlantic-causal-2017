# Takes initdata, dataframe or list containing elements: A, Y, Q, g
update <- function(initdata) {
  
  n = length(initdata$A)
  
  # Create ``clever covariate'', which is 0 for treated units. This will only be used to update Q01.
  H0W = with(initdata, - g / (1 - g))
  HAW = with(initdata, ifelse(A == 1, 0, H0W))
  
  # Fit a glm with the ``clever covariate'' moved to the weight.
  # TODO: catch error and recover if fluctuation fails.
  # Suppress any warnings here.
  fit = suppressWarnings(glm(Y ~ offset(qlogis(Q.QAW)) + 1,
            weight = -HAW,
            data = initdata,
            family = "binomial"))
  
  # Extract epsilon from the MLE.
  epsilon_hat = fit$coefficients[1]
  
  if (epsilon_hat == 0) {
    warning("Epsilon hat estimated to be zero.")
  }
  
  if (is.na(epsilon_hat)) {
    warning("Epsilon hat is NA.")
  }
  
  # Update estimated Q(0, W)
  Q0Wstar = with(initdata, plogis(qlogis(Q.Q0W) + epsilon_hat))
  
  # Update for Q(1, W) merely ensures mean(Q1W*A) = mean(Y*A), i.e., same mean as Y among treated.
  # This means we can replace with Ystar in estimator, though update may improve variance estimator
  Q1Wstar = with(initdata, Q.Q1W + mean((Y - Q.Q1W) * A)/mean(A))
  
  QAWstar = with(initdata, ifelse(A == 1, Q1Wstar, Q0Wstar))
  
  num_nas0 = sum(is.na(Q0Wstar))
  num_nas1 = sum(is.na(Q1Wstar))
  if (num_nas0 + num_nas1 > 0) {
    cat("Updated Q0Wstar num NAs is:", num_nas0, "\n","Updated Q1Wstar num NAs is:", num_nas1, "\n")
    print(summary(initdata$Q.Q0W))
    print(summary(initdata$Q.Q1W))
    print(summary(qlogis(initdata$Q.Q0W)))
    print(summary(qlogis(initdata$Q.Q1W)))
    cat("Epsilon hat:", epsilon_hat, "\n")
  }
  
  # Calculate percentage of treated units that changed after fluctuation update.
  pct_changed = mean(initdata$Q.Q0W[initdata$A == 1] != Q0Wstar[initdata$A == 1] |
                       initdata$Q.Q1W[initdata$A == 1] != Q1Wstar[initdata$A == 1])
  
  # If the treated units didn't change there is some sort of issue.
  if (pct_changed == 0) {
    warning("Fluctuation did not change any potential outcomes for treated units.")
  }
  
  psi   = with(initdata, sum((A == 1) * (Y - Q0Wstar)) / sum(A))
  
  Dstar = with(initdata, ((A == 1) - (A == 0) * g / (1 - g)) / mean(A) * (Y - QAWstar) + 
                          (A == 1) / mean(A) * (Q1Wstar - Q0Wstar - psi))

  num_nas = sum(is.na(Dstar))
  if (num_nas > 0) {
    cat("Dstar num NAs is:", num_nas, "\n")
  }
  
  # Compile results.
  results = c(psi     = psi,
              var.psi = (n-1)*var(Dstar)/n^2)
  
  return(results)
}
