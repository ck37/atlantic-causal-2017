# basic truncate function
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
  pmin(pmax(x, lower), upper)
}


# takes initdata, dataframe or list containing elements, Q0W,Y,g,A 
update <- function(initdata) {
  H = with(initdata, (A==0)*g/(mean(A)*(1-g)))
  # fit a glm with the weight
  fit = glm(Y ~ offset(qlogis(Q0W)),data=initdata,
               weights=H,family="binomial")

  # update
  Q0Wstar = with(initdata, plogis(qlogis(Q0W)+fit$coef))
  return(Q0Wstar)
}


# takes initdata, dataframe or list containing elements, Q0W,Y,g,A 
updateLS <- function(initdata) {
  H = with(initdata, (A==0)*g/(mean(A)*(1-g)))
  # fit a glm with the weight
  fit = glm(Y ~ offset(Q0W),data=initdata,
            weights=H,family="gaussian")
  
  # update yo momma
  Q0Wstar = with(initdata, truncate(Q0W+fit$coef,lower=1e-4))
  return(Q0Wstar)
}

