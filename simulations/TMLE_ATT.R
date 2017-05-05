# Susan Gruber and Mark van der Laan
# October 26, 2015
# Supplemental Materials for One-Step Targeted Minimum Loss-based Estimation
# Based on Universal Least Favorable One- Dimensional Submodels.
# R source code for
#  - iterative TMLE for ATT parameter using one epsilon
#  - one-Step TMLE for ATT parameter
#  - Simulation Study 1 from the paper
#  - Simulation Study 2 from the paper

#--------------Iterative TMLE for ATT parameter (using one epsilon) -------------------------------
calcATT.iter <- function(Y, A, Q, g1W, tol = tol, maxIter=50, gbounds, Qbounds){
	iter <- 0
	epsilon <- Inf
	calc_h2 <- function(Q){
		psi <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q)
		return(Q[,"Q1W"] - Q[,"Q0W"] - psi)
	}
	while(iter <= maxIter & (abs(epsilon) > tol)){
		iter <- iter + 1
		h1.A1 <-  rep(1, length(Y))
		h1.A0 <-  -g1W/(1-g1W)
		h1 <- h1.A1
		h1[A == 0] <- h1.A0[A==0]
		Q.off <- qlogis(Q[,"QAW"])
		h2 <- calc_h2(Q)
		g1W.off <- qlogis(g1W)
		d <- data.frame(X = c(Y, A), offset =  c(Q.off, g1W.off), h = c(h1, h2))
		epsilon <- coef(glm(X ~ -1 + offset(offset) + h, data = d, family = "binomial"))
		epsilon[is.na(epsilon)] <- 0
		Q <- .bound(plogis(qlogis(Q) + epsilon * cbind(h1, h1.A0, h1.A1)), Qbounds)
		g1W <- .bound(plogis(g1W.off + epsilon * h2), gbounds)
	}
	 psi <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q)
     IC <- ( (A - (1-A)*g1W / (1-g1W))*(Y-Q[, "QAW"]) + A *(Q[,"Q1W"]-Q[,"Q0W"] - psi)) / mean(A)
  	return(list(psi = psi, var.psi = var(IC)/length(Y), conv = abs(epsilon) <= tol, iter=iter))
}

#-----------------------------------One-Step TMLE for ATT parameter  ----------------------------------------
oneStepATT <- function(Y, A, Q, g1W, depsilon, max_iter, gbounds, Qbounds, N=NULL){
	n <- length(Y)
	if (is.null(N)) N = n
	q <- mean(A)
	calcLoss <- function(Q, g1W){
			-mean(Y * log(Q[,"QAW"]) + (1-Y) * log(1 - Q[,"QAW"]) + A * log(g1W) + (1-A) * log(1 - g1W))
	}
	psi.prev <- psi  <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q)
	H1.AW =  (A -(1-A) * g1W / (1-g1W))/q
	IC.prev <- IC.cur <- H1.AW* (Y-Q[, "QAW"]) + A*(Q[,"Q1W"]-Q[,"Q0W"] - psi)/q
	deriv <- deriv.prev <-  mean(IC.cur)
	if (deriv > 0) { depsilon <- -depsilon}
	loss.prev <- Inf
 	loss.cur <-  calcLoss(Q, g1W)
 	if(is.nan(loss.cur) | is.na(loss.cur) | is.infinite(loss.cur)) {
 		loss.cur <- Inf
 		loss.prev <- 0
 	}
 	iter <-  0
 	while (loss.prev > loss.cur & iter < max_iter){
 	  IC.prev <-  IC.cur
 	  Q.prev <- Q
 	  g1W.prev <- g1W
 	  psi.prev <- psi
 	  loss.prev <- loss.cur
 	  iter = iter+1
 	  deriv.prev=deriv
 	  if (abs(deriv.prev) < 1/N) break

 		H1 <- cbind(HAW = A/q - (1-A) * g1W / (q * (1-g1W)),
					  H0W =   - g1W/(q * (1-g1W)),
					  H1W =   1/q)
 		H2 <- (Q[,"Q1W"]-Q[,"Q0W"]-psi)/q

 		g1W <- .bound(plogis(qlogis(g1W)-depsilon*H2),gbounds)
 		Q  <- .bound(plogis(qlogis(Q) - depsilon * H1), Qbounds)

 		psi <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q)
 		loss.prev <- loss.cur
 		loss.cur <- calcLoss(Q, g1W)
 		IC.cur <- ((A - (1-A) * g1W / (1-g1W)) * (Y-Q[, "QAW"]) + A *(Q[,"Q1W"]-Q[,"Q0W"] - psi))/q
 		if(is.nan(loss.cur) | is.infinite(loss.cur) | is.na(loss.cur)) {loss.cur <- Inf}
 		deriv = mean(IC.cur)
 		print(psi.prev)
 	}
 	return(list(psi=psi.prev, epsilon = (iter-1)*depsilon,IC=IC.prev,
 	            Q=Q.prev))
}

#------------- bound function --------------
.bound <- function(x, bds){
	x[x > max(bds)] <- max(bds)
	x[x < min(bds)] <- min(bds)
	x
}

n=1000
W1 = rnorm(n)
W2 = abs(rnorm(n))
A = rbinom(n,1,plogis(rnorm(n)+W1+W2))
Y = A*rbinom(n,1,plogis(W2+W1-.4*A))+(1-A)*rbinom(n,1,plogis(W2+W1))
X = data.frame(W1=W1,W2=W2,A=A)
X0 = X
X0$A=0
X1=X
X1$A=1
newX=rbind(X,X1,X0)
W=X
W$A=NULL

library(SuperLearner)

Qkfit = SuperLearner(Y=Y, X=X, newX = newX, family = binomial(),
                     SL.library="SL.glm",method = "method.NNLS",
                     id = NULL, verbose = FALSE, control = list(), cvControl = list(),
                     obsWeights = NULL)
gkfit = SuperLearner(Y=A, X=W, newX = W, family = binomial(),
                     SL.library="SL.glm",method = "method.NNLS",
                     id = NULL, verbose = FALSE, control = list(), cvControl = list(),
                     obsWeights = NULL)
g1W = gkfit$SL.predict
QAW = Qkfit$SL.predict[1:n]
Q0W = Qkfit$SL.predict[(n+1):(2*n)]
Q1W = Qkfit$SL.predict[(2*n+1):(3*n)]
Q = cbind(QAW,Q0W,Q1W)
colnames(Q)=c("QAW","Q0W","Q1W")

res = oneStepATT(Y=Y, A=A, Q=Q, g1W=g1W, depsilon=.0001, max_iter=100000,
                 gbounds=c(1e-4,1-1e-4), Qbounds=c(1e-4,1-1e-4), N=NULL)


res$epsilon

