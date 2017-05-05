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
oneStepATT <- function(Y, A, Q, g1W, depsilon, max_iter, gbounds, Qbounds){
	n <- length(Y)
	q <- mean(A)
	calcLoss <- function(Q, g1W){
			-mean(Y * log(Q[,"QAW"]) + (1-Y) * log(1 - Q[,"QAW"]) + A * log(g1W) + (1-A) * log(1 - g1W))
	}
	psi.prev <- psi  <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q)
	H1.AW =  (A -(1-A) * g1W / (1-g1W))/q
	IC.prev <- IC.cur <- H1.AW* (Y-Q[, "QAW"]) + A*(Q[,"Q1W"]-Q[,"Q0W"] - psi)/q
	deriv <-  mean(H1.AW* (Y-Q[, "QAW"]) + A*(Q[,"Q1W"]-Q[,"Q0W"] - psi)/q )
	g1W.prev = g1W.cur = g1W
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
		g1W <- .bound(g1W,gbounds)
 		H1 <- cbind(HAW = A/q - (1-A) * g1W.prev / (q * (1-g1W.prev)),
					  H0W =   - g1W.prev/(q * (1-g1W.prev)),
					  H1W =   1/q)
 		Q  <- .bound(plogis(qlogis(Q.prev) - depsilon * H1), Qbounds)
 		psi.prev <- psi
 		psi <- mean((Q[,"Q1W"] - Q[, "Q0W"]) * g1W/q)
 		loss.prev <- loss.cur
 		loss.cur <- calcLoss(Q, g1W)
 		IC.cur <- ((A - (1-A) * g1W / (1-g1W)) * (Y-Q[, "QAW"]) + A *(Q[,"Q1W"]-Q[,"Q0W"] - psi))/q
 		if(is.nan(loss.cur) | is.infinite(loss.cur) | is.na(loss.cur)) {loss.cur <- Inf}
 		iter <- iter + 1
 		print(psi.prev)
 	}
 	 	return(list(psi = psi.prev, var.psi = var(IC.prev)/n,  epsilon = (iter-1)*depsilon))
 }

#------------- bound function --------------
.bound <- function(x, bds){
	x[x > max(bds)] <- max(bds)
	x[x < min(bds)] <- min(bds)
	x
}

#----------Simulation Study I----------------
set.seed(1)
N <- c(100, 1000)
gbounds <- c(10^-9, 1-10^-9)
depsilon <-  0.001
niter <- 1000
Qformulas<- list(Qcor = "Y~A + W1 + W2 + W3", Qmis1 = "Y~A + W1", Qmis2 = "Y ~ A")
gformulas <- list(gcor = "A~W1 + W2 + W3", gmis1 = "A ~ W1", gmis2 = "A ~ 1")
Qbounds <- matrix(c(-Inf, Inf, 10^-9, 1-10^-9), byrow=TRUE, ncol = 2)
for (n in N){
  Qcgc <- Qcgm1 <- Qc1gm2 <- Qm1gc <- Qm1gm1 <- Qm1gm2 <- Qm2gc <- Qm2gm1 <- Qm2gm2<- matrix(NA, nrow = niter, ncol = 4,
	dimnames = list(NULL, c("oneStep_unbd",  "iter_unbd", "oneStep_bd",  "iter_bd")))
   res.psi <- list(Qcgc, Qcgm1, Qc1gm2 , Qm1gc , Qm1gm1, Qm1gm2, Qm2gc, Qm2gm1, Qm2gm2)
   for (i in 1:niter){
     W1 <- rnorm(n)
	 W2 <- rnorm(n)
	 W3 <- rbinom(n, 1, .5)
	 A <- rbinom(n, 1, plogis(-0.4 - 0.2* W1 - 0.4*W2 + 0.3*W3))
	 logitY <- -1.2 -1.2*A - 0.1*W1 - 0.2*W2  - 0.1*W3
	 Y <- cbind(rbinom(n, 1, plogis(logitY)))
	 for (j in 1:length(Qformulas)){
		Qform <- Qformulas[[j]]
	  	m <- glm(Qform, family = "binomial")
	    Q <- .bound(plogis(cbind(QAW = predict(m),
					Q0W = predict(m, newdata = data.frame(A= 0, W1, W2, W3)),
					Q1W = predict(m, newdata = data.frame(A= 1, W1, W2, W3)))), c(10^-5, 1-10^-5))
		g1W <- q <- mean(A)
	    for (k in 1:length(gformulas)){
	    	for (qq in 1:nrow(Qbounds)){
	    	g1W <- .bound(predict(glm(gformulas[[k]], family = "binomial"), type = "response"), gbounds)
		 	res.oneStep <- oneStepATT(Y = Y, A = A, Q = Q, g1W = g1W,
  			  			depsilon = depsilon, max_iter = max(1000, 2/depsilon), gbounds = gbounds, Qbounds = Qbounds[qq,])
    		res.iter <- try(calcATT.iter(Y = Y, A = A, Q = Q,  g1W = g1W, tol = 10^-5, maxIter = 500, gbounds = gbounds,
    					 Qbounds = Qbounds[qq,]))
  			if(class(res.iter) == "try-error") {
  				res.iter <- NULL
  				res.iter$psi <- res.iter$var.psi <- NA
  			}
  			res.psi[[3*(j-1)+k]][i, ((qq-1)*2 + 1):(2 + (qq-1)*2)] <-  c(res.oneStep$psi,  res.iter$psi)
  		 }}
	 }}
   save(res.psi, Qformulas, gformulas, depsilon, Qbounds, file = paste0("Sim1_", n, ".Rdata"))
}

#----------Simulation Study II---------------
N <- 100
gbounds <- c(10^-9, 1-10^-9)
depsilon <-  0.001
niter <- 1000
Qformulas<- list(Qc1 = "Y ~ A + W1 + W2 + W3",
						Qc2 = "Y ~ A + W1 + W2 + W3 + W4",
						Qc3 = "Y ~ A + W1 + W2 + W3 + W4 + W5",
						Qc4 = "Y ~ A + W1 +W2 +W3 +W4 +W5 +W6",
						Qc5 = "Y ~ A + W1 + W2 +W3 +W4 +W5 + W6 + W7",
						Qc6 = "Y ~ A + W1 +W2 +W3 +W4 +W5 +W6 +W7 + W8 ",
						Qc7 = "Y ~ A + W1 +W2 +W3 +W4 +W5 +W6 +W7 + W8 + W9 ",
						Qc8 = "Y ~ A + W1 +W2 +W3 +W4 +W5 +W6 +W7 + W8 + W9 + W10",
						Qc9 = "Y ~ A + W1 +W2 +W3 +W4 +W5 +W6 +W7 + W8 + W9 + W10 + W11",
						Qc10 = "Y ~ A + W1 +W2 +W3 +W4 +W5 +W6 +W7 + W8 + W9 + W10 + W11 + W12")
gformulas <- list(gc = "A~W1 + W2 + W3")
ngform <- length(gformulas)
Qbounds <- matrix(c(-Inf, Inf, 10^-9, 1-10^-9), byrow=TRUE, ncol = 2)
for (n in N){
	Qc1gc <- Qc2gc <- Qc3gc   <- Qc4gc  <- Qc5gc  <- Qc6gc <- matrix(NA, nrow = niter, ncol = 5,
		dimnames = list(NULL, c("psi.init", "oneStep_unbd", "iter_unbd", "oneStep_bd",   "iter_bd")))
	Qc7gc  <- Qc8gc  <- Qc9gc   <- Qc10gc <- Qc6gc
	 res.psi <- list(Qc1gc,  Qc2gc , Qc3gc, Qc4gc, Qc5gc,
 						Qc6gc, Qc7gc, Qc8gc , Qc9gc,  Qc10gc)
	for (i in 1:niter){
	  W1 <- rnorm(n)
	  W2 <- rnorm(n)
	  W3 <- rbinom(n, 1, .5)
	  W4 <- rnorm(n)
	  W5 <- rnorm(n)
	  W6 <- rnorm(n)
	  W7 <- rnorm(n)
	  W8 <- rnorm(n)
	  W9 <- rnorm(n)
	  W10 <- rnorm(n)
	  W11 <- rnorm(n)
	  W12 <- rnorm(n)
	  A <- rbinom(n, 1, plogis(-0.4 - .2* W1 - .4*W2 + .3*W3))
	  logitY1 <- -2 - 1.2*A - .1*W1 + .3*W2  + .3*W3
	  Y <- rbinom(n, 1, plogis(logitY1))
	   for (j in 1:length(Qformulas)){
			Qform <- Qformulas[[j]]
		  	m <- glm(Qform, family = "binomial")
		    Q <- .bound(plogis(cbind(QAW = predict(m),
						Q0W = predict(m, newdata = data.frame(A= rep(0,n))),
						Q1W = predict(m, newdata = data.frame(A= rep(1,n))))), c(10^-5, 1-10^-5))
			g1W <- q <- mean(A)
			psi.init <- mean(Q[A==1,"Q1W"] - Q[A==1,"Q0W"])
		    for (k in 1:ngform){
		    	for (qq in 1:nrow(Qbounds)){
		    	g1W <- .bound(predict(glm(gformulas[[k]], family = "binomial"), type = "response"), gbounds)
			 	res.oneStep <- oneStepATT(Y = Y, A = A, Q = Q, g1W = g1W,
	  			  			depsilon = depsilon, max_iter = max(1000, 2/depsilon), gbounds = gbounds, Qbounds = Qbounds[qq,])
	    		res.iter <- try(calcATT.iter(Y = Y, A = A, Q = Q,  g1W = g1W, tol = 10^-5, maxIter = 500, gbounds = gbounds,
	    						Qbounds = Qbounds[qq,]))
	  			if(class(res.iter) == "try-error") {
	  				res.iter <- NULL
	  				res.iter$psi <- NA
	  			}
	  			res.psi[[ngform*(j-1)+k]][i, 1] <- psi.init
	  			res.psi[[ngform*(j-1)+k]][i, ((qq-1)*2 + 2):(3 + (qq-1)*2)] <-  c(res.oneStep$psi, res.iter$psi)
			}}
	        save(res.psi, Qformulas, gformulas, depsilon, Qbounds, file = paste0("Sim2_", n, ".Rdata"))
}}}
