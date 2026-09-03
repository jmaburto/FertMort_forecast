## Function file to fit discrete Lee-Carter with Newton-Raphson
## Modified by Giancarlo Camarda 25-11-2009 for teaching purposes
## (EDSD 220)


###################################################
## Update Alpha
Update.alpha <- function(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit){
  difD <- Dth - D.fit
  Alpha <- Alpha + difD %*% One / (D.fit %*% One)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Alpha = Alpha, D.fit = D.fit)
}

###################################################
## Update Beta
Update.beta <- function(Alpha, Beta, Kappa,
                        One, Dth, Exp, D.fit){
  difD <- Dth - D.fit
  Kappa2 <- Kappa * Kappa
  Beta <- Beta + difD %*% Kappa / (D.fit %*% Kappa2)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Beta = Beta, D.fit = D.fit)
}

###################################################
## Update Kappa
Update.kappa <- function(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit){
  difD <- Dth - D.fit
  Beta2 <- Beta * Beta
  Kappa <- Kappa + t(difD) %*% Beta / (t(D.fit) %*% Beta2)
  Kappa <- Kappa - mean(Kappa)
  Kappa <- Kappa / sqrt(sum(Kappa * Kappa))
  Kappa <- matrix(Kappa, ncol = 1)
  Eta <- Alpha %*% t(One) + Beta %*% t(Kappa)
  D.fit <- Exp * exp(Eta)
  list(Kappa = Kappa, D.fit = D.fit)
}

###################################################
##     MAIN FUNCTION for fitting the Lee-Carter
LCpoi <- function(Dth, Exp){
  # dimensions
  m <- nrow(Dth)
  n <- ncol(Dth)
  # Initialise
  One <- matrix(1, nrow=n, ncol = 1)    
  Fit.init <- log((Dth + 1)/(Exp + 2))
  Alpha <- Fit.init %*% One / n
  Beta <- matrix(1 * Alpha, ncol = 1)
  sum.Beta <- sum(Beta) 
  Beta <- Beta / sum.Beta
  Kappa <- matrix(seq(n, 1, by = -1), nrow = n, ncol = 1)
  Kappa <- Kappa - mean(Kappa)
  Kappa <- Kappa / sqrt(sum(Kappa * Kappa))
  Kappa <- Kappa * sum.Beta
  # Iteration
  D.fit <- Exp * exp(Fit.init)
  for (iter in 1:50){
    Alpha.old <- Alpha
    Beta.old <- Beta
    Kappa.old <- Kappa
    #
    temp <- Update.alpha(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit)
    D.fit <- temp$D.fit
    Alpha <- temp$Alpha
    #
    temp <- Update.beta(Alpha, Beta, Kappa,
                        One, Dth, Exp, D.fit)
    D.fit <- temp$D.fit
    Beta <- temp$Beta
    #
    temp <- Update.kappa(Alpha, Beta, Kappa,
                         One, Dth, Exp, D.fit)
    D.fit <- temp$D.fit
    Kappa <- temp$Kappa
    crit <- max(max(abs(Alpha - Alpha.old)),
                max(abs(Beta - Beta.old)),
                max(abs(Kappa - Kappa.old)))
    if(crit < 1e-04) break
  }
  # constraints
  sum.Beta <- sum(Beta) 
  Beta <- Beta / sum.Beta
  Kappa <- Kappa * sum.Beta
  # output
  out <- list(Alpha=Alpha, Beta=Beta, Kappa=Kappa)
}
