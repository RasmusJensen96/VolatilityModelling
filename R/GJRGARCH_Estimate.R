Estimate_GJRGARCH <- function(vY,Dates) {
  ## Function fitting or estimatinig a GARCH(1,1)-model on a certain vector of data
  ## Input:  vY:      Data-vector
  ## output: vPar:    parameters resulting in the maximum likeliihood given data c(omega, alpha, beta)
  ##         dLLK:    Value of the likelihood at the Maximum likelihood given data and estimated parameters
  ##         BIC:     Average BIC value at the (numerical) maximum likelihood estimate
  ##         Filter: Vector of conditional volatility

  # Initialize the values at some arbitrary values, except omega initialized at the unconditional value.
  require(Rsolnp)
  dAlpha = 0.05
  dBeta  = 0.85
  dAlphaplus = 0.05
  dOmega = var(vY) * (1.0 - dAlpha - dBeta-dAlphaplus)
  vPar = c(dOmega, dAlpha, dBeta, dAlphaplus)
  ## Optimize using the non-linear solnp, allowing for inequality constraints such as the one contained in
  ## StationarityConstr ensuring htat the estimated model is covariance (weakly) stationary
  # if (dist == "Gaussian"){
    optimizer = solnp(vPar, fun = LLH_GJRGARCH11, vY = as.vector(vY),
                      ineqfun = StationarityConstrGJR,
                      ineqLB  = 0,
                      ineqUB = 0.99999,
                      LB = c(0.00001, 0.0001, 0.0001,0.0001), UB = c(10.0, 0.999, 0.999,0.999))
    vPar = optimizer$pars
    vSigma2 = GJRGARCH11_Filter(vY, vPar[1], vPar[2], vPar[3], vPar[4])$vSigma2
  # } else if (dist == "Student-t"){
  #   dNu = 2.99
  #   optimizer = solnp(c(vPar, dNu), fun = LLH_GARCH11T, vY = vY,
  #                     ineqfun = StationarityConstrT,
  #                     ineqLB  = 0.0001,
  #                     ineqUB = 0.99999,
  #                     LB = c(0.00001, 0.0001, 0.0001, 2.1), UB = c(10.0, 0.999, 0.999, 1000))
  #   vPar = optimizer$pars
  #   vSigma2 = GARCH11T_Filter(vY, vPar[1], vPar[2], vPar[3], vPar[4])$vSigma2
  # } else {
  #   stop("Please choose either 'Gaussian' or 'Student-t'")
  # }

  ## extract the likelihood computed at its maximum
  dLLK = -tail(optimizer$values, 1)
  ## compute filtered Variance
  ## Compute the Average BIC
  ABIC = (-2 * dLLK + log(length(vY)) * length(vPar))/length(vY)
  ## return a list with estimated parameters, likelihood value and BIC
  lOut = list(vPar = vPar,
              dLLK = dLLK,
              BIC = ABIC,
              Filtered = sqrt(vSigma2),
              Observations = vY,
              Model = "GJRGARCH",
              Dist = "Gaussian",
              Dates = Dates,
              Targetfilter = "Conditional volatility")
  class(lOut) <- "FitModel"
  return(lOut)
}

