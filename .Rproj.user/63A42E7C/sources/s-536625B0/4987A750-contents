GARCH11T_Filter <- function(vY, dOmega, dAlpha, dBeta, dNu) {
  # Function filtering and calculating LLH of a GARCH(1,1) on supplied data, and GARCH parameters:
  # Inputs: vY: vector of Data
  #         dOmega: Intercept of the GARCH-Equation
  #         dAlpha: ARCH parameter
  #         dBeta: GARCH: parameter
  # Outputs: dLLK: log-likelihood given parameters and data
  #          vSigma2: Filtered Variance (squared volatility)
  # Initializing and preallocating vectors:
  iT = length(vY)
  vSigma2 = numeric(iT)
  ## Use trick to initialize the variance at an empirical estimate (using 10% of the observations)
  vSigma2[1] = var(vY[1:floor(iT * 0.1)])
  # Init. LLK
  ## Filtering and summing marginal Log likelihood contributions
  for (t in 2:iT) {
    vSigma2[t] = dOmega + dAlpha * vY[t - 1] ^ 2 + dBeta * vSigma2[t - 1]
  }

  dLLK = iT * log(gamma((dNu+1)/2)/(sqrt(pi*(dNu-2))*gamma(dNu/2)))-1/2*sum(log(vSigma2))-(dNu+1)/2*sum(log(1+(vY^2)/(vSigma2*(dNu-2))))

  lOut = list(dLLK = dLLK,
              vSigma2 = vSigma2)
  return(lOut)
}

LLH_GARCH11T <- function(vPar, vY) {
  # Auxiliary function to extract the likelihood of a given parameterization:
  # inputs: vPar: vPar[1] = omega
  #               vPar[2] = alpha
  #               vPar[3] = beta
  # Extract likelihood from the filter (objective of the maximization proc)
  dOmega = vPar[1]
  dAlpha = vPar[2]
  dBeta  = vPar[3]
  dNu    = vPar[4]
  dLLK = GARCH11T_Filter(vY, dOmega, dAlpha, dBeta, dNu)$dLLK
  return(-dLLK)
}
StationarityConstrT <- function(vPar, ...) {
  ## Auxiliary function used as a constraint in the solnp proc in estimate GARCH:
  ## The function is used to determine if a given combination of ARCH and GARCH parameters in a GARCH(1,1) equation
  ## Results in a stationary GARCH(1,1) model. That is alpha+beta < 1
  dAlpha = vPar[2]
  dBeta  = vPar[3]
  return(dAlpha + dBeta)
}
