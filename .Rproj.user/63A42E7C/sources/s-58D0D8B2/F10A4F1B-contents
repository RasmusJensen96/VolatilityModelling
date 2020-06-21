# Function to extracting the Gamma-GAS model likelihood
Extract_LLH <- function(vPar, vY) {
  ## Auxiliary function extracting the likelihood of the GARCH model.
  ## Inputs: vPar, parameters c(omega, alpha, beta, a) where a is the scale of the gamma distribution of Engle & Gallo (2006),
  ##        vY, vector of observations
  ## Output: Total log-likelihood for the model given parameters and observations

  Filter = GASGAMMA_Filter(vPar, vY)

  dNLL = -Filter[["dLLK"]]

  if (is.infinite(dNLL) | is.nan(dNLL)) {
    dNLL = 1e5
  }

  return(dNLL)
}
