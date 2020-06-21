# Gamma-GAS filter
GASGAMMA_Filter <- function(vPar, vY) {
  ## Filtering the values of the filter
  ## Input: Parameters of the model (vPar) c(omega, alpha, beta, a)
  ##        Values to filter (vY)
  ## Output: filtered Mu_t
  ##         LLK with chosen parameters
  iT = length(vY)
  vMu = numeric(iT)
  vMuTilde = numeric(iT)

  dOmega = vPar[1]
  dAlpha = vPar[2]
  dBeta  = vPar[3]
  dA     = vPar[4]
  #initialize at the unconditional value
  vMuTilde[1] = dOmega/(1 - dBeta)
  vMu[1] = exp(vMuTilde[1])
  #The first loglikelihood value
  dLLK = Density_Gamma(vY[1], vMu[1], dA,cLog = TRUE)
  for(t in 2:iT) {
    vMuTilde[t] = dOmega + dBeta * vMuTilde[t - 1] + dAlpha * ( (dA / (vMuTilde[t-1] ^2))^(-1/2) * (dA*(vY[t-1]-vMu[t-1]))/(vMu[t-1])) # Updating the conditional mean (filtering)
    vMu[t] = exp(vMuTilde[t])
    dLLK = dLLK + Density_Gamma(vY[t], vMu[t], dA,cLog = TRUE) ## Calculating the total likelihood for selected model-parameters
  }
  lOut = list("dLLK" = dLLK,
              "vMu" = vMu)
  return(lOut)
}
