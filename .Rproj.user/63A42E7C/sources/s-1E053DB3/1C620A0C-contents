Estimate_MEM <- function(vY, Dates, iMaxIter) {
  ## Function fitting a Multiplicative Error Model to data:
  ## inputs: vY: Data to which the ME-Model is fitted
  ## Outputs: vPar: Parameters of the model evaluated at the optimum
  ##                                      vPar[1] = kappa
  ##                                      vPar[2] = eta
  ##                                      vPar[3] = phi
  ##                                      vPar[4] = a, (Scale of the gamma distribution)
  ##          dLLK: Log-Likelihood evaluated at the optiimum
  ##          ABIC: Average Bayes Information Criterion evaluated at the optimum
  ##          FilteredMean: Vector of filtered means (fitted values) at the parameters yielding the optimum LLH
  require(DEoptim)
  # maximize the likelihood using a global optimizer
  optimizer = DEoptim(Extract_LLH_MEM, lower = c(0.1, 0.01,0.01,0.1),
                      upper = c(10, 0.99,0.99,300), DEoptim.control(itermax=iMaxIter), vY = vY)
  vPar = optimizer$optim$bestmem
  dLLK = -optimizer$optim$bestval
  # Retrieve the fitted values
  FilteredValues = MEM_Filter(vPar, vY)[["vMu"]]
  #compute the average BIC: having 4 parameters
  ABIC = (log(length(vY)) * 4 - 2 * dLLK)/length(vY)

  lOut = list(vPar = vPar,
              dLLK= dLLK,
              ABIC = ABIC,
              Filtered = FilteredValues,
              Observations = vY,
              Dates = Dates,
              Model = "MEM",
              Dist = "Gamma")
  class(lOut) <- "FitModel"
  return(lOut)
}

MEM_Filter <- function(vPar, vY) {
  ## Filtering the values of the filter
  ## Input: Parameters of the model (vPar)
  ##                                      vPar[1] = kappa
  ##                                      vPar[2] = eta
  ##                                      vPar[3] = phi
  ##                                      vPar[4] = a, (Scale of the gamma distribution)
  ##        Values to filter (vY)
  ##        Output: dLLK:         likelihood evaluated at the data, vY and the parameters vPar
  ##        FilteredMean: Filtered mu_[t] evaluated at the parameters vPar
  iT = length(vY)
  vMu = numeric(iT)
  dKappa = vPar[1]
  dEta   = vPar[2]
  dPhi   = vPar[3]
  dA     = vPar[4]
  #unconditional value
  vMu[1] = dKappa/(1-dEta)
  dLLK = Density_Gamma(vY[1], vMu[1], dA,cLog = TRUE)
  for(t in 2:iT) {
    vMu[t] = dKappa + dEta * vY[t-1] + dPhi * vMu[t-1]
    #The accumulated loglikelihood
    dLLK = dLLK + Density_Gamma(vY[t], vMu[t], dA,cLog = TRUE)
  }
  lOut = list("dLLK" = dLLK,
              "vMu" = vMu)
  return(lOut)
}
