Estimate_Gamma_GAS <- function(vY,Dates,iMaxIter) {
  require(DEoptim) # <- global optimizer
  ## Funcion "Estimate_Gamma_GAS":
  ## Input:   vY: Vector of observations, to which the Gamma-GAS model should be fitted
  ## Output:  vPar: estimated parameters
  #vPar: dOmega  = vPar[1]
  #dAlpha = vPar[2]
  #dBeta  = vPar[3]
  #dA     = vPar[4]
  ##         dLLK: Log Likelihood evaluated at the estimated (optimized) parameter values
  ##         ABIC: Average value of the Bayes Information Criterion, evaluated at optimum
  ##         FilteredMean: The filtered values of mu_t evaluated at optimum
  # maximize the likelihood using a global optimizer
  optimizer = DEoptim(Extract_LLH, lower = c(-0.5, 1e-3, 0.01, 0.1),
                      upper = c(0.5, 1.5, 0.999, 300), DEoptim.control(itermax=iMaxIter), vY = vY)
  vPar = optimizer$optim$bestmem
  dLLK = -optimizer$optim$bestval

  if (optimizer$optim$bestval != 1e5){
    cat(paste0("Optimizer converged: Gamma-GAS model fit complete: with LLH  ", dLLK))
  } else {
    cat(paste0("Optimizer not converged, try a different sample or series"))
  }
  # run the filter
  FilteredValues = GASGAMMA_Filter(vPar, vY)[["vMu"]]

  iT = length(vY)

  #compute the average BIC: having 4 paramaters
  ABIC = (log(iT) * 4 - 2 * dLLK)/iT

  lOut = list()
  class(lOut) <- "GASModel"

  lOut[["vPar"]] = vPar
  lOut[["dLLK"]] = dLLK
  lOut[["ABIC"]] = ABIC
  lOut[["FilteredMean"]] = FilteredValues
  lOut[["Observations"]] = vY
  lOut[["Dates"]] = Dates
  return(lOut)
}
