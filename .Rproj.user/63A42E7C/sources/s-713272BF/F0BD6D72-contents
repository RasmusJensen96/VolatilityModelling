Extract_LLH_MEM <- function(vPar, vY) {
  ## Auxiliary function:
  ## Inputs = inputs(MEM_Filter)
  ## Output: DNLL: negative log likelihood evaluated at vPar, on vY observations
  Filter = MEM_Filter(vPar, vY)

  dNLL = -Filter[["dLLK"]]

  if (is.infinite(dNLL) | is.nan(dNLL)) {
    dNLL = 1e5
  }

  return(dNLL)
}
