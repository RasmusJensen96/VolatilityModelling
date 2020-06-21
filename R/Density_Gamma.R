Density_Gamma <- function(dY, dMu, dA,cLog=T){
  ## Function returning the density (PDF) of the gamma PDF, as parameterized in the exercise:
  ## Input: Point to evaluate (dY)
  ##        Mean (dMu)
  ##        Scale (dA)
  ##        cLog <- logical indicating whether or not the output should be loggged or not
  ## Output: Density at the point.
  pdf = 1/gamma(dA) * dA^dA * dY^(dA-1) * dMu ^ (-dA) * exp( - dA * dY/dMu )
  if (cLog == T){
    pdf = log(pdf)
  }
  return(pdf)
}
