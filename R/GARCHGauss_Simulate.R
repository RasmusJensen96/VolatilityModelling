GARCHGauss_Simulate <- function(iT, dOmega, dAlpha, dBeta) {

    ## initialize the vector of simulated returns and variances
    vY = numeric(iT)
    vSigma2 = numeric(iT)

    ## initialize the variance at time t = 1 with its unconditional value
    vSigma2[1] = dOmega/(1.0 - dAlpha - dBeta)
    ## sample the first observations
    vY[1] = rnorm(1, mean = 0, sd = sqrt(vSigma2[1]))

    ##loop over iT. We start from t = 2 since t = 1 has already been sampled
    for (t in 2:iT) {
      #update the volatility
      vSigma2[t] = dOmega + dAlpha * vY[t - 1]^2 + dBeta * vSigma2[t - 1]
      #sample a new observarions
      vY[t] = rnorm(1, mean = 0, sd = sqrt(vSigma2[t]))
    }

    ## we return a list with two components: the sampled returns and the volatility
    lOut = list()
    lOut[["vY"]] = vY
    lOut[["vSigma2"]] = vSigma2

    ## output lOut
    return(lOut)
}
