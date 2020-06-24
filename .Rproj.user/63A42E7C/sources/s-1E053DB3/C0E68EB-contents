library("quantmod")
##                                Dow Jones Industrial                            ##
####################################################################################
##                                  Retrieve returns                              ##
####################################################################################
DJI = getSymbols("^DJI", auto.assign = FALSE, from = "2007-12-01", to = "2010-06-01")
rets = diff(log(DJI$DJI.Adjusted))*100
rets = as.vector(rets)[-1]
Dates = index(DJI)[-1]
rets = rets^2
####################################################################################
##                                Estimate the model                              ##
####################################################################################
fit = Estimate_Gamma_GAS(rets,Dates, 100)
fit$Name  = "DJI"
plot(fit)
##                                        VIX                                     ##
####################################################################################
##                                  Retrieve returns                              ##
####################################################################################
VIX = getSymbols("^VIX", auto.assign = FALSE, from = "2007-12-01", to = "2012-06-01")
VIX = VIX$VIX.Adjusted
Dates = index(VIX)
VIX = as.numeric(VIX)
####################################################################################
##                                Estimate the model                              ##
####################################################################################
fit = Estimate_Gamma_GAS(VIX,Dates, 100)
fit$Name  = "VIX"
plot(fit)
