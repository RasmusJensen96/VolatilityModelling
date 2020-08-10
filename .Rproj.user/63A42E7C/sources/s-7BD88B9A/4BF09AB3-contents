GARCHHansenGenT_LLH = function(vPar,vY, returnresult = F){

if (vPar[10]+vPar[12] > 1 || vPar[11]+vPar[12]>1){
    return(-1e9)
}

yplus  = max(yt,0)
yminus = max(-yt,0)

vEta = numeric(length(yt))
vLambda = vEta
vSigma2 = vEta
s
vEta[1] = vPar[13]
vLambda[1] = vPar[14]
vSigma2[1] = vPar[15]


for (t in 2:length(yt)){
  vEta[t] = vPar[1] + vPar[2] * yplus[t-1] + vPar[3] * yminus[t-1] + vPar[4] * vEta[t-1]
  vLambda[t] = vPar[5] + vPar[6] * yplus[t-1] + vPar[7] * yminus[t-1] + vPar[8] * vLambda[t-1]
  vSigma2[t] = vPar[9] + vPar[10] * yplus[t-1]^2 + vPar[11] * yminus[t-1]^2 + vPar[12]*vSigma2[t-1]
   if (vEta[t] <= 2 || abs(vLambda[t]) > 1){
     return(-1e9)
   }
}

if (returnresult == T){
  return(list(vSigma = sqrt(vSigma2),
              vLambda = vLambda,
              vEta = vEta,
              LLH = -(sum(dnorm(yt, 0, sqrt(vSigma2), log = T)) + sum(log(GenStudT(yt, vLambda, vEta))))))
}else{
  return(-(sum(dnorm(yt, 0, sqrt(vSigma2), log = T)) + sum(log(GenStudT(yt, vLambda, vEta)))))
}
}
