GenStudT <- function(x, dLambda = 0, dNu = 20){
  if (dNu <= 2){
    stop("\n Nu can't be smaller nor equal to 2")
  }

  if (abs(dLambda)>1){
    stop("Lambda must be smaller than or equal to 1 in modulus")
  }

  fooc = gamma((dNu+1)/2)/(sqrt(pi*(dNu-2))*gamma(dNu/2))
  fooa = 4*dLambda*fooc*(dNu-2)/(dNu-1)
  foob = sqrt(1+3*dLambda^2-fooa^2)

  return(ifelse(x >= -fooa/foob,     # p. 1702 Journal of Economic Dynamics & Control 27
                foob*fooc*(1+1/(dNu-2)*((foob*x+fooa)/(1+dLambda))^2)^(-(dNu+1)/2),
                foob*fooc*(1+1/(dNu-2)*((foob*x+fooa)/(1-dLambda))^2)^(-(dNu+1)/2)))
}


