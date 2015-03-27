#function to find a simple sinusoidal to fit data
FitSimpleSinusoidal <- function(ts, freq){	
#  freq -> frequency (modulus 1)
  mu <- mean(ts)
  alpha <- 0
  beta <- 0
  for (i in 1:length(ts)) {
    alpha <- alpha + ts[i]*cos(2*pi*i*omega)
    beta  <- beta  + ts[i]*sin(2*pi*i*omega)
  }
  alpha <- 2*alpha/length(ts)
  beta <- 2*beta/length(ts)
  out <- c(mu,alpha,beta)
  out
}

#function to generate data once the fit is done
GenerateTsSinusoidal <- function(freq, parm){
  ts <- numeric(N)
  for (i in 1:N) 
    ts[i] = parm[1] + parm[2]*cos(2*pi*i*freq) + parm[3]*sin(2*pi*i*freq)
  ts
}


#####################################################
# least squares method to fit a fourier polynomial  #
#####################################################

GenerateHarmMatrix <- function (ts.length, freq, n.harm) {
# auxiliary function to generate a fit of an harmonic functions to a 
# given time series #  freq is modulus 1
  v<-numeric()
  for(i in 1:ts.length){
    v[(i-1)*(2*n.harm+1)+1] =1.
    for(j in 1:n.harm){
      v[(i-1)*(2*n.harm+1)+2*j] <- cos(2*pi*freq*i*j);
      v[(i-1)*(2*n.harm+1)+2*j+1] <- sin(2*pi*freq*i*j);
    }
  }
  t(matrix(v,nrow=2*n.harm+1,ncol=ts.length))
}

GenerateHarmMatrixFast <- function(ts.length, freq, n.harm) {
# as the function above but with a faster algorithm
  v<-numeric(ts.length)
  for(i in 1:ts.length)
    v[(i-1)*(2*n.harm+1)+1] =1.;
  
  c1 <- cos(2*pi*freq);
  s1 <- sin(2*pi*freq);
  cj <- 1.;
  sj <- 0.;
  for(j in 1:n.harm){
    aux <- cj;
    cj  <- cj*c1  - sj*s1;
    sj  <- aux*s1 + sj*c1;
    cij <- 1.;
    sij <- 0.;
    for(i in 1:ts.length){
      aux <- cij;
      cij <- cij*cj - sij*sj;
      sij <- aux*sj + sij*cj;
      v[(i-1)*(2*n.harm+1)+2*j]   <- cij;
      v[(i-1)*(2*n.harm+1)+2*j+1] <- sij;
    }
  }
  t(matrix(v,nrow=2*n.harm+1,ncol=ts.length))
}

LeastSquaresFitHarm <- function(ts, freq, n.harm){
# computes the coeficients of the harmonic function that best fits ts
   m   <- GenerateHarmMatrixFast(length(ts), freq, n.harm)
   aux <- lsfit(m,ts)
   aux$coef
}

SimulateHarmFunc <- function(n.ts, n.harm, omega, coef){
  ts <- numeric(n.ts)
  for(i in 1:n.ts){
    ts[i] = coef[1];
    for(j in 1:n.harm){
      ts[i]=ts[i]+coef[2*j]*cos(2*pi*i*j*omega) + coef[2*j+1]*sin(2*pi*i*j*omega);
    }
  }
  ts
}

SubstractPeriodSignal <- function(ts, freq, n.harm){ 
# function to extract periodic oscilations of a given signal.
# freq -> vector or frequencies 
# n.harm -> number of harmonics for each frequency
# returns a list with three components 
#    1st vector with the periodic part of the original signal
#    2nd vector with the stochastic part of the original signal
#    3rd a list of vectors with the corresponding Fourrier coeficients
  if(length(freq) != length(n.harm)) 
    stop("SubstractPeriodSignal: freq and n.harm must be same length!")

  output <- list(periodic= vector(length(ts), mode="numeric"),
                 residual= vector(length(ts), mode="numeric"), 
                 f.coef  = list())
 
  
  for (i in 1:length(freq)){ 
    output$coef.fit[[i]] <- LeastSquaresFitHarm(ts - output$periodic, 
                                                freq[i], n.harm[i])
    output$periodic <- output$periodic + 
                       SimulateHarmFunc(length(ts), n.harm[i], freq[i], 
                                        output$coef.fit[[i]])
  }
  output$residual <- ts - output$periodic
  return(output)
} 
