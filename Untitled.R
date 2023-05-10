  #setwd("/Users/filo/Dropbox/Macrometrics/US-Monetary-Policy-and-the-Global-Financial-Cycle-Replication")
  library(readxl)
  library(ggplot2)
  library(rmarkdown)
  library(matrixcalc)
  library(mvtnorm)
  library(parallel)
  library(fredr)
  library(sovereign)
  library(MASS)
  #rm(list=ls())
  fredr_set_key("3330a1a80c13478e9e927a10d2c65c04")
  url <- "http://silviamirandaagrippino.com/s/DFM-Blocks.zip"
  download.file(url, "MAR(2020)")
  unzip("MAR(2020)")
  data <- read_excel("GFC_VARdata+WEB.xlsx", skip=1, .name_repair = "unique_quiet")
  inst <- read_excel("GFC_VARdata+WEB.xlsx", sheet = 2, skip=1, .name_repair = "unique_quiet")

  # Create a list of variable names for each plot
  #order industrial production that affects price deflator and 1y treasury rate, and also BIS real EER. then the global factor and global risk aversion, global real economic activity ex us, and global domestic credit. then global inflows all sectors, gz credit spread and leverage of 1)us brokers and dealers;2)leverage EU bankers;3)leverage US banks;4)leverage EU banks
  set.seed(0) 
  #monthly data
  fed <- fredr(series_id = "FEDFUNDS", 
                           observation_start = as.Date("1980-01-01"), 
                           frequency = "m",
                           aggregation_method = "average")

  fed        = ts((fed[,3]), start=c(1980,1), frequency=12) #FED federal funds rate
  instrument = ts((inst$FF4), start=c(1980,1), frequency=12) #FFF4 instrument
  pce        = ts((data$PCEPI), start=c(1980,1), frequency=12) #PCE Deflator; 
  #dgs1       = ts((data$DGS1), start=c(1980,1), frequency=12) #1 Year Treasury Rate
  #bis        = ts((data$BISREER), start=c(1980,1), frequency=12) #BIS real EER
  globalf    = ts((data$GLOBALF), start=c(1980,1), frequency=12) #Global Factor
  #globalra   = ts((data$GLOBALRA), start=c(1980,1), frequency=12) #Global Risk Aversion
  greaexus   = ts((data$GREAEXUS), start=c(1980,1), frequency=12) #Global Real Economic Activity Ex US
  indpro     = ts((data$INDPRO), start=c(1980,1), frequency=12) #industrial production
  #glbcredit    = ts((data$GLBCREDIT), start=c(1980,1), frequency=12) #Global Domestic Credit
  glbinflows    = ts((data$GLBINFLALL), start=c(1980,1), frequency=12) #Global Inflows All Sectors
  #usbdlev    = ts((data$USBDLEV), start=c(1980,1), frequency=12) #Leverage US Brokers and Dealers
  eubdlev    = ts((data$EUBDLEV), start=c(1980,1), frequency=12) #Leverage EU Global Banks
  #usbanksl   = ts((data$USBANKSL), start=c(1980,1), frequency=12) #Leverage US Banks
  #eubanksl   = ts((data$EUBANKSL), start=c(1980,1), frequency=12) #Leverage EU Banks

  #create the bigy matrix with my data
  y = cbind(instrument,
    fed,
    pce, 
    #dgs1, 
    #bis, 
    globalf, 
    #globalra, 
    greaexus, 
    #glbcredit, 
    indpro,
    glbinflows, 
    #usbdlev, 
    eubdlev 
    #usbanksl, 
    #eubanksl)
  )
  ############################################################
  #useful functions
  # setpriors <- function(N,K,kappa.3){
  #   A.prior     = matrix(0,K,N)
  # A.prior[2:(N+1),] <- kappa.3*diag(N)
  # V.prior     = (diag(c(10,1*((1:p)^(-2))%x%rep(1,N)))) #10 is kappa.2, 1 is kappa.1
  # S.prior     = diag(N)
  # nu.prior    = N+1
  # return(list("A.prior"=A.prior,
  #             "V.prior"=V.prior,
  #             "S.prior"=S.prior,
  #             "nu.prior"=nu.prior))
  # }
  ############################################################



  #start of analysis
  ############################################################
  #y = na.omit(y)
  #y <- window(y, start=c(1990,2), end=c(2010,12)) #subset in our time window: 1990,1 until 2011,1
  
  # setup function for my analysis
  basic.model <- function(bigy, p, start, end){
  
  #bigy=y
  bigy <- window(bigy, start=start, end=end)
  #bigy = yart
  #p=1
  S=5
  #kappa.1=1
  #kappa.2=10
  kappa.3=0.95
  #start=c(1990,2)
  #end=c(2010,12)
  #start=c()
  #end=c()
  ############################################################
  N       = ncol(bigy)
  K       = 1+N*p
  S.burnin= 100
  ############################################################
  Y       = bigy[(p+1):nrow(bigy),]
  X       = matrix(1,nrow(Y),1)
  for (i in 1:p){
    X     = cbind(X,bigy[(p+1):nrow(bigy)-i,])
  }
  # set the priors
  ############################################################
  # priors  <- setpriors(N,K,kappa.3)
  # A.prior  = priors$A.prior 
  # V.prior  = priors$V.prior
  # S.prior  = priors$S.prior
  # nu.prior = priors$nu.prior 
  A.prior     = matrix(0,K,N)
  A.prior[2:(N+1),] <- kappa.3*diag(N)
  V.prior     = (diag(c(10,1*((1:p)^(-2))%x%rep(1,N)))) #10 is kappa.2, 1 is kappa.1
  S.prior     = diag(N)
  nu.prior    = N+1
  
  # normal-inverse Wishart posterior parameters
  ############################################################
  V.bar.inv   = t(X)%*%X + diag(1/diag(V.prior)) #X'X+diag(V^-1)
  V.bar       = solve(V.bar.inv) #inv(X'X+diag(V^-1))
  A.bar       = V.bar%*%(t(X)%*%Y + diag(1/diag(V.prior))%*%A.prior) #V.bar(X'X+diag((V.prior)^-1(A.prior)))
  nu.bar      = nrow(Y) + nu.prior 
  S.bar       = S.prior + t(Y)%*%Y + t(A.prior)%*%diag(1/diag(V.prior))%*%A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
  #S.prior+Y'Y+(A.prior)'*diag(diag((V.prior)^-1))*A.prior-A.bar'*(X'X+diag(V^-1))*A.bar
  S.bar.inv   = solve(S.bar)
  
  # posterior draws
  ############################################################
  Sigma.posterior   = rWishart(S, df=nu.bar, Sigma=S.bar.inv)
  Sigma.posterior   = apply(Sigma.posterior,3,solve)
  Sigma.posterior   = array(Sigma.posterior,c(N,N,S))
  A.posterior       = array(rnorm(prod(c(dim(A.bar),S))),c(dim(A.bar),S)) #3 dimensional, S repetitions of rnorm
  B.posterior       = array(NA,c(N,N,S))
  #Bplus.posterior   = array(NA,c(N,K,S))
  L                 = t(chol(V.bar))
  for (s in 1:S){
    cholSigma.s     = chol(Sigma.posterior[,,s])
    B.posterior[,,s]= t(cholSigma.s)
    A.posterior[,,s]= A.bar + L%*%A.posterior[,,s]%*%cholSigma.s
}
  
  results = list("B" = B.posterior, "A" = A.posterior )
  return(results)
  }
  

  #run the function with my data
  set.seed(23)
  basic.results <- basic.model(y, 12, c(1990,2), c(2010,12))


  #artificial data with 1000 observations from a bi-variate Gaussian random walk
  #VAR model with N=2, p=1, and a constant term
  set.seed(0)   # set seed for reproducibility
  # simulate 1000 samples from the bivariate Gaussian random walk process with cumsum(rnorm())
  y1 <- ts(cumsum(rnorm(1000, 0, sd=1)))
  y2 <- ts(cumsum(rnorm(1000, 0, sd=1)))
  y.art <- cbind(y1,y2)
  plot(y.art)
  # run the function of the basic model and check the posterior parameters B0 and B1
  basic.artificial <- basic.model(y.art, 1, c(), c())

  ############################################################
  ############################################################
  #Extendend model computation
  ############################################################
  ############################################################
  
  extended.model <- function(bigy, p, S1, S2, kappa.3, start, end){
  #bigy <- window(bigy, start=start, end=end)
   #y <- window(y, start=c(1990,2), end=c(2010,12))
   #bigy=y
   #p=12
   #S1=5 #repetitions to discard
   #S2=5 #repetitions to keep
   S=S1+S2
   kappa.1=1
   kappa.2=10
   #kappa.3=0.95
  ############################################################
  
  Y       = bigy[(p+1):nrow(bigy),]
  X       = matrix(1,nrow(Y),1)
  for (i in 1:p){
    X     = cbind(X,bigy[(p+1):nrow(bigy)-i,])
  }
  
  N       = ncol(bigy)
  K       = 1+N*p
  S.burnin= 100
  # set the priors
  ############################################################
  # priors <- setpriors(N,K,kappa.3)
  # priors$A.prior  = A.prior
  # priors$V.prior  = V.prior
  # priors$S.prior  = S.prior
  A.prior     = matrix(0,K,N)
  A.prior[2:(N+1),] <- kappa.3*diag(N)
  V.prior     = (diag(c(10,1*((1:p)^(-2))%x%rep(1,N)))) #10 is kappa.2, 1 is kappa.1
  S.prior     = diag(N)
  nu.prior    = N+1
  s.sigma.prior.k = diag(N)   #for kappa.sigma and kappa.a
  nu.prior.k      = N + 1      #same, to be used in ig2 and rgamma
  #priors$nu.prior = nu.prior
  s.kappa         = 2
  s.sigma.prior   = 2
  a.sigma.prior   = 2
  
  k.sigma = matrix(0, S, 1)
  k.a     = matrix(0, S, 1)
  
  k.sigma[1] <- 1
  k.a[1] <- 1
  for (s in 1:(S1+S2)) {
  # normal-inverse Wishart, IG2, and Gamma posterior parameters
  ############################################################
  V.bar.inv   = t(X)%*%X + solve(k.a[1]*V.prior)
  V.bar       = solve(V.bar.inv)
  A.bar       = V.bar%*%(t(X)%*%Y + diag(1/diag(k.a[1]*V.prior))%*%A.prior)
  nu.bar      = nrow(Y) + nu.prior
  nu.bar.k    = nu.prior.k + N*K
  S.bar       = k.sigma[1]*diag(N) + t(Y)%*%Y + t(A.prior)%*%diag(1/diag(k.a[1]*V.prior))%*%A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
  S.bar.inv   = solve(S.bar)
  a.sigma.bar = (N*V.prior)/2+a.sigma.prior
  
  # posterior draws
  ############################################################
  Sigma.posterior     = rWishart(S, df=nu.bar, Sigma=S.bar.inv)
  Sigma.posterior     = apply(Sigma.posterior,3,solve)
  Sigma.posterior     = array(Sigma.posterior,c(N,N,S))
  sigma.post.temp     = Sigma.posterior[,,s]
  A.posterior         = array(rnorm(prod(c(dim(A.bar),S))),c(dim(A.bar),S))
  B.posterior         = array(NA,c(N,N,S))
  Bplus.posterior     = array(NA,c(N,K,S))
  L                   = t(chol(V.bar))
    cholSigma.s       = chol(Sigma.posterior[,,s])
    B.posterior[,,s]  = t(cholSigma.s)
    A.posterior[,,s]  = A.bar + L%*%A.posterior[,,s]%*%cholSigma.s
    #k.a from IG2
    s.k.a.bar     <- s.kappa + sum(diag( sigma.post.temp %*% t(A.posterior[,,s]-A.prior)%*% diag(1/diag(V.prior)) %*% (A.posterior[,,s]-A.prior)))
    k.a[s+1]      <- s.k.a.bar / rchisq(1, df=nu.bar.k)
    #k.sigma from gamma distribution
    s.sigma.bar   <- (2 + sum(diag(sigma.post.temp*s.sigma.prior))^-1+(s.sigma.prior)^-1)^-1
    k.sigma[s+1]  <- rgamma(1, shape = a.sigma.bar, scale = s.sigma.bar)
    
    
  }
  
  results = list("B" = B.posterior[,,(S1+1):S2], "A" = A.posterior[,,(S1+1):S2] )
  return(results)
  }
  
  y                   <- window(y, start=c(1990,2), end=c(2010,12))
  extended.result     <- extended.model(y, 12, 5, 5, 0.95, c(1990,2), c(2010,12))
  extended.artificial <- extended.model(y.art, 1, 500, 5000, 0.95, c(), c())
  