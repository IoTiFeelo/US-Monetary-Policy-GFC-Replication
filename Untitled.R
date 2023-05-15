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
  #library(gridExtra)
  library(patchwork)
  library(fredr)
  #rm(list=ls())
  fredr_set_key("3330a1a80c13478e9e927a10d2c65c04")
  url <- "http://silviamirandaagrippino.com/s/DFM-Blocks.zip"
  download.file(url, "MAR(2020)")
  unzip("MAR(2020)")
  data <- read_excel("GFC_VARdata+WEB.xlsx", skip=1, .name_repair = "unique_quiet")
  inst <- read_excel("GFC_VARdata+WEB.xlsx", sheet = 2, skip=1, .name_repair = "unique_quiet")
  gf <- read_excel("GFC_VARdata+WEB.xlsx", sheet= 4, skip=2, .name_repair = "unique_quiet", na = "blank")

  # Create a list of variable names for each plot
  #monthly data
  # fed <- fredr(series_id = "FEDFUNDS", 
  #                          observation_start = as.Date("1980-01-01"), 
  #                          frequency = "m",
  #                          aggregation_method = "average")
  
  
  # instrument,
  # fed,
  # pce,
  # #dgs1,
  # #bis,

  
  
  ggplot() + 
    geom_line(data = gf, aes(x = `...1`, y = `GLOBAL FACTOR 1975-2010`, linetype = "solid"), color ="blue", na.rm = TRUE, linewidth=1.2) + 
    geom_line(data = gf, aes(x = `...1`, y = `GLOBAL FACTOR 1990-2012`, linetype = "twodash"), color="blue4", na.rm = TRUE, linewidth=1.2) +
    geom_line(data = data, aes(x = `LABEL`, y = VIX, linetype = "dotted"), color="black", na.rm = TRUE, linewidth=1.2) +
    scale_linetype_manual(values=c("solid", "twodash", "dotted"), name="Global Factors", 
                          labels=c("VIX", "GF 1975-2010", "GF 1990-2012")) +
    guides(linetype = guide_legend(override.aes = list(color = c("black", "blue", "blue4")))) +
    labs(x ="", y="", title= "Global Factor for Risky Asset Prices") + theme(plot.title = element_text(size=20)) +
    theme(legend.text = element_text(size=20),legend.title = element_text(size=16, face="bold"),  legend.position = "bottom", plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"), axis.text = element_text(size = 14))
  
  ##########
  fed.gdp.plot <-  (ggplot(data=data) + 
    geom_line(aes(x = `LABEL`, y = `DGS1`, linetype = "solid"), na.rm = TRUE, color = "blue") 
    + geom_line(aes(x = `LABEL`, y = `GREAEXUS`, linetype = "dashed"), na.rm = TRUE, color = "black")+
    scale_linetype_manual(values=c("solid", "dashed"), name="", labels=c("Global Real Economic Activity Ex US", "1 Year Treasury Rate"), guide = guide_legend(override.aes = list(color = c("blue", "black")))) +
    guides(linetype = guide_legend(override.aes = list(color = c("black", "blue")))) +
    labs(x ="", y="", title= "Treasury Rate and Golbal Real Activity Index") +
      theme(legend.text = element_text(size=10),legend.title = element_text(size=16, face="bold"),  legend.position = "bottom", plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"), axis.text = element_text(size = 14)))
  #######
  pce.plot <- (ggplot(data=data) +
     geom_line(aes(x = `LABEL`, y = PCEPI, linetype = "solid"), na.rm = TRUE, color = "blue") +
     scale_linetype_manual(values=c("solid"), name="", labels=c("PCE Deflator"), guide = guide_legend(override.aes = list(color = c("blue")))) +
     guides(linetype = guide_legend(override.aes = list(color = c("blue")))) +
     labs(x ="", y="", title= "Price Deflator") +
       theme(legend.text = element_text(size=10),legend.title = element_text(size=16, face="bold"),  legend.position = "bottom", plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white"),
             plot.background = element_rect(fill = "white"), axis.text = element_text(size = 14)))
  
  
  ########
  leverage.plot <-  (ggplot(data=data) + 
    geom_line(aes(x = `LABEL`, y = EUBDLEV, linetype = "solid"), na.rm = TRUE, color = "blue", linewidth=1.1) + 
    geom_line(aes(x = `LABEL`, y = USBANKSL, linetype = "dashed"), na.rm = TRUE, color = "black", linewidth=1.1) + 
    geom_line(aes(x = `LABEL`, y = EUBANKSL, linetype = "twodash"), na.rm = TRUE, color = "purple", linewidth=1.1) +
    geom_line(aes(x = `LABEL`, y = USBDLEV, linetype = "dotted"), na.rm = TRUE, color = "violet", linewidth=1.1) +
    scale_linetype_manual(values=c("solid", "dashed", "twodash", "dotted"), name="", labels=c("Leverage EU Banks", "Leverage EU Global Banks", "Leverage US Banks", "Leverage US Brokers and Dealers"), guide = guide_legend(override.aes = list(color = c("blue", "black", "purple", "violet")))) +
    guides(linetype = guide_legend(override.aes = list(color = c("black", "blue", "purple", "violet")))) +
    labs(x ="", y="", title= "Leverage") +
      theme(legend.text = element_text(size=10),legend.title = element_text(size=16, face="bold"),  legend.position = "bottom", plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"), axis.text = element_text(size = 14)))
  ##########
  mixed.types.plot <-   (ggplot(data=data) + 
     geom_line(aes(x = `LABEL`, y = GLBINFLALL, linetype = "solid"), na.rm = TRUE, color = "blue", linewidth=1.1) + 
     geom_line(aes(x = `LABEL`, y = INDPRO, linetype = "dashed"), na.rm = TRUE, color = "black", linewidth=1.1) + 
     geom_line(aes(x = `LABEL`, y = GLBCREDIT, linetype = "twodash"), na.rm = TRUE, color = "purple", linewidth=1.1) +
       geom_line(aes(x = `LABEL`, y = BISREER, linetype = "dotted"), na.rm = TRUE, color = "violet", linewidth=1.1) +
     scale_linetype_manual(values=c("solid", "dashed", "twodash", "dotted"), name="", labels=c("Industrial Production","Global Inflows All Sectors", "Global Domestic Credit", "BIS real EER"), guide = guide_legend(override.aes = list(color = c("blue", "black", "purple", "violet")))) +
     guides(linetype = guide_legend(override.aes = list(color = c("black", "blue", "purple", "violet")))) +
     labs(x ="", y="", title= "Production, Inflows, Credit, BIS REER") +
       theme(legend.text = element_text(size=10),legend.title = element_text(size=16, face="bold"),  legend.position = "bottom", plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white"),
             plot.background = element_rect(fill = "white"), axis.text = element_text(size = 14)))
  
  ###############
  risk.types.plot <-  (ggplot(data=data) + 
    geom_line(aes(x = `LABEL`, y = GLOBALF, linetype = "solid"), na.rm = TRUE, color = "blue") + 
    geom_line(aes(x = `LABEL`, y = GLOBALRA, linetype = "dashed"), na.rm = TRUE, color = "black") +
    scale_linetype_manual(values=c("solid", "dashed"), name="", labels=c("Global Risk Aversion", "Global Factor"), guide = guide_legend(override.aes = list(color = c("black", "blue")))) +
    labs(x ="", y="", title= "Risk Aversion") +
      theme(legend.text = element_text(size=10),legend.title = element_text(size=16, face="bold"),  legend.position = "bottom", plot.title = element_text(hjust = 0.5), panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"), axis.text = element_text(size = 14)))
  ############################################################
  ############################################################
  ############################################################
  ############################################################
  
  #grid.arrange(global.credit, global.inflows, leverage.types, production.types, risk.types, nrow = 5)
  combinedplot.1 <- mixed.types.plot / leverage.plot 
    
  combinedplot.2 <- fed.gdp.plot + risk.types.plot + pce.plot
  
  combinedplot.1
  combinedplot.2
  

  #fed       = ts((fed[,3]), start=c(1980,1), frequency=12) #FED federal funds rate
  instrument = ts((inst$FF4), start=c(1980,1), frequency=12) #FFF4 instrument
  dgs1       = ts((data$DGS1), start=c(1980,1), frequency=12) #1 Year Treasury Rate
  pce        = ts((data$PCEPI), start=c(1980,1), frequency=12) #PCE Deflator; 
  bis        = ts((data$BISREER), start=c(1980,1), frequency=12) #BIS real EER
  globalf    = ts((data$GLOBALF), start=c(1980,1), frequency=12) #Global Factor
  globalra   = ts((data$GLOBALRA), start=c(1980,1), frequency=12) #Global Risk Aversion
  greaexus   = ts((data$GREAEXUS), start=c(1980,1), frequency=12) #Global Real Economic Activity Ex US
  indpro     = ts((data$INDPRO), start=c(1980,1), frequency=12) #industrial production
  glbcredit  = ts((data$GLBCREDIT), start=c(1980,1), frequency=12) #Global Domestic Credit
  glbinflows = ts((data$GLBINFLALL), start=c(1980,1), frequency=12) #Global Inflows All Sectors
  usbdlev    = ts((data$USBDLEV), start=c(1980,1), frequency=12) #Leverage US Brokers and Dealers
  eubdlev    = ts((data$EUBDLEV), start=c(1980,1), frequency=12) #Leverage EU Global Banks
  usbanksl   = ts((data$USBANKSL), start=c(1980,1), frequency=12) #Leverage US Banks
  eubanksl   = ts((data$EUBANKSL), start=c(1980,1), frequency=12) #Leverage EU Banks

  #create the bigy matrix with my data
  y = cbind(instrument,
    pce, 
    dgs1, 
    bis, 
    globalf, 
    globalra, 
    greaexus, 
    glbcredit, 
    indpro,
    glbinflows, 
    usbdlev, 
    eubdlev, 
    usbanksl, 
    eubanksl)
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
  # setup function for my analysis
  y           <- window(y, start=c(1990,2), end=c(2010,12))
  basic.model <- function(bigy, p, S, start, end){
  
  
  ############################################################
  N       = ncol(bigy)
  K       = 1+N*p
  ############################################################
  Y       = bigy[(p+1):nrow(bigy),]
  X       = matrix(1,nrow(Y),1)
  for (i in 1:p){
    X     = cbind(X,bigy[(p+1):nrow(bigy)-i,])
  }
  
  A.prior     = matrix(0,K,N)
  A.prior[2:(N+1),] <- diag(c(0,rep(1,N-1))) #0 for the instrument (stationary variable)
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
  basic.results <- basic.model(y, 12, 10000, c(1990,2), c(2010,12))


  #artificial data with 1000 observations from a bi-variate Gaussian random walk
  #VAR model with N=2, p=1, and a constant term
  set.seed(0)   # set seed for reproducibility
  # simulate 1000 samples from the bivariate Gaussian random walk process with cumsum(rnorm())
  y1 <- ts(cumsum(rnorm(1000, 0, sd=1)))
  y2 <- ts(cumsum(rnorm(1000, 0, sd=1)))
  y.art <- cbind(y1,y2)
  plot(y.art)
  # run the function of the basic model and check the posterior parameters B0 and B1
  basic.artificial <- basic.model(y.art, 1, 10000, c(), c())

  ############################################################
  ############################################################
  #Extendend model computation
  ############################################################
  ############################################################
  
  extended.model <- function(bigy, p, S1, S2){
    ############################################################
     #checks
     # bigy <- window(bigy, start=start, end=end)
     # y <- window(y, start=c(1990,2), end=c(2010,12))
     # bigy=y
     # p=12
     # S1=5 #repetitions to discard
     # S2=500 #repetitions to keep
     # kappa.3=0.95
     # s=1
  ############################################################
  
  Y       = bigy[(p+1):nrow(bigy),]
  X       = matrix(1,nrow(Y),1)
  for (i in 1:p){
    X     = cbind(X,bigy[(p+1):nrow(bigy)-i,])
  }
  
  N       = ncol(bigy)
  K       = 1+N*p
  S=S1+S2
  # set the priors
  ############################################################
  # priors <- setpriors(N,K,kappa.3)
  # priors$A.prior  = A.prior
  # priors$V.prior  = V.prior
  # priors$S.prior  = S.prior
  # priors$nu.prior = nu.prior
  #MVIW prior
  A.prior     = matrix(0,K,N)
  A.prior[2:(N+1),] <- diag(c(0,rep(1,N-1)))
  V.prior         = (diag(c(10,1*((1:p)^(-2))%x%rep(1,N)))) #10 is kappa.2, 1 is kappa.1
  nu.prior        = N+1
  #kappa.a priors
  nu.prior.a      = N + 1        #ig2
  s.a.prior       = 0.01
  s.sigma.prior   = rep(1,N)     #ig2
  #kappa.sigma priors
  s.a.prior       = 0.01           #rgamma
  a.sigma.prior   = 1           #rgamma
  
  A.posterior         = array(rnorm(prod(c(K,N,S))), c(K, N, S))
  B.posterior         = array(NA,c(N,N,S))
  Sigma.post          = array(NA,c(N,N,S))
  k.sigma             = matrix(NA, S, 1)
  k.a                 = matrix(NA, S, 1)
  k.sigma[1] <- 1
  k.a[1]     <- 1
  
  for (s in 1:S) {
    # normal-inverse Wishart, IG2, and Gamma posterior parameters
    ############################################################
    V.bar.inv   = t(X)%*%X + solve(k.a[s]*V.prior)
    V.bar       = solve(V.bar.inv)
    A.bar       = V.bar%*%(t(X)%*%Y + diag(1/diag(k.a[s]*V.prior))%*%A.prior)
    nu.bar      = nrow(Y) + nu.prior
    nu.bar.a    = nu.prior.a + N*K
    S.bar       = k.sigma[s]*diag(N) + t(Y)%*%Y + t(A.prior)%*%diag(1/diag(k.a[s]*V.prior))%*%A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
    S.bar.inv   = solve(S.bar)
    a.sigma.bar = (N*nu.prior)/2+a.sigma.prior
    # posterior draws
    ############################################################
    Sigma.posterior.inv = rWishart(1, df=nu.bar, Sigma=S.bar.inv)
    Sigma.posterior     = apply(Sigma.posterior.inv,3,solve)
    Sigma.posterior     = matrix(Sigma.posterior, N, N)
    L                   = t(chol(V.bar))
    cholSigma.s       = chol(Sigma.posterior)
    A.posterior[,,s]  = A.bar + L%*%A.posterior[,,s]%*%cholSigma.s
    
    s.a.bar     <- s.a.prior + sum(diag( Sigma.posterior.inv[,,1] %*% t(A.posterior[,,s]-A.prior)%*% diag(1/diag(V.prior)) %*% (A.posterior[,,s]-A.prior)))
    s.sigma.bar   <- (2*(sum(diag(Sigma.posterior.inv[,,1]))^-1)+((s.sigma.prior)^-1))^-1
    if (s <= S){
      #k.a from IG2
      k.a[s+1]      <- s.a.bar / rchisq(1, df=nu.bar.a)
      #k.sigma from gamma distribution
      k.sigma[s+1]  <- rgamma(1, shape = a.sigma.bar, scale = s.sigma.bar)
    }
    
    B.posterior[,,s]  = t(chol(Sigma.posterior.inv[,,]))
    Sigma.post[,,s]   = Sigma.posterior.inv
  }
  
  results = list("B" = B.posterior[,,S1+1:S2], 
                 "A" = A.posterior[,,S1+1:S2], 
                 "Sigma" = Sigma.post[,,S1+1:S2], 
                 "k.a" = k.a[S1+1:S2], 
                 "k.sigma" = k.sigma[S1+1:S2], 
                 "B.convergence" = B.posterior[,,S], 
                 "A.convergence" = A.posterior[,,S], 
                 "Sigma.convergence" = Sigma.post[,,S], 
                 "k.a.convergence" = k.a[S], 
                 "k.sigma.convergence" = k.sigma[S])
  return(results)
  }
  
  y                   <- window(y, start=c(1990,2), end=c(2010,12))
  
  extended.result     <- extended.model(y, 12, 500, 10000)
  extended.artificial <- extended.model(y.art, 1, 500, 10000)
  
  ############################################################
  ############################################################
  #plots for artificial data convergence
  
  
  ############################################################
  ############################################################  
  #presentation contents, 10 minutes
  #data
  #properties of data
  #how I obtain shocks of reference
  #present IRF
  
  ############################################################
  ############################################################  
  #learning repo
  #write a section on how to draw hyperparameters (for instance)
  #one page A4 with codes should be okay
  
  ############################################################
  ############################################################  
  #next week unconditional heteroskedasticity
  
  ####
  #adf test
  
  
  
  adf <- as.data.frame(matrix(nrow=14,ncol=3,NA))
  rownames(adf) <- colnames(y[1:14])
  colnames(adf) <- c("Dickey-Fuller","Lag-order", "p-value")
  
  for (i in 1:ncol(y)){
    adf_tmp                 <-  adf.test(y[,i])
    adf[i,"Dickey-Fuller"]  <-  round(as.numeric(adf_tmp[1]),3)
    adf[i,"Lag-order"]      <-  as.numeric(adf_tmp[2])
    adf[i,"p-value"]        <-  round(as.numeric(adf_tmp[4]),3)
  }
  
  adf.diff <- as.data.frame(matrix(nrow=14,ncol=3,NA))
  rownames(adf.diff) <- colnames(y[1:14])
  colnames(adf.diff) <- c("Dickey-Fuller diff","Lag-order", "p-value")
  
  for (i in 1: ncol(y)){
    tmp.diff                    <-  adf.test(diff(y[,i]))
    adf.diff[i,"Dickey-Fuller diff"] <-  round(as.numeric(tmp.diff[1]),3)
    adf.diff[i,"Lag-order"]     <-  as.numeric(tmp.diff[2])
    adf.diff[i,"p-value"]       <-  round(as.numeric(tmp.diff[4]),3)
  }
  
  knitr::kable(cbind(adf, adf.diff), index=TRUE)
  # USlev.asset    <- fredr_series_observations("BOGZ1FU664090005Q", observation_start = as.Date("1960-01-01"), observation_end = as.Date("2020-12-31"), aggregation_method = "eop", units="pch")
  # USlev.liab     <- fredr_series_observations("BOGZ1FL664194005Q", observation_start = as.Date("1960-01-01"), observation_end = as.Date("2020-12-31"), aggregation_method = "eop", units="pch")
  # USlev.leverage <- USlev.asset$value/USlev.liab$value
  # plot(x=USlev.leverage, y=USlev.asset$value, main="Procyclical Leverage of US Brokers and Dealers", xlab="Leverage Growth (% change)", xlim=c(-40,40), ylab="Total Asset Growth (% change)", ylim=c(-300,200))
  # 