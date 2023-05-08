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
                         observation_start = as.Date("1990-01-01"), 
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
y = cbind(#instrument,
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
y <- window(y, start=c(1990,1), end=c(2010,12)) #subset in our time window: 1990,1 until 2011,1
# setup function for my analysis
basic.model <- function(bigy, p, S, sign.restrictions, kappa.1, kappa.2, kappa.3){
############################################################
N       = ncol(bigy)
K       = 1+N*p
h       = 6 #forecast horizon
S.burnin= 100
############################################################

Y       = bigy[(p+1):nrow(bigy),]
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,bigy[(p+1):nrow(bigy)-i,])
}

# MLE
############################################################
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/nrow(Y)

# set the priors
############################################################
A.prior     = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] <- kappa.3*diag(N)
V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N)))
S.prior     = diag(diag(Sigma.hat))
nu.prior    = N+1

# normal-inverse Wishart posterior parameters
############################################################
V.bar.inv   = t(X)%*%X + diag(1/diag(V.prior))
V.bar       = solve(V.bar.inv)
A.bar       = V.bar%*%(t(X)%*%Y + diag(1/diag(V.prior))%*%A.prior)
nu.bar      = nrow(Y) + nu.prior
S.bar       = S.prior + t(Y)%*%Y + t(A.prior)%*%diag(1/diag(V.prior))%*%A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
S.bar.inv   = solve(S.bar)

# posterior draws
############################################################
Sigma.posterior   = rWishart(S, df=nu.bar, Sigma=S.bar.inv)
Sigma.posterior   = apply(Sigma.posterior,3,solve)
Sigma.posterior   = array(Sigma.posterior,c(N,N,S))
A.posterior       = array(rnorm(prod(c(dim(A.bar),S))),c(dim(A.bar),S))
B.posterior       = array(NA,c(N,N,S))
Bplus.posterior   = array(NA,c(N,K,S))
L                 = t(chol(V.bar))
for (s in 1:S){
  cholSigma.s     = chol(Sigma.posterior[,,s])
  B.posterior[,,s]= t(cholSigma.s)
  A.posterior[,,s]= A.bar + L%*%A.posterior[,,s]%*%cholSigma.s
  Bplus.posterior[,,s] = B.posterior[,,s]%*%t(A.posterior[,,s])
}

# Estimating models with sign restrictions: example
# Use Algorithm 2
############################################################
R1            = diag(sign.restrictions)
#A             = t(matrix(c(.5,.5,0,-1.25,.25,0,-1,0,.5),3,3))
#Sigma         = matrix(c(1,.5,1,.5,4.25,2.5,1,2.5,3),3,3)

#store draws results for identified parameters from desired restrictions in
i.vec = c()
Q.vec = array(NA, c(N,N,S))
B0.rot = array(NA, c(N,N,S))
B1.rot = array(NA, c(N,K,S))
A.rot = array(NA, c(K,N,S))

for (s in 1:S){
  
  B0.tilde      = B.posterior[,,s]
  B1.tilde      = Bplus.posterior[,,s]
  IR.0.tilde    = solve(B0.tilde)
  #IR.1.tilde    = solve(B0.tilde)%*%B1.tilde%*%solve(B0.tilde)
  
  sign.restrictions.do.not.hold = TRUE
  i=1
  while (sign.restrictions.do.not.hold){
    X           = matrix(rnorm(N^2),N,N)
    QR          = qr(X, tol = 1e-10)
    Q           = qr.Q(QR,complete=TRUE)
    R           = qr.R(QR,complete=TRUE)
    Q           = t(Q %*% diag(sign(diag(R))))
    B0          = Q%*%B0.tilde
    B1          = Q%*%B1.tilde
    B0.inv      = solve(B0)
    A           = t(B0.inv %*% B1)
    check       = prod(R1 %*% B0.inv %*% diag(N)[,1] > 0)
    if (check==1){sign.restrictions.do.not.hold=FALSE}
    i=i+1
  }
  
  #i.vec = c(i.vec, i)
  B0.rot[,,s] = B0
  B1.rot[,,s] = B1
  A.rot[,,s]  =  A
  #Q.vec[,,s] = Q
  #IR.0        = B0.inv
  #IR.1        = B0.inv%*%B1%*%B0.inv
  #IR.0
  #IR.1
  #R1 %*% rbind(IR.0,IR.1) %*% diag(N)[,1]
  
}

results = list("B0" = B0, "B1" = B1, "A" = A )
return(results)
}


#run the function with my data
set.seed(23)
basic.model(y, 12, 5, c(1,-1,-1,-1,-1,-1,1), (0.2)^2, 100, 1)






#artificial data with 1000 observations from a bi-variate Gaussian random walk
#VAR model with N=2, p=1, and a constant term
set.seed(0)   # set seed for reproducibility
# simulate 1000 samples from the bivariate Gaussian random walk process with cumsum(rnorm())
y1 <- ts(cumsum(rnorm(1000, 0, sd=1)))
y2 <- ts(cumsum(rnorm(1000, 0, sd=1)))
yart <- cbind(y1,y2)
plot(yart)
# run the function of the basic model and check the posterior parameters B0 and B1
basic.model(yart, 1, 5000, c(1,1), (00.2)^2, 100, 1)