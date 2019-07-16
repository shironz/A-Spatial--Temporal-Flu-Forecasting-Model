
##################################################
#  8820 Introduction to Bayesian Statistics
#  Project 1: Predicts Flu
#  The Simulation Part
#  Shirong Zhao
#  Email: shironz@clemson.edu
#################################################


library(MASS)
library(coda)
library(mvtnorm)
library(Matrix) 
library(mnormt)
library(lattice)

##########################################
# Generates data and necessary components

data.gen<-function(P1,P2){
  
  V<-P1*P2
  Vid<-array(1:V,c(P1,P2)) # Arbitrarily identifies voxels   
  
  Y<-rep(-99,V) # Arbitrary response vector                 
  Dj<-rep(-99,V)
  Nj<-matrix(-99,nrow=V,ncol=8)  # This will list the neighbors for each voxel
  
  
  for(p1 in 1:P1){
    for(p2 in 1:P2){
      
      # Create the neighborhood matrix, the design matrix, and the segmentation
      # vectors
      l1<-max(c(1,p1-1))
      l2<-max(c(1,p2-1))
      u1<-min(c(P1,p1+1))
      u2<-min(c(P2,p2+1))
      
      neighbors<-as.vector(Vid[l1:u1,l2:u2])
      neighbors<-neighbors[neighbors!=Vid[p1,p2]]
      
      Y[Vid[p1,p2]]<- 10*dnorm((p1-p2),0,10) # rnorm(1,0,0.1 is the noise data, here no noise
      
      wjv<-length(neighbors)
      Dj[Vid[p1,p2]]<-wjv
      Nj[Vid[p1,p2],1:wjv]<-neighbors           
    }
  }
  # Creating the necessary sparse matrices
  i<-NULL
  j<-NULL
  
  for(v in 1:V){
    i<-c(i,rep(v,Dj[v]))
    j<-c(j,Nj[v,1:Dj[v]])
  }
  
  W<- sparseMatrix(i, j, x = 1)  
  D<-sparseMatrix(1:V,1:V,x=Dj)
  
  return(list("Y"=Y,"D"=D,"W"=W))
}



P1 = 6
P2 = 6

data<-data.gen(P1=P1,P2=P2)
Y0<-matrix( data$Y, nrow =P1 , ncol = P2) 
#windows()
levelplot(Y0)
N=P1*P2
M=5

phi = matrix(0, nrow = N, ncol = M)
epsilon = matrix(0, nrow = N, ncol = M)

D<-data$D
W<-data$W
rho.true<-0.5
DW.true<-D-rho.true*W  # rho will be specified in advance, for example, let rho = 1
DWI.true<-solve(DW.true)
DW.true<-as.matrix(DW.true)
DWI.true<-as.matrix(DWI.true)
tau2.true<-2
theta.true <- 0.5
# phi = t(rmvnorm(1, mean=rep(0,N), sigma=tau2.true*DWI.true))

for (t in 1:M) {
  phi[,t] = t(rmvnorm(1, mean=rep(0,N), sigma=tau2.true*DWI.true))
}
epsilon[,1] = phi[,1]
for (t in 2:M) {
  epsilon[,t] = theta.true*epsilon[,t-1] + phi[,t]
}


## simulate data
epsilon.true<-as.vector(epsilon)
NM<-N*M
X<-cbind(1,rnorm(NM))
beta.true<-c(1,2)
lambda<-exp(X%*%beta.true + epsilon.true)
Y<-rpois(NM,lambda)


length(Y)

dim(X)

# First using MLE find the sd for proposal distribution of beta
fit <- glm(Y ~ X[,2], family=poisson()) # X[,1] is the intercept
summary(fit) 

#############################################################################
#############################################################################
## Simulation Part
#############################################################################
# Inputs:
# Y = response vector for current data, N*M by 1 
# N = the numner of states
# M = the number of months
# p: number of covariates
# X = design matrix for current data N*M by p
# R = prior covariance matrix for beta (i.e. beta~N(0,R))   
# a0 = prior parameter for tau2
# b0 = prior parameter for tau2, where tau2^{-1}~gamma(a0,b0)
# beta = initial value of regression coefficients
# tau2  = initial value of precission parameter
# theta: initial value for autoregressive coefficients
# rho: initial value for rho
# beta.var.prop: variance for the beta proposal distribution,(i.e., beta.p~N(beta(s-1), beta.var.prop))
# phi.var.prop: variance for the phi proposal distribution,(i.e., phi.p~N(phi(s-1), phi.var.prop))
# c: variance for the rho proposal distribution


#N = 51 
#M = 7 

p = dim(X)[2]
NM = dim(X)[1]

iter = 2e5
thin = 1e2

# a_chol<-chol(DW)
# chol2inv(a_chol)

beta.var.prop=vcov(fit)
#beta.var.prop<-diag(rep(0.00005,p), p, p) # need to specify later and it's better to use var of parameters in poisson regression
# for here I just specify var as 0.1
# we could also consider var.prop<-var(log(Y))*solve(t(X)%*%X)
phi.var.prop=diag(rep(0.00002,N), N, N) # need to consider later
delta=0.010 # used in proposal disttribution for theta,reflected random walk 
c=2 # specify the variance of proposal distribution for rho, log(rho.p/(1-rho.p)) follows normal (log(rho/(1-rho)), c)
# or rho.p/(1-rho.p) follows the lognormal(log(rho/(1-rho)), c)


# Specify the priors
R<-rep(10, p) 
a0 = 1
b0 = 1
phi = matrix(0, nrow = N, ncol = M)
epsilon = matrix(0, nrow = N, ncol = M)
epsilon<-as.vector(epsilon) # nrow=N*M
tau2 = 1
theta = 0
rho = 0.5
rho1 = rho/(1-rho)

DW<-D-rho*W  
DWI<-solve(DW)
DW<-as.matrix(DW)
DWI<-as.matrix(DWI)

# save the parameters
Beta = matrix(-99, nrow=iter/thin, ncol=p) 
beta = rep(0,p)
Beta[1, ] = beta
Tau2 = rep(-99, iter/thin)
Tau2[1] = tau2
Theta = rep(-99, iter/thin)
Theta[1] = theta
Rho = rep(-99, iter/thin)
Rho[1] = rho
Phi = matrix(-99, nrow=iter/thin, ncol=N*M)
Phi[1, ] = as.vector(phi)

acc0 = acc1 = acc2 = acc3 = 0

llik = sum(dpois(Y, exp(X%*%beta + epsilon), log = TRUE))

############################################################################   
# Burn in loop

for(i in (thin + 1):iter){
  
 
  ## update all beta simultaneously, using Metropolis Algorithm
  beta.p = t(rmvnorm(1, beta, beta.var.prop))
  beta.prior = sum(dnorm(beta, 0, R, log = TRUE))
  llik.p = sum(dpois(Y, exp(X%*%beta.p + epsilon), log = TRUE))
  beta.prior.p = sum(dnorm(beta.p, 0, R, log = TRUE))
  r = exp(llik.p -llik + beta.prior.p - beta.prior)
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    beta = beta.p
    llik = llik.p 
    acc0 = acc0 + 1 	
  }
  
  
  
  ## update epsilon and phi (the spacial random effect), using Metropolis Algorithm
  phi.p = matrix(-99, nrow = N, ncol = M) 
  epsilon.p = matrix(-99, nrow = N, ncol = M)
  phi.prior = 0
  phi.prior.p = 0
  for (t in 1:M) {
    phi.p[,t] = t(rmvnorm(1, phi[,t], phi.var.prop)) 
    phi.prior.t = dmvnorm(phi[,t], mean=rep(0,N), sigma=tau2*DWI, log = TRUE) # sigma is covariance matrix 
    phi.prior = phi.prior + phi.prior.t
    phi.prior.p.t = dmvnorm(phi.p[,t], mean=rep(0,N), sigma=tau2*DWI, log = TRUE)
    phi.prior.p = phi.prior.p + phi.prior.p.t
  }
  epsilon.p[,1] = phi.p[,1]
  for (t in 2:M) {
    epsilon.p[,t] = theta*epsilon.p[,t-1] + phi.p[,t]
  }
  epsilon.p = as.vector(epsilon.p) # change to a vector
  llik.p = sum(dpois(Y, exp(X%*%beta + epsilon.p), log = TRUE))
  r = exp(llik.p -llik + phi.prior.p - phi.prior) # need to add density of proposal distribution
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    phi = phi.p # phi is a matrix
    epsilon = epsilon.p # epsilon is a vector
    llik = llik.p 
    acc1 = acc1 + 1 	
  }
  
  
  
  ## update tau2, using Gibbs sampler
  sumt = 0
  for (t in 1:M) {
    sumt = sumt + t(phi[,t])%*%DW%*%phi[,t]/2
  }
  at = a0 + N*M/2
  bt = b0 + sumt
  tauI2 = rgamma(1,at,bt)
  tau2 = 1/tauI2
  
  
  
  ## update theta, using Metropolis-Hastings Algorithm, reflected random walk
  # prior for theta is uniform(-1,1)
  theta.p = runif(1, min=theta-delta, max=theta+delta) 
  if (theta.p < -1){
    theta.p = -2- theta.p
  } else if (theta.p > 1){
    theta.p = 2-theta.p}
  epsilon.p = matrix(-99, nrow = N, ncol = M)
  epsilon.p[,1] = phi[,1]
  for (t in 2:M) {
    epsilon.p[,t] = theta.p*epsilon.p[,t-1] + phi[,t]
  }
  epsilon.p = as.vector(epsilon.p) # change to a vector
  llik.p = sum(dpois(Y, exp(X%*%beta + epsilon.p), log = TRUE))
  r = exp(llik.p -llik)
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    theta = theta.p 
    epsilon = epsilon.p # epsilon is a vector
    llik = llik.p 
    acc2 = acc2 + 1 	
  }
  
  
  ## update rho, using Metropolis-Hastings Algorithm, symmetric random walk
  # prior for rho is uniform(0,1)
  rho1.p= rlnorm(1, meanlog = log(rho/(1-rho)), sdlog = c)
  rho.p=rho1.p/(1+rho1.p) # rho1.p=rho.p/(1-rho.p)
  DW.p<-D-rho.p*W  
  DWI.p<-solve(DW.p)
  DW.p<-as.matrix(DW.p)
  DWI.p<-as.matrix(DWI.p)  
  
  llik.rho = 0
  llik.rho.p = 0
  
  for (t in 1:M) {
    llik.rho.t = dmvnorm(phi[,t], mean=rep(0,N), sigma=tau2*DWI, log = TRUE) # sigma is covariance matrix 
    llik.rho = llik.rho + llik.rho.t
    llik.rho.p.t = dmvnorm(phi[,t], mean=rep(0,N), sigma=tau2*DWI.p, log = TRUE)
    llik.rho.p = llik.rho.p + llik.rho.p.t
  }
  
  r=exp(llik.rho.p - llik.rho 
        + dlnorm(rho1, meanlog = log(rho1.p), sdlog = c, log = TRUE)
        - dlnorm(rho1.p, meanlog = log(rho1), sdlog = c, log = TRUE))
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    rho =rho.p
    rho1 = rho1.p
    DWI = DWI.p
    acc3 = acc3 + 1 	
  }
  
  
  ## tuning the necessary parameters
   if(i %% 1000 == 0){
      beta.var.prop<-  beta.var.prop + (acc0/1000 >0.55)*0.75*beta.var.prop - (acc0/1000 < 0.35)*0.75*beta.var.prop        
      phi.var.prop<-  phi.var.prop + (acc1/1000 >0.55)*0.75*phi.var.prop - (acc1/1000 < 0.35)*0.75*phi.var.prop       
      delta<-  delta + (acc2/1000 >0.55)*0.75*delta - (acc2/1000 < 0.35)*0.75*delta
      c<-c + (acc3/1000 >0.55)*0.75*c - (acc3/1000 < 0.35)*0.75*c  
      print(c(acc0,acc1,acc2,acc3))
      acc0<-0
      acc1<-0
      acc2<-0
      acc3<-0
      
  }
}


#################################################################################################
# Sampling loop

for(i in (thin + 1):iter){
    
 
  ## update all beta simultaneously, using Metropolis Algorithm
  beta.p = t(rmvnorm(1, beta, beta.var.prop))
  beta.prior = sum(dnorm(beta, 0, R, log = TRUE))
  llik.p = sum(dpois(Y, exp(X%*%beta.p + epsilon), log = TRUE))
  beta.prior.p = sum(dnorm(beta.p, 0, R, log = TRUE))
  r = exp(llik.p -llik + beta.prior.p - beta.prior)
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    beta = beta.p
    llik = llik.p 
    acc0 = acc0 + 1 	
  }
  
  
  
  ## update epsilon and phi (the spacial random effect), using Metropolis Algorithm
  phi.p = matrix(-99, nrow = N, ncol = M) 
  epsilon.p = matrix(-99, nrow = N, ncol = M)
  phi.prior = 0
  phi.prior.p = 0
  for (t in 1:M) {
    phi.p[,t] = t(rmvnorm(1, phi[,t], phi.var.prop)) 
    phi.prior.t = dmvnorm(phi[,t], mean=rep(0,N), sigma=tau2*DWI, log = TRUE) # sigma is covariance matrix 
    phi.prior = phi.prior + phi.prior.t
    phi.prior.p.t = dmvnorm(phi.p[,t], mean=rep(0,N), sigma=tau2*DWI, log = TRUE)
    phi.prior.p = phi.prior.p + phi.prior.p.t
  }
  epsilon.p[,1] = phi.p[,1]
  for (t in 2:M) {
    epsilon.p[,t] = theta*epsilon.p[,t-1] + phi.p[,t]
  }
  epsilon.p = as.vector(epsilon.p) # change to a vector
  llik.p = sum(dpois(Y, exp(X%*%beta + epsilon.p), log = TRUE))
  r = exp(llik.p -llik + phi.prior.p - phi.prior)
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    phi = phi.p # phi is a matrix
    epsilon = epsilon.p # epsilon is a vector
    llik = llik.p 
    acc1 = acc1 + 1 	
  }
  
  
  
  ## update tau2, using Gibbs sampler
  sumt = 0
  for (t in 1:M) {
    sumt = sumt + t(phi[,t])%*%DW%*%phi[,t]/2
  }
  at = a0 + N*M/2
  bt = b0 + sumt
  tauI2 = rgamma(1,at,bt)
  tau2 = 1/tauI2
  
  
  
  ## update theta, using Metropolis-Hastings Algorithm, reflected random walk
  # prior for theta is uniform(-1,1)
  theta.p = runif(1, min=theta-delta, max=theta+delta) 
  if (theta.p < -1){
    theta.p = -2- theta.p
  } else if (theta.p > 1){
    theta.p = 2-theta.p}
  epsilon.p = matrix(-99, nrow = N, ncol = M)
  epsilon.p[,1] = phi[,1]
  for (t in 2:M) {
    epsilon.p[,t] = theta.p*epsilon.p[,t-1] + phi[,t]
  }
  epsilon.p = as.vector(epsilon.p) # change to a vector
  llik.p = sum(dpois(Y, exp(X%*%beta + epsilon.p), log = TRUE))
  r = exp(llik.p -llik)
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    theta = theta.p 
    epsilon = epsilon.p # epsilon is a vector
    llik = llik.p 
    acc2 = acc2 + 1 	
  }
  
  
  ## update rho, using Metropolis-Hastings Algorithm, symmetric random walk
  # prior for rho is uniform(0,1)
  rho1.p= rlnorm(1, meanlog = log(rho/(1-rho)), sdlog = c)
  rho.p=rho1.p/(1+rho1.p) # rho1.p=rho.p/(1-rho.p)
  DW.p<-D-rho.p*W  
  DWI.p<-solve(DW.p)
  DW.p<-as.matrix(DW.p)
  DWI.p<-as.matrix(DWI.p)  
  
  llik.rho = 0
  llik.rho.p = 0
  
  for (t in 1:M) {
    llik.rho.t = dmvnorm(phi[,t], mean=rep(0,N), sigma=tau2*DWI, log = TRUE) # sigma is covariance matrix 
    llik.rho = llik.rho + llik.rho.t
    llik.rho.p.t = dmvnorm(phi[,t], mean=rep(0,N), sigma=tau2*DWI.p, log = TRUE)
    llik.rho.p = llik.rho.p + llik.rho.p.t
  }
  
  r=exp(llik.rho.p - llik.rho 
        + dlnorm(rho1, meanlog = log(rho1.p), sdlog = c, log = TRUE)
        - dlnorm(rho1.p, meanlog = log(rho1), sdlog = c, log = TRUE))
  Z<-rbinom(1,1,min(r,1))  
  if(Z==1){
    rho =rho.p
    rho1 = rho1.p
    DWI = DWI.p
    acc3 = acc3 + 1 	
  }
  
  
  if(i %% thin == 0){
    Beta[i / thin, ] = beta
    Tau2[i / thin] = tau2
    Theta[i / thin] = theta
    Rho[i / thin] = rho
    Phi[i / thin, ] = as.vector(phi) 
    print(i)
  }
  
  
}


###############################################################################
###  Summarize the estimation outcome


acc0 / (iter-thin)
acc1 / (iter-thin)
acc2 / (iter-thin)
acc3 / (iter-thin)


Beta.mcmc = as.mcmc(Beta)
print(paste0("Estimate Mean of beta:  ", apply(Beta.mcmc, 2, mean)))
print(paste0("True beta:  ", beta.true))
HPDinterval(Beta.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Beta.mcmc)))
plot(Beta.mcmc)
autocorr.plot(Beta.mcmc)
plot(Beta[, 1], typ = 'l')
plot(Beta[, 2], typ = 'l')

Tau2.mcmc = as.mcmc(Tau2) 
print(paste0("Estimate Mean of tau2:  ", mean(Tau2.mcmc)))
print(paste0("True theta:  ", tau2.true))
HPDinterval(Tau2.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Tau2.mcmc)))
plot(Tau2.mcmc)
autocorr.plot(Tau2.mcmc)
plot(Tau2, typ = 'l')

Theta.mcmc = as.mcmc(Theta) 
print(paste0("Estimate Mean of theta:  ", mean(Theta.mcmc)))
print(paste0("True theta:  ", theta.true))
HPDinterval(Theta.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Theta.mcmc)))
plot(Theta.mcmc)
autocorr.plot(Theta.mcmc)
plot(Theta, typ = 'l')


Rho.mcmc = as.mcmc(Rho) 
print(paste0("Estimate Mean of rho:  ", mean(Rho.mcmc)))
print(paste0("True rho:  ", rho.true))
HPDinterval(Rho.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Rho.mcmc)))
plot(Rho.mcmc)
autocorr.plot(Rho.mcmc)
plot(Rho, typ = 'l')

proc.time()
