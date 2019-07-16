
##################################################
#  8820 Introduction to Bayesian Statistics
#   Project 1: Predicts Flu
#   Shirong Zhao
#################################################
library(MASS)
library(mvtnorm)
library(coda)
library(Matrix) 
library(mnormt)
library(gdata) # use inside command write.fwf
#################################################
#######################  BEGIN Import Cleaned Data ##########################\

formatInfo6<-read.csv("./formatInfo6.csv")
df6 <- read.fwf(file="./outcome6.txt", widths=formatInfo6$width + 1, skip=1, strip.white=TRUE, na.strings="n.a.") 
# V14 is the number of flu for state i and time t

######################## End Import Cleaned Data #########################

Beta=df6[,1:12]
Tau2=df6[,13]
Theta=df6[,14]
Rho=df6[,15]

Beta.mcmc = as.mcmc(Beta)
print(paste0("Estimate Mean of beta:  ", apply(Beta.mcmc, 2, mean)))
HPDinterval(Beta.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Beta.mcmc)))
plot(Beta.mcmc)
autocorr.plot(Beta.mcmc)
plot(Beta[, 1], typ = 'l')
plot(Beta[, 2], typ = 'l')


Tau2.mcmc = as.mcmc(Tau2) 
print(paste0("Estimate Mean of tau2:  ", mean(Tau2.mcmc)))
HPDinterval(Tau2.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Tau2.mcmc)))
plot(Tau2.mcmc)
autocorr.plot(Tau2.mcmc)
plot(Tau2, typ = 'l')


Theta.mcmc = as.mcmc(Theta) 
print(paste0("Estimate Mean of theta:  ", mean(Theta.mcmc)))
HPDinterval(Theta.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Theta.mcmc)))
plot(Theta.mcmc)
autocorr.plot(Theta.mcmc)
plot(Theta, typ = 'l')


Rho.mcmc = as.mcmc(Rho) 
print(paste0("Estimate Mean of rho:  ", mean(Rho.mcmc)))
HPDinterval(Rho.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(Rho.mcmc)))
plot(Rho.mcmc)
autocorr.plot(Rho.mcmc)
plot(Rho, typ = 'l')


# Here we predict the flu from Jan 2017 to Dec 2017

formatInfo<-read.csv("./formatInfoall.csv")
df <- read.fwf(file="./all.txt", widths=formatInfo$width + 1, skip=1, strip.white=TRUE, na.strings="n.a.") 
# V14 is the number of flu for state i and time t


formatInfoW<-read.csv("./formatInfoW.csv")
W<- read.fwf(file="./w.txt", widths=formatInfoW$width + 1, skip=1, strip.white=TRUE, na.strings="n.a.") 

######################## End Import Cleaned Data #########################
dim(df)
W = as.matrix(W)
d=rowSums(W[,1:47])
D=diag(d,47,47)

df$V5<-as.numeric(gsub(",", "", df$V5))

Y=df$V14
x0=1
x1=df$V4/100000000 # population in 100 million
x2=df$V5/1000 # income in 1000
x3=df$V6 # temperature
x4=df$V7 # precipation
x5=df$V8/df$V4 # share of white
x6=df$V9/df$V4 # share of African American
x7=df$V10/df$V4 # share of American Indian
x8=df$V11/df$V4 # share of Asian

df$spring = 0  # summer as a benchmark
df$autumn = 0
df$winter = 0

df$spring[which(df$V3==3 | df$V3==4 | df$V3==5)] = 1
df$autumn[which(df$V3==9 | df$V3==10 | df$V3==11)] = 1
df$winter[which(df$V3==12 | df$V3==1 | df$V3==2)] = 1

spring=df$spring
autumn=df$autumn
winter=df$winter

X=cbind(x0,x1,x2,x3,x4,x5,x6,x7,x8,spring,autumn,winter)

# here we only use the data from Oct,2010-Dec,2016 to train the data
# the data in 2017 will be used to forecast

#Y=Y[1:(4089-564)]  # 47*12=564
#X=X[1:(4089-564),]

Y=Y[(4089-564+1):4089]  # 47*12=564
X=X[(4089-564+1):4089,]
dim(X)
# popu: one hundred million
# income: thousand
# temperature: Degrees Fahrenheit,
# Temperature is in Fahrenheit and Precipitation is in Inches 

# First using MLE find the sd for proposal distribution of beta
fit <- glm(Y ~ X[,1:12]-1, family=poisson()) # X[,1] is the intercept
summary(fit) 


# plug in the estimated parameter

beta = apply(Beta.mcmc, 2, mean)
tau2 =mean(Tau2.mcmc)
theta=mean(Theta.mcmc)
rho=mean(Rho.mcmc)

DW<-D-rho*W  
DWI<-solve(DW)
DW<-as.matrix(DW)
DWI<-as.matrix(DWI)

N=47
M=87
phi = matrix(-99, nrow = N, ncol = M) 
epsilon = matrix(-99, nrow = N, ncol = M) 
for (t in 1:M) {
    phi[,t] = t(rmvnorm(1, mean=rep(0,N), sigma=tau2*DWI)) 
  }
  epsilon[,1] = phi[,1]
  for (t in 2:M) {
    epsilon[,t] = theta*epsilon[,t-1] + phi[,t]
  }

epsilon=as.vector(epsilon[,76:87])

Y.pred=exp(X%*%beta+epsilon)

mean((Y-Y.pred)^2)

Y<-matrix(Y, nrow=47, ncol=12)

Y.pred<-matrix(Y.pred, nrow=47, ncol=12)

apply(Y, 2, mean)

apply(Y.pred, 2, mean)
