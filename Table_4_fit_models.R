###--------------------------------------------------------------------------------------------###
# In this script we present the dataset generation and the estimation of the parameters.
# The dataset follows the characteristics of DOPPS data.
# The estimation steps follow the sequence outlined in the paper about dependent censoring, 
# for Weibull and PEM adjustment.
###--------------------------------------------------------------------------------------------###

rm(list=ls(all=TRUE))  

set.seed(123456789)

source("functions_Weibull.R",local=TRUE)
require(dlm)
require(survival)
require(rootSolve)

###--------------------------------------------------------------------------------------------###
# Data generated with the characteristics of the dataset provided by DOPPS.
###--------------------------------------------------------------------------------------------###

alpha.t <- 1.5
alpha.c <- 1.3
alpha.r <- 1.3
gamma.t <- 0.03
gamma.c <- 0.02
gamma.r <- 0.02
beta.t <- c(0.43, 0.01, -0.44, 0.21)
beta.c <- c(0.25, -0.29, -0.46, -0.18)
beta.r <- c(-0.72, 0.22, -0.83, -0.51)
sigma2 <- 0.115
phi.c <- 2.05 
phi.r <- 0.81

n <- 2500
ident <- c(rep(1,7), rep(2:10,20), rep(11:30,30), rep(31:50,45),rep(51:60,39), rep(61:70,30), rep(71,123))  
age <- rnorm(n,mean = 62.83, sd = 13)
gender <- rbinom(n,1,0.55)                   # 1 is male and 0 is female
race <- rbinom(n,1,0.32)                     # 1 if it is black and 0 if it is other
diabetic <- rbinom(n,1,0.47)                 # 1 if is diabetic and 0 otherwise

X <- cbind(scale(age), gender, race, diabetic)  
w_aux <- rnorm(71,0,sqrt(sigma2))
ww = w_aux[ident]

dados <- simula.frailty(alpha.c,alpha.r,alpha.t,beta.c,beta.r,beta.t,gamma.c,gamma.r,gamma.t,phi.c,phi.r, X=X, w=ww, ident=ident)

###--------------------------------------------------------------------------------------------###
# Models Adjustment
###--------------------------------------------------------------------------------------------###

Weibull_dep <- Weibull.dep(dados)
Weibull_indep <- Weibull.indep(dados)

source("functions_PEM.R",local=TRUE)
PEM_dep <- PEM.dep(dados)
PEM_indep <- PEM.indep(dados)

stat<-matrix(rep(NA,10*16),ncol=10,nrow=16)
stat<-as.data.frame(stat)

row.names(stat) <- c("$M1$", "$\\beta^{T}_{age}$","$\\beta^{T}_{gender}$","$\\beta^{T}_{race}$","$\\beta^{T}_{diabetic}$","$\\sigma^{2}$",
                     "$\\beta^{C}_{age}$","$\\beta^{C}_{gender}$","$\\beta^{C}_{race}$","$\\beta^{C}_{diabetic}$","$\\phi^{C}$",
                     "$\\beta^{R}_{age}$","$\\beta^{R}_{gender}$","$\\beta^{R}_{race}$","$\\beta^{R}_{diabetic}$","$\\phi^{R}$")

colnames(stat)<- c( "Estimate", "SE", "LI", "LS", " ", "Estimate", "SE", "LI", "LS", "RD(%)")

stat[1,1:10] <- c(" ")
stat[2:16,1] <- Weibull_dep[,1]
stat[2:16,2] <- Weibull_dep[,2] 
stat[2:16,3] <- round(Weibull_dep[,1]-1.96*Weibull_dep[,2],3)
stat[2:16,4] <- round(Weibull_dep[,1]+1.96*Weibull_dep[,2],3)
stat[2:16,5] <- c(" ") 
stat[2:16,6] <- Weibull_indep[,1] 
stat[2:16,7] <- Weibull_indep[,2]
stat[2:16,8] <- round(Weibull_indep[,1]-1.96*Weibull_indep[,2],3)
stat[2:16,9] <- round(Weibull_indep[,1]+1.96*Weibull_indep[,2],3)
stat[2:6,10] <- round(((Weibull_indep[,1]-Weibull_dep[1:5,1])/abs(Weibull_dep[1:5,1]))*100 ,3)
stat[7:16,10] <- c(" ")


stat.PE<-matrix(rep(NA,10*16),ncol=10,nrow=16)
stat.PE<-as.data.frame(stat.PE)

row.names(stat.PE) <- c("$M2$", "$\\beta^{T}_{age}$","$\\beta^{T}_{gender}$","$\\beta^{T}_{race}$","$\\beta^{T}_{diabetic}$","$\\sigma^{2}$",
                     "$\\beta^{C}_{age}$","$\\beta^{C}_{gender}$","$\\beta^{C}_{race}$","$\\beta^{C}_{diabetic}$","$\\phi^{C}$",
                     "$\\beta^{R}_{age}$","$\\beta^{R}_{gender}$","$\\beta^{R}_{race}$","$\\beta^{R}_{diabetic}$","$\\phi^{R}$")

colnames(stat.PE)<- c( "Estimate", "SE", "LI", "LS", " ", "Estimate", "SE", "LI", "LS", "RD(%)")

stat.PE[1,1:10] <- c(" ")
stat.PE[2:16,1] <- PEM_dep[,1]
stat.PE[2:16,2] <- PEM_dep[,2] 
stat.PE[2:16,3] <- round(PEM_dep[,1]-1.96*PEM_dep[,2],3)
stat.PE[2:16,4] <- round(PEM_dep[,1]+1.96*PEM_dep[,2],3)
stat.PE[2:16,5] <- c(" ") 
stat.PE[2:16,6] <- PEM_indep[,1] 
stat.PE[2:16,7] <- PEM_indep[,2]
stat.PE[2:16,8] <- round(PEM_indep[,1]-1.96*PEM_indep[,2],3)
stat.PE[2:16,9] <- round(PEM_indep[,1]+1.96*PEM_indep[,2],3)
stat.PE[2:6,10] <- round(((PEM_indep[,1]-PEM_dep[1:5,1])/abs(PEM_dep[1:5,1]))*100 ,3)
stat.PE[7:16,10] <- c(" ")

# Table 4
rbind(stat,stat.PE)

