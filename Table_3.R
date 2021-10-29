###--------------------------------------------------------------------------------------------###
# The data generation follows the characteristics of DOPPS data.
# In this script we present the generation of the dataset simulated that mimics the real data,
# and we show the construction of Table 3.
###--------------------------------------------------------------------------------------------###

rm(list=ls(all=TRUE)) 

set.seed(123456)

source("functions_Weibull.R",local=TRUE)
library(finalfit)
library(dplyr)

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
sigma2 <- 0.11
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
# Table 3: Frequency distribution between failure, withdrawal and transplant.
###--------------------------------------------------------------------------------------------###

age.dic  <- ifelse(age >= 50,1,0)
age.dic  <- factor(age.dic, label = c("<50", ">=50"), levels = 0:1); age.dic <- relevel(age.dic,ref=">=50")
gender   <- factor(gender, label = c("female", "male"), levels = 0:1); gender <- relevel(gender,ref="male")
race     <- factor(race, label = c("other", "black"), levels = 0:1); race <- relevel(race,ref="black")
diabetic <- factor(diabetic, label = c("otherwise", "diabetic"), levels = 0:1); diabetic <- relevel(diabetic,ref="diabetic")

cens <- factor(dados$cens, label = c("Failure", "Withdrawal", "Transplant", "Adm.Censoring" ), levels = 1:4)

data.tb <- data.frame(age.dic, gender, race, diabetic, cens)
                      
explanatory = c("age.dic", "gender", "race", "diabetic")
dependent = "cens" 
data.tb %>%
  summary_factorlist(dependent, explanatory,add_dependent_label=TRUE)
                      
###--------------------------------------------------------------------------------------------###
