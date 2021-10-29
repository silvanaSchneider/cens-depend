
###------------------------------------------------------------------------------------------------###
# Function to generate the fragilities of the complete log conditional used in the arms function of 
# the dlm package.
###------------------------------------------------------------------------------------------------###
# assuming dependence  
log_cond_wi <- function(w_k_aux,phi_2,phi_3,delta_T,delta_C,delta_R,X_T,X_C, X_R, Betas_T,Betas_C,Betas_R,risco_a_T,risco_a_C,risco_a_R, Sigma2){
  w_k <- log(w_k_aux/(1-w_k_aux)) 
  pred_linear_T <- exp((X_T[,]%*%Betas_T)+w_k)
  pred_linear_C <- exp((X_C[,]%*%Betas_C)+phi_2*w_k)
  pred_linear_R <- exp((X_R[,]%*%Betas_R)+phi_3*w_k)
  log_vero_w <- sum(delta_T*(w_k) - risco_a_T*pred_linear_T + delta_C*phi_2*(w_k)-risco_a_C*pred_linear_C + delta_R*phi_3*(w_k)-risco_a_R*pred_linear_R)-((w_k^2)/(2*Sigma2))-log(w_k_aux*(1-w_k_aux))   
  return( log_vero_w)
}

# assuming independence
log_cond_wi_indep <- function(w_k_aux,delta_T,X_T,Betas_T,risco_a_T,Sigma2){
  w_k <- log(w_k_aux/(1-w_k_aux)) 
  pred_linear_T <- exp((X_T[,]%*%Betas_T)+w_k)
  log_vero_w <- sum(delta_T*(w_k) - risco_a_T*pred_linear_T)-((w_k^2)/(2*Sigma2))-log(w_k_aux*(1-w_k_aux))   
  return( log_vero_w)
}
###------------------------------------------------------------------------------------------------###
# specification of support for the arms function
###------------------------------------------------------------------------------------------------###
# assuming dependence
support_wi <-  function(w_k_aux,phi_2,phi_3,delta_T,delta_C,delta_R,X_T,X_C, X_R, Betas_T,Betas_C,Betas_R,risco_a_T,risco_a_C,risco_a_R, Sigma2){(w_k_aux>0)*(w_k_aux<1)}  
# assuming independence
support_wi_indep <-  function(w_k_aux,delta_T,X_T,Betas_T,risco_a_T,Sigma2){(w_k_aux>0)*(w_k_aux<1)} 
 
###------------------------------------------------------------------------------------------------###
# Function to generate the data with the characteristics of the dataset provided by DOPPS.
###------------------------------------------------------------------------------------------------###
simula.frailty <- function(alpha.c,alpha.r,alpha.t,beta.c,beta.r,beta.t,gamma.c,gamma.r,gamma.t,phi.c,phi.r, X, w, ident)
{
  t <- rep(NA,n)
  c <- rep(NA,n)
  r <- rep(NA,n)
  e <- rep(NA,n)
  y <- rep(NA,n)
  u <- runif(n)
  v <- runif(n)
  z <- runif(n)  
  cens <- rep(NA,n)
  lambda.t <- gamma.t*exp( X%*%beta.t + w)
  lambda.c <- gamma.c*exp( X%*%beta.c + phi.c*w)  
  lambda.r <- gamma.r*exp( X%*%beta.r + phi.r*w)
  s <- runif(n,0,10) 
  for(i in 1:n)
  {
    t[i] <- (-log(u[i])/lambda.t[i])^(1/alpha.t)      
    c[i] <- (-log(v[i])/lambda.c[i])^(1/alpha.c) 
    r[i] <- (-log(z[i])/lambda.r[i])^(1/alpha.r)     
    y[i] <- min(t[i], c[i], r[i], s[i])
    e[i] <- as.numeric(t[i] < min(c[i], r[i], s[i]))
    if(t[i] == min(t[i], c[i], r[i], s[i]))
    {
      cens[i] <- 1
    }
    if(c[i] == min(t[i], c[i], r[i], s[i]))
    {
      cens[i] <- 2
    }
    if(r[i] == min(t[i], c[i], r[i], s[i]))
    {
      cens[i] <- 3
    }
    if(s[i] == min(t[i], c[i], s[i]))
    {
      cens[i] <- 4
    }
  }  
  data <- as.data.frame(cbind(u, v, z, t, c,r, y, e, X, cens, ident))
  names(data) <- c("u","v","z","t","c","r","time","event","x1","x2","x3", "x4", "cens", "ident")  
  return(data)
}

###------------------------------------------------------------------------------------------------###
# Functions for calculating the equation U.
###------------------------------------------------------------------------------------------------###

# Estimation for failure times parameters
modelo_param_T2 <-  function(par,t,delta_T, X_T, wi){
  
  pred <- X_T%*%par[3:6]
  pred_T<- exp(pred)*rowMeans(exp(wi))
  exp_T <- (t^(exp(par[1])))*(exp(par[2]))
  
  U <- -sum(delta_T*((par[1]) + (exp(par[1])-1)*log(t) + par[2] + pred + rowMeans(wi))- exp_T*pred_T)
}  

# Estimation for withdrawal times parameters
modelo_param_C2 <-  function(par,t,delta_C, X_C, wi){
  
  pred <- X_C%*%par[3:6]
  pred_C<- exp(pred)*rowMeans(exp(par[7]*wi))
  exp_C <- (t^(exp(par[1])))*(exp(par[2]))
  
  U <- -sum(delta_C*((par[1]) + (exp(par[1])-1)*log(t) + par[2] + pred + par[7]*rowMeans(wi))- exp_C*pred_C)
}    


###------------------------------------------------------------------------------------------------###
# Calculation of first order derivatives
###------------------------------------------------------------------------------------------------###
# assuming dependence
Esp.Deriv.Prim.Ordem <-  function(t,delta_T, delta_C, delta_R, X_T, X_C, X_R, beta_T, beta_C,beta_R, alpha2,alpha3,Sigma2, alpha_T, alpha_C,alpha_R, lambda_T, lambda_C,lambda_R, w_k_grupo, ident){
  wk = w_k_grupo[ident,]  
  num_param <- length(beta_T)+2+length(beta_C)+2+length(beta_R)+2+length(alpha2)+length(alpha3)+length(Sigma2)
  deriv1 <- matrix(NA,num_param,ncol(wk))   
  
  pred_T <- as.vector(exp(X_T%*%beta_T))
  pred_C <- as.vector(exp(X_C%*%beta_C))
  pred_R <- as.vector(exp(X_R%*%beta_R))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  Delta_r <- delta_R%*%t(rep(1,ncol(wk)))
  
  deriv1[1,] <- colSums(Delta_t*(alpha_T^(-1) + log(t)) - (t^(alpha_T))*log(t)*lambda_T*pred_T*exp(wk)) # alpha.t
  deriv1[2,] <- colSums(Delta_t*(lambda_T^(-1)) - (t^(alpha_T))*pred_T*exp(wk)) # lambda.t
  deriv1[3,] <- colSums(Delta_t*(X_T[,1]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,1]*exp(wk)) # beta_T_1
  deriv1[4,] <- colSums(Delta_t*(X_T[,2]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,2]*exp(wk)) # beta_T_2
  deriv1[5,] <- colSums(Delta_t*(X_T[,3]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,3]*exp(wk)) # beta_T_3
  deriv1[6,] <- colSums(Delta_t*(X_T[,4]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*exp(wk)) # beta_T_4
  deriv1[7,] <- colSums(Delta_c*(alpha_C^(-1) + log(t)) - (t^(alpha_C))*log(t)*lambda_C*pred_C*exp(alpha2*wk)) # alpha.c
  deriv1[8,] <- colSums(Delta_c*(lambda_C^(-1)) - (t^(alpha_C))*pred_C*exp(alpha2*wk)) # lambda.c
  deriv1[9,] <- colSums(Delta_c*(X_C[,1]) - (t^(alpha_C))*lambda_C*pred_C*X_C[,1]*exp(alpha2*wk))  # beta_C_1
  deriv1[10,] <- colSums(Delta_c*(X_C[,2]) - (t^(alpha_C))*lambda_C*pred_C*X_C[,2]*exp(alpha2*wk)) # beta_C_2
  deriv1[11,] <- colSums(Delta_c*(X_C[,3]) - (t^(alpha_C))*lambda_C*pred_C*X_C[,3]*exp(alpha2*wk)) # beta_C_3
  deriv1[12,] <- colSums(Delta_c*(X_C[,4]) - (t^(alpha_C))*lambda_C*pred_C*X_C[,4]*exp(alpha2*wk)) # beta_C_4
  deriv1[13,] <- colSums(Delta_c*(wk) - (t^(alpha_C))*lambda_C*pred_C*exp(alpha2*wk)*wk)  # phi.C
  deriv1[14,] <- colSums(Delta_r*(alpha_R^(-1) + log(t)) - (t^(alpha_R))*log(t)*lambda_R*pred_R*exp(alpha3*wk))  # alpha.r
  deriv1[15,] <- colSums(Delta_r*(lambda_R^(-1)) - (t^(alpha_R))*pred_R*exp(alpha3*wk)) # lambda.r
  deriv1[16,] <- colSums(Delta_r*(X_R[,1]) - (t^(alpha_R))*lambda_R*pred_R*X_R[,1]*exp(alpha3*wk))  # beta_R_1
  deriv1[17,] <- colSums(Delta_r*(X_R[,2]) - (t^(alpha_R))*lambda_R*pred_R*X_R[,2]*exp(alpha3*wk))  # beta_R_2
  deriv1[18,] <- colSums(Delta_r*(X_R[,3]) - (t^(alpha_R))*lambda_R*pred_R*X_R[,3]*exp(alpha3*wk))  # beta_R_3
  deriv1[19,] <- colSums(Delta_r*(X_R[,4]) - (t^(alpha_R))*lambda_R*pred_R*X_R[,4]*exp(alpha3*wk))  # beta_R_4
  deriv1[20,] <- colSums(Delta_r*(wk) - (t^(alpha_R))*lambda_R*pred_R*wk*exp(alpha3*wk))  # phi.R
  deriv1[21,] <- colSums(-0.5*Sigma2^(-1) + 0.5*(Sigma2^(-2))*w_k_grupo^(2))  # sigma2
  
  aux <- deriv1[,1]%*%t(deriv1[,1])
  for( i in 2:ncol(wk)){
   aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }
  return(aux/ncol(wk))
}

# assuming independence
Esp.Deriv.Prim.Ordem.indep <-  function(t,delta_T, X_T, beta_T, Sigma2, alpha_T, lambda_T, w_k_grupo, ident){
  wk = w_k_grupo[ident,]  
  num_param <- length(beta_T)+2+length(Sigma2)
  deriv1 <- matrix(NA,num_param,ncol(wk))   
  
  pred_T <- as.vector(exp(X_T%*%beta_T))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  
  deriv1[1,] <- colSums(Delta_t*(alpha_T^(-1) + log(t)) - (t^(alpha_T))*log(t)*lambda_T*pred_T*exp(wk)) # alpha.t
  deriv1[2,] <- colSums(Delta_t*(lambda_T^(-1)) - (t^(alpha_T))*pred_T*exp(wk)) # lambda.t
  deriv1[3,] <- colSums(Delta_t*(X_T[,1]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,1]*exp(wk)) # beta_T_1
  deriv1[4,] <- colSums(Delta_t*(X_T[,2]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,2]*exp(wk)) # beta_T_2
  deriv1[5,] <- colSums(Delta_t*(X_T[,3]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,3]*exp(wk)) # beta_T_3
  deriv1[6,] <- colSums(Delta_t*(X_T[,4]) - (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*exp(wk)) # beta_T_4
  deriv1[7,] <- colSums(-0.5*Sigma2^(-1) + 0.5*(Sigma2^(-2))*w_k_grupo^(2))  # sigma2
  
  aux <- deriv1[,1]%*%t(deriv1[,1])
  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }
  return(aux/ncol(wk))
}

###------------------------------------------------------------------------------------------------###
# Calculation of second order derivatives
###------------------------------------------------------------------------------------------------###
# assuming dependence
Esp.DerivParciais <-  function(t,delta_T, delta_C,delta_R, X_T, X_C,X_R, beta_T, beta_C,beta_R, alpha2,alpha3, Sigma2=sigma2, alpha_T, alpha_C,alpha_R, lambda_T, lambda_C,lambda_R, w_k_grupo, ident){
  
  num_param <- length(beta_T)+2+length(beta_C)+2+length(beta_R)+2+length(alpha2)+length(alpha3)+length(Sigma2)
  deriv2 <- matrix(0,num_param,num_param)  
  wk = w_k_grupo[ident,]   
  
  pred_T <- as.vector(exp(X_T%*%beta_T))
  pred_C <- as.vector(exp(X_C%*%beta_C))
  pred_R <- as.vector(exp(X_R%*%beta_R))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  Delta_r <- delta_R%*%t(rep(1,ncol(wk)))
  
  deriv2[1,1] <- mean(colSums(Delta_t*(-alpha_T^(-2)) - (t^(alpha_T))*(log(t)^(2))*lambda_T*pred_T*exp(wk))) # derived from alpha.t derived in relation to alpha.t
  deriv2[1,2] <- mean(colSums(- (t^(alpha_T))*(log(t))*pred_T*exp(wk)))  # derived from alpha.t derived in relation to lambda.t
  deriv2[2,1] <- deriv2[1,2]
  deriv2[1,3] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,1]*exp(wk))) # derived from alpha.t derived in relation to beta_T_1
  deriv2[3,1] <- deriv2[1,3]
  deriv2[1,4] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,2]*exp(wk))) # derived from alpha.t derived in relation to beta_T_2
  deriv2[4,1] <- deriv2[1,4]
  deriv2[1,5] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,3]*exp(wk))) # derived from alpha.t derived in relation to beta_T_3
  deriv2[5,1] <- deriv2[1,5]
  deriv2[1,6] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,4]*exp(wk))) # derived from alpha.t derived in relation to beta_T_4
  deriv2[6,1] <- deriv2[1,6]
  deriv2[2,2] <- mean(colSums( Delta_t*(-lambda_T^(-2))))  # derived from lambda.t derived in relation to lambda.t
  deriv2[2,3] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,1]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_1
  deriv2[3,2] <- deriv2[2,3]
  deriv2[2,4] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,2]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_2
  deriv2[4,2] <- deriv2[2,4]
  deriv2[2,5] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,3]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_3
  deriv2[5,2] <- deriv2[2,5]
  deriv2[2,6] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,4]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_4
  deriv2[6,2] <- deriv2[2,6]
  deriv2[3,3] <- mean(colSums( - ((t^(alpha_T))*lambda_T*pred_T*(X_T[,1]^(2))*exp(wk)))) # derived from beta_T_1 derived in relation to beta_T_1
  deriv2[3,4] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,2]*X_T[,1]*exp(wk)))  # derived from beta_T_1 derived in relation to beta_T_2
  deriv2[4,3] <- deriv2[3,4]
  deriv2[3,5] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,3]*X_T[,1]*exp(wk)))  # derived from beta_T_1 derived in relation to beta_T_3
  deriv2[5,3] <- deriv2[3,5]
  deriv2[3,6] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*X_T[,1]*exp(wk)))  # derived from beta_T_1 derived in relation to beta_T_4
  deriv2[6,3] <- deriv2[3,6]
  deriv2[4,4] <- mean(colSums(- ( (t^(alpha_T))*lambda_T*pred_T*(X_T[,2]^(2))*exp(wk)))) # derived from beta_T_2 derived in relation to beta_T_2
  deriv2[4,5] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,3]*X_T[,2]*exp(wk)))  # derived from beta_T_2 derived in relation to beta_T_3
  deriv2[5,4] <- deriv2[4,5]
  deriv2[4,6] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*X_T[,2]*exp(wk)))  # derived from beta_T_2 derived in relation to beta_T_4
  deriv2[6,4] <- deriv2[4,6]
  deriv2[5,5] <- mean(colSums(- ( (t^(alpha_T))*lambda_T*pred_T*(X_T[,3]^(2))*exp(wk)))) # derived from beta_T_3 derived in relation to beta_T_3
  deriv2[5,6] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*X_T[,3]*exp(wk)))  # derived from beta_T_3 derived in relation to beta_T_4
  deriv2[6,5] <- deriv2[5,6]
  deriv2[6,6] <- mean(colSums(- ( (t^(alpha_T))*lambda_T*pred_T*(X_T[,4]^(2))*exp(wk)))) # derived from beta_T_4 derived in relation to beta_T_4
  deriv2[7,7] <- mean(colSums(Delta_c*(-alpha_C^(-2)) - (t^(alpha_C))*(log(t)^(2))*lambda_C*pred_C*exp(alpha2*wk)))  # derived from alpha.c derived in relation to alpha.c
  deriv2[7,8] <- mean(colSums( - (t^(alpha_C))*log(t)*pred_C*exp(alpha2*wk)))  # derived from alpha.c derived in relation to lambda.c
  deriv2[8,7] <- deriv2[7,8]
  deriv2[7,9] <- mean(colSums( - (t^(alpha_C))*log(t)*lambda_C*pred_C*X_C[,1]*exp(alpha2*wk)))  # derived from alpha.c derived in relation to beta_C_1
  deriv2[9,7] <- deriv2[7,9]
  deriv2[7,10] <- mean(colSums(- (t^(alpha_C))*log(t)*lambda_C*pred_C*X_C[,2]*exp(alpha2*wk)))  # derived from alpha.c derived in relation to beta_C_2
  deriv2[10,7] <- deriv2[7,10]
  deriv2[7,11] <- mean(colSums(- (t^(alpha_C))*log(t)*lambda_C*pred_C*X_C[,3]*exp(alpha2*wk)))  # derived from alpha.c derived in relation to beta_C_3
  deriv2[11,7] <- deriv2[7,11]
  deriv2[7,12] <- mean(colSums(- (t^(alpha_C))*log(t)*lambda_C*pred_C*X_C[,4]*exp(alpha2*wk)))  # derived from alpha.c derived in relation to beta_C_4
  deriv2[12,7] <- deriv2[7,12]
  deriv2[7,13] <- mean(colSums( - (t^(alpha_C))*log(t)*lambda_C*pred_C*exp(alpha2*wk)*wk))  # derived from alpha.c derived in relation to phi.C
  deriv2[13,7] <- deriv2[7,13]
  deriv2[8,8] <- mean(colSums(Delta_c*(-lambda_C^(-2)) ))   # derived from lambda.c derived in relation to lambda.c
  deriv2[8,9] <- mean(colSums( - (t^(alpha_C))*pred_C*X_C[,1]*exp(alpha2*wk)))  # derived from lambda.c derived in relation to beta_C_1
  deriv2[9,8] <- deriv2[8,9]
  deriv2[8,10] <- mean(colSums( - (t^(alpha_C))*pred_C*X_C[,2]*exp(alpha2*wk)))  # derived from lambda.c derived in relation to beta_C_2
  deriv2[10,8] <- deriv2[8,10]
  deriv2[8,11] <- mean(colSums( - (t^(alpha_C))*pred_C*X_C[,3]*exp(alpha2*wk)))  # derived from lambda.c derived in relation to beta_C_3
  deriv2[11,8] <- deriv2[8,11]
  deriv2[8,12] <- mean(colSums( - (t^(alpha_C))*pred_C*X_C[,4]*exp(alpha2*wk)))                   # derived from lambda.c derived in relation to beta_C_4
  deriv2[12,8] <- deriv2[8,12]
  deriv2[8,13] <- mean(colSums( - (t^(alpha_C))*pred_C*exp(alpha2*wk)*wk))                        # derived from lambda.c derived in relation to phi.C
  deriv2[13,8] <- deriv2[8,13]
  deriv2[9,9] <- mean(colSums(-((t^(alpha_C))*lambda_C*pred_C*(X_C[,1]^(2))*exp(alpha2*wk))))     # derived from beta_C_1 derived in relation to beta_C_1
  deriv2[9,10] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,1]*X_C[,2]*exp(alpha2*wk)))  # derived from beta_C_1 derived in relation to beta_C_2
  deriv2[10,9] <- deriv2[9,10]
  deriv2[9,11] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,1]*X_C[,3]*exp(alpha2*wk)))  # derived from beta_C_1 derived in relation to beta_C_3
  deriv2[11,9] <- deriv2[9,11]
  deriv2[9,12] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,1]*X_C[,4]*exp(alpha2*wk)))  # derived from beta_C_1 derived in relation to beta_C_4
  deriv2[12,9] <- deriv2[9,12]
  deriv2[9,13] <- mean(colSums(- t^(alpha_C)*lambda_C*pred_C*X_C[,1]*exp(alpha2*wk)*wk))          # derived from beta_C_1 derived in relation to phi.C
  deriv2[13,9] <- deriv2[9,13]
  deriv2[10,10] <- mean(colSums( -((t^(alpha_C))*lambda_C*pred_C*(X_C[,2]^(2))*exp(alpha2*wk))))  # derived from beta_C_2 derived in relation to beta_C_2
  deriv2[10,11] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,2]*X_C[,3]*exp(alpha2*wk))) # derived from beta_C_2 derived in relation to beta_C_3
  deriv2[11,10] <- deriv2[10,11]
  deriv2[10,12] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,2]*X_C[,4]*exp(alpha2*wk))) # derived from beta_C_2 derived in relation to beta_C_4
  deriv2[12,10] <- deriv2[10,12]
  deriv2[10,13] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,2]*exp(alpha2*wk)*wk))      # derived from beta_C_2 derived in relation to phi.C
  deriv2[13,10] <- deriv2[10,13] 
  deriv2[11,11] <- mean(colSums( -((t^(alpha_C))*lambda_C*pred_C*(X_C[,3]^(2))*exp(alpha2*wk))))  # derived from beta_C_3 derived in relation to beta_C_3
  deriv2[11,12] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,3]*X_C[,4]*exp(alpha2*wk))) # derived from beta_C_3 derived in relation to beta_C_4
  deriv2[12,11] <- deriv2[11,12]
  deriv2[11,13] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,3]*exp(alpha2*wk)*wk))      # derived from beta_C_3 derived in relation to phi.C
  deriv2[13,11] <- deriv2[11,13] 
  deriv2[12,12] <- mean(colSums( -((t^(alpha_C))*lambda_C*pred_C*(X_C[,4]^(2))*exp(alpha2*wk))))  # derived from beta_C_4 derived in relation to beta_C_4
  deriv2[12,13] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*X_C[,4]*exp(alpha2*wk)*wk))      # derived from beta_C_4 derived in relation to phi.C
  deriv2[13,12] <- deriv2[12,13] 
  deriv2[13,13] <- mean(colSums( - (t^(alpha_C))*lambda_C*pred_C*exp(alpha2*wk)*wk^(2)))          # derived from phi.C derived in relation to phi.C
  deriv2[14,14] <- mean(colSums(Delta_r*(-alpha_R^(-2)) - (t^(alpha_R))*(log(t)^(2))*lambda_R*pred_R*exp(alpha3*wk)))  # derived from alpha.r derived in relation to alpha.r
  deriv2[14,15] <- mean(colSums( - (t^(alpha_R))*log(t)*pred_R*exp(alpha3*wk)))   # derived from alpha.r derived in relation to lambda.r
  deriv2[15,14] <- deriv2[14,15]
  deriv2[14,16] <- mean(colSums( - (t^(alpha_R))*log(t)*lambda_R*pred_R*X_R[,1]*exp(alpha3*wk)))  # derived from alpha.r derived in relation to beta_R_1
  deriv2[16,14] <- deriv2[14,16]
  deriv2[14,17] <- mean(colSums(- (t^(alpha_R))*log(t)*lambda_R*pred_R*X_R[,2]*exp(alpha3*wk)))   # derived from alpha.r derived in relation to beta_R_2
  deriv2[17,14] <- deriv2[14,17]
  deriv2[14,18] <- mean(colSums(- (t^(alpha_R))*log(t)*lambda_R*pred_R*X_R[,3]*exp(alpha3*wk)))   # derived from alpha.r derived in relation to beta_R_3
  deriv2[18,14] <- deriv2[14,18]
  deriv2[14,19] <- mean(colSums(- (t^(alpha_R))*log(t)*lambda_R*pred_R*X_R[,4]*exp(alpha3*wk)))   # derived from alpha.r derived in relation to beta_R_4
  deriv2[19,14] <- deriv2[14,19]
  deriv2[14,20] <- mean(colSums(- (t^(alpha_R))*log(t)*lambda_R*pred_R*wk*exp(alpha3*wk)))        # derived from alpha.r derived in relation to phi.R
  deriv2[20,14] <- deriv2[14,20]
  deriv2[15,15] <- mean(colSums( Delta_r*(-lambda_R^(-2))))   # derived from lambda.r derived in relation to lambda.r
  deriv2[15,16] <- mean(colSums( - (t^(alpha_R))*pred_R*X_R[,1]*exp(alpha3*wk)))  # derived from lambda.r derived in relation to beta_R_1
  deriv2[16,15] <- deriv2[15,16]
  deriv2[15,17] <- mean(colSums( - (t^(alpha_R))*pred_R*X_R[,2]*exp(alpha3*wk)))  # derived from lambda.r derived in relation to beta_R_2
  deriv2[17,15] <- deriv2[15,17]
  deriv2[15,18] <- mean(colSums( - (t^(alpha_R))*pred_R*X_R[,3]*exp(alpha3*wk)))  # derived from lambda.r derived in relation to beta_R_3
  deriv2[18,15] <- deriv2[15,18]
  deriv2[15,19] <- mean(colSums( - (t^(alpha_R))*pred_R*X_R[,4]*exp(alpha3*wk)))  # derived from lambda.r derived in relation to beta_R_4
  deriv2[19,15] <- deriv2[15,19]
  deriv2[15,20] <- mean(colSums( - (t^(alpha_R))*pred_R*wk*exp(alpha3*wk)))       # derived from lambda.r derived in relation to phi.R
  deriv2[20,15] <- deriv2[15,20]
  deriv2[16,16] <- mean(colSums(-((t^(alpha_R))*lambda_R*pred_R*(X_R[,1]^(2))*exp(alpha3*wk))))   # derived from beta_R_1 derived in relation to beta_R_1
  deriv2[16,17] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,1]*X_R[,2]*exp(alpha3*wk))) # derived from beta_R_1 derived in relation to beta_R_2
  deriv2[17,16] <- deriv2[16,17]
  deriv2[16,18] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,1]*X_R[,3]*exp(alpha3*wk))) # derived from beta_R_1 derived in relation to beta_R_3
  deriv2[18,16] <- deriv2[16,18]
  deriv2[16,19] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,1]*X_R[,4]*exp(alpha3*wk))) # derived from beta_R_1 derived in relation to beta_R_4
  deriv2[19,16] <- deriv2[16,19]
  deriv2[16,20] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,1]*wk*exp(alpha3*wk)))      # derived from beta_R_1 derived in relation to phi.R
  deriv2[20,16] <- deriv2[16,20]
  deriv2[17,17] <- mean(colSums(-((t^(alpha_R))*lambda_R*pred_R*(X_R[,2]^(2))*exp(alpha3*wk))))   # derived from beta_R_2 derived in relation to beta_R_2
  deriv2[17,18] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,2]*X_R[,3]*exp(alpha3*wk))) # derived from beta_R_2 derived in relation to beta_R_3
  deriv2[18,17] <- deriv2[17,18]
  deriv2[17,19] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,2]*X_R[,4]*exp(alpha3*wk))) # derived from beta_R_2 derived in relation to beta_R_4
  deriv2[19,17] <- deriv2[17,19]
  deriv2[17,20] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,2]*wk*exp(alpha3*wk)))      # derived from beta_R_2 derived in relation to phi.R
  deriv2[20,17] <- deriv2[17,20]
  deriv2[18,18] <- mean(colSums(-((t^(alpha_R))*lambda_R*pred_R*(X_R[,3]^(2))*exp(alpha3*wk))))   # derived from beta_R_3 derived in relation to beta_R_3
  deriv2[18,19] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,3]*X_R[,4]*exp(alpha3*wk))) # derived from beta_R_3 derived in relation to beta_R_4
  deriv2[19,18] <- deriv2[18,19]
  deriv2[18,20] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,3]*wk*exp(alpha3*wk)))      # derived from beta_R_3 derived in relation to phi.R
  deriv2[20,18] <- deriv2[18,20]
  deriv2[19,19] <- mean(colSums(-((t^(alpha_R))*lambda_R*pred_R*(X_R[,4]^(2))*exp(alpha3*wk))))   # derived from beta_R_4 derived in relation to beta_R_4
  deriv2[19,20] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*X_R[,4]*wk*exp(alpha3*wk)))      # derived from beta_R_4 derived in relation to phi.R
  deriv2[20,19] <- deriv2[19,20]
  deriv2[20,20] <- mean(colSums( - (t^(alpha_R))*lambda_R*pred_R*exp(alpha3*wk)*wk^(2)))          # derived from phi.R derived in relation to phi.R
  deriv2[21,21] <- mean(colSums( 0.5*Sigma2^(-2) - (Sigma2^(-3))*w_k_grupo^(2)))                  # derived from sigma2 derived in relation to sigma2
  
  return((as.matrix(deriv2)))
}

# assuming independence
Esp.DerivParciais.indep <-function(t,delta_T, X_T, beta_T, Sigma2=sigma2, alpha_T, lambda_T, w_k_grupo, ident){
  
  num_param <- length(beta_T)+2+length(Sigma2)
  deriv2 <- matrix(0,num_param,num_param)  
  wk = w_k_grupo[ident,]   
  
  pred_T <- as.vector(exp(X_T%*%beta_T))
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  
  deriv2[1,1] <- mean(colSums(Delta_t*(-alpha_T^(-2)) - (t^(alpha_T))*(log(t)^(2))*lambda_T*pred_T*exp(wk))) # derived from alpha.t derived in relation to alpha.t
  deriv2[1,2] <- mean(colSums(- (t^(alpha_T))*(log(t))*pred_T*exp(wk)))  # derived from alpha.t derived in relation to lambda.t
  deriv2[2,1] <- deriv2[1,2]
  deriv2[1,3] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,1]*exp(wk))) # derived from alpha.t derived in relation to beta_T_1
  deriv2[3,1] <- deriv2[1,3]
  deriv2[1,4] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,2]*exp(wk))) # derived from alpha.t derived in relation to beta_T_2
  deriv2[4,1] <- deriv2[1,4]
  deriv2[1,5] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,3]*exp(wk))) # derived from alpha.t derived in relation to beta_T_3
  deriv2[5,1] <- deriv2[1,5]
  deriv2[1,6] <- mean(colSums(- (t^(alpha_T))*(log(t))*lambda_T*pred_T*X_T[,4]*exp(wk))) # derived from alpha.t derived in relation to beta_T_4
  deriv2[6,1] <- deriv2[1,6]
  deriv2[2,2] <- mean(colSums( Delta_t*(-lambda_T^(-2))))  # derived from lambda.t derived in relation to lambda.t
  deriv2[2,3] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,1]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_1
  deriv2[3,2] <- deriv2[2,3]
  deriv2[2,4] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,2]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_2
  deriv2[4,2] <- deriv2[2,4]
  deriv2[2,5] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,3]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_3
  deriv2[5,2] <- deriv2[2,5]
  deriv2[2,6] <- mean(colSums( - (t^(alpha_T))*pred_T*X_T[,4]*exp(wk)))  # derived from lambda.t derived in relation to beta_T_4
  deriv2[6,2] <- deriv2[2,6]
  deriv2[3,3] <- mean(colSums( - ((t^(alpha_T))*lambda_T*pred_T*(X_T[,1]^(2))*exp(wk)))) # derived from beta_T_1 derived in relation to beta_T_1
  deriv2[3,4] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,2]*X_T[,1]*exp(wk)))  # derived from beta_T_1 derived in relation to beta_T_2
  deriv2[4,3] <- deriv2[3,4]
  deriv2[3,5] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,3]*X_T[,1]*exp(wk)))  # derived from beta_T_1 derived in relation to beta_T_3
  deriv2[5,3] <- deriv2[3,5]
  deriv2[3,6] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*X_T[,1]*exp(wk)))  # derived from beta_T_1 derived in relation to beta_T_4
  deriv2[6,3] <- deriv2[3,6]
  deriv2[4,4] <- mean(colSums(- ( (t^(alpha_T))*lambda_T*pred_T*(X_T[,2]^(2))*exp(wk)))) # derived from beta_T_2 derived in relation to beta_T_2
  deriv2[4,5] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,3]*X_T[,2]*exp(wk)))  # derived from beta_T_2 derived in relation to beta_T_3
  deriv2[5,4] <- deriv2[4,5]
  deriv2[4,6] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*X_T[,2]*exp(wk)))  # derived from beta_T_2 derived in relation to beta_T_4
  deriv2[6,4] <- deriv2[4,6]
  deriv2[5,5] <- mean(colSums(- ( (t^(alpha_T))*lambda_T*pred_T*(X_T[,3]^(2))*exp(wk)))) # derived from beta_T_3 derived in relation to beta_T_3
  deriv2[5,6] <- mean(colSums(- (t^(alpha_T))*lambda_T*pred_T*X_T[,4]*X_T[,3]*exp(wk)))  # derived from beta_T_3 derived in relation to beta_T_4
  deriv2[6,5] <- deriv2[5,6]
  deriv2[6,6] <- mean(colSums(- ( (t^(alpha_T))*lambda_T*pred_T*(X_T[,4]^(2))*exp(wk)))) # derived from beta_T_4 derived in relation to beta_T_4
  deriv2[7,7] <- mean(colSums( 0.5*Sigma2^(-2) - (Sigma2^(-3))*w_k_grupo^(2)))           # derived from sigma2 derived in relation to sigma2
  
  return((as.matrix(deriv2)))
}
###------------------------------------------------------------------------------------------------###
# Weibull Proposed model
###------------------------------------------------------------------------------------------------###

Weibull.dep  <- function(dados){
  n <- nrow(dados)
  o <- order(dados$time)               # sort the data
  dados <- dados[o,]   
  ident <- dados$ident
  m <- max(ident)                      # number of clusters
  age <- dados$x1
  gender <- dados$x2                   # 1 is male and 0 is female
  race <- dados$x3                     # 1 if it is black and 0 if it is other
  diabetic <- dados$x4                 # 1 if is diabetic and 0 otherwise
  time <- dados$time                   # time in years
  
  delta_t <- ifelse(dados$cens==1,1,0) # if it is failure
  delta_c <- ifelse(dados$cens==2,1,0) # if it is withdrawal
  delta_r <- ifelse(dados$cens==3,1,0) # if it is transplant
  delta_s <- ifelse(dados$cens==4,1,0) # if it is administrative censoring 
  
  X <- cbind(scale(as.numeric(age)), as.numeric(gender),as.numeric(race), as.numeric(diabetic))
  p <- ncol(X)
  ###--------------------------------------------------------------------------------------------###
  # Specification of the initial values for the parameters.
  ###--------------------------------------------------------------------------------------------###
  
  ajuste_coxph_T <- coxph(Surv(time, delta_t) ~  X[,1] + X[,2] + X[,3] + X[,4], method="breslow")
  risco_a_T      <- basehaz(ajuste_coxph_T, centered=FALSE)$hazard      
  beta_T         <- ajuste_coxph_T$coef
  beta_T0        <- beta_T
  
  ajuste_coxph_C <- coxph(Surv(time, delta_c) ~  X[,1] + X[,2] + X[,3] + X[,4], method="breslow")
  risco_a_C      <- basehaz(ajuste_coxph_C, centered=FALSE)$hazard  
  beta_C         <- ajuste_coxph_C$coef
  beta_C0        <- beta_C
  
  ajuste_coxph_R <- coxph(Surv(time, delta_r) ~  X[,1] + X[,2] + X[,3] + X[,4], method="breslow")
  risco_a_R      <- basehaz(ajuste_coxph_R, centered=FALSE)$hazard  
  beta_R         <- ajuste_coxph_R$coef
  beta_R0        <- beta_R
  
  phi.C <- 0; phi.R <- 0; Sigma2 <- 1; alpha_T <- 1; lambda_T <- 1; alpha_C <- 1; lambda_C <- 1;
  alpha_R <- 1;  lambda_R <- 1
  param <- c(beta_T,beta_C,phi.C,beta_R,phi.R,Sigma2) 
  ###--------------------------------------------------------------------------------------------###
  # EMMC algorithm specifications
  ###--------------------------------------------------------------------------------------------###
  maxit <- 100                             # maximum number of iterations
  eps1= rep(1e-3, length(param))
  eps2= rep(1e-4, length(param))
  n_intMCs = c(rep(10,20),rep(25,30), rep(50,40), rep(75,10))
  out =  matrix(NA,maxit+1,length(param))     
  dif =  matrix(NA,maxit+1,length(param))   
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE
  X_T <- X; N_T <- length(X_T[,1])
  X_C <- X; N_C <- length(X_C[,1])
  X_R <- X; N_R <- length(X_R[,1])
  ###--------------------------------------------------------------------------------------------###
  # Initiation of the EMMC algorithm 
  ###--------------------------------------------------------------------------------------------###
  while (continue == TRUE) {
    count = rep(0,maxit+1)
    out[s,] = c(beta_T,beta_C,phi.C,beta_R,phi.R,Sigma2) 
    n_intMC = n_intMCs[s]
    w_chapeu_grupo <- matrix(NA, m, n_intMC)
    ###------------------------------------------------------------------------------------
    # calculation of frailty for each cluster 
    for ( k in 1:m){  # k refers to the k-th cluster                                                    
      w_trans <- arms(0.5, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,phi_2=phi.C,phi_3=phi.R,delta_T=delta_t[ident==k],delta_C=delta_c[ident==k],delta_R=delta_r[ident==k], X_T=X_T[ident==k,],X_C=X_C[ident==k,],X_R=X_R[ident==k,],Betas_T=beta_T,Betas_C=beta_C,Betas_R=beta_R,risco_a_T=risco_a_T[ident==k],risco_a_C=risco_a_C[ident==k],risco_a_R=risco_a_R[ident==k],Sigma2=Sigma2) 
      w_auxi <- log(w_trans)-log(1-w_trans) 
      w_chapeu_grupo[k,] <- w_auxi
    }
    wi = w_chapeu_grupo[ident,]
    Sigma2 <- mean(w_chapeu_grupo^2)  
    ###------------------------------------------------------------------------------------
    # Estimate for failure time	parameters
    param_T   <- optim(par = c(0.01,0.01,beta_T0), fn = modelo_param_T2,NULL, control = list(maxit = 5000), method="BFGS",t=time,delta_T=delta_t, X_T=X_T, wi=wi)
    par_T     <- param_T$par   
    alpha_T   <- exp(par_T[1])
    lambda_T  <- exp(par_T[2])
    beta_T    <- par_T[3:6]
    risco_a_T <- (time^alpha_T)*(lambda_T)
    ###-----------------------------------------------------------------------------------  
    # Estimate for withdrawal time	parameters 
    param_C   <- optim(par = c(0.01,0.01,beta_C0,0.01), fn = modelo_param_C2,NULL, control = list(maxit = 5000), method = "BFGS",t=time,delta_C=delta_c, X_C=X_C, wi=wi)
    par_C     <- param_C$par   
    alpha_C   <- exp(par_C[1])
    lambda_C  <- exp(par_C[2])
    beta_C    <- par_C[3:6]
    phi.C     <- par_C[7]
    risco_a_C <- (time^alpha_C)*(lambda_C)
    ###-----------------------------------------------------------------------------------  
    # Estimate for transplant time	parameters 
    # The estimation for the transplant parameters has the same calculations as the withdrawal estimation. So, can use the same function
    # of withdrawal, but with the transplant variables.  
    param_R   <- optim(par = c(0.01,0.01,beta_R0,0.01), fn = modelo_param_C2,gr=NULL, control = list(maxit = 5000), method = "BFGS",t=time,delta_C=delta_r, X_C=X_R, wi=wi)
    par_R     <- param_R$par   
    alpha_R   <- exp(par_R[1])
    lambda_R  <- exp(par_R[2])
    beta_R    <- par_R[3:6]
    phi.R     <- par_R[7]
    risco_a_R <- (time^alpha_R)*(lambda_R)
    ###------------------------------------------------------------------------------------
    # stop criterion specification
    out[s+1,]=  c(beta_T,beta_C,phi.C,beta_R,phi.R,Sigma2) 
    dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)
    for (p in 1: length(param)){
      if (dif[s,p]<eps2[p]) {
        count[s] = count[s] + 1
      }
    }
    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param))) 
    }
  } 
  ###--------------------------------------------------------------------------------------------###
  # End of the EMMC algorithm
  ###--------------------------------------------------------------------------------------------###
  # final estimative for the parameters
  beta_T <- apply(out[(s-2):s,1:4], 2, mean); 
  Sigma2 <- mean(out[(s-2):s,15]); 
  beta_C <- apply(out[(s-2):s,5:8], 2, mean); 
  phi.C  <- mean(out[(s-2):s,9]); 
  beta_R <- apply(out[(s-2):s,10:13], 2, mean);
  phi.R  <- mean(out[(s-2):s,14]); 
  Estimate <- c(beta_T,Sigma2,beta_C,phi.C,beta_R,phi.R)
  ###--------------------------------------------------------------------------------------------###
  # Standard errors calculation  
  ###--------------------------------------------------------------------------------------------###
  Esp_deriv_ordem1 <- Esp.Deriv.Prim.Ordem(t=time,delta_T=delta_t, delta_C=delta_c,delta_R=delta_r, X_T=X_T, X_C=X_C,X_R=X_R,beta_T=beta_T,beta_C=beta_C,beta_R=beta_R,
                                           alpha2=phi.C, alpha3=phi.R, Sigma2=Sigma2, alpha_T=alpha_T, alpha_C=alpha_C,alpha_R=alpha_R, 
                                           lambda_T=lambda_T, lambda_C=lambda_C,lambda_R=lambda_R, w_k_grupo=w_chapeu_grupo, ident=ident)
  
  Esp_deriv_ordem2 <- Esp.DerivParciais(t=time,delta_T=delta_t, delta_C=delta_c, delta_R=delta_r, X_T=X_T, X_C=X_C,X_R=X_R,beta_T=beta_T,beta_C=beta_C,beta_R=beta_R,
                                        alpha2=phi.C, alpha3=phi.R, Sigma2=Sigma2, alpha_T=alpha_T, alpha_C=alpha_C,alpha_R=alpha_R,
                                        lambda_T=lambda_T, lambda_C=lambda_C,lambda_R=lambda_R, w_k_grupo=w_chapeu_grupo, ident=ident)
  
  
  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  EP_louis <- sqrt(diag(Var)) 
  #---------------------------------------------------------------------------------------------------------------
  # Standard errors for the parameters of interest
  p <- ncol(X)
  SE.beta_T <- EP_louis[(2+1):(2+p)];  
  SE.Sigma2 <- EP_louis[2+p+2+p+1+2+p+1+1];
  SE.beta_C <- EP_louis[(2+p+2+1):(2+p+2+p)]; 
  SE.phi.C  <- EP_louis[2+p+2+p+1]; 
  SE.beta_R <- EP_louis[(2+p+2+p+1+2+1):(2+p+2+p+1+2+p)]; 
  SE.phi.R  <- EP_louis[2+p+2+p+1+2+p+1]; 
  SE <- c(SE.beta_T,SE.Sigma2,SE.beta_C,SE.phi.C,SE.beta_R,SE.phi.R)
  
  return(cbind(round(Estimate,3),round(SE,3)))
}
  
###------------------------------------------------------------------------------------------------###
# Weibull Standard model
###------------------------------------------------------------------------------------------------###

Weibull.indep  <- function(dados){
  n <- nrow(dados)
  o <- order(dados$time)               # sort the data
  dados <- dados[o,]   
  ident <- dados$ident
  m <- max(ident)                      # number of clusters
  age <- dados$x1
  gender <- dados$x2                   # 1 is male and 0 is female
  race <- dados$x3                     # 1 if it is black and 0 if it is other
  diabetic <- dados$x4                 # 1 if is diabetic and 0 otherwise
  time <- dados$time                   # time in years
  
  delta_t <- ifelse(dados$cens==1,1,0) # if it is failure
  delta_c <- ifelse(dados$cens==2,1,0) # if it is withdrawal
  delta_r <- ifelse(dados$cens==3,1,0) # if it is transplant
  delta_s <- ifelse(dados$cens==4,1,0) # if it is administrative censoring 
  
  X <- cbind(scale(as.numeric(age)), as.numeric(gender),as.numeric(race), as.numeric(diabetic))
  p <- ncol(X)
  ###--------------------------------------------------------------------------------------------###
  # Specification of the initial values for the parameters.
  ###--------------------------------------------------------------------------------------------###
  
  ajuste_coxph_T <- coxph(Surv(time, delta_t) ~  X[,1] + X[,2] + X[,3] + X[,4], method="breslow")
  risco_a_T      <- basehaz(ajuste_coxph_T, centered=FALSE)$hazard      
  beta_T         <- ajuste_coxph_T$coef
  beta_T0        <- beta_T
  
  Sigma2 <- 1; alpha_T <- 1; lambda_T <- 1; 
  param <- c(beta_T,Sigma2) 
  ###--------------------------------------------------------------------------------------------###
  # EMMC algorithm specifications
  ###--------------------------------------------------------------------------------------------###
  maxit <- 100                             # maximum number of iterations
  eps1= rep(1e-3, length(param))
  eps2= rep(1e-4, length(param))
  n_intMCs = c(rep(10,20),rep(25,30), rep(50,40), rep(75,10))
  out =  matrix(NA,maxit+1,length(param))     
  dif =  matrix(NA,maxit+1,length(param))   
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE
  X_T <- X; N_T <- length(X_T[,1])
  ###--------------------------------------------------------------------------------------------###
  # Initiation of the EMMC algorithm 
  ###--------------------------------------------------------------------------------------------###
  while (continue == TRUE) {
    count = rep(0,maxit+1)
    out[s,] = c(beta_T,Sigma2) 
    n_intMC = n_intMCs[s]
    w_chapeu_grupo <- matrix(NA, m, n_intMC)
    ###------------------------------------------------------------------------------------
    # calculation of frailty for each cluster 
    for ( k in 1:m){  # k refers to the k-th cluster                                                    
      w_trans <- arms(0.5, myldens=log_cond_wi_indep, indFunc=support_wi_indep,n.sample=n_intMC,delta_T=delta_t[ident==k],X_T=X_T[ident==k,],Betas_T=beta_T,risco_a_T=risco_a_T[ident==k],Sigma2=Sigma2) 
      w_auxi <- log(w_trans)-log(1-w_trans) 
      w_chapeu_grupo[k,] <- w_auxi
    }
    wi = w_chapeu_grupo[ident,]
    Sigma2 <- mean(w_chapeu_grupo^2)  
    ###------------------------------------------------------------------------------------
    # Estimate for failure time	parameters
    param_T   <- optim(par = c(0.01,0.01,beta_T0), fn = modelo_param_T2,NULL, control = list(maxit = 5000), method="BFGS",t=time,delta_T=delta_t, X_T=X_T, wi=wi)
    par_T     <- param_T$par   
    alpha_T   <- exp(par_T[1])
    lambda_T  <- exp(par_T[2])
    beta_T    <- par_T[3:6]
    risco_a_T <- (time^alpha_T)*(lambda_T)
    ###-----------------------------------------------------------------------------------  
    # stop criterion specification
    out[s+1,]=  c(beta_T,Sigma2) 
    dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)
    for (p in 1: length(param)){
      if (dif[s,p]<eps2[p]) {
        count[s] = count[s] + 1
      }
    }
    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param))) 
    }
  } 
  ###--------------------------------------------------------------------------------------------###
  # End of the EMMC algorithm
  ###--------------------------------------------------------------------------------------------###
  # final estimative for the parameters
  beta_T <- apply(out[(s-2):s,1:4], 2, mean); 
  Sigma2 <- mean(out[(s-2):s,5]); 
  Estimate <- c(beta_T,Sigma2)
  ###--------------------------------------------------------------------------------------------###
  # Standard errors calculation  
  ###--------------------------------------------------------------------------------------------###
  Esp_deriv_ordem1 <- Esp.Deriv.Prim.Ordem.indep(t=time,delta_T=delta_t,X_T=X_T,beta_T=beta_T,Sigma2=Sigma2,
                                                 alpha_T=alpha_T,lambda_T=lambda_T, w_k_grupo=w_chapeu_grupo, ident=ident)
  
  Esp_deriv_ordem2 <- Esp.DerivParciais.indep(t=time,delta_T=delta_t,X_T=X_T,beta_T=beta_T, Sigma2=Sigma2, 
                                              alpha_T=alpha_T,lambda_T=lambda_T,w_k_grupo=w_chapeu_grupo, ident=ident)
  
  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  EP_louis <- sqrt(diag(Var)) 
  #---------------------------------------------------------------------------------------------------------------
  # Standard errors for the parameters of interest
  p <- ncol(X)
  SE.beta_T <- EP_louis[(2+1):(2+p)];  
  SE.Sigma2 <- EP_louis[2+p+1];
  SE <- c(SE.beta_T,SE.Sigma2)
  
  return(cbind(round(Estimate,3),round(SE,3)))
}

  