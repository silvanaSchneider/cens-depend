
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

###-------------------------------------------------------------------------------------------------###
# specification of support for the arms function
###-------------------------------------------------------------------------------------------------###
# assuming dependence
support_wi <-  function(w_k_aux,phi_2,phi_3,delta_T,delta_C,delta_R,X_T,X_C, X_R, Betas_T,Betas_C,Betas_R,risco_a_T,risco_a_C,risco_a_R, Sigma2){(w_k_aux>0)*(w_k_aux<1)}  

# assuming independence
support_wi_indep <-  function(w_k_aux,delta_T,X_T,Betas_T,risco_a_T,Sigma2){(w_k_aux>0)*(w_k_aux<1)}  

###------------------------------------------------------------------------------------------------###
# Function that calculates the grid times.
###------------------------------------------------------------------------------------------------###
time.grid <- function(time, event, n.int=NULL)
{
  o <- order(time)  
  time <- time[o]    
  event <- event[o]
  time.aux <- unique(time[event==1])
  if(is.null(n.int))
  {
    n.int <- length(time.aux)
  }
  m <- length(time.aux)
  if(n.int > m)
  {
    a <- c(0,unique(time[event==1]))
    a[length(a)] <- Inf
  }
  else
  {
    b <- min(m,n.int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    a_inf <- c(0,time.aux[idf])
    a_inf[length(a_inf)] <- Inf
    a_s_inf  <- c(0,time.aux[idf])
  }
  saida <- list(a_inf,a_s_inf)
  return(saida)
}

###------------------------------------------------------------------------------------------------###
# Function that calculates the differences (t-a[j]) for each individual
###------------------------------------------------------------------------------------------------###
deltas <- function(a_inf,t,IDE,b,n){
  tij <- matrix(0,n,b)
  deltaij <- matrix(0,n,b)
  a <- a_inf
  for(i in 1:n){
    for(j in 1:b){
      deltaij[i,j] <- (min(t[i], a[j+1])-a[j])*((t[i]-a[j])>0)
    }
  }
  return(deltaij)
}

###------------------------------------------------------------------------------------------------###
# function that calculates the cumulative failure rate function for the MEP.
###------------------------------------------------------------------------------------------------###
H0ti <- function(a_inf,t,deltaij,lambda,b,n){
  Hoti <- matrix(0,n,b)
  H0ti <- NULL
  for(i in 1:n){
    for(j in 1:b){
      Hoti[i,j] <- (deltaij[i,j]*lambda[j] )
    }
  }
  for(i in 1:n){
    H0ti[i]<- sum(Hoti[i,])
  }
  return(H0ti)
}

###------------------------------------------------------------------------------------------------###
# Functions for calculating the equation U.
###------------------------------------------------------------------------------------------------###
cumsumrev = function(x){rev(cumsum(rev(x)))}

# Estimation of coefficients for failure times
modelo_T_MEP <-  function(beta_T, X_T, delta.t,risco_a_T,wi,n_intMC){
  w_kl_beta_T <- risco_a_T*exp((X_T[,]%*%beta_T))*(rowSums(exp(wi[,]))/n_intMC)      
  w_kl_beta_T_num <- cbind(w_kl_beta_T, w_kl_beta_T,w_kl_beta_T,w_kl_beta_T)*X_T  
  U_T_1 <- colSums((X_T*delta.t - w_kl_beta_T_num))
  c(U_T_1 = U_T_1)
}

# Estimation of coefficients for withdrawal times
modelo_C_MEP <-  function(beta_C, X_C,delta.c,risco_a_C,wi,n_intMC){
  w_kl_beta_C <- risco_a_C*exp((X_C[,]%*%beta_C[1:4]))*(rowSums(exp(beta_C[5]*wi[,]))/n_intMC)   
  w_kl_beta_C_num <- cbind(w_kl_beta_C, w_kl_beta_C,w_kl_beta_C,w_kl_beta_C)*X_C  
  U_C_1 <- colSums((X_C*delta.c - w_kl_beta_C_num))
  
  w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]%*%beta_C[1:4])*(rowSums(wi[,]*exp(beta_C[5]*wi[,]))/n_intMC)
  U_C_alpha <- sum(delta.c*(rowSums(wi[,])/n_intMC) - w_kl_beta_C_alpha_num)
  c(U_C_1 = U_C_1,U_C_alpha=U_C_alpha)
}

###------------------------------------------------------------------------------------------------###
# Calculation of first order derivatives
###------------------------------------------------------------------------------------------------###
# assuming dependence
Esp.Deriv.Prim.Ordem <-  function( X_T, X_C,X_R,delta_T,delta_C,delta_R, beta_T, beta_C,beta_R, alpha2,alpha3, Sigma2,nu_ind_j_T,nu_ind_j_C,nu_ind_j_R,
                                   lambda_T_j,lambda_C_j,lambda_R_j, risco_a_T, risco_a_C,risco_a_R,deltaij_T,deltaij_C,deltaij_R,w_k_grupo, ident){
n <- nrow(X_T)
wk = w_k_grupo[ident,] 
num_param <- length(beta_T)+length(beta_C)+length(beta_R)+length(alpha2)+length(alpha3)+length(Sigma2)
deriv1 <- matrix(NA,(num_param+length(lambda_T_j)+length(lambda_C_j)+length(lambda_R_j)),ncol(wk))  
  
pred_T <- as.vector(exp(X_T%*%beta_T))*exp(wk)
pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha2*wk)
pred_R <- as.vector(exp(X_R%*%beta_R))*exp(alpha3*wk)
Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
Delta_r <- delta_R%*%t(rep(1,ncol(wk)))
lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
lambda_c <- t(lambda_C_j%*%t(rep(1,n)))
lambda_r <- t(lambda_R_j%*%t(rep(1,n)))

deriv1[1,] <- colSums(X_T[,1]*(Delta_t - risco_a_T*pred_T)) # beta_T_1
deriv1[2,] <- colSums(X_T[,2]*(Delta_t - risco_a_T*pred_T)) # beta_T_2
deriv1[3,] <- colSums(X_T[,3]*(Delta_t - risco_a_T*pred_T)) # beta_T_3
deriv1[4,] <- colSums(X_T[,4]*(Delta_t - risco_a_T*pred_T)) # beta_T_4
  
for ( j in 1:length(lambda_T_j)){ 
  deriv1[(j+4),] <- nu_ind_j_T[j]*(lambda_T_j[j]^(-1)) - colSums(deltaij_T[,j]*pred_T)   # vector of lambda_T_j
}
  
deriv1[(length(lambda_T_j)+5),] <- colSums(X_C[,1]*(Delta_c - risco_a_C*pred_C)) # beta_C_1
deriv1[(length(lambda_T_j)+6),] <- colSums(X_C[,2]*(Delta_c - risco_a_C*pred_C)) # beta_C_2
deriv1[(length(lambda_T_j)+7),] <- colSums(X_C[,3]*(Delta_c - risco_a_C*pred_C)) # beta_C_3
deriv1[(length(lambda_T_j)+8),] <- colSums(X_C[,4]*(Delta_c - risco_a_C*pred_C)) # beta_C_4
deriv1[(length(lambda_T_j)+9),] <- colSums(wk*(Delta_c - risco_a_C*pred_C))      # phi_C
  
for ( j in 1:length(lambda_C_j)){     
  deriv1[(j+(length(lambda_T_j)+9)),] <- nu_ind_j_C[j]*(lambda_C_j[j]^(-1)) - colSums(deltaij_C[,j]*pred_C) # vector of lambda_C_j
}
  
deriv1[(length(lambda_T_j)+length(lambda_C_j)+10),] <- colSums(X_R[,1]*(Delta_r - risco_a_R*pred_R)) # beta_R_1
deriv1[(length(lambda_T_j)+length(lambda_C_j)+11),] <- colSums(X_R[,2]*(Delta_r - risco_a_R*pred_R)) # beta_R_2
deriv1[(length(lambda_T_j)+length(lambda_C_j)+12),] <- colSums(X_R[,3]*(Delta_r - risco_a_R*pred_R)) # beta_R_3
deriv1[(length(lambda_T_j)+length(lambda_C_j)+13),] <- colSums(X_R[,4]*(Delta_r - risco_a_R*pred_R)) # beta_R_4
deriv1[(length(lambda_T_j)+length(lambda_C_j)+14),] <- colSums(wk*(Delta_r - risco_a_R*pred_R))      # phi_R

for ( j in 1:length(lambda_R_j)){     
  deriv1[(j+(length(lambda_T_j)+length(lambda_C_j)+14)),] <- nu_ind_j_R[j]*(lambda_R_j[j]^(-1)) - colSums(deltaij_R[,j]*pred_R) # vector of lambda_R_j
}
  
deriv1[length(lambda_C_j)+length(lambda_T_j)+length(lambda_R_j)+15,] <- colSums(-0.5*Sigma2^(-1) + 0.5*(Sigma2^(-2))*w_k_grupo^(2)) # sigma2
  
aux <- deriv1[,1]%*%t(deriv1[,1])
for( i in 2:ncol(wk)){
  aux <- aux + deriv1[,i]%*%t(deriv1[,i])
}
return(aux/ncol(wk))
}  

# assuming independence
Esp.Deriv.Prim.Ordem.indep <-  function(X_T, delta_T,beta_T, Sigma2,nu_ind_j_T,lambda_T_j,risco_a_T,deltaij_T,w_k_grupo, ident){
  n <- nrow(X_T)
  wk = w_k_grupo[ident,] 
  num_param <- length(beta_T)+length(Sigma2)
  deriv1 <- matrix(NA,(num_param+length(lambda_T_j)),ncol(wk))  
  
  pred_T <- as.vector(exp(X_T%*%beta_T))*exp(wk)
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  
  deriv1[1,] <- colSums(X_T[,1]*(Delta_t - risco_a_T*pred_T)) # beta_T_1
  deriv1[2,] <- colSums(X_T[,2]*(Delta_t - risco_a_T*pred_T)) # beta_T_2
  deriv1[3,] <- colSums(X_T[,3]*(Delta_t - risco_a_T*pred_T)) # beta_T_3
  deriv1[4,] <- colSums(X_T[,4]*(Delta_t - risco_a_T*pred_T)) # beta_T_4
  
  for ( j in 1:length(lambda_T_j)){ 
    deriv1[(j+4),] <- nu_ind_j_T[j]*(lambda_T_j[j]^(-1)) - colSums(deltaij_T[,j]*pred_T)   # vector of lambda_T_j
  }
  
  deriv1[length(lambda_T_j)+5,] <- colSums(-0.5*Sigma2^(-1) + 0.5*(Sigma2^(-2))*w_k_grupo^(2)) # sigma2
  
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
Esp.DerivParciais <-  function( X_T, X_C,X_R,delta_T,delta_C,delta_R, beta_T, beta_C,beta_R, alpha2,alpha3, Sigma2, nu_ind_j_T,nu_ind_j_C,nu_ind_j_R, 
                                lambda_T_j,lambda_C_j,lambda_R_j, risco_a_T, risco_a_C,risco_a_R,deltaij_T,deltaij_C,deltaij_R, w_k_grupo, ident){
n <- nrow(X_T)
num_param <- length(beta_T)+length(beta_C)+length(beta_R)+length(alpha2)+length(alpha3)+length(Sigma2)
deriv2 <- matrix(0,(num_param+length(lambda_T_j)+length(lambda_C_j)+length(lambda_R_j)),(num_param+length(lambda_T_j)+length(lambda_C_j)+length(lambda_R_j)))  #vetor que com as derivadas de primeira ordem, 
wk = w_k_grupo[ident,]   
  
pred_T <- as.vector(exp(X_T%*%beta_T))*exp(wk)
pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha2*wk)
pred_R <- as.vector(exp(X_R%*%beta_R))*exp(alpha3*wk)
Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
Delta_r <- delta_R%*%t(rep(1,ncol(wk)))
lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
lambda_c <- t(lambda_C_j%*%t(rep(1,n)))
lambda_r <- t(lambda_R_j%*%t(rep(1,n)))
  
deriv2[1,1] <- - mean(colSums(X_T[,1]*X_T[,1]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_1
deriv2[1,2] <- - mean(colSums(X_T[,1]*X_T[,2]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_2
deriv2[2,1] <- deriv2[1,2]                       
deriv2[1,3] <- - mean(colSums(X_T[,1]*X_T[,3]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_3
deriv2[3,1] <- deriv2[1,3]   
deriv2[1,4] <- - mean(colSums(X_T[,1]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_4
deriv2[4,1] <- deriv2[1,4]
deriv2[2,2] <- - mean(colSums(X_T[,2]*X_T[,2]*risco_a_T*pred_T)) # derived from beta_T_2's derived in relation to beta_T_2
deriv2[2,3] <- - mean(colSums(X_T[,2]*X_T[,3]*risco_a_T*pred_T)) # derived from beta_T_2's derived in relation to beta_T_3
deriv2[3,2] <- deriv2[2,3]
deriv2[2,4] <- - mean(colSums(X_T[,2]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_2's derived in relation to beta_T_4
deriv2[4,2] <- deriv2[2,4]
deriv2[3,3] <- - mean(colSums(X_T[,3]*X_T[,3]*risco_a_T*pred_T)) # derived from beta_T_3's derived in relation to beta_T_3
deriv2[3,4] <- - mean(colSums(X_T[,3]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_3's derived in relation to beta_T_4
deriv2[4,3] <- deriv2[3,4]
deriv2[4,4] <- - mean(colSums(X_T[,4]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_4's derived in relation to beta_T_4

for (j in 1:length(lambda_T_j)){    
  deriv2[(j+4),(j+4)] <- - nu_ind_j_T[j]*(lambda_T_j[j]^(-2))    # derived from lambda_T_j's derived in relation to lambda_T_j
}
for (j in 1:length(lambda_T_j)){
  deriv2[1,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,1])) # derived from lambda_T_j's derived in relation to beta_T_1
  deriv2[(j+4),1] <- deriv2[1,(j+4)]
}
for (j in 1:length(lambda_T_j)){
  deriv2[2,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,2])) # derived from lambda_T_j's derived in relation to beta_T_2
  deriv2[(j+4),2] <- deriv2[2,(j+4)]
}
for (j in 1:length(lambda_T_j)){
  deriv2[3,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,3])) # derived from lambda_T_j's derived in relation to beta_T_3
  deriv2[(j+4),3] <- deriv2[3,(j+4)]
}
for (j in 1:length(lambda_T_j)){
  deriv2[4,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,4])) # derived from lambda_T_j's derived in relation to beta_T_4
  deriv2[(j+4),4] <- deriv2[4,(j+4)]
}

deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+5)] <- - mean(colSums(X_C[,1]*X_C[,1]*risco_a_C*pred_C)) # derived from beta_C_1 derived in relation to beta_C_1
deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+6)] <- - mean(colSums(X_C[,1]*X_C[,2]*risco_a_C*pred_C)) # derived from beta_C_1 derived in relation to beta_C_2
deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+5)] <- deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+6)]
deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+7)] <- - mean(colSums(X_C[,1]*X_C[,3]*risco_a_C*pred_C)) # derived from beta_C_1 derived in relation to beta_C_3
deriv2[(length(lambda_T_j)+7),(length(lambda_T_j)+5)] <- deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+7)]
deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+8)] <- - mean(colSums(X_C[,1]*X_C[,4]*risco_a_C*pred_C)) # derived from beta_C_1 derived in relation to beta_C_4
deriv2[(length(lambda_T_j)+8),(length(lambda_T_j)+5)] <- deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+8)]
deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+9)] <- - mean(colSums(X_C[,1]*wk*risco_a_C*pred_C)) # derived from beta_C_1 derived in relation to phi.C
deriv2[(length(lambda_T_j)+9),(length(lambda_T_j)+5)] <- deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+9)]
deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+6)] <- - mean(colSums(X_C[,2]*X_C[,2]*risco_a_C*pred_C)) # derived from beta_C_2 derived in relation to beta_C_2
deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+7)] <- - mean(colSums(X_C[,2]*X_C[,3]*risco_a_C*pred_C)) # derived from beta_C_2 derived in relation to beta_C_3
deriv2[(length(lambda_T_j)+7),(length(lambda_T_j)+6)] <- deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+7)]
deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+8)] <- - mean(colSums(X_C[,2]*X_C[,4]*risco_a_C*pred_C)) # derived from beta_C_2 derived in relation to beta_C_4
deriv2[(length(lambda_T_j)+8),(length(lambda_T_j)+6)] <- deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+8)]
deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+9)] <- - mean(colSums(X_C[,2]*wk*risco_a_C*pred_C)) # derived from beta_C_2 derived in relation to phi.C
deriv2[(length(lambda_T_j)+9),(length(lambda_T_j)+6)] <- deriv2[(length(lambda_T_j)+6),(length(lambda_T_j)+9)]
deriv2[(length(lambda_T_j)+7),(length(lambda_T_j)+7)] <- - mean(colSums(X_C[,3]*X_C[,3]*risco_a_C*pred_C)) # derived from beta_C_3 derived in relation to beta_C_3
deriv2[(length(lambda_T_j)+7),(length(lambda_T_j)+8)] <- - mean(colSums(X_C[,3]*X_C[,4]*risco_a_C*pred_C)) # derived from beta_C_3 derived in relation to beta_C_4
deriv2[(length(lambda_T_j)+8),(length(lambda_T_j)+7)] <- deriv2[(length(lambda_T_j)+7),(length(lambda_T_j)+8)]
deriv2[(length(lambda_T_j)+7),(length(lambda_T_j)+9)] <- - mean(colSums(X_C[,3]*wk*risco_a_C*pred_C)) # derived from beta_C_3 derived in relation to phi.C
deriv2[(length(lambda_T_j)+9),(length(lambda_T_j)+7)] <- deriv2[(length(lambda_T_j)+7),(length(lambda_T_j)+9)]
deriv2[(length(lambda_T_j)+8),(length(lambda_T_j)+8)] <- - mean(colSums(X_C[,4]*X_C[,4]*risco_a_C*pred_C)) # derived from beta_C_4 derived in relation to beta_C_4
deriv2[(length(lambda_T_j)+8),(length(lambda_T_j)+9)] <- - mean(colSums(X_C[,4]*wk*risco_a_C*pred_C)) # derived from beta_C_4 derived in relation to phi.C
deriv2[(length(lambda_T_j)+9),(length(lambda_T_j)+8)] <- deriv2[(length(lambda_T_j)+8),(length(lambda_T_j)+9)]
deriv2[(length(lambda_T_j)+9),(length(lambda_T_j)+9)] <- - mean(colSums(wk*wk*risco_a_C*pred_C)) # derived from phi.C derived in relation to phi.C

for (j in 1:length(lambda_C_j)){
  deriv2[(j+length(lambda_T_j)+9),(j+length(lambda_T_j)+9)] <- - nu_ind_j_C[j]*(lambda_C_j[j]^(-2)) # derived from lambda_C_j's derived in relation to lambda_C_j
}
for (j in 1:length(lambda_C_j)){
  deriv2[(length(lambda_T_j)+5),(j+length(lambda_T_j)+9)] <- -mean(colSums(deltaij_C[,j]*pred_C*X_C[,1])) # derived from lambda_C_j's derived in relation to beta_C_1
  deriv2[(j+length(lambda_T_j)+9),(length(lambda_T_j)+5)] <- deriv2[(length(lambda_T_j)+5),(j+length(lambda_T_j)+9)]
}
for (j in 1:length(lambda_C_j)){
  deriv2[(length(lambda_T_j)+6),(j+length(lambda_T_j)+9)] <- -mean(colSums(deltaij_C[,j]*pred_C*X_C[,2])) # derived from lambda_C_j's derived in relation to beta_C_2
  deriv2[(j+length(lambda_T_j)+9),(length(lambda_T_j)+6)] <- deriv2[(length(lambda_T_j)+6),(j+length(lambda_T_j)+9)]
}
for (j in 1:length(lambda_C_j)){
  deriv2[(length(lambda_T_j)+7),(j+length(lambda_T_j)+9)] <- -mean(colSums(deltaij_C[,j]*pred_C*X_C[,3])) # derived from lambda_C_j's derived in relation to beta_C_3
  deriv2[(j+length(lambda_T_j)+9),(length(lambda_T_j)+7)] <- deriv2[(length(lambda_T_j)+7),(j+length(lambda_T_j)+9)]
}
for (j in 1:length(lambda_C_j)){
  deriv2[(length(lambda_T_j)+8),(j+length(lambda_T_j)+9)] <- -mean(colSums(deltaij_C[,j]*pred_C*X_C[,4])) # derived from lambda_C_j's derived in relation to beta_C_4
  deriv2[(j+length(lambda_T_j)+9),(length(lambda_T_j)+8)] <- deriv2[(length(lambda_T_j)+8),(j+length(lambda_T_j)+9)]
}
for (j in 1:length(lambda_C_j)){
  deriv2[(length(lambda_T_j)+9),(j+length(lambda_T_j)+9)] <- -mean(colSums(deltaij_C[,j]*pred_C*wk)) # derived from lambda_C_j's derived in relation to phi.C
  deriv2[(j+length(lambda_T_j)+9),(length(lambda_T_j)+9)] <- deriv2[(length(lambda_T_j)+9),(j+length(lambda_T_j)+9)]
}
 
deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+10)] <- - mean(colSums(X_R[,1]*X_R[,1]*risco_a_R*pred_R)) # derived from beta_R_1 derived in relation to beta_R_1
deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+11)] <- - mean(colSums(X_R[,1]*X_R[,2]*risco_a_R*pred_R)) # derived from beta_R_1 derived in relation to beta_R_2
deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+10)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+11)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+12)] <- - mean(colSums(X_R[,1]*X_R[,3]*risco_a_R*pred_R)) # derived from beta_R_1 derived in relation to beta_R_3
deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(length(lambda_T_j)+length(lambda_C_j)+10)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+12)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+13)] <- - mean(colSums(X_R[,1]*X_R[,4]*risco_a_R*pred_R)) # derived from beta_R_1 derived in relation to beta_R_4
deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(length(lambda_T_j)+length(lambda_C_j)+10)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+13)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+14)] <- - mean(colSums(X_R[,1]*wk*risco_a_R*pred_R)) # derived from beta_R_1 derived in relation to phi.R
deriv2[(length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+10)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(length(lambda_T_j)+length(lambda_C_j)+14)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+11)] <- - mean(colSums(X_R[,2]*X_R[,2]*risco_a_R*pred_R)) # derived from beta_R_2 derived in relation to beta_R_2
deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+12)] <- - mean(colSums(X_R[,2]*X_R[,3]*risco_a_R*pred_R)) # derived from beta_R_2 derived in relation to beta_R_3
deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(length(lambda_T_j)+length(lambda_C_j)+11)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+12)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+13)] <- - mean(colSums(X_R[,2]*X_R[,4]*risco_a_R*pred_R)) # derived from beta_R_2 derived in relation to beta_R_4
deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(length(lambda_T_j)+length(lambda_C_j)+11)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+13)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+14)] <- - mean(colSums(X_R[,2]*wk*risco_a_R*pred_R)) # derived from beta_R_2 derived in relation to phi.R
deriv2[(length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+11)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(length(lambda_T_j)+length(lambda_C_j)+14)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(length(lambda_T_j)+length(lambda_C_j)+12)] <- - mean(colSums(X_R[,3]*X_R[,3]*risco_a_R*pred_R)) # derived from beta_R_3 derived in relation to beta_R_3
deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(length(lambda_T_j)+length(lambda_C_j)+13)] <- - mean(colSums(X_R[,3]*X_R[,4]*risco_a_R*pred_R)) # derived from beta_R_3 derived in relation to beta_R_4
deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(length(lambda_T_j)+length(lambda_C_j)+12)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(length(lambda_T_j)+length(lambda_C_j)+13)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(length(lambda_T_j)+length(lambda_C_j)+14)] <- - mean(colSums(X_R[,3]*wk*risco_a_R*pred_R)) # derived from beta_R_3 derived in relation to phi.R
deriv2[(length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+12)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(length(lambda_T_j)+length(lambda_C_j)+14)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(length(lambda_T_j)+length(lambda_C_j)+13)] <- - mean(colSums(X_R[,4]*X_R[,4]*risco_a_R*pred_R)) # derived from beta_R_4 derived in relation to beta_R_4
deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(length(lambda_T_j)+length(lambda_C_j)+14)] <- - mean(colSums(X_R[,4]*wk*risco_a_R*pred_R)) # derived from beta_R_4 derived in relation to phi.R
deriv2[(length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+13)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(length(lambda_T_j)+length(lambda_C_j)+14)]
deriv2[(length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+14)] <- - mean(colSums(wk*wk*risco_a_R*pred_R)) # derived from phi.R derived in relation to phi.R
  
for (j in 1:length(lambda_R_j)){
  deriv2[(j+length(lambda_T_j)+length(lambda_C_j)+14),(j+length(lambda_T_j)+length(lambda_C_j)+14)] <- - nu_ind_j_R[j]*(lambda_R_j[j]^(-2)) # derived from lambda_R_j's derived in relation to lambda_R_j
}
for (j in 1:length(lambda_R_j)){
  deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(j+length(lambda_T_j)+length(lambda_C_j)+14)] <- -mean(colSums(deltaij_R[,j]*pred_R*X_R[,1])) # derived from lambda_R_j's derived in relation to beta_R_1
  deriv2[(j+length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+10)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+10),(j+length(lambda_T_j)+length(lambda_C_j)+14)]
}
for (j in 1:length(lambda_R_j)){
  deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(j+length(lambda_T_j)+length(lambda_C_j)+14)] <- -mean(colSums(deltaij_R[,j]*pred_R*X_R[,2])) # derived from lambda_R_j's derived in relation to beta_R_2
  deriv2[(j+length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+11)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+11),(j+length(lambda_T_j)+length(lambda_C_j)+14)]
}
for (j in 1:length(lambda_R_j)){
  deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(j+length(lambda_T_j)+length(lambda_C_j)+14)] <- -mean(colSums(deltaij_R[,j]*pred_R*X_R[,3])) # derived from lambda_R_j's derived in relation to beta_R_3
  deriv2[(j+length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+12)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+12),(j+length(lambda_T_j)+length(lambda_C_j)+14)]
}
for (j in 1:length(lambda_R_j)){
  deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(j+length(lambda_T_j)+length(lambda_C_j)+14)] <- -mean(colSums(deltaij_R[,j]*pred_R*X_R[,4])) # derived from lambda_R_j's derived in relation to beta_R_4
  deriv2[(j+length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+13)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+13),(j+length(lambda_T_j)+length(lambda_C_j)+14)]
}
for (j in 1:length(lambda_R_j)){
  deriv2[(length(lambda_T_j)+length(lambda_C_j)+14),(j+length(lambda_T_j)+length(lambda_C_j)+14)] <- -mean(colSums(deltaij_R[,j]*pred_R*wk)) # derived from lambda_R_j's derived in relation to phi.R
  deriv2[(j+length(lambda_T_j)+length(lambda_C_j)+14),(length(lambda_T_j)+length(lambda_C_j)+14)] <- deriv2[(length(lambda_T_j)+length(lambda_C_j)+14),(j+length(lambda_T_j)+length(lambda_C_j)+14)]
}

deriv2[(length(lambda_T_j)+length(lambda_C_j)+length(lambda_R_j)+15),(length(lambda_T_j)+length(lambda_C_j)+length(lambda_R_j)+15)] <- mean(colSums( 0.5*Sigma2^(-2) - (Sigma2^(-3))*w_k_grupo^(2))) # derived from sigma2 derived in relation to sigma2
return((as.matrix(deriv2)))
}

# assuming independence
Esp.DerivParciais.indep <-  function(X_T,delta_T,beta_T, Sigma2, nu_ind_j_T,lambda_T_j,risco_a_T,deltaij_T,w_k_grupo, ident){
  n <- nrow(X_T)
  num_param <- length(beta_T)+length(Sigma2)
  deriv2 <- matrix(0,(num_param+length(lambda_T_j)),(num_param+length(lambda_T_j)))  #vetor que com as derivadas de primeira ordem, 
  wk = w_k_grupo[ident,]   
  
  pred_T <- as.vector(exp(X_T%*%beta_T))*exp(wk)
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  
  deriv2[1,1] <- - mean(colSums(X_T[,1]*X_T[,1]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_1
  deriv2[1,2] <- - mean(colSums(X_T[,1]*X_T[,2]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_2
  deriv2[2,1] <- deriv2[1,2]                       
  deriv2[1,3] <- - mean(colSums(X_T[,1]*X_T[,3]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_3
  deriv2[3,1] <- deriv2[1,3]   
  deriv2[1,4] <- - mean(colSums(X_T[,1]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_1's derived in relation to beta_T_4
  deriv2[4,1] <- deriv2[1,4]
  deriv2[2,2] <- - mean(colSums(X_T[,2]*X_T[,2]*risco_a_T*pred_T)) # derived from beta_T_2's derived in relation to beta_T_2
  deriv2[2,3] <- - mean(colSums(X_T[,2]*X_T[,3]*risco_a_T*pred_T)) # derived from beta_T_2's derived in relation to beta_T_3
  deriv2[3,2] <- deriv2[2,3]
  deriv2[2,4] <- - mean(colSums(X_T[,2]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_2's derived in relation to beta_T_4
  deriv2[4,2] <- deriv2[2,4]
  deriv2[3,3] <- - mean(colSums(X_T[,3]*X_T[,3]*risco_a_T*pred_T)) # derived from beta_T_3's derived in relation to beta_T_3
  deriv2[3,4] <- - mean(colSums(X_T[,3]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_3's derived in relation to beta_T_4
  deriv2[4,3] <- deriv2[3,4]
  deriv2[4,4] <- - mean(colSums(X_T[,4]*X_T[,4]*risco_a_T*pred_T)) # derived from beta_T_4's derived in relation to beta_T_4
  
  for (j in 1:length(lambda_T_j)){    
    deriv2[(j+4),(j+4)] <- - nu_ind_j_T[j]*(lambda_T_j[j]^(-2))    # derived from lambda_T_j's derived in relation to lambda_T_j
  }
  for (j in 1:length(lambda_T_j)){
    deriv2[1,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,1])) # derived from lambda_T_j's derived in relation to beta_T_1
    deriv2[(j+4),1] <- deriv2[1,(j+4)]
  }
  for (j in 1:length(lambda_T_j)){
    deriv2[2,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,2])) # derived from lambda_T_j's derived in relation to beta_T_2
    deriv2[(j+4),2] <- deriv2[2,(j+4)]
  }
  for (j in 1:length(lambda_T_j)){
    deriv2[3,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,3])) # derived from lambda_T_j's derived in relation to beta_T_3
    deriv2[(j+4),3] <- deriv2[3,(j+4)]
  }
  for (j in 1:length(lambda_T_j)){
    deriv2[4,(j+4)] <- -mean(colSums(deltaij_T[,j]*pred_T*X_T[,4])) # derived from lambda_T_j's derived in relation to beta_T_4
    deriv2[(j+4),4] <- deriv2[4,(j+4)]
  }
  
  deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+5)] <- mean(colSums( 0.5*Sigma2^(-2) - (Sigma2^(-3))*w_k_grupo^(2))) # derived from sigma2 derived in relation to sigma2
  return((as.matrix(deriv2)))
}

###------------------------------------------------------------------------------------------------###
# PEM Proposed model
###------------------------------------------------------------------------------------------------###

PEM.dep  <- function(dados){
  
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
risco_a_T      <- basehaz(ajuste_coxph_T, centered=FALSE)$hazard  #cumulative hazard
beta_T         <- ajuste_coxph_T$coef
beta_T0        <- beta_T

ajuste_coxph_C <- coxph(Surv(time, delta_c) ~  X[,1] + X[,2] + X[,3] + X[,4], method="breslow")
risco_a_C      <- basehaz(ajuste_coxph_C, centered=FALSE)$hazard  #cumulative hazard
beta_C         <- ajuste_coxph_C$coef
beta_C0        <- beta_C

ajuste_coxph_R <- coxph(Surv(time, delta_r) ~  X[,1] + X[,2] + X[,3] + X[,4], method="breslow")
risco_a_R      <- basehaz(ajuste_coxph_R, centered=FALSE)$hazard  #cumulative hazard
beta_R         <- ajuste_coxph_R$coef
beta_R0        <- beta_R

phi.C <- 0
phi.R <- 0
Sigma2 <- 1
param <- c(beta_T,beta_C,phi.C,beta_R, phi.R, Sigma2) 

###--------------------------------------------------------------------------------------------###
# EMMC algorithm specifications
###--------------------------------------------------------------------------------------------###
maxit <- 100                             # maximum number of iterations
eps1= rep(1e-3, length(param))
eps2= rep(1e-4, length(param))
n_intMCs = c(rep(10,20),rep(25,30), rep(50,40), rep(75,10))
out =  matrix(NA,maxit+1,length(param))  # start the objects used
dif =  matrix(NA,maxit+1,length(param))   
final = length(param)
count = rep(0,maxit+1)
s=1
continue=TRUE
X_T <- X; N_T <- length(X_T[,1])
X_C <- X; N_C <- length(X_C[,1])
X_R <- X; N_R <- length(X_R[,1])
bmax <- NULL  

###--------------------------------------------------------------------------------------------###
# Initiation of the EMMC algorithm 
###--------------------------------------------------------------------------------------------###

while (continue == TRUE) {
  count = rep(0,maxit+1)
  out[s,] = c(beta_T,beta_C,phi.C,beta_R, phi.R, Sigma2)
  n_intMC = n_intMCs[s]
  w_chapeu_grupo = matrix(NA, m, n_intMC)
  
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
  # estimate for failure rate for failure time	
  a_T <- time.grid(time=time, event=delta_t, n.int=bmax)  
  num_int_T  <- (length(a_T[[1]])-1)
  id_T <- as.numeric(cut(time,a_T[[1]]))  
  U_T <- rep(0,length(unique(id_T))-1)
  nu_T <- tapply(delta_t,id_T ,sum)
  
  pred_linear_T <- exp(X_T%*%beta_T)*rowMeans(exp(wi[,]))
  A_T <- a_T[[2]]  
  xi <- NULL
  for(j in 1:num_int_T){
    y <- time
    y[id_T<j] <- A_T[j]
    y[id_T>j] <- A_T[j+1]
    xi[j] <- sum((y-A_T[j])*pred_linear_T)
  }
  sumX_T <- xi
  lambda_T_j <- nu_T/sumX_T
  ###------------------------------------------------------------------------------------
  # estimate for failure rate for withdrawal time
  a_C <- time.grid(time=time, event=delta_c, n.int=bmax)
  num_int_C  <- (length(a_C[[1]])-1)
  id_C <- as.numeric(cut(time,a_C[[1]]))  
  U_C <- rep(0,length(unique(id_C))-1)
  nu_C <- tapply(delta_c,id_C ,sum)
  
  pred_linear_C <- exp(X_C%*%beta_C)*rowMeans(exp(phi.C*wi[,]))
  A_C <- a_C[[2]]
  xi <- NULL
  for(j in 1:num_int_C){
    y <- time
    y[id_C<j] <- A_C[j]
    y[id_C>j] <- A_C[j+1]
    xi[j] <- sum((y-A_C[j])*pred_linear_C)
  }
  sumX_C <- xi
  lambda_C_j <- nu_C/sumX_C
  ###------------------------------------------------------------------------------------
  # estimate for failure rate for transplant time
  a_R <- time.grid(time=time, event=delta_r, n.int=bmax)
  num_int_R  <- (length(a_R[[1]])-1)
  id_R <- as.numeric(cut(time,a_R[[1]]))
  U_R <- rep(0,length(unique(id_R))-1)
  nu_R <- tapply(delta_r,id_R ,sum)
  
  pred_linear_R <- exp(X_R%*%beta_R)*rowMeans(exp(phi.R*wi[,]))
  A_R <- a_R[[2]]
  xi <- NULL
  for(j in 1:num_int_R){
    y <- time
    y[id_R<j] <- A_R[j]
    y[id_R>j] <- A_R[j+1]
    xi[j] <- sum((y-A_R[j])*pred_linear_R)
  }
  sumX_R <- xi
  lambda_R_j <- nu_R/sumX_R
  ###------------------------------------------------------------------------------------
  # calculation of the cumulative failure rate function   
  deltaij_T <- deltas(A_T,time,id_T,num_int_T,n)
  risco_a_T <- H0ti(a_T[[1]],time,deltaij_T,lambda_T_j,num_int_T,n)
  deltaij_C <- deltas(A_C,time,id_C,num_int_C,n)
  risco_a_C <- H0ti(a_C[[1]],time,deltaij_C,lambda_C_j,num_int_C,n)
  deltaij_R <- deltas(A_R,time,id_R,num_int_R,n)
  risco_a_R <- H0ti(a_R[[1]],time,deltaij_R,lambda_R_j,num_int_R,n)
  ###------------------------------------------------------------------------------------
  # estimation of the coefficients, by the U equations 
  S_T <- multiroot(f = modelo_T_MEP, start = beta_T0,  X_T=X_T, delta.t=delta_t,risco_a_T=risco_a_T,wi=wi,n_intMC=n_intMC)
  beta_T <- S_T$root
  S_C <- multiroot(f = modelo_C_MEP, start = c(beta_C0, 0), X_C=X_C,delta.c=delta_c,risco_a_C=risco_a_C,wi=wi,n_intMC=n_intMC)
  beta_C <- S_C$root[1:4]
  phi.C <- S_C$root[5]
  #The estimation for the transplant parameters has the same calculations as the withdrawal estimation. So, can use the same function
  #of withdrawal, but with the transplant variables.  
  S_R <- multiroot(f = modelo_C_MEP, start = c(beta_R0, 0), X_C=X_R,delta.c=delta_r,risco_a_C=risco_a_R,wi=wi,n_intMC=n_intMC)
  beta_R <- S_R$root[1:4]
  phi.R <- S_R$root[5]
  ###------------------------------------------------------------------------------------
  # stop criterion specification
  out[s+1,]= c(beta_T,beta_C,phi.C,beta_R, phi.R, Sigma2)
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
beta_C <- apply(out[(s-2):s,5:8], 2, mean); 
phi.C  <- mean(out[(s-2):s,9]); 
beta_R <- apply(out[(s-2):s,10:13], 2, mean); 
phi.R  <- mean(out[(s-2):s,14]); 
Sigma2 <- mean(out[(s-2):s,15]); 
Estimate <- c(beta_T,Sigma2,beta_C,phi.C,beta_R,phi.R)
###--------------------------------------------------------------------------------------------###
# Standard errors calculation  
###--------------------------------------------------------------------------------------------###

Esp_deriv_ordem1 <- Esp.Deriv.Prim.Ordem( X_T=X_T,X_C=X_C,X_R=X_R,delta_T=delta_t,delta_C=delta_c,delta_R=delta_r, beta_T=beta_T,beta_C=beta_C,beta_R=beta_R,alpha2=phi.C,alpha3=phi.R,Sigma2=Sigma2,
                                          nu_ind_j_T=nu_T,nu_ind_j_C=nu_C,nu_ind_j_R=nu_R, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j,lambda_R_j=lambda_R_j, risco_a_T=risco_a_T,risco_a_C=risco_a_C,
                                          risco_a_R=risco_a_R,deltaij_T=deltaij_T,deltaij_C=deltaij_C,deltaij_R=deltaij_R,w_k_grupo=w_chapeu_grupo,ident=ident)

Esp_deriv_ordem2 <- Esp.DerivParciais(X_T=X_T,X_C=X_C,X_R=X_R,delta_T=delta_t,delta_C=delta_c,delta_R=delta_r, beta_T=beta_T,beta_C=beta_C,beta_R=beta_R,alpha2=phi.C,alpha3=phi.R,Sigma2=Sigma2,
                                      nu_ind_j_T=nu_T,nu_ind_j_C=nu_C,nu_ind_j_R=nu_R, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j,lambda_R_j=lambda_R_j, risco_a_T=risco_a_T,risco_a_C=risco_a_C,
                                      risco_a_R=risco_a_R,deltaij_T=deltaij_T,deltaij_C=deltaij_C,deltaij_R=deltaij_R,w_k_grupo=w_chapeu_grupo,ident=ident)

InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
Var <- solve(InfFisher)
EP_louis <- sqrt(diag(Var)) 
#---------------------------------------------------------------------------------------------------------------
# Standard errors for the parameters of interest
p <- ncol(X)
SE.beta_T <- EP_louis[1:p];    
SE.beta_C <- EP_louis[(length(lambda_T_j)+p+1):(length(lambda_T_j)+p+p)]; 
SE.phi.C  <- EP_louis[length(lambda_T_j)+p+p+1]; 
SE.beta_R <- EP_louis[(length(lambda_T_j)+p+p+1+length(lambda_C_j)+1):(length(lambda_T_j)+p+p+1+length(lambda_C_j)+p)];
SE.phi.R  <- EP_louis[length(lambda_T_j)+p+p+1+length(lambda_C_j)+p+1]; 
SE.Sigma2 <- EP_louis[length(lambda_T_j)+p+p+1+length(lambda_C_j)+p+1+length(lambda_R_j)+1]; 

SE <- c(SE.beta_T,SE.Sigma2,SE.beta_C,SE.phi.C,SE.beta_R,SE.phi.R)

return(cbind(round(Estimate,3),round(SE,3)))

}

###------------------------------------------------------------------------------------------------###
# PEM Standard model
###------------------------------------------------------------------------------------------------###

PEM.indep  <- function(dados){
  
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
  risco_a_T      <- basehaz(ajuste_coxph_T, centered=FALSE)$hazard  #cumulative hazard
  beta_T         <- ajuste_coxph_T$coef
  beta_T0        <- beta_T
  Sigma2 <- 1
  param <- c(beta_T,Sigma2) 
  ###--------------------------------------------------------------------------------------------###
  # EMMC algorithm specifications
  ###--------------------------------------------------------------------------------------------###
  maxit <- 100                             # maximum number of iterations
  eps1= rep(1e-3, length(param))
  eps2= rep(1e-4, length(param))
  n_intMCs = c(rep(10,20),rep(25,30), rep(50,40), rep(75,10))
  out =  matrix(NA,maxit+1,length(param))  # start the objects used
  dif =  matrix(NA,maxit+1,length(param))   
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE
  X_T <- X; N_T <- length(X_T[,1])
  bmax <- NULL  
  
  ###--------------------------------------------------------------------------------------------###
  # Initiation of the EMMC algorithm 
  ###--------------------------------------------------------------------------------------------###
  
  while (continue == TRUE) {
    count = rep(0,maxit+1)
    out[s,] = c(beta_T,Sigma2)
    n_intMC = n_intMCs[s]
    w_chapeu_grupo = matrix(NA, m, n_intMC)
    
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
    # estimate for failure rate for failure time	
    a_T <- time.grid(time=time, event=delta_t, n.int=bmax)  
    num_int_T  <- (length(a_T[[1]])-1)
    id_T <- as.numeric(cut(time,a_T[[1]]))  
    U_T <- rep(0,length(unique(id_T))-1)
    nu_T <- tapply(delta_t,id_T ,sum)
    
    pred_linear_T <- exp(X_T%*%beta_T)*rowMeans(exp(wi[,]))
    A_T <- a_T[[2]]  
    xi <- NULL
    for(j in 1:num_int_T){
      y <- time
      y[id_T<j] <- A_T[j]
      y[id_T>j] <- A_T[j+1]
      xi[j] <- sum((y-A_T[j])*pred_linear_T)
    }
    sumX_T <- xi
    lambda_T_j <- nu_T/sumX_T
    ###------------------------------------------------------------------------------------
    # calculation of the cumulative failure rate function   
    deltaij_T <- deltas(A_T,time,id_T,num_int_T,n)
    risco_a_T <- H0ti(a_T[[1]],time,deltaij_T,lambda_T_j,num_int_T,n)
    ###------------------------------------------------------------------------------------
    # estimation of the coefficients, by the U equations 
    S_T <- multiroot(f = modelo_T_MEP, start = beta_T0,  X_T=X_T, delta.t=delta_t,risco_a_T=risco_a_T,wi=wi,n_intMC=n_intMC)
    beta_T <- S_T$root
    ###------------------------------------------------------------------------------------
    # stop criterion specification
    out[s+1,]= c(beta_T,Sigma2)
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
  
  Esp_deriv_ordem1 <- Esp.Deriv.Prim.Ordem.indep(X_T=X_T,delta_T=delta_t,beta_T=beta_T,Sigma2=Sigma2,
                                            nu_ind_j_T=nu_T,lambda_T_j=lambda_T_j,risco_a_T=risco_a_T,
                                            deltaij_T=deltaij_T,w_k_grupo=w_chapeu_grupo,ident=ident)
  
  Esp_deriv_ordem2 <- Esp.DerivParciais.indep(X_T=X_T,delta_T=delta_t,beta_T=beta_T,Sigma2=Sigma2,
                                        nu_ind_j_T=nu_T,lambda_T_j=lambda_T_j,risco_a_T=risco_a_T,
                                        deltaij_T=deltaij_T,w_k_grupo=w_chapeu_grupo,ident=ident)
  
  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  EP_louis <- sqrt(diag(Var)) 
  #---------------------------------------------------------------------------------------------------------------
  # Standard errors for the parameters of interest
  p <- ncol(X)
  SE.beta_T <- EP_louis[1:p];    
  SE.Sigma2 <- EP_louis[length(lambda_T_j)+p+1]; 
  
  SE <- c(SE.beta_T,SE.Sigma2)
  
  return(cbind(round(Estimate,3),round(SE,3)))
  
}
