# Functions for Spatial Equilibrium Model
# Annotation and dependencies in the solver functions

# MODEL INVERSION (inversionModel)
# 4 Steps. 
# Step 1 wages_inversion is iterative (algorithmic). The rest are mechanical based on model equations
# Function to invert model, so amenities, wages, productivities, and development density are chosen to match model to data.
# Fundamentals solving for: 
# b_i residential amenities
# a_i productivity
# phi_i (varphi) density of land development 

# The assumptions of how the economy works (e.g. utility function, firm production function, etc.) make much of the model mechanical.

# Step 1: wages_inversion
# Key step for inverting the model
# Computes equilibrium wages that make the model labor in every location in equal to the observed data. 
# It finds the w's such that equation (3.2) holds.
# model Labor supply = Residence * commuting share from i to j.
# Commuting share is a function of wages in each location and commuting costs

# Steps:
# a) Start with initial vector of wages. 
# b) construct commuting shares 
# c) construct total employment 
# d) Create new vector of wages 
# e) update new wage vector based on contraction mapping 
# f) Normalize the vector of wages 
# g) repeat algorithm until difference in wages is less than some tolerance factor.

# Input functions
# commuting_matrix -> creates (d_ij) commuting costs (aka trade costs) from travel times

# Output from wages_inversion
# w_i = w_tr -> wages
# W_i -> Market access measure in each location
# lambda_ij_i -> Probability of individuals in each location of working in each location

# Output function from wages_inversion
# Step 1b - compute average income y_bar
# function: av_income_simple
# Probability of individuals in each location of working in each location * wage in location

# Step 2: density_development
# Computes residential and commercial floorspace supply and equilibrium prices
# MECHANICAL using model equations
# Based on parameters, residence, employment, average income and wages

# Output 
# average floorspace price
# normalized floorspace price
# floorspace for firms
# floorspace for residential (people)
# total floorspace (firms + residents)
# varphi - development density

# Step 3: productivity
#' Computes productivity levels in each location
#' MECHANICAL using model equations
# FOC from firm cost min

# Output
# a = productivity

# Step 4: amenities	
#' Function to estimate amenity parameters of locations where users live.
#' MECHANICAL using model equations

# Output 
# b = amenities

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Re-running counterfactuals
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# inversionModel recovers fundamentals (productivity (a), amenities (b), density of land development(varphi))
# counterfactuals would be policies that impact productivity, amenities, density of land development

# This creates model defined parameters for economic fundamentals
# a = local productivity
# b = local amenities
# varph = local development density
# w_eq = local wages
# u_eq = local utilities
# Q_eq = local floorspace prices (housing prices)
# ttheta_eq = local share floorspace for commercial


# solveModel
# starts with EXOgenous fundamentals (and model defined parameters) recovered by inversion
  # productivities, amenities, commuting costs, available floor space

# Estimate new level of ENDogenous variables
  # ttheta_i        Share commercial floorspace
  # Q_i             Floorspace prices (i.e. housing prices)
  # w_i             wages
  # L_Ri            Resident population
  # L_Mi            Employment
#----------------------------------------------------
# Average income in each location
#----------------------------------------------------
#' Computes average income in each location, which is the weighted average of
#' the income of the people living in the location.
#'
#' @param lambda_ij_i NxN matrix - Probability of individuals in each
#'     location of working in each location.
#' @param w_tr NxS - Wages in each location in each sector.
av_income_simple = function(lambda_ij_i,
                            w_tr){
  
  y_bar = (sumDims2(array_operator(lambda_ij_i, w_tr, '*'), 2));
  return(list(y_bar=y_bar))
}

#------------------------------------------------------
# Travel times -> iceberg commuting costs
#-------------------------------------------------------
#' Function to transform travel times into iceberg commuting costs
# @param t_ij NxN matrix - Travel time matrix across locations
# @param epsilon Float - Parameter that transforms travel times to commuting costs
#' 
# @return A NxN matrix of commuting costs
commuting_matrix = function(t_ij,
                            epsilon){
  tau = exp(epsilon*t_ij)
  return(list(tau=tau))
}

#------------------------------------
# Wage Inversion
#------------------------------------

#' Function to compute equilibrium wages that make the model labor in every
#' location in equal to the observed data. It finds the w's
#' such that equation (3.2) holds. Maybe means C.2
#'
#' @param N Integer - Number of locations.
#' @param w_init Initial vector of wages.
#' @param theta Float - Commuting elasticity.
#' @param tau NxN matrix - Commuting cost matrix across all locations.
#' @param L_i Nx1 matrix - Number of residents in each location.
#' @param L_j Nx1 matrix - Number of workers in each location.
#' @param nu_init Float - Convergence parameter to update wages.
#'     Default nu=0.01.
#' @param tol Float - Maximum tolerable error for estimating total labor.
#'     Default tol=10^-10.
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=10000.
#'
#' @return A list with equilibrium wages and probability of workers in each
#'     location working in every other location.
wages_inversion = function(N,
                           w_init,
                           theta,
                           tau,
                           L_i,
                           L_j,
                           nu_init=0.05,
                           tol=10^-10,
                           maxiter=10000){
  
  # Settings
  outerdiff = Inf
  w = w_init
  iter = 0
  nu = nu_init
  
  cat("Inverting wages...\n")
  while(outerdiff>tol & iter<maxiter){
    # 1) Labor supply
    # Indirect utility
    w_tr = aperm(array(w, dim=c(N,1)), c(2,1));
    rep_w_tr = kronecker(w_tr^theta, array(1, dim=c(N, 1)));
    # Constructing emp` loyment shares
    w_tr_tau = array_operator(w_tr^theta, tau^(-theta), '*');
    lambda_ij_i = array_operator(w_tr_tau, sumDims2(w_tr_tau,2), '/');
    W_i = (sumDims2(w_tr_tau,2))^(1/theta);
    
    # Labor is equal to probabilities * total number of residents * proportion of workers in each sector.
    L_ij = array_operator(L_i, lambda_ij_i, '*')
    L_j_tr = sumDims2(L_ij, 1)
    #    L_j_model = aperm(L_j_tr, c(2, 1));
    Ratio_supply = array_operator(L_j_tr, w, "/");
    w_prime = array_operator(L_j, Ratio_supply, "/");
    
    z_L = array_operator(w, w_prime, '-');
    w = array_operator(w*(1-nu), w_prime*nu, '+');
    w_mean = exp(mean(log(w)))
    w = w/w_mean;
    outerdiff = max(abs(z_L))
    
    iter = iter+1;
    
    if(iter %% 10 == 0){
      cat(paste0("Iteration: ", iter, ", error: ", round(outerdiff, 10), ".\n"))
    }
  }
  if(outerdiff<=tol){
    cat(paste0("Converged after ", iter, " iterations. Error=", round(outerdiff, 10), ".\n"))
  } else{
    cat(paste0("Reached maximum number of iterations (", iter, "). Error=", round(outerdiff, 10), ".\n"))
  }
  
  return(list(w=w, w_tr=w_tr, W_i=W_i, lambda_ij_i=lambda_ij_i))
}


#------------------------------------------------------
# Density-Development
#-------------------------------------------------------

#' Computes residential and commercial floorspace supply and equilibrium prices.
#' MECHANICAL
#'    Based on parameters, residence, employment and wages 
#' @param Q Nx1 array - Floorspaces prices.
#' @param K Nx1 array - Land supply.
#' @param w NxS - Wages in each location in each sector.
#' @param L_j Nx1 matrix - Number of workers in each location.
#' @param y_bar - Average income in each location.
#' @param L_i Nx1 matrix - Number of residents in each location.
#' @param beta Float - Cobb-Douglas parameter output elasticity wrt labor.
#' @param alpha Float - Utility parameter that determines preferences for
#'     consumption.
#' @param mu Float - Floorspace prod function: output elast wrt capita, 1-mu wrt land.     
density_development = function(Q,
                               K,
                               w,
                               L_j,
                               y_bar,
                               L_i,
                               beta,
                               alpha,
                               mu){
  
  Q_mean = exp(mean(log(Q)));
  Q_norm = Q/Q_mean;
  FS_f = ((1-beta)/beta)*(array_operator(array_operator(w, L_j, '*'), Q_norm, '/'));
  FS_r = (1-alpha)*(array_operator(array_operator(y_bar, L_i, '*'),Q_norm,'/'));
  FS = FS_f+FS_r;
  varphi = array_operator(FS, K^(1-mu), '/'); 
  return(list(Q_mean = Q_mean, Q_norm = Q_norm, FS_f = FS_f, FS_r = FS_r, FS = FS, varphi=varphi))
}
# output 
  # average floorspace price
  # normalized floorspace price
  # floorspace firms
  # floorspace residential
  # total floorspace (firms + residents)
  # varphi - development density

#-----------------------------------------
# Productivity
#-----------------------------------------

#' Computes productivity levels in each location
#' MECHANICAL
#' 
#' @param N Float - Number of locations.
#' @param Q Nx1 matrix - Floorspace prices in each location.
#' @param w Nx1 matrix - wages in each location.
#' @param L_j Nx1 matrix - Employment in each location.
#' @param K Nx1 matrix - Land in each location.
#' @param t_ij NxN matrix - Travel times matrix.
#' @param delta Float - decay parameter agglomeration.
#' @param lambda Float - agglomeration force.
#' @param beta Float - Output elasticity wrt labor
productivity = function(N,
                        Q,
                        w,
                        L_j,
                        K,
                        t_ij,
                        delta,
                        lambda,
                        beta){
  
  Q_mean = exp(mean(log(Q)));
  Q_norm = Q/Q_mean;
  beta_tilde = ((1-beta)^(1-beta))*(beta^beta); 
  A = (1/beta_tilde)*(array_operator(Q_norm^(1-beta), w^beta, '*'));
  L_j_dens = (array_operator(L_j, K, '/'));
  L_j_dens_per = aperm(array(L_j_dens, dim=c(N,1)), c(2,1));
  L_j_dens_rep = kronecker(L_j_dens_per, array(1, dim=c(N, 1)));
  Upsilon = sumDims2(array_operator(exp(-delta*t_ij), L_j_dens_rep, '*'), 2);
  a = array_operator(A, Upsilon^lambda, "/");
  a = a*(L_j>0)
  return(list(A = A, a = a))
}


#------------------------------------------------
# Amenity parameters
#------------------------------------------------

#' Function to estimate amenity parameters of locations where users live.
#' MECHANICAL
#'  
#' @param theta Float - Parameter that governs the reallocation of workers across
#'     locations in the city. This parameter measures how sensible are migration
#'     flows within the city to changes in real income.
#' @param N Integer - Number of locations.
#' @param L_i Nx1 matrix - Total residents.
#' @param W_i Nx1 matrix - Market access measure in each location.
#' @param Q Nx1 matrix - Floor space prices.
#' @param K Nx1 matrix - Land area
#' @param alpha Float - Para     
#' @param t_ij NxN matrix - Travel times across locations.
#' @param rho Float - decay parameter for amenities.
#' @param eta Float - congestion force
#'
#' @return Matrix with the amenity distribution of living in each location.
living_amenities_simple = function(theta,
                                   N,
                                   L_i,
                                   W_i,
                                   Q,
                                   K,
                                   alpha,
                                   t_ij,
                                   rho,
                                   eta){
  Q_mean = exp(mean(log(Q)));
  Q_norm = Q/Q_mean;
  L_i_mean = exp(mean(log(L_i)));
  L_i_norm = L_i/L_i_mean;
  W_i_mean = exp(mean(log(W_i)));
  W_i_norm = W_i/W_i_mean;
  B = array_operator(array_operator(L_i_norm^(1/theta), Q_norm^(1-alpha), '*'), W_i_norm^((-1)), '*');
  L_i_dens = (array_operator(L_i, K, '/'));
  L_i_dens_per = aperm(array(L_i_dens, dim=c(N,1)), c(2,1));
  L_i_dens_rep = kronecker(L_i_dens_per, array(1, dim=c(N, 1)));
  Omega = sumDims2(array_operator(exp(-rho*t_ij), L_i_dens_rep, '*'), 2);  
  b = array_operator(B, Omega^(-eta), "/");
  b = b*(L_i>0);
  return(list(B = B, b = b))
}



#----------------------------------------------
# Model inversion
#----------------------------------------------

#' Function to invert model, so amenities, wages, productivities, and development density
#'  are chosen to match model to data.
#'
#' @param N Integer - Number of locations.
#' @param L_i Nx1 matrix - Number of residents in each location.
#' @param L_j Nx1 matrix - Number of workers in each location. 
#' @param Q Nx1 matrix - Floorspace prices
#' @param K Nx1 matrix - Land area
#' @param t_ij NxN matrix - Travel times across all possible locations.
#' @param alpha Float - Utility parameter that determines preferences for
#'     consumption.
#' @param beta Float - Output elasticity wrt labor
#' @param theta Float - Commuting elasticity and migration elasticity.
#' @param delta Float - Decay parameter agglomeration
#' @param rho Float - Decay parameter congestion
#' @param lambda Float - Agglomeration force
#' @param epsilon Float - Parameter that transforms travel times to commuting costs
#' @param mu Float - Floorspace prod function: output elast wrt capital, 1-mu wrt land.     
#' @param eta Float - Congestion force
#' @param nu_init Float - Convergence parameter to update wages.
#'     Default nu=0.01.
#' @param tol Int - tolerance factor
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=5000.
#'
#' @return Equilibrium values.
#' @export
#'
#' @examples
#' N=5
#' L_i = c(63, 261, 213, 182, 113)
#' L_j = c(86, 278, 189, 180, 99)
#' Q = c(2123, 1576, 1371, 1931, 1637)
#' K = c(0.44, 1.45, 1.15, 0.87, 0.58)
#' t_ij = rbind(c(0.0, 6.6, 5.5, 5.6, 6.4),
#'              c(6.7, 0.0, 3.9, 4.6, 4.4),
#'              c(5.5, 3.9, 0.0, 2.8, 3.0),
#'              c(5.6, 4.6, 2.8, 0.0, 2.7),
#'              c(6.4, 4.4, 3.0, 2.7, 0.0))
#' 
#' inversionModel(N=N,
#'                L_i=L_i,
#'                L_j=L_j,
#'                Q=Q,
#'                K=K,
#'                t_ij=t_ij)
#'        
#'        
# Parameters
alpha=0.7
beta=0.7
theta=7
delta=0.3585
rho=0.9094
lambda=0.01
epsilon=0.01
mu=0.3
eta=0.1548
nu_init=0.005
tol=10^-10
maxiter=5000        
inversionModel = function(N,
                          L_i,
                          L_j,
                          Q,
                          K,
                          t_ij,
                          alpha=0.7,
                          beta=0.7,
                          theta=7,
                          delta=0.3585,
                          rho=0.9094,
                          lambda=0.01,
                          epsilon=0.01,
                          mu=0.3,
                          eta=0.1548,
                          nu_init=0.005,
                          tol=10^-10,
                          maxiter=5000){
  
  # Formatting of input data making them all arrays
  # also helps with column and row naming
  if(is.data.frame(L_i)){
    L_i = array(unlist(L_i), dim(L_i))
  } else if(is.null(dim(L_i))){
    L_i = array(L_i, dim=c(N,1))
  }
  
  if(is.data.frame(L_j)){
    L_j = array(unlist(L_j), dim(L_j))
  } else if(is.null(dim(L_j))){
    L_j = array(L_j, dim=c(N,1))
  }
  if(is.data.frame(K)){
    K = array(unlist(K), dim(K))  
  } else if(is.null(dim(K))){
    K = array(K, dim=c(N,1))
  }
  if(is.data.frame(Q)){
    Q = array(unlist(Q), dim(Q))
  } else if(is.null(dim(Q))){
    Q = array(Q, dim=c(N,1))
  }
  t_ij = array(unlist(t_ij), dim(t_ij))  
  
  # Normalize L_i to have the same size as L_j
  L_i=L_i*sum(L_j)/sum(L_i)
  
  # Initialization
  w_init=array(1, dim=c(N,1))
  
  # Transformation of travel times to trade costs
  D = commuting_matrix(t_ij=t_ij, 
                       epsilon=epsilon)
  tau = D$tau
  
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Step 1 - Wages
  # Start with predicted number of workers
  # compute worker commuting probabilities
  # then solve adjusted wages consistent with commuting clearing\
  # i.e. wages that match model predicted Lm_j with observed Lm_j
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# tau - commuting costs
# theta - commuting elasticity
  # Finding the wages that match the data
  WI = wages_inversion(N=N,
                       w_init=w_init,
                       theta=theta,
                       tau=tau,
                       L_i=L_i,
                       L_j=L_j,
                       nu_init=nu_init,
                       tol=tol,
                       maxiter=maxiter)
  
# Output of wages_inversion
  # Equilibrium wages
  w = WI$w # wages
  w_tr = WI$w_tr # wages
  W_i = WI$W_i # Market access measure in each location.
  lambda_ij_i = WI$lambda_ij_i # Probability of individuals in each location of working in each location.

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
 # Step 1b - compute average income y_bar
  # Probability of individuals in each location of working in each location * wage in location
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#  
  
# Average income

Inc = av_income_simple(lambda_ij_i=lambda_ij_i,
                         w_tr = w_tr
  )
  y_bar = Inc$y_bar
  
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#Density of development
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
    DensD = density_development(Q=Q,
                              K=K,
                              w=w,
                              L_j=L_j,
                              y_bar=y_bar,
                              L_i=L_i,
                              beta=beta,
                              alpha=alpha,
                              mu=mu
  )
  Q_mean = DensD$Q_mean
  Q_norm = DensD$Q_norm
  FS_f = DensD$FS_f
  FS_r = DensD$FS_r
  FS = DensD$FS
  varphi = DensD$varphi
  ttheta = FS_f/FS
  
  #Productivities
  Prod = productivity(N=N,
                      Q=Q,
                      w=w,
                      L_j=L_j,
                      K=K,
                      t_ij = t_ij,
                      delta=delta,
                      lambda=lambda,
                      beta=beta
  )
  A = Prod$A
  a = Prod$a
  
  # Amenities
  AM = living_amenities_simple(theta=theta,
                               N=N,
                               L_i=L_i,
                               W_i=W_i,
                               Q=Q,
                               K=K,
                               alpha=alpha,
                               t_ij=t_ij,
                               rho=rho,
                               eta=eta
  )
  
  B = AM$B
  b = AM$b  
  
  # Save and export
  Q_alpha = Q_norm^(1-alpha)
  u = array_operator(array_operator(W_i,Q_alpha,'/'),B,'*')
  U = (sumDims(u,1))^(1/theta)
  
  return(list(A=A, a=a, u=u, B=B, b=b, w=w, varphi=varphi, U=U, Q_norm=Q_norm, ttheta=ttheta))
}




#-------------------------------------------------------
# Solving Model
#-------------------------------------------------------

#' Function to solve counterfactuals.
#'
#' @param N Integer - Number of locations.
  #' @param L_i Nx1 array - Number of residents in each location
#' @param L_j Nx1 array - Number of workers in each location
#' @param K Nx1 array - Land supply
#' @param t_ij NxN matrix - Travel times across locations
#' @param a Nx1 array - Total Factor Productivity in each location
#' @param b Nx1 array - Vector of amenities in each location
#' @param varphi Nx1 array - Density of development
#' @param w_eq Nx1 array - Initial vector of wages
#' @param u_eq Nx1 array - Initial vector of welfare
#' @param Q_eq Nx1 array - Initial price for floorspace
#' @param ttheta_eq Nx1 array - Share of floorspace used commercially 
#' @param alpha Float - Exp. share in consumption, 1-alpha exp. share in housing
#' @param beta Float - Output elasticity with respect to labor
#' @param theta Float - Commuting and migration elasticity.
#' @param mu Float - Floorspace prod function: output elasticity wrt capital
#' @param delta Float - Decay parameter agglomeration force
#' @param lambda Float - agglomeration externality
#' @param rho Float - decay parameter for amenities
#' @param eta Float - amenity externality
#' @param epsilon Float - Parameter that transforms travel times to commuting costs
#' @param zeta Float - convergence parameter
#' @param tol Int - tolerance factor
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=5000.
#' 
#' @return Counterfactual values.
#' @export
#'
#' @examples
#' N=5
#' L_i = c(63, 261, 213, 182, 113)
#' L_j = c(86, 278, 189, 180, 99)
#' Q = c(2123, 1576, 1371, 1931, 1637)
#' K = c(0.44, 1.45, 1.15, 0.87, 0.58)
#' t_ij = rbind(c(0.0, 6.6, 5.5, 5.6, 6.4),
#'              c(6.7, 0.0, 3.9, 4.6, 4.4),
#'              c(5.5, 3.9, 0.0, 2.8, 3.0),
#'              c(5.6, 4.6, 2.8, 0.0, 2.7),
#'              c(6.4, 4.4, 3.0, 2.7, 0.0))
#' 
#' a = c(1.7, 1.7, 1.6, 1.8, 1.6)
#' b = c(2.2, 2.5, 2.4, 2.6, 2.3)
#' varphi = c(95, 219, 215, 167, 148)
#' w_eq = c(0.9, 1.0, 1.0, 1.0, 0.9)
#' u_eq = c(1.0, 1.3, 1.2, 1.2, 1.1)
#' Q_eq = c(1.2, 0.9, 0.8, 1.1, 0.9)
#' ttheta_eq = c(0.5, 0.4, 0.4, 0.4, 0.4)
#' solveModel(N=N,
#'            L_i=L_i,
#'            L_j=L_j,
#'            K=K,
#'            t_ij=t_ij,
#'            a=a,
#'            b=b,
#'            varphi=varphi,
#'            w_eq=w_eq,
#'            u_eq=u_eq,
#'            Q_eq=Q_eq,
#'            ttheta_eq=ttheta_eq)
#'            
solveModel = function(N,
                      L_i,
                      L_j,
                      K,
                      t_ij,
                      a,
                      b,
                      varphi,
                      w_eq,
                      u_eq,
                      Q_eq,
                      ttheta_eq,
                      alpha=0.7,
                      beta=0.7,
                      theta=7,
                      mu=0.3,
                      delta=0.3585,
                      lambda=0.01,
                      rho=0.9094,
                      eta=0.1548,
                      epsilon=0.01,
                      zeta=0.95,
                      tol=10^-10,
                      maxiter=5000){
  
  # Formatting of input data
  if(is.data.frame(L_i)){
    L_i = array(unlist(L_i), dim(L_i))
  } else if(is.null(dim(L_i))){
    L_i = array(L_i, dim=c(N,1))
  }
  
  if(is.data.frame(L_j)){
    L_j = array(unlist(L_j), dim(L_j))
  } else if(is.null(dim(L_j))){
    L_j = array(L_j, dim=c(N,1))
  }
  if(is.data.frame(K)){
    K = array(unlist(K), dim(K))  
  } else if(is.null(dim(K))){
    K = array(K, dim=c(N,1))
  }
  
  t_ij = array(unlist(t_ij), dim(t_ij))  
  
  if(is.null(dim(a))){
    a = array(a, dim=c(N,1))
  }
  if(is.null(dim(b))){
    b = array(b, dim=c(N,1))
  }
  if(is.null(dim(varphi))){
    varphi = array(varphi, dim=c(N,1))
  }
  if(is.null(dim(w_eq))){
    w_eq = array(w_eq, dim=c(N,1))
  }
  if(is.null(dim(u_eq))){
    u_eq = array(u_eq, dim=c(N,1))
  }
  if(is.null(dim(Q_eq))){
    Q_eq = array(Q_eq, dim=c(N,1))
  }
  if(is.null(dim(ttheta_eq))){
    ttheta_eq = array(ttheta_eq, dim=c(N,1))
  }
  
  # Normalize L_i to have the same size as L_j
  L_i=L_i*sum(L_j)/sum(L_i)
  
  D = commuting_matrix(t_ij=t_ij, epsilon = epsilon)
  tau = D$tau
  L_i = array(unlist(L_i),dim(L_i))
  L_j = array(unlist(L_j),dim(L_j))
  K = array(unlist(K), dim(K))
  
  # Settings
  outerdiff = Inf;
  w = w_eq;
  u = u_eq;
  Q = Q_eq;
  ttheta = ttheta_eq
  iter = 0;
  zeta_init = zeta;
  
  cat("Solving model...\n")
  while(outerdiff>tol & iter < maxiter){
# 1) Labor supply equation
    w_tr = aperm(array(w, dim=c(N,1)), c(2,1));
    rep_w_tr = kronecker(w_tr^theta, array(1, dim=c(N, 1)));
    # Constructing employment shares
    w_tr_tau = array_operator(w_tr^theta, tau^(-theta), '*');
    lambda_ij_i = array_operator(w_tr_tau, sumDims2(w_tr_tau,2), '/');
    W_i = (sumDims2(w_tr_tau,2))^(1/theta);
    # Labor is equal to probabilities * total number of residents * proportion of workers in each sector.
    L_ij = array_operator(L_i, lambda_ij_i, '*')
    L_j = sumDims2(L_ij, 1)
    L = sum(L_i)
    lambda_i = L_i/L
    
# 2 average income
    av_income = av_income_simple(lambda_ij_i=lambda_ij_i,w_tr = w_tr)
    ybar = av_income$y_bar
    
# 3 Total floorspace
    FS = array_operator(varphi,K^(1-mu),"*")
    
# 4 Agglomeration externalities
    L_j_dens = (array_operator(L_j, K, '/'));
    L_j_dens_per = aperm(array(L_j_dens, dim=c(N,1)), c(2,1));
    L_j_dens_rep = kronecker(L_j_dens_per, array(1, dim=c(N, 1)));
    Upsilon = sumDims2(array_operator(exp(-delta*t_ij), L_j_dens_rep, '*'), 2);    
    A = array_operator(a, Upsilon^lambda, '*')
    
# 5 Amenities
    L_i_dens = (array_operator(L_i, K, '/'));
    L_i_dens_per = aperm(array(L_i_dens, dim=c(N,1)), c(2,1));
    L_i_dens_rep = kronecker(L_i_dens_per, array(1, dim=c(N, 1)));
    Omega = sumDims2(array_operator(exp(-rho*t_ij), L_i_dens_rep, '*'), 2);
    B = array_operator(b, Omega^(-eta),'*')
    
# 6 Residents, probabilities, and welfare
    u =  array_operator(array_operator(W_i, Q^(1-alpha), '/'), B, '*')
    U = sum(u^theta)
    lambda_i_upd = (u^theta)/U
    U = U^(1/theta)
    
# 7 Total output by location
    FS_f = array_operator(ttheta,array_operator(varphi, K^(1-mu), '*'), '*')
    Y = array_operator(A, array_operator(L_j^beta, FS_f^(1-beta), '*'), '*')
    Q_upd1 = (1-beta)*array_operator(Y,FS_f, '/')
    w_upd = beta*array_operator(Y, L_j, '/')
    
# 8 Housing prices
    FS_r = array_operator((1-ttheta), array_operator(varphi, K^(1-mu), '*'), '*')
    X = array_operator(ybar, L_i, '*')
    Q_upd2 = (1-alpha)*array_operator(X, FS_r, '/')
    Q_upd = Q_upd1*(a>0) + Q_upd2*(a==0 & b>0)
    
# 9 Share of commercial floorspace
    LP = array_operator(Q_upd1, array_operator(varphi, K^(1-mu), '*'), '*')
    ttheta_upd = (1-beta)*array_operator(Y, LP, '/')
    ttheta_upd = (b==0)+ttheta_upd*(b>0)
    
# 10 Calculating the main differences
    z_w = array_operator(w, w_upd, '-')
    z_L = array_operator(lambda_i, lambda_i_upd, '-')
    z_Q = array_operator(Q, Q_upd, '-')
    z_theta = array_operator(ttheta, ttheta_upd, '-')
    outerdiff = max(c(max(abs(z_w)), max(abs(z_L)), max(abs(z_Q)), max(abs(z_theta))))
    iter = iter+1
    
# 11 New vector of variables
lambda_i = zeta*lambda_i + (1-zeta)*lambda_i_upd
    Q = zeta*Q + (1-zeta)*Q_upd
    w = zeta*w + (1-zeta)*w_upd
    ttheta = zeta*ttheta + (1-zeta)*ttheta_upd
    L_i = lambda_i*L
    if(iter %% 10 == 0){
      cat(paste0("Iteration: ", iter, ", error: ", round(outerdiff, 10), ".\n"))
    }
  }
  if(outerdiff<=tol){
    cat(paste0("Converged after ", iter, " iterations. Error=", round(outerdiff, 10), ".\n"))
  } else{
    cat(paste0("Reached maximum number of iterations (", iter, "). Error=", round(outerdiff, 10), ".\n"))
  }
  
  return(list(w=w, W_i=W_i, B=B, A=A, Q=Q, lambda_ij_i=lambda_ij_i, L_i=L_i, L_j=L_j,
              ybar=ybar, lambda_i=lambda_i, ttheta=ttheta, u=u, U=U))
}



