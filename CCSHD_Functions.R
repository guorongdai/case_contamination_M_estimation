suppressMessages({
  if (!require(np)) install.packages('np', repos =
                                       'https://cran.revolutionanalytics.com/')
  
  if (!require(glmnet)) install.packages('glmnet', repos =
                                           'https://cran.revolutionanalytics.com/')
  
  if (!require(randomForest)) install.packages('randomForest', repos =
                                                 'https://cran.revolutionanalytics.com/')
  
  if (!require(xgboost)) install.packages('xgboost', repos =
                                         'https://cran.revolutionanalytics.com/')
  
  if (!require(KRLS)) install.packages('KRLS', repos =
                                         'https://cran.revolutionanalytics.com/')
  
  library(np)
  library(glmnet)
  library(randomForest)
  library(xgboost)
  library(KRLS)
  
})

# Main function calculating our estimator thetahat
# tau is the disease ratio in the case-control sample.
# eta is the disease ratio in the population.
# D is an n-dimensional vector of statuses.
# muhatZ is an N1-dimensional vector of phenotyping model estimators.
# thetahatIN is an initial estimator of theta, which is a p-dimensional vector.
# Omegahat is a Jocobian estimator, which is a p * p matrix.
# The value of 'analysis can be 'primary' or 'secondary', indicating
# the primary (logistic regression of D on X) or secondary (linear 
# regression of Y on X) should be conducted.
# The output has four components:
# $ est is our estimator thetahat.
# $ cov.new is the covariance matrix estimate of thetahat calculated with thetahat.
# $ cov is the covariance matrix estimate of thetahat calculated with thetahatIN.
CCS.M.Est = function(tau, eta, D, Y = NULL, S, X, muhatZ, thetahatIN, Omegahat, psi, analysis)
# CCS.M.Est = function(tau, eta, D, muhatZ, thetahatIN, Omegahat, psi, analysis)
{
  
  n = length(D)
  N1 = length(muhatZ)
  N = nrow(psi)
  deltan = n / N1
  alpha = ifelse(analysis == "primary", tau, eta)
  
  ########################################
  # calculate the estimator thetahat
  case.part = apply( rmultiply(psi[1 : N1, ], muhatZ), 2, mean )
  control.part = apply( psi[-(1 : N1), ], 2, mean )
  validation.part = apply( rmultiply(psi[1 : n, ], D - muhatZ[1 : n]), 2, mean )
  score = alpha * case.part + (1 - alpha) * control.part + alpha * validation.part
  
  thetahat = thetahatIN + as.vector(Omegahat %*% score) # the estimator
  ########################################
  
  ########################################
  # calculate the covariance matrix estimate of thetahat
  # Here we use thetahatIN in the function psi and the Jacobian estimate Omegahat. 
  # Changing it to thetahathat may improve the covariance matrix estimate.
  X.thetahat = as.vector(X %*% thetahat)
  if(analysis == "primary") psi.new = rmultiply( X, S - h.function(X.thetahat) )
  if(analysis == "secondary") psi.new = rmultiply(X, Y - X.thetahat)
  cov.validation = cov( rmultiply(psi.new[1 : n, ], D - (1 - deltan) * muhatZ[1 : n]) )
  cov.case = cov( rmultiply(psi.new[1 : N1, ], muhatZ) )
  cov.control = cov(psi.new[-(1 : N1), ])
  An = alpha ^ 2 * cov.validation + alpha ^ 2 * deltan * (1 - deltan) * cov.case + 
       (1 - alpha) ^ 2 * tau * deltan * cov.control / (1 - tau)
  cov.thetahat.new = Omegahat %*% An %*% Omegahat / n
  
  cov.validation = cov( rmultiply(psi[1 : n, ], D - (1 - deltan) * muhatZ[1 : n]) )
  cov.case = cov( rmultiply(psi[1 : N1, ], muhatZ) )
  cov.control = cov(psi[-(1 : N1), ])
  
  An = alpha ^ 2 * cov.validation + alpha ^ 2 * deltan * (1 - deltan) * cov.case + 
       (1 - alpha) ^ 2 * tau * deltan * cov.control / (1 - tau)
  cov.thetahat = Omegahat %*% An %*% Omegahat / n
  ########################################
  
  result = list("est" = thetahat, "cov.new" = cov.thetahat.new, "cov" = cov.thetahat)
  
  return(result)
  
}

# preliminary function calculating the initial estimator along with its covariance matrix estimate, 
# the Jacobian estimator and the values of the psi function at the initial estimator.
# tau is the disease ratio in the case-control sample.
# eta is the disease ratio in the population.
# D is an n-dimensional vector of statuses.
# Y is an N-dimensional vector of responses.
# S is an N-dimensional vector, indicating being a candidate case or a control.
# X is an N * p matrix of predictors. The first column of X should be a vector of ones to capture the intercept.
# The value of 'analysis can be 'primary' or 'secondary', indicating
# the primary (logistic regression of D on X) or secondary (linear 
# regression of Y on X) should be conducted.
# The output has four components:
# $ initial.est is the initial estimator thetahatIN.
# $ initial.cov is the covariance matrix estimate of thetahatIN.
# $ Jacobian is the Jacobian estimator Omegahat
# $ psi is the values of the psi function at the initial estimator.
prelim = function(tau, eta, D, Y = NULL, S, X, analysis)
{
  
  n = length(D)
  N1 = sum(S == 1)
  N = length(S)
  N0 = N - N1
  deltan = n / N1
  
  if(analysis == "primary") 
  {
    
    alpha = tau
    regression = "binomial"
    response = S
    
  }
  
  if(analysis == "secondary") 
  {
    
    alpha = eta 
    regression = "gaussian"
    response = Y
    
  }
  
  ########################################
  # calculate the initial estimator
  complete = c(1 : n, (N1 + 1) : N) 
  w = numeric(n + N0)
  w[1 : n] = alpha * D / n
  w[-(1 : n)] = (1 - alpha) / N0
  thetahatIN = coef( suppressWarnings( glm(response[complete] ~ -1 + X[complete, ], family = regression, weights = w) ) )
  ########################################
  
  X.thetahatIN = as.vector(X %*% thetahatIN)
  
  if(analysis == "primary") 
  {
    
    psi = rmultiply(X, response - h.function(X.thetahatIN))
    
    w.Omega.V = D * hprime.function(X.thetahatIN[1 : n])
    w.Omega.C = hprime.function(X.thetahatIN[-(1 : N1)])
    
  }
  
  if(analysis == "secondary") 
  {

    psi = rmultiply(X, response - X.thetahatIN)
    
    w.Omega.V = D
    w.Omega.C = 1
    
  }
  
  ########################################
  # calculate the Jacobian estimator
  X.V = X[1 : n, ]
  X.C = X[-(1 : N1), ]
  Omegahat = solve( alpha * crossprod(X.V, rmultiply(X.V, w.Omega.V)) / n + 
                    (1 - alpha) * crossprod(X.C, rmultiply(X.C, w.Omega.C)) / N0 )
  ########################################
  
  ########################################
  # calculate the covariance matrix estimate of thetahatIN
  covIN.validation = cov( rmultiply(psi[1 : n, ], D) )
  covIN.control = cov(psi[-(1 : N1), ])
  An.In = alpha ^ 2 * covIN.validation + (1 - alpha) ^ 2 * tau * deltan * covIN.control / (1 - tau)
  cov.thetahatIN = Omegahat %*% An.In %*% Omegahat / n
  ########################################
  
  result = list("initial.est" = thetahatIN, "Jacobian" = Omegahat, "initial.cov" = cov.thetahatIN, "psi" = psi)
  
  return(result)
  
}

# calculate the phenotyping model estimator
# D is an n-dimensional vector of individuals' status.
# Z.V is an n * d matrix of the predictors Z in the validation set.
# Z.N is an (N1 - n) * d matrix of the predictors Z in the nonvalidated case pool.
# Z.V and Z.N do not need to have a column of ones to caputure intercepts.
# M is the number of folds in cross fitting.
# "LR" is logistic regression; "KS" is kernel smoothing on the direction selected by logistic regression;
# "RF" is random forest; "XGB" is XGBoost; "KMR" is kernel machine regression.
# The output is an N1-dimensional vector of the estimated nuisance function at (Z.V ^ T, Z.N ^ T) ^ T.
muhat.function = function(D, Z.V, Z.N, M, method)
{
  
  n = length(D)
  nM = floor(n / M)
  
  result.V = numeric(n)
  result.N = numeric( nrow(Z.N) )
  
  for(m in 1 : (M - 1))
  {
    
    ind = (m - 1) * nM + (1 : nM)
    which.fit = 1 : length(ind)
    
    D.fit = D[-ind]
    Z.fit = Z.V[-ind, ]
    Z.predict = rbind(Z.V[ind, ], Z.N)
    
    NF.fit = NF.Est(D.fit, Z.fit, Z.predict, method)
    
    result.V[ind] = NF.fit[which.fit]
    result.N = result.N + NF.fit[-which.fit] / M
    
  }
  
  ind = ( ((M - 1) * nM) + 1 ) : n
  which.fit = 1 : length(ind)
  
  D.fit = D[-ind]
  Z.fit = Z.V[-ind, ]
  Z.predict = rbind(Z.V[ind, ], Z.N)
  
  NF.fit = NF.Est(D.fit, Z.fit, Z.predict, method)
  
  result.V[ind] = NF.fit[which.fit]
  result.N = result.N + NF.fit[-which.fit] / M
  
  return(c(result.V, result.N))
  
}

# method specifies the method used to fit the phenotyping model:
# "LR" is logistic regression; "KS" is kernel smoothing on the direction selected by logistic regression;
# "RF" is random forest; "XGB" is XGBoost; "KMR" is kernel machine regression.
NF.Est = function(D.fit, Z.fit, Z.predict, method)
{
  
  if(method == "LR")
  {
    
    model = cv.glmnet(Z.fit, D.fit, family = "binomial")
    result = predict(model, newx =  Z.predict, s = "lambda.min", type = "response")
    
  }
    
  if(method == "KS")
  {
    
    DR.model = cv.glmnet(Z.fit, D.fit, family = "binomial")
    
    trans.fit = predict(DR.model, newx = Z.fit, s = "lambda.min", type = "link")
    
    if(length(unique(trans.fit)) == 1)
    {
      
      result = rep(mean(D.fit), nrow(Z.predict))
      
    } else {
      
      trans.predict = predict(DR.model, newx = Z.predict, s = "lambda.min", type = "link")
      
      sink("aux")
      h = npcdensbw(xdat = trans.fit, ydat = as.factor(D.fit), nmulti = 1) $ xbw
      result = npreg(txdat = trans.fit, tydat = D.fit, bws = h, exdat = trans.predict) $ mean
      sink(NULL)
      
    }
    
  }
    
  if(method == "RF")
  {
    
    model = randomForest(x = Z.fit, y = as.factor(D.fit))
    result = predict(model, Z.predict, type = "prob")[, 2]
    
  }
  
  if(method == "XGB")
  {
    
    sink("aux")
    model = xgboost(data = Z.fit, label = D.fit, nrounds = 5, objective = "binary:logistic", verbose = 0)
    sink(NULL)
    result = predict(model, Z.predict)
  }
  
  if(method == "KMR")
  {
    
    model = suppressWarnings(krls(X = Z.fit, y = D.fit, derivative = F, vcov = F, print.level = 0))
    result = as.vector(predict(model, Z.predict) $ fit)
    
  }
  
  return(result)  
  
}

truth = function(tau, eta, D.true, Y.true = NULL, S.true, X.true, analysis)
{
  
  N1 = sum(S.true == 1)
  N = length(S.true)
  N0 = N - N1
  
  if(analysis == "primary") 
  {
    
    alpha = tau
    regression = "binomial"
    response = S.true
    
  }
  
  if(analysis == "secondary") 
  {
    
    alpha = eta 
    regression = "gaussian"
    response = Y.true
    
  }
  
  ########################################
  # calculate the initial estimator
  w = numeric(N)
  ind.one = which(S.true == 1)
  ind.zero = which(S.true == 0)
  w[ind.one] = alpha * D.true[ind.one] / N1
  w[ind.zero] = (1 - alpha) / N0
  theta0 = coef( suppressWarnings( glm(response ~ -1 + X.true, family = regression, weights = w) ) )
  ########################################
  
  return(theta0)
  
}

# calculate E(D | Z, S = 1)
D.prob = function(Z, q, phenotyping)
{
  
  Zq = Z[, 1 : q]
  
  if(phenotyping == 1) result = rep(7 / 10, nrow(Z)) # constant model
  
  if(phenotyping == 2) result = h.function(5 * rowSums(Zq) / sqrt(q)) # linear model
  
  if(phenotyping == 3) result = h.function(3 * rowSums(Zq) / sqrt(q) + (3 / 2) * rowSums(Zq) ^ 2 / q) # single-index model
  
  # if(phenotyping == 3) result = h.function(3 * rowSums(Zq) / sqrt(q) + (3 / 2) * (rowSums(Zq) ^ 2 + rowSums(Zq) ^ 3) / q)
  
  # if(phenotyping == 3) result = h.function(3 * rowSums(Zq) / sqrt(q) + (3 / 2) * (rowSums(Zq) ^ 2 / q + rowSums(Zq) ^ 3 / sqrt(q)))
  
  if(phenotyping == 4) # multiple-index model
  {
    
    q.half = ceiling(q / 2)
    # omega = numeric(q)
    kappa = numeric(q)
    # omega[1 : q.half] = 1 # omega[1 : q.half] = 1. Other components of omega are zero.
    suppressWarnings( {kappa[1 : q] = c(1, 0)} ) # kappa[1 : q] = c(1, 0, 1, 0, ...). Other components of kappa are zero
   
    result = h.function( 3 * rowSums(Zq) / sqrt(q) * ( 1  + rowSums(Zq[, 1 : q.half] / sqrt(q)) ) - 
                         3 * as.vector((Zq %*% kappa)) ^ 2 / q ) 
    
  }
  
  if(phenotyping == 5) result = h.function(3 * rowSums(Zq) / sqrt(q) + (3 / 2) * rowSums(Zq ^ 2) / q) # additive model
  
  # if(phenotyping == 5) result = h.function(3 * rowSums(Zq) / sqrt(q) + (3 / 2) * ( rowSums(Zq ^ 2) + rowSums(Zq ^ 3) ) / q) # additive model
  
  # if(phenotyping == 5) result = pnorm(3 * rowSums(Zq) / sqrt(q) + (3 / 2) * ( rowSums(Zq ^ 2) / q + rowSums(Zq ^ 3) / sqrt(q) ))
  
  return(result)
  
}

# return a vector whose jth component indicates whether theta0[j] is in the interval (CI[j, 1], CI[j, 2])
cover = function(theta0, CI) { (theta0 >= CI[, 1]) & (theta0 <= CI[, 2]) }

# calculate diag(v) %*% M
rmultiply = function(M, v) { sweep(M, 1, v, "*") }

# add a column of ones to the left of a matrix
add.one = function(X) { cbind(rep(1, nrow(X)), X) }

# standard logistic function
h.function = function(x) { 1 / ( 1 + exp(-x) ) }

# derivative of the standard logistic function
hprime.function = function(x) { exp(x) / (1 + exp(x)) ^ 2 }
