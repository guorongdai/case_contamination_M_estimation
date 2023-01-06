# This is the code to generate simulation results.

suppressMessages({
  
  if (!require(MASS)) install.packages('MASS', repos =
                                         'https://cran.revolutionanalytics.com/')
  if (!require(doParallel)) install.packages('parallel', repos =
                                               'https://cran.revolutionanalytics.com/')
  
  library(MASS)
  library(doParallel)
  
})

registerDoParallel(48) # set multicores

setwd("/Users/daiguorong/Dropbox/CCSHD/Uploaded_Files")
source("CCSHD_Functions.R")

n = 200
N = 25000
d = 500
p = 11 # dimension of X_{-1} (If the intercept term is included, X is (p + 1)-dimensional.)
# q = ceiling(sqrt(d)) # sparsity level of the phenotyping model
q = n
analysis = "secondary"
method = "RF"
phenotyping = 2

M = 10
tau = 2 / 5
N1 = N * tau
N0 = N - N1
deltan = n / N1

N.raw = 2 * N
N.true = 50000

#####################################################
# set the mean and covariance matrix of Z
mean.Z = numeric(d)

a = 0.5
cov.Z = matrix(0, d, d)
for(i in 1 : d) for(j in 1 : d) cov.Z[i, j] = a ^ (abs(i - j))
##################################################### 

##################################################### 
# approximate the true value of theta0
set.seed(618)
Z.true = mvrnorm(N.true, mean.Z, cov.Z)
S.true = rbinom( N.true, 1, h.function(2 * rowSums(Z.true) / sqrt(d)) )
X.true = cbind(rep(1, N.true), Z.true[, 1 : p])
Y.true = Z.true[, p + 1]
eta = mean(S.true)

set.seed(618)
muZ = D.prob(Z.true, q, phenotyping)
D.true = rbinom(N.true, 1, muZ)

theta0 = truth(tau, eta, D.true, Y.true, S.true, X.true, analysis)
##################################################### 

##################################################### 
# approximate the optimal efficiency relative to the initial estimator
N1.true = sum(S.true == 1)
N0.true = N.true - N1.true

X.theta0 = as.vector(X.true %*% theta0)

ind = which(S.true == 1)

if(analysis == "primary") 
{
  
  alpha = tau
  psi = rmultiply(X.true, S.true - h.function(X.theta0))
  
  w.Omega.V = D.true[ind] * hprime.function(X.theta0[ind])
  w.Omega.C = hprime.function(X.theta0[-ind])
  
}
  
  
if(analysis == "secondary") 
{
  
  alpha = eta
  psi = rmultiply(X.true, Y.true - X.theta0)
  
  w.Omega.V = D.true[ind]
  w.Omega.C = 1
  
}

cov.validation = cov( rmultiply(psi[ind, ], D.true[ind] - (1 - deltan) * muZ[ind]) )
cov.case = cov( rmultiply(psi[ind, ], muZ[ind]) )
cov.control = cov(psi[-ind, ])

Omega = solve( alpha * crossprod( X.true[ind, ], rmultiply(X.true[ind, ], w.Omega.V) ) / N1.true + 
               (1 - alpha) * crossprod( X.true[-ind, ], rmultiply(X.true[-ind, ], w.Omega.C) ) / N0.true )

A = alpha ^ 2 * cov.validation + alpha ^ 2 * deltan * (1 - deltan) * cov.case + 
    (1 - alpha) ^ 2 * tau * deltan * cov.control / (1 - tau)
cov0 = Omega %*% A %*% Omega / n

covIN.validation = cov( rmultiply(psi[ind, ], D.true[ind]) )
A.IN = alpha ^ 2 * covIN.validation + (1 - alpha) ^ 2 * tau * deltan * cov.control / (1 - tau)
cov.initial = Omega %*% A.IN %*% Omega / n

opt.re = sum(diag(cov.initial)) / sum(diag(cov0))
# opt.re
##################################################### 

r = 500

output = foreach(i = 1 : r, .combine = rbind) %dopar%
{
  
  set.seed(512 * i)
  
  Z.raw = mvrnorm(N.raw, mean.Z, cov.Z)
  S.prob = h.function(2 * rowSums(Z.raw) / sqrt(d))
  S.raw = rbinom(N.raw, 1, S.prob)
  
  ind.case = which(S.raw == 1)[1 : N1]
  ind.control = which(S.raw == 0)[1 : N0]
  Z.case = Z.raw[ind.case, ]
  Z.C = Z.raw[ind.control, ]
  Z.V = Z.case[1 : n, ]
  Z.N = Z.case[-(1 : n), ]
  
  # D = rbinom( n, 1, h.function(5 * rowSums(Z.V)) )
  D = rbinom( n, 1, D.prob(Z.V, q, phenotyping) ) # The larger the variance of E(D | Z, S = 1), the better the result.
  S = c(rep(1, N1), rep(0, N0))
  Z = rbind(Z.V, Z.N, Z.C)
  X = cbind(rep(1, N), Z[, 1 : p])
  Y = Z[, p + 1]
  
  initial.fit = prelim(tau, eta, D, Y, S, X, analysis)
  thetahatIN = initial.fit $ initial.est
  Omegahat = initial.fit $ Jacobian
  psi = initial.fit $ psi
  cov.thetahatIN = initial.fit $ initial.cov
  
  std.initial = sqrt(diag(cov.thetahatIN))
  CI.initial = cbind(thetahatIN - 1.96 * std.initial, thetahatIN + 1.96 * std.initial) 
  # 95% confidence intervals of the coefficients based on thetahatIN, recorded in a p * 2 matrix
  ACIL.initial = mean(CI.initial[, 2] - CI.initial[, 1])
  # average length
  coverage.initial = cover(theta0, CI.initial)
  # whether theta0 is in the confidence interval
  
  muhatZ = muhat.function(D, Z.V, Z.N, M, method)
  
  result = CCS.M.Est(tau, eta, D, Y, S, X, muhatZ, thetahatIN, Omegahat, psi, analysis)
  thetahat = result $ est
  cov.thetahat.new = result $ cov.new
  cov.thetahat = result $ cov
  
  std.new = sqrt(diag(cov.thetahat.new))
  CI.new = cbind(thetahat - 1.96 * std.new, thetahat + 1.96 * std.new) 
  # 95% confidence intervals of the coefficients based on thetahat, recorded in a p * 2 matrix
  ACIL.new = mean(CI.new[, 2] - CI.new[, 1])
  # average length
  coverage.new = cover(theta0, CI.new)
  # whether theta0 is in the confidence interval
  
  std = sqrt(diag(cov.thetahat))
  CI = cbind(thetahat - 1.96 * std, thetahat + 1.96 * std) 
  # 95% confidence intervals of the coefficients based on thetahat, recorded in a p * 2 matrix
  ACIL = mean(CI[, 2] - CI[, 1])
  # average length
  coverage = cover(theta0, CI)
  # whether theta0 is in the confidence interval
  
  MSE.thetahat = sum((thetahat - theta0) ^ 2) # MSE of thetahat
  MSE.thetahatIN = sum((thetahatIN - theta0) ^ 2) # MSE of thetahatIN
  
  c(coverage.new, coverage, coverage.initial, ACIL.new, ACIL, ACIL.initial, MSE.thetahat, MSE.thetahatIN)
  
}

registerDoSEQ()

result = apply(output, 2, mean)

CR.new = sum( abs(result[1 : (p + 1)] - 0.95) ) / (p + 1)

CR = sum( abs(result[(p + 2) : (2 * p + 2)] - 0.95) ) / (p + 1)

CR.initial = sum( abs(result[(2 * p + 3) : (3 * p + 3)] - 0.95) ) / (p + 1)

ACIL.new = result[3 * p + 4]
ACIL = result[3 * p + 5]
ACIL.initial = result[3 * p + 6]

re = result[3 * p + 8] / result[3 * p + 7]

outcome = c(CR.new, ACIL.new, CR, ACIL, CR.initial, ACIL.initial, re, opt.re)

names(outcome) = c("CR.thetahat.new", "ACIL.thetahat.new", 
                   "CR.thetahat", "ACIL.thetahat",
                   "CR.thetahatIN", "ACIL.thetahatIN", 
                   "RE", "opt.RE")

outcome





